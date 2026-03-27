from __future__ import annotations

import streamlit as st
import pandas as pd
import json
import sys
import time
from pathlib import Path
from datetime import datetime, timezone

# Add local src and root src to sys.path
sys.path.append(str(Path(__file__).parent / "src"))
sys.path.append(str(Path(__file__).parent.parent / "src"))

from iwa_rnaseq_counter.legacy.config import get_default_session_state
from iwa_rnaseq_counter.legacy.fastq_discovery import collect_fastq_metadata, discover_fastq_files
from iwa_rnaseq_counter.legacy.sample_parser import (
    build_sample_table,
    detect_lane_groups,
    group_fastq_by_sample,
    infer_sample_layout,
    parse_sample_sheet,
)
from iwa_rnaseq_counter.legacy.strandedness import infer_strandedness
from iwa_rnaseq_counter.legacy.validators import (
    validate_analysis_name,
    validate_input_directory,
    validate_output_directory,
    validate_run_conditions,
    validate_salmon_index,
    validate_tx2gene_file,
)
from ui.sections import (
    render_analysis_section,
    render_app_header,
    render_fastq_section,
    render_reference_section,
    render_result_section,
    render_run_section,
    render_sample_section,
)
from iwa_job_runner.models.job_summary import discover_jobs, load_job_summary
from iwa_job_runner.core.run_artifacts import setup_run_directory
from iwa_job_runner.core.executor import LocalJobExecutor
from iwa_job_runner.models.job_request import JobRequestSpec
from iwa_job_runner.core.monitor import JobMonitor

def init_session_state() -> None:
    defaults = get_default_session_state()
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value
            
    # Migration: If existing session has "results", force it to "output" for v0.4.2+
    if st.session_state.get("output_dir") == "results":
        st.session_state.output_dir = "output"
    if st.session_state.get("runs_root") == "results":
        st.session_state.runs_root = "output"

def run_app() -> None:
    st.set_page_config(page_title="iwa-rnaseq-counter", layout="wide")
    render_app_header()

    # 1. Background Job Discovery & Recovery
    with st.sidebar:
        st.markdown("### 📂 Data Location")
        root_input = st.text_input("Runs Root (History)", value=st.session_state.runs_root, help="Directory to scan for previous analysis runs.")
        if root_input != st.session_state.runs_root:
            st.session_state.runs_root = root_input
            st.rerun()
            
        if st.button("🔄 Refresh Jobs", use_container_width=True):
            st.rerun()

    output_base = Path(st.session_state.runs_root)
    if not output_base.exists():
        st.sidebar.error(f"Root not found: {st.session_state.runs_root}")
        recent_summaries = []
    else:
        recent_summaries = discover_jobs(output_base)
    
    # Session Recovery: If a job is running and nothing is selected, auto-select it
    if st.session_state.selected_job_id is None:
        running_jobs = [s for s in recent_summaries if s.status in ("queued", "running")]
        if running_jobs:
            st.session_state.selected_job_id = running_jobs[0].job_id

    # 2. Sidebar Job List
    with st.sidebar:
        st.markdown("---")
        if st.button("➕ New Analysis", use_container_width=True):
            st.session_state.selected_job_id = None
            # Clear run status to avoid showing old results in New Analysis mode
            st.session_state.run_status = "idle"
            st.rerun()

        st.markdown("### 🕒 Job History")
        if not recent_summaries:
            st.caption("No recent jobs found.")
        else:
            for summary in recent_summaries:
                status_icon = "⏳" if summary.status in ("queued", "running") else "✅" if summary.status == "completed" else "❌"
                label = f"{status_icon} {summary.run_name}"
                
                # Highlight if selected
                is_selected = (st.session_state.selected_job_id == summary.job_id)
                if st.button(label, key=f"job_{summary.job_id}", help=f"ID: {summary.job_id}", use_container_width=True, type="secondary" if not is_selected else "primary"):
                    st.session_state.selected_job_id = summary.job_id
                    st.rerun()

        st.markdown("---")
        st.markdown("### メンテナンス")
        if st.button("❌ Session State Clear"):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.rerun()

    # Shared Validation results (initially empty)
    name_validation = input_validation = output_validation = index_validation = tx2gene_validation = run_validation = {}

    # 3. Main View Mode Selection
    if st.session_state.selected_job_id is None:
        # ==========================================
        # INPUT MODE
        # ==========================================
        analysis_values = render_analysis_section()
        st.session_state.analysis_name = analysis_values["analysis_name"]
        st.session_state.input_dir = analysis_values["input_dir"]
        st.session_state.output_dir = analysis_values["output_dir"]

        name_validation = validate_analysis_name(st.session_state.analysis_name)
        input_validation = validate_input_directory(st.session_state.input_dir)
        output_validation = validate_output_directory(st.session_state.output_dir)

        # FASTQ Discovery logic
        if analysis_values.get("scan_requested") and input_validation["is_valid"]:
            fastq_files = discover_fastq_files(st.session_state.input_dir)
            st.session_state.scanned_fastq_files = fastq_files
            st.session_state.fastq_df = collect_fastq_metadata(fastq_files)

            sample_groups = group_fastq_by_sample(fastq_files)
            st.session_state.detected_samples = sample_groups
            layout_map = infer_sample_layout(sample_groups)
            lane_map = detect_lane_groups(sample_groups)
            st.session_state.sample_df = build_sample_table(sample_groups, layout_map, lane_map)
            
            st.session_state.fastq_scan_input_dir = st.session_state.input_dir
            st.session_state.fastq_scan_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            st.session_state.fastq_scan_file_count = len(fastq_files)
            st.session_state.fastq_scan_completed = True
            st.success("FASTQ をスキャンしました。")

        # CSV Load logic
        if analysis_values.get("load_csv_requested") and analysis_values.get("csv_path"):
            try:
                csv_df = parse_sample_sheet(analysis_values["csv_path"], st.session_state.input_dir)
                st.session_state.sample_df = csv_df
                st.success(f"CSV から {len(csv_df)} サンプルを読み込みました。")
            except Exception as e:
                st.error(f"CSV 読み込みエラー: {e}")

        # Render Sections
        render_fastq_section(st.session_state.get("fastq_df"))

        st.session_state.sample_df = render_sample_section(
            st.session_state.sample_df if st.session_state.sample_df is not None else pd.DataFrame()
        )

        reference_values = render_reference_section()
        st.session_state.salmon_index_path = reference_values["salmon_index_path"]
        st.session_state.tx2gene_path = reference_values["tx2gene_path"]

        index_validation = validate_salmon_index(st.session_state.salmon_index_path)
        tx2gene_validation = validate_tx2gene_file(st.session_state.tx2gene_path)

        if st.session_state.get("last_salmon_index_path") != st.session_state.salmon_index_path:
            st.session_state.strandedness_prediction = None
            st.session_state.last_salmon_index_path = st.session_state.salmon_index_path

        if reference_values.get("estimate_strandedness"):
            if st.session_state.get("sample_df") is not None and not st.session_state.sample_df.empty and index_validation["is_valid"]:
                with st.spinner("Strandedness を推定中..."):
                    strandedness_result = infer_strandedness(
                        sample_df=st.session_state.sample_df,
                        input_dir=st.session_state.input_dir,
                        salmon_index_path=st.session_state.salmon_index_path,
                    )
                    st.session_state.strandedness_prediction = strandedness_result
            else:
                st.warning("サンプル表が空、または Salmon Index パスが無効なため、推定できません。")

        run_validation = validate_run_conditions(
            input_dir=st.session_state.input_dir,
            output_dir=st.session_state.output_dir,
            sample_df=st.session_state.get("sample_df"),
            salmon_index_path=st.session_state.salmon_index_path,
            tx2gene_path=st.session_state.tx2gene_path,
            strandedness_mode=st.session_state.get("strandedness_mode", "Auto-detect"),
            strandedness_result=st.session_state.get("strandedness_prediction"),
        )
        st.session_state.run_validation_result = run_validation

        run_values = render_run_section(
            strandedness_result=st.session_state.get("strandedness_prediction"),
            run_validation=st.session_state.get("run_validation_result"),
        )
        st.session_state.strandedness_mode = run_values["strandedness_mode"]
        st.session_state.threads = run_values["threads"]

        # Handle Run Execution
        if run_values["run_requested"]:
            started_at_iso = datetime.now(timezone.utc).astimezone().isoformat()
            job_id = f"RUN_GUI_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            
            dirs = setup_run_directory(Path(st.session_state.output_dir), job_id)
            
            config_data = {
                "analysis_name": st.session_state.analysis_name,
                "input_dir": str(Path(st.session_state.input_dir).resolve()),
                "output_dir": str(Path(st.session_state.output_dir).resolve()),
                "salmon_index_path": str(Path(st.session_state.salmon_index_path).resolve()),
                "tx2gene_path": str(Path(st.session_state.tx2gene_path).resolve()),
                "strandedness_mode": st.session_state.strandedness_mode,
                "threads": st.session_state.threads,
                "strandedness_prediction": st.session_state.strandedness_prediction
            }
            config_path = (dirs.inputs / "run_config.json").resolve()
            with open(config_path, "w") as f:
                json.dump(config_data, f, indent=2)
                
            def resolve_path_list(p_list):
                if not isinstance(p_list, (list, tuple)):
                    if isinstance(p_list, str) and p_list.strip():
                        return str(Path(p_list).resolve())
                    return p_list
                return [str(Path(p).resolve()) for p in p_list]

            backend_df = st.session_state.sample_df.copy()
            for col in ["r1_paths", "r2_paths", "all_paths"]:
                if col in backend_df.columns:
                    backend_df[col] = backend_df[col].apply(resolve_path_list)
            
            sample_sheet_path = (dirs.inputs / "sample_sheet.csv").resolve()
            backend_df.to_csv(sample_sheet_path, index=False)
            
            command = [
                "pixi", "run", "python", str(Path(__file__).parent.resolve() / "cli.py"),
                "run-gui-backend",
                "--config", str(config_path),
                "--sample-sheet", str(sample_sheet_path),
                "--outdir", str(dirs.root.resolve()),
                "--started-at", started_at_iso
            ]
            
            request = JobRequestSpec(
                job_id=job_id,
                requested_app="iwa_rnaseq_counter_gui_backend",
                command=command,
                parameters={
                    "config": str(config_path),
                    "sample_sheet": str(sample_sheet_path),
                    "started_at": started_at_iso
                },
                resources={"threads": st.session_state.threads}
            )
            
            executor = LocalJobExecutor(Path(st.session_state.output_dir))
            executor.submit(request)
            
            st.session_state.selected_job_id = job_id
            st.rerun()

    else:
        # ==========================================
        # JOB VIEW MODE
        # ==========================================
        # ==========================================
        # JOB VIEW MODE
        # ==========================================
        current_job = next((s for s in recent_summaries if s.job_id == st.session_state.selected_job_id), None)
        
        if not current_job:
            st.error(f"Selected job {st.session_state.selected_job_id} not found in {st.session_state.runs_root}")
            if st.button("Back to New Analysis"):
                st.session_state.selected_job_id = None
                st.rerun()
            return

        col1, col2, col3 = st.columns([2, 1, 1])
        with col1:
            st.subheader(f"Job: {current_job.run_name}")
            st.caption(f"ID: `{current_job.job_id}`")
        with col2:
            status_color = "green" if current_job.status == "completed" else "orange" if current_job.status in ("queued", "running") else "red"
            st.markdown(f"**Status:** :{status_color}[{current_job.status.upper()}]")
            st.markdown(f"**Started:** {current_job.started_at or 'N/A'}")
        with col3:
            st.markdown(f"**Samples:** {current_job.sample_count}")
            st.markdown(f"**Elapsed:** {current_job.elapsed_seconds / 60:.1f} min")
        
        st.markdown("---")

        # Monitoring uses the parent of the job's output_dir as the base
        monitor = JobMonitor(Path(current_job.output_dir).parent)
        spec = monitor.get_status(current_job.job_id)
        
        if spec:
            if spec.status in ("queued", "running"):
                with st.status(f"Job '{spec.run_id}' is running...", expanded=True):
                    st.write(f"Status: `{spec.status}`")
                    st.write("You can safely leave this page or check other jobs. The analysis continues in the background.")
                time.sleep(2)
                st.rerun()
                
            elif spec.status == "completed":
                run_dir = Path(current_job.output_dir)
                
                # Load results locally for rendering
                job_run_summary = None
                job_output_files = None
                
                try:
                    with open(run_dir / "results" / "run_summary.json", "r") as f:
                        job_run_summary = json.load(f)
                except Exception:
                    pass
                    
                from iwa_rnaseq_counter.legacy.run_artifacts import build_output_manifest
                try:
                    # Use run.log which is created by cli.py
                    job_output_files = build_output_manifest(
                        run_dir, {}, None, None, run_dir / "logs" / "run.log"
                    )
                except Exception:
                    pass

                render_result_section(
                    run_status="success",
                    output_files=job_output_files,
                    log_summary=None,
                    run_summary=job_run_summary
                )

            elif spec.status == "failed":
                st.error(f"Job {spec.run_id} failed!")
                st.info(f"Check logs at: `{Path(current_job.output_dir) / 'logs' / 'run.log'}`")
        else:
            st.warning("Could not retrieve job status.")

    # Shared Footer
    with st.expander("Validation details"):
        st.write({
            "name_validation": name_validation,
            "input_validation": input_validation,
            "output_validation": output_validation,
            "index_validation": index_validation,
            "tx2gene_validation": tx2gene_validation,
            "run_validation": run_validation,
        })

def main() -> None:
    init_session_state()
    run_app()

if __name__ == "__main__":
    main()

from __future__ import annotations

import streamlit as st
import pandas as pd
import json
import sys
import time
from pathlib import Path
from datetime import datetime, timezone

# ローカルのsrcとルートのsrcをsys.pathに追加
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
# [v0.6.0 C-04]
# validation 層が validate_salmon_index / validate_tx2gene_file という
# Salmon 前提の関数名・契約を UI へ直接露出している。
# v0.6.0 では app.py が backend 固有 validator 名を直接知らず、
# reference/index validator の抽象名へ寄せたい。
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
    render_job_history_panel,
    render_reference_section,
    render_result_section,
    render_run_section,
    render_sample_section,
    format_elapsed_time,
)
from iwa_job_runner.models.run_discovery import resolve_runs_root, discover_runs
from iwa_job_runner.models.run_selection import pick_active_run, sort_runs_for_sidebar
from iwa_job_runner.core.run_artifacts import setup_run_directory
from iwa_job_runner.core.executor import LocalJobExecutor
from iwa_job_runner.models.job_request import JobRequestSpec
from iwa_job_runner.core.monitor import JobMonitor

def init_session_state() -> None:
    defaults = get_default_session_state()
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value
            
    # Migration: 既存のセッションに「結果」がある場合は、v0.4.2以降では強制的に「出力」にします
    if st.session_state.get("output_dir") == "results":
        st.session_state.output_dir = "output"
    if st.session_state.get("runs_root") == "results":
        st.session_state.runs_root = "output"
    
    # Backward compatibility: 旧セッションで recovery フラグが無ければ初回復元を有効化
    if "needs_initial_recovery" not in st.session_state:
        st.session_state.needs_initial_recovery = True

def run_app() -> None:
    st.set_page_config(page_title="iwa-rnaseq-counter", layout="wide")
    render_app_header()
    #
    # 1. Runs Root & Job Discovery
    #
    with st.sidebar:
        st.markdown("### 📂 Data Location")
        root_default = st.session_state.get("runs_root", "output")
        root_input = st.text_input("Runs Root (History)", value=root_default, help="Directory to scan for previous analysis runs.")
        
        if root_input != st.session_state.get("runs_root"):
            st.session_state.runs_root = root_input
            st.session_state.selected_job_id = None 
            st.session_state.needs_initial_recovery = True # Trigger re-recovery on root change
            st.rerun()
            
        if st.button("🔄 Refresh Jobs", use_container_width=True):
            st.rerun()

    # Resolve and discover
    runs_root_path = resolve_runs_root(st.session_state.runs_root)
    if not runs_root_path.exists():
        st.sidebar.error(f"Root not found: {st.session_state.runs_root}")
        recent_runs = []
    else:
        recent_runs = discover_runs(runs_root_path)
    #
    # 2. Stable Session Recovery / Selection
    # 真のソースは Runs Root 配下の run artifact
    #
    if st.session_state.needs_initial_recovery:
        # 初回ロード回復またはルート変更回復
        suggestion = pick_active_run(recent_runs, None)
        st.session_state.selected_job_id = suggestion.run_dir.name if suggestion else None
        st.session_state.needs_initial_recovery = False
        
    elif st.session_state.selected_job_id is not None:
        # ユーザーは既にジョブビューモードになっています。
        # 現在の選択内容を解決してください（バックグラウンドプロセスによって更新されている可能性があります）。
        hint_id = st.session_state.selected_job_id
        hint_path = str(runs_root_path / hint_id)
        
        suggestion = pick_active_run(recent_runs, hint_path)
        if suggestion:
            new_id = suggestion.run_dir.name
            if st.session_state.selected_job_id != new_id:
                st.session_state.selected_job_id = new_id
        else:
            # stale selection (e.g. current job deleted or moved out of root)
            # 候補が取れないときは、自動的に New Analysis に落とす
            st.session_state.selected_job_id = None
    else:
        # Explicit None (New Analysis) -> 新規分析モードを維持する
        pass
    #
    # 3. View Mode Management
    # Display transient notice if set
    #
    if st.session_state.ui_notice:
        st.toast(st.session_state.ui_notice, icon="✅")
        st.session_state.ui_notice = None

    # Top-level Mode Selector
    modes = {
        "run": "🚀 Run / Input",
        "history": "🕒 Job History"
    }
    selected_mode_label = st.segmented_control(
        "View Mode",
        options=list(modes.values()),
        default=modes[st.session_state.view_mode],
        label_visibility="collapsed"
    )
    # Sync label back to key
    new_view_mode_list = [k for k, v in modes.items() if v == selected_mode_label]
    if new_view_mode_list:
        new_view_mode = new_view_mode_list[0]
        if new_view_mode != st.session_state.view_mode:
            st.session_state.view_mode = new_view_mode
            st.rerun()

    if st.session_state.view_mode == "history":
        selected_id = render_job_history_panel(sort_runs_for_sidebar(recent_runs), st.session_state.selected_job_id)
        if selected_id != st.session_state.selected_job_id:
            st.session_state.selected_job_id = selected_id
            st.session_state.needs_initial_recovery = False
            st.session_state.view_mode = "run"  # Auto-switch on selection
            st.session_state.ui_notice = f"Job {selected_id} をロードしました"
            st.rerun()
    
    else: 
        # view_mode == "run"
        # 共有検証結果（初期値は空）
        name_validation = input_validation = output_validation = index_validation = tx2gene_validation = run_validation = {}

        # Run/Input Header
        col_title, col_new = st.columns([3, 1])
        with col_title:
            if st.session_state.selected_job_id:
                current_job = next((s for s in recent_runs if s.run_dir.name == st.session_state.selected_job_id), None)
                if current_job:
                    st.info(f"**Current Job:** `{current_job.run_name}` | **Status:** `{current_job.status}` | **Elapsed:** {format_elapsed_time(current_job.elapsed_seconds, current_job.status)}")
            else:
                st.info("解析の条件を入力して実行してください。既存のジョブを確認するには Job History を開いてください。")
        
        with col_new:
            if st.button("➕ New Analysis", use_container_width=True, type="secondary"):
                st.session_state.selected_job_id = None
                st.session_state.needs_initial_recovery = False
                st.session_state.ui_notice = "新規解析モードに切り替えました"
                st.session_state.view_mode = "run"
                st.rerun()
        #
        # 4. Shared Maintainance
        #
        with st.sidebar:
            st.markdown("---")
            st.markdown("### メンテナンス")
            if st.button("❌ Session State Clear"):
                for key in list(st.session_state.keys()):
                    del st.session_state[key]
                st.rerun()
        #
        # 5. Main View Mode Selection
        #
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

            # render Sections
            render_fastq_section(st.session_state.get("fastq_df"))

            st.session_state.sample_df = render_sample_section(
                st.session_state.sample_df if st.session_state.sample_df is not None else pd.DataFrame()
            )

            reference_values = render_reference_section()

            #
            # 内部だけ quantifier を持っておく
            #
            if "quantifier" not in st.session_state:
                st.session_state.quantifier = "salmon"

            if "quantifier_version" not in st.session_state:
                st.session_state.quantifier_version = "1.10.1"

            # [v0.6.0 C-03 / C-08]
            # session_state が salmon_index_path / tx2gene_path を直接保持している。
            # つまり GUI の入力契約自体が Salmon 前提になっている。
            # v0.6.0 では backend 非依存な reference/index/resource 名へ寄せる
            st.session_state.salmon_index_path = reference_values["salmon_index_path"]
            st.session_state.tx2gene_path = reference_values["tx2gene_path"]

            # [v0.6.0 C-04]
            # UI から取得した参照情報を Salmon 専用 validator で直接検証している。
            # validator の責務分離前に名称だけ変えると危険なので、
            # まずは app.py -> generic validator facade の配線を作るのが安全。
            index_validation = validate_salmon_index(st.session_state.salmon_index_path)
            tx2gene_validation = validate_tx2gene_file(st.session_state.tx2gene_path)

            if st.session_state.get("last_salmon_index_path") != st.session_state.salmon_index_path:

                # [v0.6.0 C-03]
                # last_salmon_index_path という session key 名自体が Salmon 固有。
                # 依存除去の際は状態キー名の移行も必要。
                st.session_state.strandedness_prediction = None
                st.session_state.last_salmon_index_path = st.session_state.salmon_index_path

            if reference_values.get("estimate_strandedness"):
                if st.session_state.get("sample_df") is not None and not st.session_state.sample_df.empty and index_validation["is_valid"]:
                    with st.spinner("Strandedness を推定中..."):

                        # [v0.6.0 C-03 / C-08]
                        # strandedness 推定が salmon_index_path を直接要求している。
                        # ここは将来的に backend 別推定へ広げるか、
                        # まずは Salmon 専用機能として adapter に隔離する必要あり。
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

            # ---
            # Handle Run Execution
            # ---
            if run_values["run_requested"]:
                started_at_iso = datetime.now(timezone.utc).astimezone().isoformat()
                job_id = f"RUN_GUI_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                
                # 現在の出力ディレクトリを新しいジョブのルートとして使用
                dirs = setup_run_directory(Path(st.session_state.output_dir), job_id)
                
                # [v0.6.0 C-03 / C-08 / C-09]
                # run_config.json の内容が salmon_index_path / tx2gene_path を前提に固定されている。
                # app.py が GUI 入力を backend 固有 config へ直接落としているため、
                # v0.6.0 では app.py -> generic run request config -> backend adapterの段階を分けたい。
                config_data = {
                    "analysis_name": st.session_state.analysis_name,
                    "input_dir": str(Path(st.session_state.input_dir).resolve()),
                    "output_dir": str(Path(st.session_state.output_dir).resolve()),
                    "salmon_index_path": str(Path(st.session_state.salmon_index_path).resolve()),
                    "tx2gene_path": str(Path(st.session_state.tx2gene_path).resolve()),
                    "strandedness_mode": st.session_state.strandedness_mode,
                    "threads": st.session_state.threads,
                    "strandedness_prediction": st.session_state.strandedness_prediction,

                    # v0.6.0 quantifier追加
                    "quantifier": st.session_state.quantifier,
                    "quantifier_version": st.session_state.quantifier_version,
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
                
                # 追跡可能性のためにGUIステータスを書き込み (v0.5.0 requirement)
                # [v0.6.0 C-03]
                # gui_input_status.csv にも salmon_index_path / tx2gene_path を直接書いている。
                # トレーサビリティ要件自体は維持しつつ、
                # backend 固有語彙の露出をどの層まで許容するか整理が必要。
                status_df = pd.DataFrame([
                    {"parameter": "analysis_name", "value": st.session_state.analysis_name},
                    {"parameter": "input_dir", "value": str(Path(st.session_state.input_dir).resolve())},
                    {"parameter": "output_dir", "value": str(Path(st.session_state.output_dir).resolve())},
                    {"parameter": "salmon_index_path", "value": str(Path(st.session_state.salmon_index_path).resolve())},
                    {"parameter": "tx2gene_path", "value": str(Path(st.session_state.tx2gene_path).resolve())},
                    {"parameter": "strandedness_mode", "value": st.session_state.strandedness_mode},
                    {"parameter": "threads", "value": st.session_state.threads},
                    # v0.6.0 追加
                    {"parameter": "quantifier", "value": st.session_state.quantifier},
                    {"parameter": "quantifier_version", "value": st.session_state.quantifier_version},

                    {"parameter": "started_at", "value": started_at_iso},
                    {"parameter": "job_id", "value": job_id},
                ])
                status_path = (dirs.inputs / "gui_input_status.csv").resolve()
                status_df.to_csv(status_path, index=False)
                
                # [v0.6.0 C-03 / C-09]
                # app.py は最終的に run-gui-backend へ config file を渡すだけだが、
                # その config schema がすでに Salmon 固有である点が本質。
                # ここは CLI コマンド自体より、前段の config 契約を抽象化するのが先。
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
                        "started_at": started_at_iso,
                        # v0.6.0 追加
                        "quantifier": st.session_state.quantifier,
                        "quantifier_version": st.session_state.quantifier_version,
                    },
                    resources={"threads": st.session_state.threads}
                )
                
                executor = LocalJobExecutor(Path(st.session_state.output_dir))
                executor.submit(request)
                
                st.session_state.selected_job_id = job_id
                st.session_state.needs_initial_recovery = False
                st.rerun()

        else:
            # ==========================================
            # JOB VIEW MODE
            # ==========================================
            current_job = next((s for s in recent_runs if s.run_dir.name == st.session_state.selected_job_id), None)
            
            if not current_job:
                # Stale selection detected inside Job View mode
                # Try once to pick the best active run, or fall back to New Analysis
                suggestion = pick_active_run(recent_runs, None)
                st.session_state.selected_job_id = suggestion.run_dir.name if suggestion else None
                st.session_state.needs_initial_recovery = False
                st.rerun()
                return

            col1, col2, col3 = st.columns([2, 1, 1])
            with col1:
                st.subheader(f"Job: {current_job.run_name}")
                st.caption(f"ID: `{current_job.run_dir.name}`")
            with col2:
                status_color = "green" if current_job.status == "completed" else "orange" if current_job.status in ("queued", "running") else "red"
                st.markdown(f"**Status:** :{status_color}[{current_job.status.upper()}]")
                st.markdown(f"**Started:** {current_job.started_at or 'N/A'}")
            with col3:
                st.markdown(f"**Samples:** {current_job.sample_count}")
                st.markdown(f"**Elapsed:** {format_elapsed_time(current_job.elapsed_seconds, current_job.status)}")
            
            st.markdown("---")

            # Monitoring uses the parent of the job's run_dir as the base
            monitor = JobMonitor(current_job.run_dir.parent)
            spec = monitor.get_status(current_job.run_dir.name)
            
            if spec:
                if spec.status in ("queued", "running"):
                    with st.status(f"Job '{spec.run_id}' is running...", expanded=True):
                        st.write(f"Status: `{spec.status}`")
                        st.write("You can safely leave this page or check other jobs. The analysis continues in the background.")
                    time.sleep(2)
                    st.rerun()
                    
                elif spec.status == "completed":
                    run_dir = current_job.run_dir
                    
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
                    st.info(f"Check logs at: `{current_job.run_dir / 'logs' / 'run.log'}`")
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

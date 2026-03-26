from __future__ import annotations

import streamlit as st
import pandas as pd
import json
import sys
from pathlib import Path

# Add src to sys.path to allow importing from iwa_rnaseq_counter package
sys.path.append(str(Path(__file__).parent / "src"))

from iwa_rnaseq_counter.legacy.config import get_default_session_state
import time
from datetime import datetime
from iwa_rnaseq_counter.legacy.fastq_discovery import collect_fastq_metadata, discover_fastq_files
from iwa_rnaseq_counter.legacy.sample_parser import (
    apply_sample_table_edits,
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
from iwa_rnaseq_counter.legacy.salmon_runner import run_salmon_quant
from iwa_rnaseq_counter.legacy.gene_aggregator import (
    load_tx2gene_map,
    aggregate_transcript_to_gene,
    build_transcript_quant_table,
    save_quant_tables,
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
import shutil
import os


def init_session_state() -> None:
    defaults = get_default_session_state()
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def run_app() -> None:
    st.set_page_config(page_title="iwa-rnaseq-counter", layout="wide")
    render_app_header()

    with st.sidebar:
        st.markdown("### メンテナンス")
        if st.button("❌ Session State Clear"):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.rerun()

    analysis_values = render_analysis_section()
    st.session_state.analysis_name = analysis_values["analysis_name"]
    st.session_state.input_dir = analysis_values["input_dir"]
    st.session_state.output_dir = analysis_values["output_dir"]

    name_validation = validate_analysis_name(st.session_state.analysis_name)
    input_validation = validate_input_directory(st.session_state.input_dir)
    output_validation = validate_output_directory(st.session_state.output_dir)

    # Discovery logic
    if analysis_values.get("scan_requested") and input_validation["is_valid"]:
        fastq_files = discover_fastq_files(st.session_state.input_dir)
        st.session_state.scanned_fastq_files = fastq_files
        st.session_state.fastq_df = collect_fastq_metadata(fastq_files)

        sample_groups = group_fastq_by_sample(fastq_files)
        st.session_state.detected_samples = sample_groups
        layout_map = infer_sample_layout(sample_groups)
        lane_map = detect_lane_groups(sample_groups)
        st.session_state.sample_df = build_sample_table(sample_groups, layout_map, lane_map)
        st.success("FASTQ をスキャンしました。")

    if analysis_values.get("load_csv_requested") and analysis_values.get("csv_path"):
        try:
            csv_df = parse_sample_sheet(analysis_values["csv_path"], st.session_state.input_dir)
            st.session_state.sample_df = csv_df
            st.success(f"CSV から {len(csv_df)} サンプルを読み込みました。")
        except Exception as e:
            st.error(f"CSV 読み込みエラー: {e}")

    # UI Components
    render_fastq_section(st.session_state.get("fastq_df"))

    st.session_state.sample_df = render_sample_section(
        st.session_state.sample_df if st.session_state.sample_df is not None else pd.DataFrame()
    )

    reference_values = render_reference_section()
    st.session_state.salmon_index_path = reference_values["salmon_index_path"]
    st.session_state.tx2gene_path = reference_values["tx2gene_path"]

    index_validation = validate_salmon_index(st.session_state.salmon_index_path)
    tx2gene_validation = validate_tx2gene_file(st.session_state.tx2gene_path)

    strandedness_result = None
    if st.session_state.get("sample_df") is not None and not st.session_state.sample_df.empty and index_validation["is_valid"]:
        strandedness_result = infer_strandedness(
            sample_df=st.session_state.sample_df,
            input_dir=st.session_state.input_dir,
            salmon_index_path=st.session_state.salmon_index_path,
        )
        st.session_state.strandedness_prediction = strandedness_result

    # Final validation for RUN button state
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

    if run_values["run_requested"]:
        from iwa_rnaseq_counter.legacy.run_artifacts import setup_run_directory, save_run_config, save_sample_sheet, build_output_manifest
        
        st.session_state.run_status = "running"
        start_time = time.time()
        try:
            run_dir = setup_run_directory(st.session_state.output_dir, st.session_state.analysis_name)
            master_log_path = run_dir / "logs" / "run.log"

            if st.session_state.strandedness_prediction and "probe_cmd" in st.session_state.strandedness_prediction:
                master_log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(master_log_path, "a") as f:
                    f.write(f"\n[STRANDEDNESS PROBE COMMAND]\n{' '.join(st.session_state.strandedness_prediction['probe_cmd'])}\n")

            with st.status("Running Salmon quantification...", expanded=True) as status:
                st.write("Step 1: Running Salmon for all samples...")
                run_result = run_salmon_quant(
                    sample_df=st.session_state.sample_df,
                    salmon_index_path=st.session_state.salmon_index_path,
                    run_output_dir=str(run_dir),
                    strandedness_mode=st.session_state.strandedness_mode,
                    threads=st.session_state.threads,
                )
                st.session_state.latest_log_summary = run_result["log_summary"]
                
                outputs = run_result["outputs"]
                success_outputs = [o for o in outputs if o.get("is_success")]
                failure_count = len(outputs) - len(success_outputs)

                if success_outputs:
                    # 4. Save Artifacts (Config, Initial Sample Sheet)
                    config_data = {
                        "analysis_name": st.session_state.analysis_name,
                        "input_dir": st.session_state.input_dir,
                        "output_dir": st.session_state.output_dir,
                        "salmon_index_path": st.session_state.salmon_index_path,
                        "tx2gene_path": st.session_state.tx2gene_path,
                        "strandedness_mode": st.session_state.strandedness_mode,
                        "threads": st.session_state.threads,
                        "strandedness_prediction": st.session_state.strandedness_prediction
                    }
                    config_path = save_run_config(run_dir, config_data)
                    sample_sheet_path = save_sample_sheet(run_dir, st.session_state.sample_df)

                    # 5. Aggregate
                    status.update(label=f"Aggregating {len(success_outputs)} samples...", expanded=True)
                    st.write("Step 2: Aggregating successful results...")
                    
                    tx2gene_df = load_tx2gene_map(st.session_state.tx2gene_path)
                    t_tpm_df = build_transcript_quant_table(success_outputs, value_type="TPM")
                    t_nr_df = build_transcript_quant_table(success_outputs, value_type="NumReads")
                    g_tpm_df = aggregate_transcript_to_gene(t_tpm_df, tx2gene_df)
                    g_nr_df = aggregate_transcript_to_gene(t_nr_df, tx2gene_df)
                    
                    input_source = "unknown"
                    if "input_source" in st.session_state.sample_df.columns and not st.session_state.sample_df.empty:
                        input_source = str(st.session_state.sample_df["input_source"].iloc[0])
                        
                    from iwa_rnaseq_counter.legacy.sample_parser import METADATA_COLUMNS
                    sample_metadata_columns = [c for c in METADATA_COLUMNS if c in st.session_state.sample_df.columns]
                    
                    sample_metadata_columns_nonempty = []
                    for col in sample_metadata_columns:
                        if col == "exclude":
                            if st.session_state.sample_df[col].fillna(False).astype(bool).any():
                                sample_metadata_columns_nonempty.append(col)
                        else:
                            if st.session_state.sample_df[col].fillna("").astype(str).str.strip().replace("nan", "").ne("").any():
                                sample_metadata_columns_nonempty.append(col)
                    
                    sample_ids_all = st.session_state.sample_df["sample_id"].tolist()
                    sample_ids_success = [o["sample_id"] for o in success_outputs]
                    sample_ids_failed = [o["sample_id"] for o in outputs if not o.get("is_success")]
                    sample_ids_aggregated = sample_ids_success # v0.1.7: aggregation uses success samples

                    # Normalize output paths to be relative to run_dir
                    rel_outputs = []
                    for o in outputs:
                        rel_o = o.copy()
                        for key in ["quant_path", "aux_info_dir", "log_path"]:
                            if rel_o.get(key):
                                try:
                                    rel_o[key] = str(Path(rel_o[key]).relative_to(run_dir))
                                except ValueError:
                                    pass
                        rel_outputs.append(rel_o)

                    run_summary = {
                        "analysis_name": st.session_state.analysis_name,
                        "run_name": st.session_state.analysis_name,
                        "sample_count": len(st.session_state.sample_df),
                        "success_count": len(success_outputs),
                        "failure_count": failure_count,
                        "sample_ids_all": sample_ids_all,
                        "sample_ids_success": sample_ids_success,
                        "sample_ids_failed": sample_ids_failed,
                        "sample_ids_aggregated": sample_ids_aggregated,
                        "input_source": input_source,
                        "sample_metadata_columns": sample_metadata_columns,
                        "sample_metadata_columns_nonempty": sample_metadata_columns_nonempty,
                        "transcript_rows": len(t_tpm_df),
                        "gene_rows": len(g_tpm_df),
                        "elapsed_seconds": time.time() - start_time,
                        "outputs": outputs,
                        "save_path": str(run_dir),
                        "quantifier": "salmon",
                        "quantifier_version": "1.10.1",
                        "salmon_index_path": st.session_state.salmon_index_path,
                        "tx2gene_path": st.session_state.tx2gene_path,
                        "strandedness": st.session_state.strandedness_prediction,
                        "threads": st.session_state.threads,
                    }
                    
                    # Create a version for disk with relative paths
                    disk_summary = run_summary.copy()
                    disk_summary["save_path"] = run_dir.name
                    disk_summary["outputs"] = rel_outputs

                    matrices = {
                        "transcript_tpm": t_tpm_df, "transcript_numreads": t_nr_df,
                        "gene_tpm": g_tpm_df, "gene_numreads": g_nr_df
                    }
                    
                    from iwa_rnaseq_counter.legacy.gene_aggregator import save_quant_tables
                    output_paths = save_quant_tables(
                        matrices=matrices,
                        sample_df=st.session_state.sample_df,
                        run_summary=disk_summary, # Save the relative version
                        run_output_dir=run_dir
                    )

                    # Generate dataset_manifest.json
                    manifest_data = {
                        "manifest_version": "1.0",
                        "generated_at": datetime.now().astimezone().isoformat(),
                        "app_name": "iwa-rnaseq-counter",
                        "app_version": "v0.1.7",
                        "run_name": st.session_state.analysis_name,
                        "analysis_name": st.session_state.analysis_name,
                        "input_source": input_source,
                        "quantifier": "salmon",
                        "quantifier_version": "1.10.1",
                        "sample_count_total": len(sample_ids_all),
                        "sample_count_success": len(sample_ids_success),
                        "sample_count_failed": len(sample_ids_failed),
                        "sample_ids_all": sample_ids_all,
                        "sample_ids_success": sample_ids_success,
                        "sample_ids_failed": sample_ids_failed,
                        "sample_ids_aggregated": sample_ids_aggregated,
                        "results_dir": "results",
                        "files": {
                            "sample_metadata": "results/sample_metadata.csv",
                            "sample_qc_summary": "results/sample_qc_summary.csv",
                            "transcript_tpm": "results/transcript_tpm.csv",
                            "transcript_numreads": "results/transcript_numreads.csv",
                            "gene_tpm": "results/gene_tpm.csv",
                            "gene_numreads": "results/gene_numreads.csv",
                            "run_summary": "results/run_summary.json",
                            "sample_sheet": "sample_sheet.csv",
                            "run_config": "run_config.json",
                            "run_log": "logs/run.log"
                        }
                    }
                    from iwa_rnaseq_counter.legacy.run_artifacts import save_dataset_manifest
                    save_dataset_manifest(run_dir, manifest_data)
                    
                    st.session_state.run_summary = run_summary
                    st.session_state.output_files = build_output_manifest(
                        run_dir, output_paths, config_path, sample_sheet_path, master_log_path
                    )
                    st.session_state.run_status = "success"
                    
                    if failure_count > 0:
                        status.update(label=f"Completed with {failure_count} failures.", state="complete", expanded=True)
                        st.warning(f"{failure_count} サンプルの解析に失敗しましたが、成功したサンプルのみで集約しました。")
                    else:
                        status.update(label="Analysis completed successfully!", state="complete", expanded=False)
                else:
                    st.session_state.run_status = "failed"
                    st.error("すべてのサンプルの解析に失敗しました。")
                    if run_result["errors"]:
                        st.error("\n".join(run_result["errors"]))
                    status.update(label="All samples failed.", state="error", expanded=True)
        except Exception as e:
            st.session_state.run_status = "failed"
            st.error(f"Execution Error: {e}")

    render_result_section(
        run_status=st.session_state.run_status,
        output_files=st.session_state.output_files,
        log_summary=st.session_state.latest_log_summary,
        run_summary=st.session_state.get("run_summary")
    )

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

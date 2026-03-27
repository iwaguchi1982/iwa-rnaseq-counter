import time
import pandas as pd
from pathlib import Path
from datetime import datetime, timezone
import logging

from iwa_rnaseq_counter.legacy.salmon_runner import run_salmon_quant
from iwa_rnaseq_counter.legacy.gene_aggregator import build_transcript_quant_table, aggregate_transcript_to_gene, load_tx2gene_map, save_quant_tables
from iwa_rnaseq_counter.legacy.run_artifacts import save_dataset_manifest
from iwa_rnaseq_counter.builders.gui_artifact_export import write_gui_supporting_inputs

logger = logging.getLogger(__name__)

def run_gui_backend_pipeline(run_dir: Path, config_data: dict, sample_df: pd.DataFrame, started_at_iso: str):
    """
    Executes the exact pipeline that the GUI used to run synchronously.
    Outputs are written to run_dir.
    """
    start_time = time.time()
    
    salmon_index_path = config_data.get("salmon_index_path")
    tx2gene_path = config_data.get("tx2gene_path")
    strandedness_mode = config_data.get("strandedness_mode", "Auto-detect")
    threads = config_data.get("threads", 4)
    analysis_name = config_data.get("analysis_name", "GUI_Run")
    
    logger.info("Step 1: Running Salmon for all samples...")
    run_result = run_salmon_quant(
        sample_df=sample_df,
        salmon_index_path=salmon_index_path,
        run_output_dir=str(run_dir),
        strandedness_mode=strandedness_mode,
        threads=threads,
    )
    
    outputs = run_result["outputs"]
    success_outputs = [o for o in outputs if o.get("is_success")]
    failure_count = len(outputs) - len(success_outputs)
    
    if not success_outputs:
        raise RuntimeError("All samples failed Salmon quantification.")
        
    logger.info(f"Step 2: Aggregating {len(success_outputs)} successful results...")
    tx2gene_df = load_tx2gene_map(tx2gene_path)
    t_tpm_df = build_transcript_quant_table(success_outputs, value_type="TPM")
    t_nr_df = build_transcript_quant_table(success_outputs, value_type="NumReads")
    g_tpm_df = aggregate_transcript_to_gene(t_tpm_df, tx2gene_df)
    g_nr_df = aggregate_transcript_to_gene(t_nr_df, tx2gene_df)
    
    input_source = "unknown"
    if "input_source" in sample_df.columns and not sample_df.empty:
        input_source = str(sample_df["input_source"].iloc[0])
        
    from iwa_rnaseq_counter.legacy.sample_parser import METADATA_COLUMNS
    sample_metadata_columns = [c for c in METADATA_COLUMNS if c in sample_df.columns]
    
    sample_metadata_columns_nonempty = []
    for col in sample_metadata_columns:
        if col == "exclude":
            if sample_df[col].fillna(False).astype(bool).any():
                sample_metadata_columns_nonempty.append(col)
        else:
            if sample_df[col].fillna("").astype(str).str.strip().replace("nan", "").ne("").any():
                sample_metadata_columns_nonempty.append(col)
                
    sample_ids_all = sample_df["sample_id"].tolist()
    sample_ids_success = [o["sample_id"] for o in success_outputs]
    sample_ids_failed = [o["sample_id"] for o in outputs if not o.get("is_success")]
    sample_ids_aggregated = sample_ids_success
    
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
        "analysis_name": analysis_name,
        "run_name": analysis_name,
        "sample_count": len(sample_df),
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
        "salmon_index_path": salmon_index_path,
        "tx2gene_path": tx2gene_path,
        "strandedness": config_data.get("strandedness_prediction"),
        "threads": threads,
        "log_summary": run_result["log_summary"]
    }
    
    disk_summary = run_summary.copy()
    disk_summary["save_path"] = run_dir.name
    disk_summary["outputs"] = rel_outputs

    matrices = {
        "transcript_tpm": t_tpm_df, "transcript_numreads": t_nr_df,
        "gene_tpm": g_tpm_df, "gene_numreads": g_nr_df
    }
    
    output_paths = save_quant_tables(
        matrices=matrices,
        sample_df=sample_df,
        run_summary=disk_summary,
        run_output_dir=run_dir
    )
    
    import shutil
    try:
        shutil.copy(tx2gene_path, run_dir / "results" / "feature_annotation.tsv")
    except Exception as e:
        logger.warning(f"Failed to copy feature_annotation.tsv: {e}")
    
    manifest_data = {
        "manifest_version": "1.0",
        "generated_at": datetime.now().astimezone().isoformat(),
        "app_name": "iwa-rnaseq-counter",
        "app_version": "0.3.5",
        "run_name": analysis_name,
        "analysis_name": analysis_name,
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
    
    save_dataset_manifest(run_dir, manifest_data)
    
    try:
        write_gui_supporting_inputs(
            run_dir=run_dir,
            run_summary=run_summary,
            matrix_rel_path="results/gene_numreads.csv",
            log_rel_path="logs/run.log",
            started_at=started_at_iso
        )
    except Exception as spec_err:
        logger.warning(f"Failed to generate pipeline specs: {spec_err}")
        
    logger.info("Pipeline completed successfully.")
    return True

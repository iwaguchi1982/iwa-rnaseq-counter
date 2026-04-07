import time
import pandas as pd
from pathlib import Path
from datetime import datetime, timezone
import logging

from iwa_rnaseq_counter.pipeline.quantifiers.registry import get_quantifier
from iwa_rnaseq_counter.legacy.gene_aggregator import build_transcript_quant_table, aggregate_transcript_to_gene, load_tx2gene_map, save_quant_tables
from iwa_rnaseq_counter.legacy.run_artifacts import save_dataset_manifest
from iwa_rnaseq_counter.builders.gui_artifact_export import write_gui_supporting_inputs

logger = logging.getLogger(__name__)
 
def _load_gene_counts_table(gene_counts_path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(gene_counts_path, sep="\t")

    if {"feature_id", "count"}.issubset(df.columns):
        out = df[["feature_id", "count"]].copy()
    elif df.shape[1] >= 2:
        out = df.iloc[:, :2].copy()
        out.columns = ["feature_id", "count"]
    else:
        raise ValueError(f"Unsupported gene_counts file shape: {gene_counts_path}")

    out["feature_id"] = out["feature_id"].astype(str)
    out["count"] = pd.to_numeric(out["count"], errors="coerce").fillna(0).astype(int)
    return out


def _build_gene_counts_matrix_from_outputs(success_outputs: list[dict]) -> pd.DataFrame:
    series_list: list[pd.Series] = []

    for output in success_outputs:
        gene_counts_path = output.get("gene_counts_path")
        sample_id = output.get("sample_id")

        if not gene_counts_path:
            raise ValueError(f"gene_counts_path is missing for sample {sample_id}")

        counts_df = _load_gene_counts_table(gene_counts_path)
        s = counts_df.set_index("feature_id")["count"]
        s.name = str(sample_id)
        series_list.append(s)

    if not series_list:
        raise RuntimeError("No gene-level count outputs were available.")

    matrix_df = pd.concat(series_list, axis=1).fillna(0).astype(int)
    matrix_df.index.name = "feature_id"
    return matrix_df


def _build_gui_matrices_from_run_result(run_result: dict, tx2gene_path: str):
    aggregation_input_kind = run_result.get("aggregation_input_kind", "transcript_quant")
    outputs = run_result.get("outputs", [])
    success_outputs = [o for o in outputs if o.get("is_success")]

    if not success_outputs:
        raise RuntimeError("No successful quantifier outputs found.")

    if aggregation_input_kind == "transcript_quant":
        tx2gene_df = load_tx2gene_map(tx2gene_path)
        t_tpm_df = build_transcript_quant_table(success_outputs, value_type="TPM")
        t_nr_df = build_transcript_quant_table(success_outputs, value_type="NumReads")
        g_tpm_df = aggregate_transcript_to_gene(t_tpm_df, tx2gene_df)
        g_nr_df = aggregate_transcript_to_gene(t_nr_df, tx2gene_df)
        return success_outputs, t_tpm_df, t_nr_df, g_tpm_df, g_nr_df

    if aggregation_input_kind == "gene_counts":
        g_nr_df = _build_gene_counts_matrix_from_outputs(success_outputs)

        # v0.7.1 最小版:
        # direct gene-count backend では transcript-level / TPM をまだ持たない
        t_tpm_df = pd.DataFrame()
        t_nr_df = pd.DataFrame()
        g_tpm_df = pd.DataFrame()

        return success_outputs, t_tpm_df, t_nr_df, g_tpm_df, g_nr_df

    raise ValueError(f"Unsupported aggregation_input_kind: {aggregation_input_kind!r}")


def _summarize_output_capabilities(outputs: list[dict]) -> dict:
    success_outputs = [o for o in outputs if o.get("is_success")]

    has_mapping_metrics = any(
        any(o.get(k) is not None for k in ["num_processed", "num_mapped", "mapping_rate"])
        for o in success_outputs
    )
    has_transcript_quant = any(
        bool(o.get("transcript_quant_path") or o.get("quant_path"))
        for o in success_outputs
    )
    has_gene_counts = any(bool(o.get("gene_counts_path")) for o in success_outputs)

    return {
        "has_mapping_metrics": has_mapping_metrics,
        "has_transcript_quant": has_transcript_quant,
        "has_gene_counts": has_gene_counts,
    }


def run_gui_backend_pipeline(run_dir: Path, config_data: dict, sample_df: pd.DataFrame, started_at_iso: str):
    """
    Executes the exact pipeline that the GUI used to run synchronously.
    Outputs are written to run_dir.
    """
    start_time = time.time()
    
    quantifier_name = config_data.get("quantifier", "salmon")
    quantifier_index_path = (
        config_data.get("quantifier_index_path")
        or config_data.get("salmon_index_path")
    )
    tx2gene_path = config_data.get("tx2gene_path")
    strandedness_mode = config_data.get("strandedness_mode", "Auto-detect")
    threads = config_data.get("threads", 4)
    analysis_name = config_data.get("analysis_name", "GUI_Run")
    quantifier_version = config_data.get("quantifier_version")
    annotation_gtf_path = config_data.get("annotation_gtf_path")
    
    logger.info(f"Step 1: Running {quantifier_name} for all samples...")
    quant = get_quantifier(quantifier_name)
    run_result = quant.run_quant(
        sample_df=sample_df,
        run_output_dir=run_dir,
        threads=threads,
        strandedness_mode=strandedness_mode,
        reference_config={
            "quantifier_index": quantifier_index_path,
            "tx2gene_path": tx2gene_path,
            "annotation_gtf_path": annotation_gtf_path,
        },
    )

    quantifier_name = run_result.get("quantifier", quantifier_name)
    quantifier_version = (
        run_result.get("quantifier_version")
        or quantifier_version
        or "unknown"
    )
    
    outputs = run_result["outputs"]
    failure_count = len(outputs) - len([o for o in outputs if o.get("is_success")])
    
    success_outputs, t_tpm_df, t_nr_df, g_tpm_df, g_nr_df = _build_gui_matrices_from_run_result(
        run_result,
        tx2gene_path,
    )

    logger.info(f"Step 2: Aggregating {len(success_outputs)} successful results...")
    
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
        for key in [
            "quant_path",
            "transcript_quant_path",
            "gene_counts_path",
            "aux_info_dir",
            "meta_info_json",
            "log_path",
            "output_dir",
        ]:
            if rel_o.get(key):
                try:
                    rel_o[key] = str(Path(rel_o[key]).relative_to(run_dir))
                except ValueError:
                    pass
        rel_outputs.append(rel_o)
        
    capability_summary = _summarize_output_capabilities(outputs)

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
        "gene_rows": len(g_nr_df),
        "elapsed_seconds": time.time() - start_time,
        "has_mapping_metrics": capability_summary["has_mapping_metrics"],
        "has_transcript_quant": capability_summary["has_transcript_quant"],
        "has_gene_counts": capability_summary["has_gene_counts"],
        "outputs": outputs,
        "save_path": str(run_dir),
        "quantifier": quantifier_name,
        "quantifier_version": quantifier_version,
        "quantifier_index_path": quantifier_index_path,
        # v0.7.0 compatibility alias:
        # app.py / 過去run表示が salmon_index_path をまだ参照する可能性があるため残す
        "salmon_index_path": quantifier_index_path,
        "tx2gene_path": tx2gene_path,
        "annotation_gtf_path": annotation_gtf_path,
        "aggregation_input_kind": run_result.get("aggregation_input_kind", "transcript_quant"),
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
    
    # Step 3: Prepare feature_annotation.tsv (v0.5.0 Contract)
    from iwa_rnaseq_counter.legacy.annotation_helper import prepare_feature_annotation, get_standard_annotation_path
    annotation_out = get_standard_annotation_path(run_dir)
    has_annotation = prepare_feature_annotation(tx2gene_path, annotation_out)
    
    if has_annotation:
        logger.info(f"Feature annotation prepared at {annotation_out}")
    else:
        logger.warning("Could not prepare feature_annotation.tsv. Reporter will fall back to feature_id.")

    manifest_data = {
        "manifest_version": "1.0",
        "generated_at": datetime.now().astimezone().isoformat(),
        "app_name": "iwa-rnaseq-counter",
        "app_version": "0.3.5",
        "run_name": analysis_name,
        "analysis_name": analysis_name,
        "input_source": input_source,
        "quantifier": quantifier_name,
        "quantifier_version": quantifier_version,
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
    
    if has_annotation:
        manifest_data["files"]["feature_annotation"] = "results/feature_annotation.tsv"
    
    save_dataset_manifest(run_dir, manifest_data)
    
    try:
        write_gui_supporting_inputs(
            run_dir=run_dir,
            run_summary=run_summary,
            matrix_rel_path="results/gene_numreads.csv",
            log_rel_path="logs/run.log",
            started_at=started_at_iso,
            feature_annotation_path=str(annotation_out.resolve()) if has_annotation else None,
            feature_annotation_available=has_annotation
        )
    except Exception as spec_err:
        logger.warning(f"Failed to generate pipeline specs: {spec_err}")
        
    logger.info("Pipeline completed successfully.")
    return True

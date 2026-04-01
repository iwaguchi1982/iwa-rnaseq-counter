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

def run_gui_backend_pipeline(run_dir: Path, config_data: dict, sample_df: pd.DataFrame, started_at_iso: str):
    """
    GUIが同期実行に使用していたパイプラインと全く同じものを実行します。
    出力はrun_dirに書き込まれます。
    """
    start_time = time.time()
    
    # [v0.6.0 C-03 / C-08]
    # GUI backend が config_data から salmon_index_path / tx2gene_path を直接読んでいる。
    # 入口契約そのものが Salmon 前提になっているため、
    # v0.6.0 では backend 非依存な reference/index/resource 契約へ寄せたい。
    salmon_index_path = config_data.get("salmon_index_path")
    tx2gene_path = config_data.get("tx2gene_path")
    strandedness_mode = config_data.get("strandedness_mode", "Auto-detect")
    threads = config_data.get("threads", 4)
    analysis_name = config_data.get("analysis_name", "GUI_Run")
    quantifier_name = config_data.get("quantifier", "salmon")
    
    logger.info(f"Step 1: Running {quantifier_name} for all samples...")

    # [v0.6.0 C-01]
    # GUI backend が Salmon 実装関数 run_salmon_quant() を直接呼んでいる。
    # ここは将来的に quantifier adapter の共通 API
    # 例: quantifier_impl.run_quant(...)
    # へ置き換え、GUI backend は backend 差分を知らない層にしたい。
    quant = get_quantifier(quantifier_name)

    run_result = quant.run_quant(
        sample_df=sample_df,
        run_output_dir=run_dir,
        threads=threads,
        strandedness_mode=strandedness_mode,
        reference_config={
            "quantifier_index": salmon_index_path,
            "tx2gene_path": tx2gene_path,
        },
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
        # [v0.6.0 C-05 / C-09]
        # run_summary に quantifier 名と version が固定値で埋め込まれている。
        # GUI backend が "salmon" / "1.10.1" を直書きしており、
        # backend 情報の出所が抽象層ではなく GUI backend 自身になっている。
        # v0.6.0 では backend 実装から受け取る metadata に寄せたい。
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
    
    # Step 3: feature_annotation.tsv を準備する (v0.5.0 Contract)
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
        # [v0.6.0 C-05 / C-09]
        # dataset manifest にも quantifier 名と version が固定値で埋め込まれている。
        # backend 差分を manifest writer / adapter でなく GUI backend 本体が握っている状態。
        # v0.6.0 では backend 由来 metadata の注入点を整理したい。
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
    
    if has_annotation:
        manifest_data["files"]["feature_annotation"] = "results/feature_annotation.tsv"
    
    save_dataset_manifest(run_dir, manifest_data)
    
    try:
        # [v0.6.0 C-09]
        # GUI adapter spec 生成へ feature_annotation_path / availability を渡している。
        # ここ自体は v0.5.0 契約上妥当だが、backend 情報の責務分離後に
        # gui_backend -> spec writer の入力契約を見直す余地あり。
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

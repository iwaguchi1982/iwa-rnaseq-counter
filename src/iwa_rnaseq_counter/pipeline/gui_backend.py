# gui_backend.py
# 全体の役割：GUIから投入されたジョブのバックエンド実行を担う中心的なパイプライン。データの定量、集計、成果物（Artifact）の生成を一括で行う。
# 処理の流れ：設定読込 → Quantifier実行 → 成功結果の集約 → 行列・サマリー保存 → アノテーション準備 → マニフェスト保存 → Spec出力
# ---
# 主要な処理ステップ
# 1. 実行構成の読み込みと解決:
# GUIから渡された config_data からリファレンスパス、スレッド数、鎖性設定などの実行条件を抽出。
# 選択された Quantifier（現在は主に Salmon）を registry から取得して初期化。
# 
# 2. 定量計算 (Quantification):
# 配置されたサンプルごとに定量ツールを実行。
# 全サンプルのうち成功したものと失敗したものを分離・フィルタリングし、下流処理の安全性を確保。
# 
# 3. 遺伝子レベルへの集計 (Aggregation):
# 転写産物（Transcript）レベルの定量結果を、tx2geneマップを用いて遺伝子（Gene）レベルへ集計。
# TPM および NumReads の両方の指標について、Transcript/Gene 各層の行列を構築。
# 
# 4. 実行サマリーと統計の生成:
# 解析名、サンプル数、成功数、処理時間、Quantifier情報など、実行のすべてを記録した run_summary を生成。
# 
# 5. 成果物 (Artifact) のエクスポート:
# 行列データ（CSV）、実行サマリー（JSON）、アノテーション表（TSV）を所定のディレクトリ構造に従って保存。
# 下流の解析ツール（Reporterなど）がデータを正しく読み込めるよう、dataset_manifest.json および標準規格の Spec ファイルを出力。
#

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
    # --- GUI backend pipeline ---
    """
    この関数は、GUIで生成された実行構成とサンプルシートデータを受け取り、
    定量化を実行し、出力を集計し、レポート作成者向けの成果物を作成します。
    GUIが同期実行に使用していたパイプラインと全く同じものを実行します。
    出力はrun_dirに書き込まれます。
    """
    start_time = time.time()
    
    # [v0.6.0 C-03 / C-08]
    # GUI backend は config_data から reference path と実行条件を読み込む。
    # 現時点では reference 側の入力契約に salmon_index_path / tx2gene_path という
    # Salmon 語彙が残っているため、将来的には backend 非依存な契約へ寄せたい。
    salmon_index_path = config_data.get("salmon_index_path")
    # --- GUI config ingestion ---
    # run_config.jsonから参照/リソースパスと実行オプションを読み込みます。
    # この段階では、参照側のキーにはSalmon固有の命名規則がまだ含まれています。
    tx2gene_path = config_data.get("tx2gene_path")
    strandedness_mode = config_data.get("strandedness_mode", "Auto-detect")
    threads = config_data.get("threads", 4)
    analysis_name = config_data.get("analysis_name", "GUI_Run")
    quantifier_name = config_data.get("quantifier", "salmon")
    quantifier_version = config_data.get(
        "quantifier_version",
        "1.10.1" if quantifier_name == "salmon" else "unknown",
    )
    
    logger.info(f"Step 1: Running {quantifier_name} for all samples...")

    # [v0.6.0 C-01]
    # GUI backend は quantifier registry 経由で backend 実装を解決して実行する。
    # これにより、この層は特定 backend 実装の直 import / 直呼びから切り離された。
    # 今後の課題は、reference 入力契約や metadata 注入点も同様に整理すること。
    # --- Quantifier execution ---
    # 選択された量指定子の実装を解決し、
    # バックエンド固有のランナーを直接呼び出す代わりに、共通の量指定子インターフェースを介して実行します。
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
    # --- Successful output filtering ---
    # サンプルごとの定量化結果のうち、成功例と失敗例を分離します。
    # 下流の集計処理およびアーティファクトのエクスポート処理の前に分離します。
    outputs = run_result["outputs"]
    success_outputs = [o for o in outputs if o.get("is_success")]
    failure_count = len(outputs) - len(success_outputs)
    
    if not success_outputs:
        raise RuntimeError("All samples failed Salmon quantification.")
        
    logger.info(f"Step 2: Aggregating {len(success_outputs)} successful results...")
    # --- Gene-level aggregation ---
    # 転写産物レベルの定量出力を転写産物テーブルに変換し、
    # tx2geneマッピングを使用して遺伝子レベルのマトリックスに集計します。
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
    
    # --- Run summary assembly ---
    # 後続の処理結果レンダリング、
    # トレーサビリティ、および仕様書/マニフェスト生成のために、この実行のメモリ内サマリーを作成します。
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
        "quantifier": quantifier_name,
        "quantifier_version": quantifier_version,
        # [v0.6.0 C-05 / C-09]
        # run_summary には quantifier 名と version を記録している。
        # 現在は config_data / 実行文脈ベースで注入しているが、
        # 最終的には backend 実装または正規化された execution context から
        # 一貫して受け取る形へ寄せたい。
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
    # --- Result table export ---
    # トランスクリプトレベルおよび遺伝子レベルの結果テーブルとサマリーアーティファクトを、
    # 後で読み込みやレポート作成のために実行ディレクトリに保存します。
    output_paths = save_quant_tables(
        matrices=matrices,
        sample_df=sample_df,
        run_summary=disk_summary,
        run_output_dir=run_dir
    )
    
    # Step 3: feature_annotation.tsv を準備する (v0.5.0 Contract)
    # --- Feature annotation preparation ---
    # 利用可能なマッピングリソースから、レポーター向けの feature_annotation.tsv を準備します。
    # アノテーションが準備できない場合は、レポーターは feature_id にフォールバックできます。
    from iwa_rnaseq_counter.legacy.annotation_helper import prepare_feature_annotation, get_standard_annotation_path
    annotation_out = get_standard_annotation_path(run_dir)
    has_annotation = prepare_feature_annotation(tx2gene_path, annotation_out)
    
    if has_annotation:
        logger.info(f"Feature annotation prepared at {annotation_out}")
    else:
        logger.warning("Could not prepare feature_annotation.tsv. Reporter will fall back to feature_id.")

    # --- Dataset manifest export ---
    # 生成されたファイル、サンプル数、
    #および下流のレポーターの読み込みのための実行コンテキストを記述したデータセットレベルのマニフェストを作成します。
    manifest_data = {
        "manifest_version": "1.0",
        "generated_at": datetime.now().astimezone().isoformat(),
        "app_name": "iwa-rnaseq-counter",
        "app_version": "0.3.5",
        "run_name": analysis_name,
        "analysis_name": analysis_name,
        "input_source": input_source,
        # [v0.6.0 C-05 / C-09]
        # dataset manifest にも quantifier 名と version を記録している。
        # 現在は GUI backend で実行文脈から注入しているが、
        # 将来的には backend metadata の注入点をさらに整理したい。
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
    
    # --- Spec export ---
    # GUI の実行を、下流コンポーネントで使用される標準コントラクトを通じて利用できるように、
    # MatrixSpec / ExecutionRunSpec アダプタを生成します。
    try:
        # [v0.6.0 C-09]
        # GUI adapter spec 生成へ feature_annotation_path / availability を渡している。
        # ここ自体は v0.5.0 契約上妥当だが、backend 情報や execution context の
        # 注入責務をどこに持たせるかは今後さらに整理したい。
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

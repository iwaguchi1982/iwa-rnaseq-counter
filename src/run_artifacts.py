from __future__ import annotations

import json
from pathlib import Path
import pandas as pd
from datetime import datetime

from src.config import get_default_output_filenames


def setup_run_directory(output_dir: str, analysis_name: str) -> Path:
    """
    解析ごとのディレクトリを作成する。
    output_dir / analysis_name
    """
    base = Path(output_dir)
    run_dir = base / analysis_name
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "logs").mkdir(exist_ok=True)
    return run_dir


def save_run_config(run_output_dir: Path, config_data: dict) -> Path:
    """
    実行時の設定を JSON として保存する (再現性のため)。
    """
    names = get_default_output_filenames()
    path = run_output_dir / names["run_config"]
    
    # タイムスタンプを追加
    config_data["timestamp"] = datetime.now().isoformat()
    config_data["app_version"] = "v0.1.7"
    
    with open(path, "w", encoding="utf-8") as f:
        json.dump(config_data, f, indent=4, ensure_ascii=False)
    return path


def save_dataset_manifest(run_output_dir: Path, manifest_data: dict) -> Path:
    """
    Reporter 向けの入口 manifest を保存する。
    """
    path = run_output_dir / "results" / "dataset_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest_data, f, indent=4, ensure_ascii=False)
    return path


def save_sample_sheet(run_output_dir: Path, sample_df: pd.DataFrame) -> Path:
    """
    実行に使用した全サンプル情報を保存する (入力パス等を含む完全版)。
    """
    names = get_default_output_filenames()
    path = run_output_dir / names["sample_sheet"]
    sample_df.to_csv(path, index=False)
    return path


def save_sample_metadata_csv(run_output_dir: Path, sample_df: pd.DataFrame) -> Path:
    """
    Reporter 向けの軽量なメタデータ表を保存する。
    FASTQ パス等を除き、生物学的な条件設定のみを抽出する。
    """
    from src.sample_parser import METADATA_COLUMNS
    cols = ["sample_id"] + [c for c in METADATA_COLUMNS if c in sample_df.columns]
    
    path = run_output_dir / "sample_metadata.csv"
    sample_df[cols].to_csv(path, index=False)
    return path


def build_output_manifest(
    run_output_dir: Path,
    quant_results: dict,
    config_path: Path,
    sample_sheet_path: Path,
    master_log_path: Path
) -> list[dict]:
    """
    UI での表示や一括DL用の成果物一覧 (Manifest) を作成する。
    設計書に基づき {label, path} 構造とする。
    """
    # quant_results のキーを src/gene_aggregator.py:save_quant_tables の戻り値に合わせる
    manifest = [
        {"label": "Gene TPM Matrix", "path": quant_results.get("gene_tpm")},
        {"label": "Gene NumReads Matrix", "path": quant_results.get("gene_numreads")},
        {"label": "Transcript TPM Matrix", "path": quant_results.get("transcript_tpm")},
        {"label": "Transcript NumReads Matrix", "path": quant_results.get("transcript_numreads")},
        {"label": "Sample Metadata (CSV)", "path": quant_results.get("sample_metadata")},
        {"label": "Sample QC Summary (CSV)", "path": quant_results.get("sample_qc_summary")},
        {"label": "Dataset Manifest (JSON)", "path": str(run_output_dir / "results" / "dataset_manifest.json")},
        {"label": "Run Summary (JSON)", "path": quant_results.get("run_summary")},
        {"label": "Run Config (JSON)", "path": str(config_path)},
        {"label": "Sample Sheet (CSV)", "path": str(sample_sheet_path)},
        {"label": "Execution Log (Log)", "path": str(master_log_path)},
    ]
    # 実在するパスのみを返す
    return [m for m in manifest if m["path"] and Path(m["path"]).exists()]

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
    config_data["app_version"] = "v0.1.6"
    
    with open(path, "w", encoding="utf-8") as f:
        json.dump(config_data, f, indent=4, ensure_ascii=False)
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
    manifest = [
        {"label": "Transcript TPM", "path": quant_results.get("transcript_tpm_csv")},
        {"label": "Transcript NumReads", "path": quant_results.get("transcript_numreads_csv")},
        {"label": "Gene TPM", "path": quant_results.get("gene_tpm_csv")},
        {"label": "Gene NumReads", "path": quant_results.get("gene_numreads_csv")},
        {"label": "Run Summary (JSON)", "path": quant_results.get("run_summary_json")},
        {"label": "Run Config (JSON)", "path": str(config_path)},
        {"label": "Sample Sheet (CSV)", "path": str(sample_sheet_path)},
        {"label": "Sample Metadata (CSV)", "path": quant_results.get("sample_metadata_csv")},
        {"label": "Sample QC Summary (CSV)", "path": quant_results.get("qc_summary_csv")},
        {"label": "Execution Log (Log)", "path": str(master_log_path)},
    ]
    return [m for m in manifest if m["path"]]

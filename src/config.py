import json
from pathlib import Path
from typing import Any


def get_default_session_state() -> dict[str, Any]:
    """
    アプリケーション起動時のデフォルト状態を返す。
    config.json が存在する場合は、その値を優先する。
    """
    # 共通のデフォルト
    defaults = {
        "analysis_name": "yeast_run",
        "input_dir": "",
        "output_dir": "results",
        "scanned_fastq_files": None,
        "fastq_df": None,
        "detected_samples": None,
        "sample_df": None,
        "salmon_index_path": "",
        "tx2gene_path": "",
        "threads": 4,
        "run_status": "idle",
        "latest_log_summary": None,
        "output_files": None,
        "strandedness_mode": "Auto-detect",
        "strandedness_prediction": None,
        "run_validation_result": None,
        "salmon_version": "1.10.1",
        "discovery_mode": "auto",
        "last_csv_path": "",
    }

    # config.json からの上書き
    config_path = Path("config.json")
    if config_path.exists():
        try:
            with open(config_path, "r", encoding="utf-8") as f:
                external_config = json.load(f)
                # トップレベルのキーを上書き
                for k, v in external_config.items():
                    if k in defaults:
                        defaults[k] = v
                
                # 種別のデフォルト設定がある場合はそれも考慮 (将来用)
                species = external_config.get("default_species", "yeast")
                species_defaults = external_config.get("species_defaults", {}).get(species, {})
                for k, v in species_defaults.items():
                    if k in defaults:
                        defaults[k] = v
        except Exception:
            pass

    return defaults


def get_supported_fastq_extensions() -> tuple[str, ...]:
    return (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def get_strandedness_options() -> list[str]:
    return ["Auto-detect", "unstranded", "forward", "reverse"]


def get_default_output_filenames() -> dict[str, str]:
    return {
        "transcript_tpm": "transcript_TPM.csv",
        "transcript_numreads": "transcript_NumReads.csv",
        "gene_tpm": "gene_TPM.csv",
        "gene_numreads": "gene_NumReads.csv",
        "run_log": "run.log",
        "run_config": "run_config.json",
        "run_summary": "run_summary.json",
        "sample_sheet": "sample_sheet.csv",
        "sample_qc_summary": "sample_qc_summary.csv",
    }

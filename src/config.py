from __future__ import annotations

from typing import Any


def get_default_session_state() -> dict[str, Any]:
    return {
        "analysis_name": "run_001",
        "input_dir": "",
        "output_dir": "",
        "scanned_fastq_files": [],
        "fastq_df": None,
        "sample_df": None,
        "sample_table_edits": None,
        "salmon_index_path": "",
        "tx2gene_path": "",
        "strandedness_mode": "Auto-detect",
        "strandedness_prediction": None,
        "run_validation_result": None,
        "run_status": "idle",
        "output_files": [],
        "latest_log_summary": "",
    }


def get_supported_fastq_extensions() -> tuple[str, ...]:
    return (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def get_strandedness_options() -> list[str]:
    return ["Auto-detect", "unstranded", "forward", "reverse"]


def get_default_output_filenames() -> dict[str, str]:
    return {
        "transcript_quant": "transcript_quant.csv",
        "gene_quant": "gene_quant.csv",
        "run_log": "run.log",
        "run_config": "run_config.json",
    }

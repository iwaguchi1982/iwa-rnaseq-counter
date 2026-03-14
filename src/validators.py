from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.strandedness import validate_strandedness_selection


def validate_input_directory(input_dir: str) -> dict:
    if not input_dir:
        return _invalid("入力ディレクトリを指定してください。")
    path = Path(input_dir)
    if not path.exists():
        return _invalid("入力ディレクトリが存在しません。")
    if not path.is_dir():
        return _invalid("入力パスがディレクトリではありません。")
    return _valid()


def validate_output_directory(output_dir: str) -> dict:
    if not output_dir:
        return _invalid("出力ディレクトリを指定してください。")
    path = Path(output_dir)
    parent = path if path.exists() else path.parent
    if not parent.exists():
        return _invalid("出力先の親ディレクトリが存在しません。")
    return _valid()


def validate_salmon_index(salmon_index_path: str) -> dict:
    if not salmon_index_path:
        return _invalid("Salmon index が指定されていません。")
    path = Path(salmon_index_path)
    if not path.exists() or not path.is_dir():
        return _invalid("Salmon index パスが有効ではありません。")
    return _valid()


def validate_tx2gene_file(tx2gene_path: str) -> dict:
    if not tx2gene_path:
        return _invalid("tx2gene ファイルが指定されていません。")
    path = Path(tx2gene_path)
    if not path.exists() or not path.is_file():
        return _invalid("tx2gene ファイルが存在しません。")
    return _valid()


def validate_sample_structure(sample_df: pd.DataFrame | None) -> dict:
    if sample_df is None or sample_df.empty:
        return _invalid("サンプル情報がありません。")
    if "status" in sample_df.columns and (sample_df["status"] == "error").any():
        return _invalid("サンプル構造に不整合があります。")
    return _valid()


def validate_run_conditions(
    input_dir: str,
    output_dir: str,
    sample_df: pd.DataFrame | None,
    salmon_index_path: str,
    tx2gene_path: str,
    strandedness_mode: str,
    strandedness_result: dict | None,
) -> dict:
    checks = {
        "input_dir": validate_input_directory(input_dir),
        "output_dir": validate_output_directory(output_dir),
        "sample_structure": validate_sample_structure(sample_df),
        "salmon_index": validate_salmon_index(salmon_index_path),
        "tx2gene": validate_tx2gene_file(tx2gene_path),
        "strandedness": validate_strandedness_selection(strandedness_mode, strandedness_result),
    }

    errors = []
    warnings = []
    for result in checks.values():
        errors.extend(result.get("errors", []))
        warnings.extend(result.get("warnings", []))

    return {
        "is_valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "checks": {name: result.get("is_valid", False) for name, result in checks.items()},
    }


def _valid() -> dict:
    return {"is_valid": True, "errors": [], "warnings": []}


def _invalid(message: str) -> dict:
    return {"is_valid": False, "errors": [message], "warnings": []}

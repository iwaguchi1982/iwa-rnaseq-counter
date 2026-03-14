from __future__ import annotations

import pandas as pd


def infer_strandedness(sample_df: pd.DataFrame, input_dir: str, salmon_index_path: str) -> dict:
    if sample_df is None or sample_df.empty:
        return {
            "mode": "unknown",
            "confidence": "low",
            "reason": "No sample information available",
            "is_valid": False,
            "warnings": ["サンプル情報がないため strandedness を推定できません。"],
        }

    return {
        "mode": "unknown",
        "confidence": "low",
        "reason": "Strandedness inference is not implemented yet in scaffold",
        "is_valid": True,
        "warnings": ["現在は雛形のため、実際の自動推定は未実装です。"],
    }


def validate_strandedness_selection(strandedness_mode: str, strandedness_result: dict | None) -> dict:
    if strandedness_mode == "Auto-detect" and not strandedness_result:
        return {
            "is_valid": False,
            "errors": ["Auto-detect が選択されていますが、推定結果がありません。"],
            "warnings": [],
        }
    return {
        "is_valid": True,
        "errors": [],
        "warnings": strandedness_result.get("warnings", []) if strandedness_result else [],
    }

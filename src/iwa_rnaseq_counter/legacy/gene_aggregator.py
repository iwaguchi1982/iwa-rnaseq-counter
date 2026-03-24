from __future__ import annotations

from pathlib import Path
import json
from typing import Any

import pandas as pd

from .config import get_default_output_filenames
from .qc import evaluate_sample_qc
from .run_artifacts import save_sample_metadata_csv


def load_tx2gene_map(tx2gene_path: str) -> pd.DataFrame:
    sep = "," if tx2gene_path.endswith(".csv") else "\t"
    df = pd.read_csv(tx2gene_path, sep=sep)
    if len(df.columns) >= 2:
        df = df.iloc[:, :2]
        df.columns = ["transcript_id", "gene_id"]
    return df


def aggregate_transcript_to_gene(transcript_df: pd.DataFrame, tx2gene_df: pd.DataFrame) -> pd.DataFrame:
    if transcript_df.empty:
        return pd.DataFrame()
    
    # transcript_id は index に入っている想定
    merged = transcript_df.reset_index()
    if "transcript_id" not in merged.columns and "Name" in merged.columns:
        merged = merged.rename(columns={"Name": "transcript_id"})
        
    merged = merged.merge(tx2gene_df, on="transcript_id", how="left")
    
    # aggregation
    sample_cols = [c for c in merged.columns if c not in {"transcript_id", "gene_id"}]
    gene_df = merged.groupby("gene_id", dropna=False)[sample_cols].sum(numeric_only=True).fillna(0)
    return gene_df


def build_transcript_quant_table(sample_quant_paths: list[dict], value_type: str = "TPM") -> pd.DataFrame:
    """
    value_type: "TPM" or "NumReads"
    """
    if not sample_quant_paths:
        return pd.DataFrame()

    tables = []
    for item in sample_quant_paths:
        quant_path = item["quant_path"]
        sample_id = item["sample_id"]
        if not quant_path or not Path(quant_path).exists():
            continue
        q = pd.read_csv(quant_path, sep="\t")
        if "Name" not in q.columns:
            continue
        
        if value_type not in q.columns:
            # Fallback if preferred column is missing
            alt = "NumReads" if value_type == "TPM" else "TPM"
            if alt in q.columns:
                value_col = alt
            else:
                continue
        else:
            value_col = value_type

        q = q[["Name", value_col]].rename(columns={"Name": "transcript_id", value_col: sample_id})
        tables.append(q)

    if not tables:
        return pd.DataFrame()

    merged = tables[0]
    for table in tables[1:]:
        merged = merged.merge(table, on="transcript_id", how="outer")
    return merged.set_index("transcript_id").fillna(0)


def save_quant_tables(
    matrices: dict[str, pd.DataFrame], 
    sample_df: pd.DataFrame,
    run_summary: dict,
    run_output_dir: str
) -> dict:
    """
    matrices: {"transcript_tpm": df, "transcript_numreads": df, "gene_tpm": df, "gene_numreads": df}
    """
    result_dir = Path(run_output_dir) / "results"
    result_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {}
    
    # Standardized Matrix Filenames mapping
    matrix_mapping = {
        "transcript_tpm": "transcript_tpm.csv",
        "transcript_numreads": "transcript_numreads.csv",
        "gene_tpm": "gene_tpm.csv",
        "gene_numreads": "gene_numreads.csv"
    }

    # Save Matrices
    for key, filename in matrix_mapping.items():
        if key in matrices:
            path = result_dir / filename
            matrices[key].to_csv(path)
            output_paths[key] = str(path)

    # Save sample sheet and metadata
    sample_sheet_path = result_dir / "sample_sheet.csv"
    sample_df.to_csv(sample_sheet_path, index=False)
    output_paths["sample_sheet"] = str(sample_sheet_path)

    metadata_path = save_sample_metadata_csv(result_dir, sample_df)
    output_paths["sample_metadata"] = str(metadata_path)

    # Standardized QC Summary
    qc_data = []
    for o in run_summary.get("outputs", []):
        if not o.get("is_success"): continue
        qc = evaluate_sample_qc(
            sample_id=o["sample_id"],
            num_processed=o.get("num_processed", 0),
            num_mapped=o.get("num_mapped", 0),
            num_decoy=o.get("num_decoy", 0),
            num_filter=o.get("num_filter", 0)
        )
        # MEMO.txt spec: sample_id, qc_status, mapped_fragments, mapping_rate, decoy_fragments, decoy_rate, low_score_fragments, low_score_rate, alerts
        qc_data.append({
            "sample_id": o["sample_id"],
            "qc_status": qc["status"],
            "mapped_fragments": o.get("num_mapped", 0),
            "mapping_rate": qc["mapping_rate"],
            "decoy_fragments": o.get("num_decoy", 0),
            "decoy_rate": qc["decoy_rate"],
            "low_score_fragments": o.get("num_filter", 0),
            "low_score_rate": qc["filter_rate"],
            "alerts": "; ".join(qc["alerts"]) if qc["alerts"] else "OK"
        })
    
    if qc_data:
        qc_summary_path = result_dir / "sample_qc_summary.csv"
        pd.DataFrame(qc_data).to_csv(qc_summary_path, index=False)
        output_paths["sample_qc_summary"] = str(qc_summary_path)

    # Save Run Summary JSON (naming fixed to run_summary.json)
    summary_path = result_dir / "run_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(run_summary, f, indent=2, ensure_ascii=False)
    output_paths["run_summary"] = str(summary_path)

    return output_paths

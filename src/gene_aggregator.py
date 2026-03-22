from __future__ import annotations

from pathlib import Path
import json
from typing import Any

import pandas as pd

from src.config import get_default_output_filenames
from src.qc import evaluate_sample_qc
from src.run_artifacts import save_sample_metadata_csv


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
    names = get_default_output_filenames()
    result_dir = Path(run_output_dir) / "results"
    result_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {}
    
    # Save Matrices
    for key, df in matrices.items():
        if key in names:
            path = result_dir / names[key]
            df.to_csv(path)
            output_paths[f"{key}_csv"] = str(path)

    # 4. Save sample sheet and metadata
    sample_sheet_path = result_dir / "sample_sheet.csv"
    sample_df.to_csv(sample_sheet_path, index=False)
    output_paths["sample_sheet_csv"] = str(sample_sheet_path)

    metadata_path = save_sample_metadata_csv(result_dir, sample_df)
    output_paths["sample_metadata_csv"] = str(metadata_path)

    # QC Summary
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
        qc_data.append({
            "QC": qc["icon"],
            "Sample": o["sample_id"],
            "Group": sample_df[sample_df["sample_id"] == o["sample_id"]]["group"].values[0] if "group" in sample_df.columns else "",
            "Condition": sample_df[sample_df["sample_id"] == o["sample_id"]]["condition"].values[0] if "condition" in sample_df.columns else "",
            "Mapped%": f"{qc['mapping_rate']:.2f}%",
            "Decoy%": f"{qc['decoy_rate']:.1f}%",
            "Filt%": f"{qc['filter_rate']:.1f}%",
            "Mapped Reads": o.get("num_mapped", 0),
            "Alerts": "; ".join(qc["alerts"]) if qc["alerts"] else "OK"
        })
    
    if qc_data:
        qc_summary_path = result_dir / "sample_qc_summary.csv"
        pd.DataFrame(qc_data).to_csv(qc_summary_path, index=False)
        output_paths["qc_summary_csv"] = str(qc_summary_path)

    # Save Run Summary JSON
    summary_path = result_dir / names["run_summary"]
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(run_summary, f, indent=2, ensure_ascii=False)
    output_paths["run_summary_json"] = str(summary_path)

    return output_paths

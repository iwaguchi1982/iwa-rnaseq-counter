from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.config import get_default_output_filenames


def load_tx2gene_map(tx2gene_path: str) -> pd.DataFrame:
    sep = "," if tx2gene_path.endswith(".csv") else "\t"
    df = pd.read_csv(tx2gene_path, sep=sep)
    if len(df.columns) >= 2:
        df = df.iloc[:, :2]
        df.columns = ["transcript_id", "gene_id"]
    return df


def aggregate_transcript_to_gene(sample_quant_paths: list[dict], tx2gene_df: pd.DataFrame) -> pd.DataFrame:
    merged = build_transcript_quant_table(sample_quant_paths)
    if merged.empty:
        return pd.DataFrame()
    merged = merged.reset_index().rename(columns={"index": "transcript_id"})
    merged = merged.merge(tx2gene_df, on="transcript_id", how="left")
    sample_cols = [c for c in merged.columns if c not in {"transcript_id", "gene_id"}]
    gene_df = merged.groupby("gene_id", dropna=False)[sample_cols].sum(numeric_only=True)
    return gene_df


def build_transcript_quant_table(sample_quant_paths: list[dict]) -> pd.DataFrame:
    if not sample_quant_paths:
        return pd.DataFrame()

    tables = []
    for item in sample_quant_paths:
        quant_path = item["quant_path"]
        sample_id = item["sample_id"]
        if not Path(quant_path).exists():
            continue
        q = pd.read_csv(quant_path, sep="\t")
        if "Name" not in q.columns:
            continue
        value_col = "TPM" if "TPM" in q.columns else ("NumReads" if "NumReads" in q.columns else None)
        if value_col is None:
            continue
        q = q[["Name", value_col]].rename(columns={"Name": "transcript_id", value_col: sample_id})
        tables.append(q)

    if not tables:
        return pd.DataFrame()

    merged = tables[0]
    for table in tables[1:]:
        merged = merged.merge(table, on="transcript_id", how="outer")
    return merged.set_index("transcript_id")


def save_quant_tables(transcript_df: pd.DataFrame, gene_df: pd.DataFrame, run_output_dir: str) -> dict:
    names = get_default_output_filenames()
    result_dir = Path(run_output_dir) / "results"
    result_dir.mkdir(parents=True, exist_ok=True)

    transcript_path = result_dir / names["transcript_quant"]
    gene_path = result_dir / names["gene_quant"]
    transcript_df.to_csv(transcript_path)
    gene_df.to_csv(gene_path)
    return {
        "transcript_quant_csv": str(transcript_path),
        "gene_quant_csv": str(gene_path),
    }

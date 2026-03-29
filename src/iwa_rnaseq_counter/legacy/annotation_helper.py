import pandas as pd
from pathlib import Path
import logging

# Conceptual Roles:
# 1. tx2gene: 集約用資源 (Aggregation Resource)
#    - 管轄: Counter (Salmon の結果を遺伝子単位に纏めるために必須)
# 2. feature_annotation.tsv: 表示用契約 (Reporter Display Contract)
#    - 管轄: Reporter (シンボル表示を行うためのオプション契約)
#
# 同一の由来を持つ場合でも、Counter と Reporter の間では意味を区別して扱う。
"""
Annotation Helper for iwa-rnaseq-counter (v0.5.1)
"""

logger = logging.getLogger(__name__)

def get_standard_annotation_path(outdir: Path) -> Path:
    """Return the standardized path for the reporter-ready annotation file."""
    return outdir / "results" / "feature_annotation.tsv"

def prepare_feature_annotation(tx2gene_path: str | Path, out_path: str | Path) -> bool:
    """
    Attempt to extract a v0.5.0-compliant feature_annotation.tsv from a tx2gene file.
    Standardized columns: feature_id, gene_symbol
    """
    try:
        tx2gene_path = Path(tx2gene_path)
        out_path = Path(out_path)

        if not tx2gene_path.exists():
            logger.warning(f"tx2gene file not found: {tx2gene_path}")
            return False

        # Load tx2gene
        # tx2gene is typically TAB separated but let's be flexible
        try:
            df = pd.read_csv(tx2gene_path, sep="\t")
            if len(df.columns) < 2:
                df = pd.read_csv(tx2gene_path, sep=",")
        except Exception:
            logger.warning(f"Failed to read tx2gene as TSV or CSV: {tx2gene_path}")
            return False

        if df.empty:
            logger.warning("tx2gene file is empty.")
            return False

        # Column priority for feature_id: gene_id, feature_id, ensembl_gene_id, gene
        feature_id_col = None
        for col in ["gene_id", "feature_id", "ensembl_gene_id", "gene"]:
            if col in df.columns:
                feature_id_col = col
                break
        
        if not feature_id_col:
            # Try to see if any column name contains 'gene' and 'id'
            for col in df.columns:
                if "gene" in col.lower() and "id" in col.lower():
                    feature_id_col = col
                    break
            if not feature_id_col and len(df.columns) >= 2:
                # Blind fallback: use second column if first is transcript
                feature_id_col = df.columns[1]

        if not feature_id_col:
            logger.warning("Could not identify gene/feature ID column in tx2gene.")
            return False

        # Column priority for gene_symbol: gene_symbol, symbol, gene_name
        symbol_col = None
        for col in ["gene_symbol", "symbol", "gene_name", "gene_symbol_y"]:
            if col in df.columns:
                symbol_col = col
                break

        if not symbol_col:
            logger.info("No symbol column found in tx2gene. Skipping feature_annotation.tsv generation.")
            return False

        # Unique mappings only
        annotation_df = df[[feature_id_col, symbol_col]].drop_duplicates()
        annotation_df.columns = ["feature_id", "gene_symbol"]
        
        # Ensure string IDs and clean up
        annotation_df["feature_id"] = annotation_df["feature_id"].astype(str)
        annotation_df["gene_symbol"] = annotation_df["gene_symbol"].fillna("").astype(str)

        # Final check: are there ANY non-empty symbols?
        if annotation_df["gene_symbol"].str.strip().replace("", "nan").ne("nan").any():
            out_path.parent.mkdir(parents=True, exist_ok=True)
            annotation_df.to_csv(out_path, sep="\t", index=False)
            logger.info(f"Generated feature_annotation.tsv at {out_path}")
            return True
        else:
            logger.info("All extracted gene symbols are empty. Skipping feature_annotation.tsv generation.")
            return False

    except Exception as e:
        logger.error(f"Failed to prepare feature annotation: {e}")
        return False

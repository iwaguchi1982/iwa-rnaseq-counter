from __future__ import annotations

from pathlib import Path
import pandas as pd

from .base import BaseQuantifier
from iwa_rnaseq_counter.legacy.salmon_runner import run_salmon_quant


class SalmonQuantifier(BaseQuantifier):
    name = "salmon"

    def run_quant(
        self,
        *,
        sample_df: pd.DataFrame,
        run_output_dir: str | Path,
        threads: int,
        strandedness_mode: str,
        reference_config: dict,
    ) -> dict:
        salmon_index_path = reference_config["quantifier_index"]

        return run_salmon_quant(
            sample_df=sample_df,
            salmon_index_path=salmon_index_path,
            run_output_dir=str(run_output_dir),
            strandedness_mode=strandedness_mode,
            threads=threads,
        )
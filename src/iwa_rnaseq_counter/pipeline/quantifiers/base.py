from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
import pandas as pd


class BaseQuantifier(ABC):
    name: str

    @abstractmethod
    def run_quant(
        self,
        *,
        sample_df: pd.DataFrame,
        run_output_dir: str | Path,
        threads: int,
        strandedness_mode: str,
        reference_config: dict,
    ) -> dict:
        """Return normalized raw run result."""
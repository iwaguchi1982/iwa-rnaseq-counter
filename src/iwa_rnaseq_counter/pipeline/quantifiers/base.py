from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, TypedDict

import pandas as pd


AggregationInputKind = Literal["transcript_quant", "gene_counts"]

@dataclass
class BackendCapabilities:
    aggregation_input_kind: AggregationInputKind
    has_transcript_quant: bool
    has_gene_counts: bool
    has_mapping_metrics: bool
    requires_tx2gene: bool
    requires_annotation_gtf: bool

class QuantifierOutput(TypedDict, total=False):
    # --- required-ish common keys ---
    sample_id: str
    backend: str
    is_success: bool

    # --- common artifact paths ---
    output_dir: str | None
    log_path: str | None
    transcript_quant_path: str | None
    gene_counts_path: str | None

    # --- v0.7.x compatibility alias for legacy aggregator ---
    quant_path: str | None

    # --- optional helper artifacts ---
    aux_info_dir: str | None
    meta_info_json: str | None

    # --- normalized metrics ---
    metrics: dict[str, Any]
    backend_artifacts: dict[str, Any]

    # --- legacy compatibility fields (keep for now) ---
    mapping_rate: float
    num_mapped: int
    num_processed: int
    num_decoy: int
    num_filter: int


class QuantifierRunResult(TypedDict, total=False):
    is_success: bool
    quantifier: str
    quantifier_version: str | None

    # how runner should enter aggregation
    aggregation_input_kind: AggregationInputKind

    # normalized execution/reference context
    reference_context: dict[str, Any]

    # execution status
    errors: list[str]
    warnings: list[str]
    log_summary: str
    master_log_path: str | None

    # per-sample results
    outputs: list[QuantifierOutput]


class BaseQuantifier(ABC):
    """
    Counter 内で backend 差分を吸収するための minimal public contract。

    v0.7.0 では interface を増やしすぎず、
    runner / gui_backend が必要とする最小 surface だけを固定する。
    """

    name: str = "unknown"

    @abstractmethod
    def get_capabilities(self) -> BackendCapabilities:
        """
        backend の処理能力（出力形式や必須入力要件など）を返す。
        """
        raise NotImplementedError

    @abstractmethod
    def run_quant(
        self,
        *,
        sample_df: pd.DataFrame,
        run_output_dir: str | Path,
        threads: int,
        strandedness_mode: str,
        reference_config: dict[str, Any],
    ) -> QuantifierRunResult:
        """
        backend を実行し、共通 shape の run_result を返す。
        """
        raise NotImplementedError

    def resolve_version(self) -> str | None:
        """
        backend binary の version 文字列を返す。
        v0.7.0 では失敗時 None を許容する。
        """
        return None

    def validate_environment(
        self,
        *,
        reference_config: dict[str, Any],
    ) -> list[str]:
        """
        preflight error message の list を返す。
        v0.7.0 では optional 実装でよい。
        """
        return []

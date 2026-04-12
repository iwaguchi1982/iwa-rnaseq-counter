from __future__ import annotations

from pathlib import Path
from typing import Any
import pandas as pd

from iwa_rnaseq_counter.legacy.salmon_runner import run_salmon_quant
from .base import BaseQuantifier, QuantifierRunResult, QuantifierOutput, BackendCapabilities, BackendReferenceRequirements


class SalmonQuantifier(BaseQuantifier):
    """
    既存の legacy.salmon_runner を BaseQuantifier 規格へ適合させるためのラッパー。
    v0.7.0 では実装の移行を最小限に留め、interface の正規化を優先する。
    """

    name: str = "salmon"

    def get_capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            aggregation_input_kind="transcript_quant",
            has_transcript_quant=True,
            has_gene_counts=False,
            has_mapping_metrics=True,
            reference_requirements=BackendReferenceRequirements(
                quantifier_index="required",
                tx2gene="required",
                annotation_gtf="unused"
            )
        )

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
        legacy runner を呼び出し、結果を v0.7.x 規格へ変換して返す。
        """
        # reference_config から quantifier index を取得する。
        # v0.7.0 では互換のため salmon_index_path も受け入れるが、
        # public contract 上は quantifier_index を正とする。
        quantifier_index = (
            reference_config.get("quantifier_index")
            or reference_config.get("salmon_index_path")
        )
        
        if not quantifier_index:
            return {
                "is_success": False,
                "quantifier": self.name,
                "quantifier_version": self.resolve_version(),
                "aggregation_input_kind": "transcript_quant",
                "reference_context": {},
                "errors": ["quantifier_index is required in reference_config."],
                "warnings": [],
                "log_summary": "",
                "master_log_path": None,
                "outputs": [],
            }

        # 既存ロジックの実行
        legacy_result = run_salmon_quant(
            sample_df=sample_df,
            salmon_index_path=str(quantifier_index),
            run_output_dir=str(run_output_dir),
            strandedness_mode=strandedness_mode,
            threads=threads,
        )

        # v0.7.x 形式への変換
        outputs: list[QuantifierOutput] = []
        for o in legacy_result.get("outputs", []):
            q_out: QuantifierOutput = {
                "sample_id": o["sample_id"],
                "backend": self.name,
                "is_success": o["is_success"],
                "output_dir": str(Path(o["quant_path"]).parent) if o.get("quant_path") else None,
                "log_path": o.get("log_path"),
                "transcript_quant_path": o.get("quant_path"),  # Salmon は transcript-level
                "quant_path": o.get("quant_path"),             # 互換用
                "aux_info_dir": o.get("aux_info_dir"),
                "meta_info_json": str(Path(o["aux_info_dir"]) / "meta_info.json") if o.get("aux_info_dir") else None,
                "mapping_rate": o.get("mapping_rate", 0.0),
                "num_mapped": o.get("num_mapped", 0),
                "num_processed": o.get("num_processed", 0),
                "num_decoy": o.get("num_decoy", 0),
                "num_filter": o.get("num_filter", 0),
                "metrics": {
                    "mapping_rate": o.get("mapping_rate"),
                },
                "backend_artifacts": {}
            }
            outputs.append(q_out)

        return {
            "is_success": legacy_result["is_success"],
            "quantifier": self.name,
            "quantifier_version": self.resolve_version(),
            "aggregation_input_kind": "transcript_quant",
            "reference_context": {
                "quantifier_index_path": str(quantifier_index),
                "tx2gene_path": (
                    str(reference_config.get("tx2gene_path"))
                    if reference_config.get("tx2gene_path")
                    else None
                ),
            },
            "errors": legacy_result.get("errors", []),
            "warnings": legacy_result.get("warnings", []),
            "log_summary": legacy_result.get("log_summary", ""),
            "master_log_path": legacy_result.get("master_log_path"),
            "outputs": outputs,
        }

    def resolve_version(self) -> str | None:
        # v0.7.0 では固定値を返すが、将来的に binary から取得するように拡張可能
        return "1.10.1"

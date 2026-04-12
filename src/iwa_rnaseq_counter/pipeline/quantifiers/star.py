from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from .base import BaseQuantifier, QuantifierOutput, QuantifierRunResult, BackendCapabilities, BackendReferenceRequirements


class StarQuantifier(BaseQuantifier):
    """
    STAR backend adapter.

    v0.7.1 では STAR を gene-level backend として扱い、
    ReadsPerGene.out.tab を backend 内で正規化して Counter 共通契約へ載せる。
    """

    name: str = "star"

    def get_capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            aggregation_input_kind="gene_counts",
            has_transcript_quant=False,
            has_gene_counts=True,
            has_mapping_metrics=False,
            reference_requirements=BackendReferenceRequirements(
                quantifier_index="required",
                tx2gene="unused",
                annotation_gtf="unused"
            )
        )

    def resolve_version(self) -> str | None:
        try:
            completed = subprocess.run(
                ["STAR", "--version"],
                capture_output=True,
                text=True,
                check=False,
            )
            raw = (completed.stdout or completed.stderr or "").strip()
            if not raw:
                return None

            first_line = raw.splitlines()[0].strip()
            m = re.search(r"(\d+\.\d+(?:\.\d+)?)", first_line)
            return m.group(1) if m else first_line
        except Exception:
            return None

    def validate_environment(
        self,
        *,
        reference_config: dict[str, Any],
    ) -> list[str]:
        errors: list[str] = []

        if shutil.which("STAR") is None:
            errors.append("STAR binary was not found in PATH.")

        quantifier_index = reference_config.get("quantifier_index")
        if not quantifier_index:
            errors.append("quantifier_index is required for STAR.")
        elif not Path(str(quantifier_index)).exists():
            errors.append(f"STAR genome index does not exist: {quantifier_index}")

        return errors

    def _resolve_read_files(self, row: pd.Series) -> list[str]:
        layout_final = str(row.get("layout_final", "single-end")).strip().lower()

        if layout_final == "paired-end":
            r1_paths = row.get("r1_paths", []) or []
            r2_paths = row.get("r2_paths", []) or []
            if not r1_paths or not r2_paths:
                raise ValueError("paired-end sample requires both r1_paths and r2_paths")
            return [str(r1_paths[0]), str(r2_paths[0])]

        all_paths = row.get("all_paths", []) or []
        if all_paths:
            return [str(all_paths[0])]

        r1_paths = row.get("r1_paths", []) or []
        if r1_paths:
            return [str(r1_paths[0])]

        raise ValueError("single-end sample requires all_paths or r1_paths")

    def _needs_zcat(self, read_files: list[str]) -> bool:
        return any(str(p).endswith(".gz") for p in read_files)

    def _select_star_count_column(self, strandedness_mode: str) -> int:
        """
        ReadsPerGene.out.tab:
          col1 = gene_id
          col2 = unstranded
          col3 = 1st-strand
          col4 = 2nd-strand

        v0.7.1 の最小版では次の簡易対応にする。
        必要なら後続で library convention に合わせて調整する。
        """
        mode = str(strandedness_mode or "").strip().lower()

        if mode in {"", "auto-detect", "auto", "unknown", "unstranded", "none"}:
            return 1
        if mode in {"forward", "fr", "fr-secondstrand", "secondstrand"}:
            return 2
        if mode in {"reverse", "rf", "fr-firststrand", "firststrand"}:
            return 3

        return 1

    def _normalize_reads_per_gene(
        self,
        raw_gene_counts_path: Path,
        normalized_gene_counts_path: Path,
        strandedness_mode: str,
    ) -> None:
        col_idx = self._select_star_count_column(strandedness_mode)

        df = pd.read_csv(
            raw_gene_counts_path,
            sep="\t",
            header=None,
            comment=None,
        )

        # STAR の先頭特殊行を除外
        df = df[~df[0].astype(str).str.startswith("N_")].copy()

        norm_df = pd.DataFrame(
            {
                "feature_id": df[0].astype(str),
                "count": pd.to_numeric(df[col_idx], errors="coerce").fillna(0).astype(int),
            }
        )
        norm_df.to_csv(normalized_gene_counts_path, sep="\t", index=False)

    def run_quant(
        self,
        *,
        sample_df: pd.DataFrame,
        run_output_dir: str | Path,
        threads: int,
        strandedness_mode: str,
        reference_config: dict[str, Any],
    ) -> QuantifierRunResult:
        quantifier_index = reference_config.get("quantifier_index")
        tx2gene_path = reference_config.get("tx2gene_path")

        preflight_errors = self.validate_environment(reference_config=reference_config)
        if preflight_errors:
            return {
                "is_success": False,
                "quantifier": self.name,
                "quantifier_version": self.resolve_version(),
                "aggregation_input_kind": "gene_counts",
                "reference_context": {},
                "errors": preflight_errors,
                "warnings": [],
                "log_summary": "",
                "master_log_path": None,
                "outputs": [],
            }

        run_output_dir = Path(run_output_dir)
        outputs: list[QuantifierOutput] = []
        errors: list[str] = []
        warnings: list[str] = []

        for _, row in sample_df.iterrows():
            sample_id = str(row["sample_id"])
            sample_outdir = run_output_dir / "star" / sample_id
            sample_outdir.mkdir(parents=True, exist_ok=True)

            log_path = sample_outdir / "star.log"
            raw_gene_counts_path = sample_outdir / "ReadsPerGene.out.tab"
            normalized_gene_counts_path = sample_outdir / "gene_counts.tsv"

            try:
                read_files = self._resolve_read_files(row)

                cmd = [
                    "STAR",
                    "--runThreadN",
                    str(threads),
                    "--genomeDir",
                    str(quantifier_index),
                    "--readFilesIn",
                    *read_files,
                    "--quantMode",
                    "GeneCounts",
                    "--outSAMtype",
                    "None",
                    "--outFileNamePrefix",
                    str(sample_outdir) + "/",
                ]

                if self._needs_zcat(read_files):
                    cmd.extend(["--readFilesCommand", "zcat"])

                completed = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                )

                log_text = ""
                if completed.stdout:
                    log_text += completed.stdout
                if completed.stderr:
                    if log_text:
                        log_text += "\n"
                    log_text += completed.stderr
                log_path.write_text(log_text, encoding="utf-8")

                if completed.returncode != 0:
                    outputs.append(
                        {
                            "sample_id": sample_id,
                            "backend": self.name,
                            "is_success": False,
                            "output_dir": str(sample_outdir),
                            "log_path": str(log_path),
                            "gene_counts_path": None,
                            "transcript_quant_path": None,
                            "metrics": {},
                            "backend_artifacts": {
                                "star_returncode": completed.returncode,
                            },
                        }
                    )
                    errors.append(f"STAR failed for sample {sample_id}")
                    continue

                if not raw_gene_counts_path.exists():
                    outputs.append(
                        {
                            "sample_id": sample_id,
                            "backend": self.name,
                            "is_success": False,
                            "output_dir": str(sample_outdir),
                            "log_path": str(log_path),
                            "gene_counts_path": None,
                            "transcript_quant_path": None,
                            "metrics": {},
                            "backend_artifacts": {},
                        }
                    )
                    errors.append(
                        f"STAR completed but ReadsPerGene.out.tab was not found for sample {sample_id}"
                    )
                    continue

                self._normalize_reads_per_gene(
                    raw_gene_counts_path=raw_gene_counts_path,
                    normalized_gene_counts_path=normalized_gene_counts_path,
                    strandedness_mode=strandedness_mode,
                )

                outputs.append(
                    {
                        "sample_id": sample_id,
                        "backend": self.name,
                        "is_success": True,
                        "output_dir": str(sample_outdir),
                        "log_path": str(log_path),
                        "gene_counts_path": str(normalized_gene_counts_path),
                        "transcript_quant_path": None,
                        "metrics": {},
                        "backend_artifacts": {
                            "raw_gene_counts_path": str(raw_gene_counts_path),
                        },
                    }
                )

            except Exception as e:
                try:
                    log_path.write_text(str(e), encoding="utf-8")
                except Exception:
                    pass

                outputs.append(
                    {
                        "sample_id": sample_id,
                        "backend": self.name,
                        "is_success": False,
                        "output_dir": str(sample_outdir),
                        "log_path": str(log_path),
                        "gene_counts_path": None,
                        "transcript_quant_path": None,
                        "metrics": {},
                        "backend_artifacts": {},
                    }
                )
                errors.append(f"STAR exception for sample {sample_id}: {e}")

        is_success = any(o.get("is_success") for o in outputs)

        return {
            "is_success": is_success,
            "quantifier": self.name,
            "quantifier_version": self.resolve_version(),
            "aggregation_input_kind": "gene_counts",
            "reference_context": {
                "quantifier_index_path": str(quantifier_index),
                "tx2gene_path": str(tx2gene_path) if tx2gene_path else None,
            },
            "errors": errors,
            "warnings": warnings,
            "log_summary": f"STAR processed {len(outputs)} sample(s), success={sum(bool(o.get('is_success')) for o in outputs)}",
            "master_log_path": None,
            "outputs": outputs,
        }

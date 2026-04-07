from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from .base import BaseQuantifier, QuantifierOutput, QuantifierRunResult


class KallistoQuantifier(BaseQuantifier):
    """
    kallisto backend adapter.

    v0.7.3 では kallisto を transcript-level backend として扱い、
    abundance.tsv を quant.sf 互換の TSV に正規化して Counter 共通契約へ載せる。
    """

    name: str = "kallisto"

    def resolve_version(self) -> str | None:
        try:
            completed = subprocess.run(
                ["kallisto", "version"],
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

        if shutil.which("kallisto") is None:
            errors.append("kallisto binary was not found in PATH.")

        quantifier_index = reference_config.get("quantifier_index")
        if not quantifier_index:
            errors.append("quantifier_index is required for kallisto.")
        else:
            index_path = Path(str(quantifier_index))
            if not index_path.exists() or not index_path.is_file():
                errors.append(f"kallisto index file does not exist: {quantifier_index}")

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

    def _normalize_abundance_tsv(
        self,
        abundance_tsv_path: Path,
        quant_sf_compatible_path: Path,
    ) -> None:
        df = pd.read_csv(abundance_tsv_path, sep="\t")

        required = {"target_id", "length", "eff_length", "est_counts", "tpm"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"kallisto abundance.tsv is missing required columns: {', '.join(sorted(missing))}"
            )

        out = pd.DataFrame(
            {
                "Name": df["target_id"].astype(str),
                "Length": pd.to_numeric(df["length"], errors="coerce").fillna(0),
                "EffectiveLength": pd.to_numeric(df["eff_length"], errors="coerce").fillna(0),
                "TPM": pd.to_numeric(df["tpm"], errors="coerce").fillna(0),
                "NumReads": pd.to_numeric(df["est_counts"], errors="coerce").fillna(0),
            }
        )
        out.to_csv(quant_sf_compatible_path, sep="\t", index=False)

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
                "aggregation_input_kind": "transcript_quant",
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

        single_fragment_length = int(reference_config.get("single_fragment_length", 200))
        single_fragment_sd = int(reference_config.get("single_fragment_sd", 20))

        for _, row in sample_df.iterrows():
            sample_id = str(row["sample_id"])
            layout_final = str(row.get("layout_final", "single-end")).strip().lower()
            sample_outdir = run_output_dir / "kallisto" / sample_id
            sample_outdir.mkdir(parents=True, exist_ok=True)

            log_path = sample_outdir / "kallisto.log"
            abundance_tsv_path = sample_outdir / "abundance.tsv"
            quant_sf_path = sample_outdir / "quant.sf"

            try:
                read_files = self._resolve_read_files(row)

                cmd = [
                    "kallisto",
                    "quant",
                    "-i",
                    str(quantifier_index),
                    "-o",
                    str(sample_outdir),
                    "-t",
                    str(threads),
                    "--plaintext",
                ]

                if layout_final == "paired-end":
                    cmd.extend(read_files)
                else:
                    warnings.append(
                        f"kallisto single-end defaults were used for sample {sample_id}: "
                        f"fragment_length={single_fragment_length}, fragment_sd={single_fragment_sd}"
                    )
                    cmd.extend(
                        [
                            "--single",
                            "-l",
                            str(single_fragment_length),
                            "-s",
                            str(single_fragment_sd),
                            read_files[0],
                        ]
                    )

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
                            "transcript_quant_path": None,
                            "quant_path": None,
                            "gene_counts_path": None,
                            "metrics": {},
                            "backend_artifacts": {
                                "kallisto_returncode": completed.returncode,
                            },
                        }
                    )
                    errors.append(f"kallisto failed for sample {sample_id}")
                    continue

                if not abundance_tsv_path.exists():
                    outputs.append(
                        {
                            "sample_id": sample_id,
                            "backend": self.name,
                            "is_success": False,
                            "output_dir": str(sample_outdir),
                            "log_path": str(log_path),
                            "transcript_quant_path": None,
                            "quant_path": None,
                            "gene_counts_path": None,
                            "metrics": {},
                            "backend_artifacts": {},
                        }
                    )
                    errors.append(
                        f"kallisto completed but abundance.tsv was not found for sample {sample_id}"
                    )
                    continue

                self._normalize_abundance_tsv(
                    abundance_tsv_path=abundance_tsv_path,
                    quant_sf_compatible_path=quant_sf_path,
                )

                outputs.append(
                    {
                        "sample_id": sample_id,
                        "backend": self.name,
                        "is_success": True,
                        "output_dir": str(sample_outdir),
                        "log_path": str(log_path),
                        "transcript_quant_path": str(quant_sf_path),
                        "quant_path": str(quant_sf_path),  # legacy/current aggregator compatibility
                        "gene_counts_path": None,
                        "metrics": {},
                        "backend_artifacts": {
                            "abundance_tsv_path": str(abundance_tsv_path),
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
                        "transcript_quant_path": None,
                        "quant_path": None,
                        "gene_counts_path": None,
                        "metrics": {},
                        "backend_artifacts": {},
                    }
                )
                errors.append(f"kallisto exception for sample {sample_id}: {e}")

        is_success = any(o.get("is_success") for o in outputs)

        return {
            "is_success": is_success,
            "quantifier": self.name,
            "quantifier_version": self.resolve_version(),
            "aggregation_input_kind": "transcript_quant",
            "reference_context": {
                "quantifier_index_path": str(quantifier_index),
                "tx2gene_path": str(tx2gene_path) if tx2gene_path else None,
                "single_fragment_length": single_fragment_length,
                "single_fragment_sd": single_fragment_sd,
            },
            "errors": errors,
            "warnings": warnings,
            "log_summary": (
                f"kallisto processed {len(outputs)} sample(s), "
                f"success={sum(bool(o.get('is_success')) for o in outputs)}"
            ),
            "master_log_path": None,
            "outputs": outputs,
        }

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from .base import BaseQuantifier, QuantifierOutput, QuantifierRunResult, BackendCapabilities, BackendReferenceRequirements


class Hisat2Quantifier(BaseQuantifier):
    """
    HISAT2 backend adapter.

    v0.7.2 最小版では:
      1) hisat2 で alignment
      2) samtools sort で BAM 化
      3) featureCounts で gene-level count 作成
      4) gene_counts.tsv を Counter 共通契約へ載せる
    """

    name: str = "hisat2"

    def get_capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            aggregation_input_kind="gene_counts",
            has_transcript_quant=False,
            has_gene_counts=True,
            has_mapping_metrics=False,
            reference_requirements=BackendReferenceRequirements(
                quantifier_index="required",
                tx2gene="unused",
                annotation_gtf="required"
            )
        )

    def resolve_version(self) -> str | None:
        try:
            completed = subprocess.run(
                ["hisat2", "--version"],
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

        if shutil.which("hisat2") is None:
            errors.append("hisat2 binary was not found in PATH.")
        if shutil.which("samtools") is None:
            errors.append("samtools binary was not found in PATH.")
        if shutil.which("featureCounts") is None:
            errors.append("featureCounts binary was not found in PATH.")

        quantifier_index = reference_config.get("quantifier_index")
        annotation_gtf_path = reference_config.get("annotation_gtf_path")

        if not quantifier_index:
            errors.append("quantifier_index is required for HISAT2.")
        else:
            # hisat2 index prefix を最低限確認
            prefix = Path(str(quantifier_index))
            parent = prefix.parent if prefix.parent != Path("") else Path(".")
            stem = prefix.name
            candidates = list(parent.glob(f"{stem}*.ht2")) + list(parent.glob(f"{stem}*.ht2l"))
            if not candidates:
                errors.append(f"HISAT2 index prefix looks invalid: {quantifier_index}")

        if not annotation_gtf_path:
            errors.append("annotation_gtf_path is required for HISAT2 gene-level counting.")
        elif not Path(str(annotation_gtf_path)).exists():
            errors.append(f"annotation_gtf_path does not exist: {annotation_gtf_path}")

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

    def _featurecounts_strand_arg(self, strandedness_mode: str) -> str:
        """
        featureCounts -s:
          0 = unstranded
          1 = stranded
          2 = reversely stranded

        v0.7.2 最小版の簡易対応。
        必要なら後続で library convention を厳密化する。
        """
        mode = str(strandedness_mode or "").strip().lower()

        if mode in {"", "auto-detect", "auto", "unknown", "unstranded", "none"}:
            return "0"
        if mode in {"forward", "fr", "fr-secondstrand", "secondstrand"}:
            return "1"
        if mode in {"reverse", "rf", "fr-firststrand", "firststrand"}:
            return "2"

        return "0"

    def _normalize_featurecounts_output(
        self,
        raw_counts_path: Path,
        normalized_counts_path: Path,
    ) -> None:
        df = pd.read_csv(raw_counts_path, sep="\t", comment="#")

        # featureCounts 標準出力:
        # Geneid Chr Start End Strand Length sample.bam
        if "Geneid" not in df.columns:
            raise ValueError(f"featureCounts output does not contain 'Geneid': {raw_counts_path}")

        sample_count_columns = [
            c for c in df.columns
            if c not in {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
        ]
        if not sample_count_columns:
            raise ValueError(f"No count column found in featureCounts output: {raw_counts_path}")

        count_col = sample_count_columns[-1]

        norm_df = pd.DataFrame(
            {
                "feature_id": df["Geneid"].astype(str),
                "count": pd.to_numeric(df[count_col], errors="coerce").fillna(0).astype(int),
            }
        )
        norm_df.to_csv(normalized_counts_path, sep="\t", index=False)

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
        annotation_gtf_path = reference_config.get("annotation_gtf_path")

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
            sample_outdir = run_output_dir / "hisat2" / sample_id
            sample_outdir.mkdir(parents=True, exist_ok=True)

            sam_path = sample_outdir / "aligned.sam"
            bam_path = sample_outdir / "aligned.sorted.bam"
            raw_counts_path = sample_outdir / "featurecounts.txt"
            gene_counts_path = sample_outdir / "gene_counts.tsv"
            log_path = sample_outdir / "hisat2.log"

            try:
                read_files = self._resolve_read_files(row)
                layout_final = str(row.get("layout_final", "single-end")).strip().lower()

                hisat2_cmd = [
                    "hisat2",
                    "-x",
                    str(quantifier_index),
                    "-p",
                    str(threads),
                    "-S",
                    str(sam_path),
                ]
                if layout_final == "paired-end":
                    hisat2_cmd.extend(["-1", read_files[0], "-2", read_files[1]])
                else:
                    hisat2_cmd.extend(["-U", read_files[0]])

                hisat2_completed = subprocess.run(
                    hisat2_cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                )

                log_text = ""
                if hisat2_completed.stdout:
                    log_text += hisat2_completed.stdout
                if hisat2_completed.stderr:
                    if log_text:
                        log_text += "\n"
                    log_text += hisat2_completed.stderr

                if hisat2_completed.returncode != 0:
                    log_path.write_text(log_text, encoding="utf-8")
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
                                "hisat2_returncode": hisat2_completed.returncode,
                            },
                        }
                    )
                    errors.append(f"HISAT2 failed for sample {sample_id}")
                    continue

                samtools_cmd = [
                    "samtools",
                    "sort",
                    "-@",
                    str(threads),
                    "-o",
                    str(bam_path),
                    str(sam_path),
                ]
                samtools_completed = subprocess.run(
                    samtools_cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                )

                if samtools_completed.stdout:
                    log_text += "\n" + samtools_completed.stdout
                if samtools_completed.stderr:
                    log_text += "\n" + samtools_completed.stderr

                if samtools_completed.returncode != 0:
                    log_path.write_text(log_text, encoding="utf-8")
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
                                "samtools_returncode": samtools_completed.returncode,
                            },
                        }
                    )
                    errors.append(f"samtools sort failed for sample {sample_id}")
                    continue

                featurecounts_cmd = [
                    "featureCounts",
                    "-T",
                    str(threads),
                    "-a",
                    str(annotation_gtf_path),
                    "-o",
                    str(raw_counts_path),
                    "-s",
                    self._featurecounts_strand_arg(strandedness_mode),
                    str(bam_path),
                ]
                fc_completed = subprocess.run(
                    featurecounts_cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                )

                if fc_completed.stdout:
                    log_text += "\n" + fc_completed.stdout
                if fc_completed.stderr:
                    log_text += "\n" + fc_completed.stderr

                log_path.write_text(log_text, encoding="utf-8")

                if fc_completed.returncode != 0:
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
                                "featurecounts_returncode": fc_completed.returncode,
                            },
                        }
                    )
                    errors.append(f"featureCounts failed for sample {sample_id}")
                    continue

                if not raw_counts_path.exists():
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
                    errors.append(f"featureCounts output was not found for sample {sample_id}")
                    continue

                self._normalize_featurecounts_output(
                    raw_counts_path=raw_counts_path,
                    normalized_counts_path=gene_counts_path,
                )

                outputs.append(
                    {
                        "sample_id": sample_id,
                        "backend": self.name,
                        "is_success": True,
                        "output_dir": str(sample_outdir),
                        "log_path": str(log_path),
                        "gene_counts_path": str(gene_counts_path),
                        "transcript_quant_path": None,
                        "metrics": {},
                        "backend_artifacts": {
                            "sam_path": str(sam_path),
                            "bam_path": str(bam_path),
                            "featurecounts_raw_path": str(raw_counts_path),
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
                errors.append(f"HISAT2 exception for sample {sample_id}: {e}")

        is_success = any(o.get("is_success") for o in outputs)

        return {
            "is_success": is_success,
            "quantifier": self.name,
            "quantifier_version": self.resolve_version(),
            "aggregation_input_kind": "gene_counts",
            "reference_context": {
                "quantifier_index_path": str(quantifier_index),
                "tx2gene_path": str(tx2gene_path) if tx2gene_path else None,
                "annotation_gtf_path": str(annotation_gtf_path) if annotation_gtf_path else None,
            },
            "errors": errors,
            "warnings": warnings,
            "log_summary": f"HISAT2 processed {len(outputs)} sample(s), success={sum(bool(o.get('is_success')) for o in outputs)}",
            "master_log_path": None,
            "outputs": outputs,
        }

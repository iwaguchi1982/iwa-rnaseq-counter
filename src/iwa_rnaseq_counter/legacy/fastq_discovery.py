from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import get_supported_fastq_extensions


def discover_fastq_files(input_dir: str) -> list[str]:
    if not input_dir:
        return []
    base = Path(input_dir)
    if not base.exists() or not base.is_dir():
        return []

    extensions = get_supported_fastq_extensions()
    files = [str(p) for p in base.rglob("*") if p.is_file() and p.name.endswith(extensions)]
    return sorted(files)


def collect_fastq_metadata(fastq_files: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for file_path in fastq_files:
        p = Path(file_path)
        rows.append(
            {
                "file_name": p.name,
                "path": str(p),
                "size_bytes": p.stat().st_size if p.exists() else 0,
                "extension": _detect_extension(p.name),
                "is_gzipped": p.name.endswith(".gz"),
                "detected_sample_token": _detect_sample_token(p.name),
                "detected_read_token": _detect_read_token(p.name),
                "detected_lane_token": _detect_lane_token(p.name),
            }
        )
    return pd.DataFrame(rows)


def _detect_extension(filename: str) -> str:
    for ext in get_supported_fastq_extensions():
        if filename.endswith(ext):
            return ext
    return ""


def _detect_read_token(filename: str) -> str | None:
    for token in ("_R1", "_R2", "_1", "_2"):
        if token in filename:
            return token.replace("_", "")
    return None


def _detect_lane_token(filename: str) -> str | None:
    for i in range(1, 9):
        token = f"L00{i}"
        if token in filename:
            return token
    return None


def _detect_sample_token(filename: str) -> str:
    name = filename
    for token in ("_R1", "_R2", "_1", "_2", ".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        name = name.replace(token, "")
    for i in range(1, 9):
        name = name.replace(f"_L00{i}", "")
    return name.strip("_-")

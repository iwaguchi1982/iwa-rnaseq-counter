from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from ..io.read_matrix_spec import read_matrix_spec


@dataclass
class AnalysisBundlePaths:
    manifest_path: str
    matrix_spec_path: str | None
    execution_run_spec_path: str | None
    merged_matrix_path: str | None
    sample_metadata_input_path: str | None
    aligned_sample_metadata_path: str | None
    analysis_merge_summary_path: str | None
    log_path: str | None
    feature_annotation_path: str | None


@dataclass
class AnalysisBundle:
    manifest: dict[str, Any]
    paths: AnalysisBundlePaths
    matrix_spec: Any | None
    execution_run_spec: dict[str, Any] | None
    analysis_merge_summary: dict[str, Any] | None
    aligned_sample_metadata: pd.DataFrame | None
    merged_matrix: pd.DataFrame | None


def _read_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _resolve_manifest_path(path_or_dir: str | Path) -> Path:
    """
    analysis_bundle_manifest.json 自体、または bundle root directory を受け取る。
    """
    p = Path(path_or_dir)

    if p.is_file():
        return p.resolve()

    manifest_path = p / "results" / "analysis_bundle_manifest.json"
    if manifest_path.exists():
        return manifest_path.resolve()

    raise FileNotFoundError(
        "analysis bundle manifest not found. "
        f"Expected file or bundle root containing results/analysis_bundle_manifest.json: {p}"
    )


def _resolve_path_value(value: str | None, manifest_path: Path) -> str | None:
    if value is None:
        return None

    raw = str(value).strip()
    if not raw:
        return None

    p = Path(raw)
    if p.is_absolute():
        return str(p.resolve())

    # relative path は manifest file 基準で解決
    return str((manifest_path.parent / p).resolve())


def _resolve_bundle_paths(manifest: dict[str, Any], manifest_path: Path) -> AnalysisBundlePaths:
    paths = manifest.get("paths", {})

    return AnalysisBundlePaths(
        manifest_path=str(manifest_path.resolve()),
        matrix_spec_path=_resolve_path_value(paths.get("matrix_spec"), manifest_path),
        execution_run_spec_path=_resolve_path_value(paths.get("execution_run_spec"), manifest_path),
        merged_matrix_path=_resolve_path_value(paths.get("merged_matrix"), manifest_path),
        sample_metadata_input_path=_resolve_path_value(paths.get("sample_metadata_input"), manifest_path),
        aligned_sample_metadata_path=_resolve_path_value(paths.get("aligned_sample_metadata"), manifest_path),
        analysis_merge_summary_path=_resolve_path_value(paths.get("analysis_merge_summary"), manifest_path),
        log_path=_resolve_path_value(paths.get("log"), manifest_path),
        feature_annotation_path=_resolve_path_value(paths.get("feature_annotation"), manifest_path),
    )


def _read_execution_run_spec_dict(path_value: str | None) -> dict[str, Any] | None:
    if not path_value:
        return None

    path = Path(path_value)
    if not path.exists():
        return None

    return _read_json(path)


def _read_analysis_merge_summary_dict(path_value: str | None) -> dict[str, Any] | None:
    if not path_value:
        return None

    path = Path(path_value)
    if not path.exists():
        return None

    return _read_json(path)


def _read_optional_table(path_value: str | None, sep: str = "\t") -> pd.DataFrame | None:
    if not path_value:
        return None

    path = Path(path_value)
    if not path.exists():
        return None

    return pd.read_csv(path, sep=sep)


def read_analysis_bundle(
    path_or_dir: str | Path,
    *,
    load_matrix_spec: bool = True,
    load_execution_run_spec: bool = True,
    load_analysis_merge_summary: bool = True,
    load_aligned_sample_metadata: bool = True,
    load_merged_matrix: bool = False,
) -> AnalysisBundle:
    """
    v0.9.1-2:
    analysis bundle manifest から主要成果物をまとめて読む。

    Parameters
    ----------
    path_or_dir:
        analysis_bundle_manifest.json そのもの、
        または bundle root directory.
    load_matrix_spec:
        matrix.spec.json を読むか
    load_execution_run_spec:
        execution-run.spec.json を読むか
    load_analysis_merge_summary:
        analysis_merge_summary.json を読むか
    load_aligned_sample_metadata:
        aligned_sample_metadata.tsv を読むか
    load_merged_matrix:
        merged_gene_numreads.tsv を読むか（重いので default False）
    """
    manifest_path = _resolve_manifest_path(path_or_dir)
    manifest = _read_json(manifest_path)
    paths = _resolve_bundle_paths(manifest, manifest_path)

    matrix_spec_obj = None
    if load_matrix_spec and paths.matrix_spec_path:
        matrix_spec_path = Path(paths.matrix_spec_path)
        if matrix_spec_path.exists():
            matrix_spec_obj = read_matrix_spec(matrix_spec_path)

    execution_run_spec_obj = None
    if load_execution_run_spec:
        execution_run_spec_obj = _read_execution_run_spec_dict(paths.execution_run_spec_path)

    analysis_merge_summary_obj = None
    if load_analysis_merge_summary:
        analysis_merge_summary_obj = _read_analysis_merge_summary_dict(paths.analysis_merge_summary_path)

    aligned_sample_metadata_df = None
    if load_aligned_sample_metadata:
        aligned_sample_metadata_df = _read_optional_table(paths.aligned_sample_metadata_path, sep="\t")

    merged_matrix_df = None
    if load_merged_matrix:
        merged_matrix_df = _read_optional_table(paths.merged_matrix_path, sep="\t")

    return AnalysisBundle(
        manifest=manifest,
        paths=paths,
        matrix_spec=matrix_spec_obj,
        execution_run_spec=execution_run_spec_obj,
        analysis_merge_summary=analysis_merge_summary_obj,
        aligned_sample_metadata=aligned_sample_metadata_df,
        merged_matrix=merged_matrix_df,
    )

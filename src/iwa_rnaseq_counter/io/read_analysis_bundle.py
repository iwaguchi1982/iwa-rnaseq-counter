import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from ..io.read_matrix_spec import read_matrix_spec
from ..models.execution_run import ExecutionRunSpec
from ..models.matrix import MatrixSpec


@dataclass(frozen=True)
class AnalysisBundleContractInfo:
    contract_name: str
    contract_version: str
    bundle_kind: str
    producer: str | None
    producer_version: str | None
    is_supported: bool
    compatibility_status: str


@dataclass
class AnalysisBundlePaths:
    manifest_path: Path
    bundle_root: Path
    matrix_spec_path: Path
    execution_run_spec_path: Path
    merged_matrix_path: Path
    aligned_sample_metadata_path: Path
    analysis_merge_summary_path: Path
    build_analysis_matrix_log_path: Path
    feature_annotation_path: Path | None = None


@dataclass
class AnalysisBundle:
    manifest: dict[str, Any]
    paths: AnalysisBundlePaths
    matrix_spec: MatrixSpec
    execution_run_spec: ExecutionRunSpec
    analysis_merge_summary: dict[str, Any]
    aligned_sample_metadata: pd.DataFrame
    contract_info: AnalysisBundleContractInfo
    merged_matrix: pd.DataFrame | None = None


def _read_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _resolve_manifest_path(path_like: str | Path) -> Path:
    """
    受け口は 2 通り許容する:
      1. results/analysis_bundle_manifest.json
      2. analysis bundle の outdir
    """
    path = Path(path_like)

    if path.is_dir():
        manifest_path = path / "results" / "analysis_bundle_manifest.json"
    else:
        manifest_path = path

    if not manifest_path.exists() or not manifest_path.is_file():
        raise FileNotFoundError(f"analysis bundle manifest not found: {manifest_path}")

    return manifest_path.resolve()


def _resolve_bundle_root(manifest: dict[str, Any], manifest_path: Path) -> Path:
    """
    manifest に bundle_root があればそれを優先。
    無ければ results/analysis_bundle_manifest.json の 2 つ上を bundle root とみなす。
    """
    bundle_root_raw = manifest.get("bundle_root")
    if bundle_root_raw:
        return Path(str(bundle_root_raw)).resolve()

    return manifest_path.parent.parent.resolve()


def _extract_artifact_relpath(
    manifest: dict[str, Any],
    artifact_name: str,
) -> str | None:
    """
    manifest["artifacts"][name] は以下のどちらでも受ける:
      - {"path": "...", ...}
      - "..."
    """
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict):
        raise ValueError("analysis bundle manifest must contain dict field 'artifacts'")

    entry = artifacts.get(artifact_name)
    if entry is None:
        return None

    if isinstance(entry, str):
        value = entry.strip()
        return value or None

    if isinstance(entry, dict):
        raw = entry.get("path")
        if raw is None:
            return None
        value = str(raw).strip()
        return value or None

    raise TypeError(
        f"artifact entry for {artifact_name!r} must be dict or str, got {type(entry).__name__}"
    )


def _resolve_artifact_path(
    manifest: dict[str, Any],
    *,
    artifact_name: str,
    bundle_root: Path,
    required: bool = True,
) -> Path | None:
    raw_path = _extract_artifact_relpath(manifest, artifact_name)

    if raw_path is None:
        if required:
            raise ValueError(
                f"required artifact {artifact_name!r} is missing in analysis bundle manifest"
            )
        return None

    path = Path(raw_path)
    if not path.is_absolute():
        path = bundle_root / path

    return path.resolve()


def _read_execution_run_spec(file_path: str | Path) -> ExecutionRunSpec:
    """
    現状の repo では ExecutionRunSpec に from_dict() がないため、
    bundle loader 側で最小の reader を持つ。
    """
    data = _read_json(Path(file_path))

    schema_name = data.get("$schema_name", data.get("schema_name"))
    if schema_name != "ExecutionRunSpec":
        raise ValueError(f"Expected ExecutionRunSpec, got {schema_name!r}")

    return ExecutionRunSpec(
        schema_name=data.get("$schema_name", data.get("schema_name", "")),
        schema_version=data.get("$schema_version", data.get("schema_version", "")),
        run_id=data.get("run_id", ""),
        app_name=data.get("app_name", ""),
        app_version=data.get("app_version", ""),
        started_at=data.get("started_at", ""),
        input_refs=data.get("input_refs", []) or [],
        output_refs=data.get("output_refs", []) or [],
        parameters=data.get("parameters", {}) or {},
        execution_backend=data.get("execution_backend"),
        finished_at=data.get("finished_at"),
        status=data.get("status", "pending"),
        log_path=data.get("log_path", ""),
        metadata=data.get("metadata", {}) or {},
        overlay=data.get("overlay", {}) or {},
    )


def _read_tabular_file(path: Path) -> pd.DataFrame:
    ...

def _parse_major_version(version_text: str) -> int:
    first = str(version_text).split(".", 1)[0]
    return int(first)


def evaluate_analysis_bundle_contract(
    manifest: dict[str, Any],
) -> AnalysisBundleContractInfo:
    contract_name = str(manifest.get("contract_name") or "")
    contract_version = str(manifest.get("contract_version") or "")
    bundle_kind = str(manifest.get("bundle_kind") or "")
    producer = manifest.get("producer")
    producer_version = manifest.get("producer_version")

    if contract_name != "analysis_bundle":
        return AnalysisBundleContractInfo(
            contract_name=contract_name,
            contract_version=contract_version,
            bundle_kind=bundle_kind,
            producer=producer,
            producer_version=producer_version,
            is_supported=False,
            compatibility_status="invalid_contract_name",
        )

    if bundle_kind != "rna_seq_analysis_bundle":
        return AnalysisBundleContractInfo(
            contract_name=contract_name,
            contract_version=contract_version,
            bundle_kind=bundle_kind,
            producer=producer,
            producer_version=producer_version,
            is_supported=False,
            compatibility_status="unsupported_bundle_kind",
        )

    if not contract_version:
        return AnalysisBundleContractInfo(
            contract_name=contract_name,
            contract_version=contract_version,
            bundle_kind=bundle_kind,
            producer=producer,
            producer_version=producer_version,
            is_supported=False,
            compatibility_status="missing_contract_version",
        )

    try:
        major = _parse_major_version(contract_version)
    except Exception:
        return AnalysisBundleContractInfo(
            contract_name=contract_name,
            contract_version=contract_version,
            bundle_kind=bundle_kind,
            producer=producer,
            producer_version=producer_version,
            is_supported=False,
            compatibility_status="invalid_contract_version",
        )

    if major != 1:
        return AnalysisBundleContractInfo(
            contract_name=contract_name,
            contract_version=contract_version,
            bundle_kind=bundle_kind,
            producer=producer,
            producer_version=producer_version,
            is_supported=False,
            compatibility_status="unsupported_contract_major",
        )

    return AnalysisBundleContractInfo(
        contract_name=contract_name,
        contract_version=contract_version,
        bundle_kind=bundle_kind,
        producer=producer,
        producer_version=producer_version,
        is_supported=True,
        compatibility_status="supported",
    )


def _raise_if_unsupported_analysis_bundle_contract(
    contract_info: AnalysisBundleContractInfo,
) -> None:
    if contract_info.is_supported:
        return

    raise ValueError(
        "unsupported analysis bundle contract: "
        f"status={contract_info.compatibility_status}, "
        f"contract_name={contract_info.contract_name!r}, "
        f"contract_version={contract_info.contract_version!r}, "
        f"bundle_kind={contract_info.bundle_kind!r}"
    )
    suffix = path.suffix.lower()

    if suffix in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t")
        if len(df.columns) == 1:
            return pd.read_csv(path)
        return df

    if suffix == ".csv":
        df = pd.read_csv(path)
        if len(df.columns) == 1:
            return pd.read_csv(path, sep="\t")
        return df

    try:
        df = pd.read_csv(path)
        if len(df.columns) > 1:
            return df
    except Exception:
        pass

    return pd.read_csv(path, sep="\t")


def resolve_analysis_bundle_paths(
    manifest_path: str | Path,
) -> AnalysisBundlePaths:
    manifest_path_obj = _resolve_manifest_path(manifest_path)
    manifest = _read_json(manifest_path_obj)

    schema_name = manifest.get("schema_name")
    if schema_name and schema_name != "AnalysisBundleManifest":
        raise ValueError(f"Expected AnalysisBundleManifest, got {schema_name!r}")

    bundle_root = _resolve_bundle_root(manifest, manifest_path_obj)

    return AnalysisBundlePaths(
        manifest_path=manifest_path_obj,
        bundle_root=bundle_root,
        matrix_spec_path=_resolve_artifact_path(
            manifest, artifact_name="matrix_spec", bundle_root=bundle_root, required=True
        ),
        execution_run_spec_path=_resolve_artifact_path(
            manifest,
            artifact_name="execution_run_spec",
            bundle_root=bundle_root,
            required=True,
        ),
        merged_matrix_path=_resolve_artifact_path(
            manifest, artifact_name="merged_matrix", bundle_root=bundle_root, required=True
        ),
        aligned_sample_metadata_path=_resolve_artifact_path(
            manifest,
            artifact_name="aligned_sample_metadata",
            bundle_root=bundle_root,
            required=True,
        ),
        analysis_merge_summary_path=_resolve_artifact_path(
            manifest,
            artifact_name="analysis_merge_summary",
            bundle_root=bundle_root,
            required=True,
        ),
        build_analysis_matrix_log_path=_resolve_artifact_path(
            manifest,
            artifact_name="build_analysis_matrix_log",
            bundle_root=bundle_root,
            required=True,
        ),
        feature_annotation_path=_resolve_artifact_path(
            manifest,
            artifact_name="feature_annotation",
            bundle_root=bundle_root,
            required=False,
        ),
    )


def read_analysis_bundle(
    manifest_path: str | Path,
    *,
    load_merged_matrix: bool = False,
) -> AnalysisBundle:
    """
    analysis_bundle_manifest.json を入口に、
    downstream が必要とする analysis handoff artifact 群をまとめて読む。

    Parameters
    ----------
    manifest_path:
        以下のどちらでもよい:
        - results/analysis_bundle_manifest.json
        - analysis bundle の outdir
    load_merged_matrix:
        True の時だけ merged matrix を読み込む。
        count matrix は重くなりうるため default=False。

    Returns
    -------
    AnalysisBundle
    """
    paths = resolve_analysis_bundle_paths(manifest_path)
    manifest = _read_json(paths.manifest_path)

    contract_info = evaluate_analysis_bundle_contract(manifest)
    _raise_if_unsupported_analysis_bundle_contract(contract_info)

    matrix_spec = read_matrix_spec(paths.matrix_spec_path)
    execution_run_spec = _read_execution_run_spec(paths.execution_run_spec_path)
    analysis_merge_summary = _read_json(paths.analysis_merge_summary_path)
    aligned_sample_metadata = _read_tabular_file(paths.aligned_sample_metadata_path)

    merged_matrix = None
    if load_merged_matrix:
        merged_matrix = pd.read_csv(paths.merged_matrix_path, sep="\t", index_col=0)

    return AnalysisBundle(
        manifest=manifest,
        paths=paths,
        matrix_spec=matrix_spec,
        execution_run_spec=execution_run_spec,
        analysis_merge_summary=analysis_merge_summary,
        aligned_sample_metadata=aligned_sample_metadata,
        contract_info=contract_info,
        merged_matrix=merged_matrix,
    )

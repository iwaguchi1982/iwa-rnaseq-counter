import json
import logging
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from importlib.metadata import PackageNotFoundError, version

import pandas as pd

from ..models.execution_run import ExecutionRunSpec
from ..models.matrix import MatrixSpec

logger = logging.getLogger(__name__)


def _canonical_source_context(spec: MatrixSpec) -> dict:
    """
    source MatrixSpec から downstream が読みやすい provenance を抽出する。
    """
    metadata = spec.metadata or {}

    return {
        "matrix_id": spec.matrix_id,
        "matrix_path": spec.matrix_path,
        "matrix_scope": spec.matrix_scope,
        "matrix_kind": spec.matrix_kind,
        "feature_type": spec.feature_type,
        "value_type": spec.value_type,
        "normalization": spec.normalization,
        "sample_axis": spec.sample_axis,
        "feature_id_system": spec.feature_id_system,
        "feature_annotation_path": spec.feature_annotation_path,
        "source_assay_ids": list(spec.source_assay_ids or []),
        "source_specimen_ids": list(spec.source_specimen_ids or []),
        "source_subject_ids": list(spec.source_subject_ids or []),
        "producer_app": metadata.get("producer_app"),
        "producer_version": metadata.get("producer_version"),
        "quantifier": metadata.get("quantifier"),
        "quantifier_version": metadata.get("quantifier_version"),
        "aggregation_input_kind": metadata.get("aggregation_input_kind"),
        "quantifier_index_path": metadata.get("quantifier_index_path"),
        "tx2gene_path": metadata.get("tx2gene_path"),
        "annotation_gtf_path": metadata.get("annotation_gtf_path"),
        "reference_context": metadata.get("reference_context", {}),
        "feature_annotation_available": metadata.get("feature_annotation_available"),
    }


def _unique_preserve_order(values: list) -> list:
    seen = set()
    out = []
    for v in values:
        key = str(v)
        if key in seen:
            continue
        seen.add(key)
        out.append(v)
    return out


RECOMMENDED_SAMPLE_METADATA_COLUMNS = (
    "subject_id",
    "condition",
    "group",
    "batch",
)

ANALYSIS_BUNDLE_SCHEMA_NAME = "AnalysisBundleManifest"
ANALYSIS_BUNDLE_SCHEMA_VERSION = "0.2.0"

ANALYSIS_BUNDLE_CONTRACT_NAME = "analysis_bundle"
ANALYSIS_BUNDLE_CONTRACT_VERSION = "1.0.0"
ANALYSIS_BUNDLE_KIND = "rna_seq_analysis_bundle"
ANALYSIS_BUNDLE_PRODUCER = "iwa_rnaseq_counter"


def _get_counter_producer_version() -> str:
    try:
        return version("iwa-rnaseq-counter")
    except PackageNotFoundError:
        return "unknown"
    except Exception:
        return "unknown"


def _build_analysis_bundle_contract_block(*, producer_version: str) -> dict[str, Any]:
    return {
        "schema_name": ANALYSIS_BUNDLE_SCHEMA_NAME,
        "schema_version": ANALYSIS_BUNDLE_SCHEMA_VERSION,
        "contract_name": ANALYSIS_BUNDLE_CONTRACT_NAME,
        "contract_version": ANALYSIS_BUNDLE_CONTRACT_VERSION,
        "bundle_kind": ANALYSIS_BUNDLE_KIND,
        "producer": ANALYSIS_BUNDLE_PRODUCER,
        "producer_version": producer_version,
    }


def _inspect_recommended_sample_metadata_columns(
    columns: Sequence[str] | pd.Index | None,
) -> dict[str, Any]:
    """
    sample metadata の推奨列充足状況を軽く点検する。
    欠けていても error にはせず、warning / summary 用の情報として返す。
    """
    recommended = list(RECOMMENDED_SAMPLE_METADATA_COLUMNS)

    if columns is None:
        return {
            "status": "not_available",
            "recommended": recommended,
            "present": [],
            "missing": recommended,
            "available_columns": [],
        }

    available_columns = [str(col) for col in list(columns)]
    present = [col for col in recommended if col in available_columns]
    missing = [col for col in recommended if col not in available_columns]

    if not missing:
        status = "complete"
    elif present:
        status = "partial"
    else:
        status = "missing"

    return {
        "status": status,
        "recommended": recommended,
        "present": present,
        "missing": missing,
        "available_columns": available_columns,
    }


def _normalize_optional_path(value: str | None) -> str | None:
    if value is None:
        return None
    s = str(value).strip()
    return s if s else None


def _inspect_feature_annotation_file(path_value: str | None) -> dict[str, Any]:
    """
    feature_annotation.tsv を analysis handoff 用 artifact として使ってよいかを検査する。
    標準列: feature_id, gene_symbol
    """
    normalized = _normalize_optional_path(path_value)
    if not normalized:
        return {
            "path": None,
            "exists": False,
            "readable": False,
            "is_usable": False,
            "columns": [],
            "row_count": 0,
            "status": "missing",
            "reason": "path_missing",
        }

    path = Path(normalized)
    if not path.exists() or not path.is_file():
        return {
            "path": str(path),
            "exists": False,
            "readable": False,
            "is_usable": False,
            "columns": [],
            "row_count": 0,
            "status": "missing",
            "reason": "file_not_found",
        }

    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        return {
            "path": str(path),
            "exists": True,
            "readable": False,
            "is_usable": False,
            "columns": [],
            "row_count": 0,
            "status": "invalid",
            "reason": f"read_error:{e}",
        }

    columns = [str(c) for c in df.columns]
    required_columns = {"feature_id", "gene_symbol"}
    has_required_columns = required_columns.issubset(set(columns))

    return {
        "path": str(path),
        "exists": True,
        "readable": True,
        "is_usable": has_required_columns,
        "columns": columns,
        "row_count": int(len(df)),
        "status": "usable" if has_required_columns else "invalid",
        "reason": "ok" if has_required_columns else "missing_required_columns",
    }


def _validate_mergeable_matrix_specs(matrix_specs: list[MatrixSpec]) -> None:
    """
    analysis merge 前に、列結合してよい最低限の共通性を確認する。
    """
    if not matrix_specs:
        raise ValueError("matrix_specs must not be empty")

    first = matrix_specs[0]
    expected = {
        "matrix_kind": first.matrix_kind,
        "feature_type": first.feature_type,
        "value_type": first.value_type,
        "normalization": first.normalization,
        "sample_axis": first.sample_axis,
        "feature_id_system": first.feature_id_system,
    }

    if first.sample_axis != "specimen":
        raise ValueError(
            f"build_analysis_matrix currently requires sample_axis='specimen', got {first.sample_axis!r}"
        )

    specimen_ids_all: list[str] = []

    for spec in matrix_specs:
        observed = {
            "matrix_kind": spec.matrix_kind,
            "feature_type": spec.feature_type,
            "value_type": spec.value_type,
            "normalization": spec.normalization,
            "sample_axis": spec.sample_axis,
            "feature_id_system": spec.feature_id_system,
        }

        mismatches = {
            k: (expected[k], observed[k])
            for k in expected
            if expected[k] != observed[k]
        }
        if mismatches:
            raise ValueError(
                f"MatrixSpec {spec.matrix_id} is not merge-compatible. mismatches={mismatches}"
            )

        if not spec.source_specimen_ids:
            raise ValueError(f"MatrixSpec {spec.matrix_id} has no source_specimen_ids")

        specimen_ids_all.extend(spec.source_specimen_ids)

    duplicated = sorted({sid for sid in specimen_ids_all if specimen_ids_all.count(sid) > 1})
    if duplicated:
        raise ValueError(
            f"Duplicate specimen IDs detected across source matrices: {duplicated}"
        )


def _build_merge_provenance(matrix_specs: list[MatrixSpec]) -> dict:
    source_contexts = [_canonical_source_context(spec) for spec in matrix_specs]

    quantifiers = _unique_preserve_order(
        [ctx["quantifier"] for ctx in source_contexts if ctx.get("quantifier")]
    )
    quantifier_versions = _unique_preserve_order(
        [ctx["quantifier_version"] for ctx in source_contexts if ctx.get("quantifier_version")]
    )
    aggregation_input_kinds = _unique_preserve_order(
        [ctx["aggregation_input_kind"] for ctx in source_contexts if ctx.get("aggregation_input_kind")]
    )
    quantifier_index_paths = _unique_preserve_order(
        [
            _normalize_optional_path(ctx.get("quantifier_index_path"))
            for ctx in source_contexts
            if _normalize_optional_path(ctx.get("quantifier_index_path"))
        ]
    )
    tx2gene_paths = _unique_preserve_order(
        [
            _normalize_optional_path(ctx.get("tx2gene_path"))
            for ctx in source_contexts
            if _normalize_optional_path(ctx.get("tx2gene_path"))
        ]
    )
    annotation_gtf_paths = _unique_preserve_order(
        [
            _normalize_optional_path(ctx.get("annotation_gtf_path"))
            for ctx in source_contexts
            if _normalize_optional_path(ctx.get("annotation_gtf_path"))
        ]
    )

    annotation_inspections = [
        _inspect_feature_annotation_file(ctx.get("feature_annotation_path"))
        for ctx in source_contexts
    ]

    existing_annotation_paths = _unique_preserve_order(
        [
            item["path"]
            for item in annotation_inspections
            if item.get("path") and item.get("exists")
        ]
    )
    usable_annotation_paths = _unique_preserve_order(
        [
            item["path"]
            for item in annotation_inspections
            if item.get("path") and item.get("is_usable")
        ]
    )
    usable_column_signatures = _unique_preserve_order(
        [
            tuple(item.get("columns", []))
            for item in annotation_inspections
            if item.get("is_usable")
        ]
    )

    missing_count = sum(1 for item in annotation_inspections if item.get("status") == "missing")
    invalid_count = sum(1 for item in annotation_inspections if item.get("status") == "invalid")
    usable_count = sum(1 for item in annotation_inspections if item.get("status") == "usable")

    if usable_count == 0:
        feature_annotation_consensus_status = "missing_or_invalid"
        feature_annotation_consensus_path = None
        feature_annotation_is_usable = False
    elif len(usable_annotation_paths) == 1 and len(usable_column_signatures) == 1:
        feature_annotation_consensus_status = "consistent"
        feature_annotation_consensus_path = usable_annotation_paths[0]
        feature_annotation_is_usable = True
    elif len(usable_column_signatures) == 1:
        feature_annotation_consensus_status = "path_inconsistent_but_schema_consistent"
        feature_annotation_consensus_path = usable_annotation_paths[0]
        feature_annotation_is_usable = True
    else:
        feature_annotation_consensus_status = "schema_inconsistent"
        feature_annotation_consensus_path = None
        feature_annotation_is_usable = False

    feature_annotation_consensus_reason = {
        "missing_count": missing_count,
        "invalid_count": invalid_count,
        "usable_count": usable_count,
        "existing_paths": existing_annotation_paths,
        "usable_paths": usable_annotation_paths,
        "usable_column_signatures": [list(sig) for sig in usable_column_signatures],
    }

    return {
        "source_matrix_count": len(matrix_specs),
        "source_quantifiers": quantifiers,
        "source_quantifier_versions": quantifier_versions,
        "source_aggregation_input_kinds": aggregation_input_kinds,
        "source_feature_annotation_paths": existing_annotation_paths,
        "source_feature_annotation_inspections": annotation_inspections,
        "source_quantifier_index_paths": quantifier_index_paths,
        "source_tx2gene_paths": tx2gene_paths,
        "source_annotation_gtf_paths": annotation_gtf_paths,
        "feature_annotation_consensus_status": feature_annotation_consensus_status,
        "feature_annotation_consensus_path": feature_annotation_consensus_path,
        "feature_annotation_consensus_reason": feature_annotation_consensus_reason,
        "feature_annotation_is_usable": feature_annotation_is_usable,
        "source_matrix_contexts": source_contexts,
    }


def _read_sample_metadata_table(sample_metadata_path: Path) -> pd.DataFrame:
    """
    sample metadata を TSV / CSV のどちらでも読み込めるようにする。
    """
    if not sample_metadata_path.exists() or not sample_metadata_path.is_file():
        raise FileNotFoundError(f"sample metadata file not found: {sample_metadata_path}")

    suffix = sample_metadata_path.suffix.lower()

    if suffix in {".tsv", ".txt"}:
        df = pd.read_csv(sample_metadata_path, sep="\t")
        if len(df.columns) == 1:
            df = pd.read_csv(sample_metadata_path)
        return df

    if suffix == ".csv":
        df = pd.read_csv(sample_metadata_path)
        if len(df.columns) == 1:
            df = pd.read_csv(sample_metadata_path, sep="\t")
        return df

    try:
        df = pd.read_csv(sample_metadata_path)
        if len(df.columns) > 1:
            return df
    except Exception:
        pass

    return pd.read_csv(sample_metadata_path, sep="\t")


def _resolve_sample_metadata_id_column(df: pd.DataFrame) -> str:
    """
    specimen axis に整列するための ID 列を決める。
    """
    priority = ["specimen_id", "sample_id", "sample"]
    for col in priority:
        if col in df.columns:
            return col

    raise ValueError(
        "sample metadata must include one of the following ID columns: "
        "specimen_id, sample_id, sample"
    )


def _validate_and_align_sample_metadata(
    sample_metadata_path: Path,
    specimen_ids_in_order: list[str],
    metadata_dir: Path | None = None,
) -> dict[str, Any]:
    """
    merged matrix の specimen 列順に sample metadata を整列する。
    metadata_dir が与えられた場合のみ aligned artifact を保存する。

    error:
      - ID 列が無い
      - duplicate ID がある
      - merged matrix 側 specimen_id が metadata に存在しない

    allow:
      - metadata 側に余分な row がある
    """
    df = _read_sample_metadata_table(sample_metadata_path)
    if df.empty:
        raise ValueError("sample metadata is empty")

    recommended_columns = _inspect_recommended_sample_metadata_columns(df.columns)
    id_col = _resolve_sample_metadata_id_column(df)

    df = df.copy()
    df[id_col] = df[id_col].astype(str).str.strip()

    duplicate_ids = sorted(df.loc[df[id_col].duplicated(), id_col].unique().tolist())
    if duplicate_ids:
        raise ValueError(
            f"sample metadata contains duplicate IDs in column '{id_col}': {duplicate_ids}"
        )

    metadata_ids = df[id_col].tolist()
    metadata_id_set = set(metadata_ids)
    specimen_id_set = set(specimen_ids_in_order)

    missing_specimen_ids = [sid for sid in specimen_ids_in_order if sid not in metadata_id_set]
    extra_metadata_ids = [sid for sid in metadata_ids if sid not in specimen_id_set]

    if missing_specimen_ids:
        raise ValueError(
            "sample metadata is missing rows for merged specimen IDs: "
            f"{missing_specimen_ids}"
        )

    aligned_df = df.set_index(id_col).loc[specimen_ids_in_order].reset_index()

    if id_col != "specimen_id":
        aligned_df.insert(0, "specimen_id", aligned_df[id_col].astype(str))

    aligned_path = None
    if metadata_dir is not None:
        metadata_dir.mkdir(parents=True, exist_ok=True)
        out_path = metadata_dir / "aligned_sample_metadata.tsv"
        aligned_df.to_csv(out_path, sep="\t", index=False)
        aligned_path = str(out_path.resolve())

    return {
        "status": "aligned",
        "id_column": id_col,
        "aligned_sample_metadata_path": aligned_path,
        "row_count_input": int(len(df)),
        "row_count_aligned": int(len(aligned_df)),
        "missing_specimen_ids": missing_specimen_ids,
        "extra_metadata_ids": extra_metadata_ids,
        "duplicate_ids": duplicate_ids,
        "aligned_columns": [str(c) for c in aligned_df.columns],
        "recommended_columns": recommended_columns,
    }


def _load_matrix_for_merge(spec: MatrixSpec) -> tuple[pd.DataFrame, str]:
    """
    単一 assay matrix を specimen_id 列へ正規化して返す。
    """
    df = pd.read_csv(spec.matrix_path, sep="\t", index_col=0)

    if not spec.source_specimen_ids:
        raise ValueError(f"MatrixSpec {spec.matrix_id} has no source_specimen_ids")

    specimen_id = str(spec.source_specimen_ids[0])

    if df.shape[1] != 1:
        if specimen_id not in df.columns:
            raise ValueError(
                f"MatrixSpec {spec.matrix_id} expected single-sample matrix or column named {specimen_id}, "
                f"got columns={list(df.columns)}"
            )
        df = df[[specimen_id]]
    else:
        df.columns = [specimen_id]

    return df, specimen_id


def _collect_analysis_merge_warnings(
    merge_provenance: dict[str, Any],
    sample_metadata_alignment: dict[str, Any],
) -> list[str]:
    warnings: list[str] = []

    annotation_status = merge_provenance.get("feature_annotation_consensus_status")
    if annotation_status and annotation_status != "consistent":
        warnings.append(
            f"feature_annotation_consensus_status={annotation_status}"
        )

    extra_metadata_ids = sample_metadata_alignment.get("extra_metadata_ids") or []
    if extra_metadata_ids:
        warnings.append(
            f"sample metadata contains extra rows not used in merged matrix: {extra_metadata_ids}"
        )

    source_quantifiers = merge_provenance.get("source_quantifiers") or []
    if len(source_quantifiers) > 1:
        warnings.append(
            f"multiple quantifiers are mixed in analysis merge: {source_quantifiers}"
        )

    recommended_columns = sample_metadata_alignment.get("recommended_columns", {}) or {}
    missing_recommended_columns = recommended_columns.get("missing", []) or []
    recommended_status = recommended_columns.get("status")

    if recommended_status != "not_available" and missing_recommended_columns:
        warnings.append(
            "sample metadata is missing recommended columns: "
            + ", ".join(missing_recommended_columns)
        )

    return warnings


def preview_build_analysis_matrix(
    matrix_specs: list[MatrixSpec],
    sample_metadata_path: Path,
    outdir: Path | None = None,
    matrix_id: str | None = None,
    run_id: str | None = None,
) -> dict[str, Any]:
    """
    v0.9.1-4:
    build-analysis-matrix の dry-run 用 preview。
    書き込みは行わず、merge readiness に加えて
    analysis bundle の planned manifest shape も返す。
    """
    _validate_mergeable_matrix_specs(matrix_specs)

    specimen_ids_in_order: list[str] = []
    for spec in matrix_specs:
        if not spec.source_specimen_ids:
            raise ValueError(f"MatrixSpec {spec.matrix_id} has no source_specimen_ids")
        specimen_ids_in_order.append(str(spec.source_specimen_ids[0]))

    merge_provenance = _build_merge_provenance(matrix_specs)
    sample_metadata_alignment = _validate_and_align_sample_metadata(
        sample_metadata_path=sample_metadata_path,
        specimen_ids_in_order=specimen_ids_in_order,
        metadata_dir=None,
    )
    warnings = _collect_analysis_merge_warnings(
        merge_provenance=merge_provenance,
        sample_metadata_alignment=sample_metadata_alignment,
    )

    preview: dict[str, Any] = {
        "status": "preview_ok",
        "source_matrix_count": len(matrix_specs),
        "source_specimen_ids": specimen_ids_in_order,
        "source_quantifiers": merge_provenance["source_quantifiers"],
        "source_quantifier_versions": merge_provenance["source_quantifier_versions"],
        "source_aggregation_input_kinds": merge_provenance["source_aggregation_input_kinds"],
        "feature_annotation_consensus_status": merge_provenance["feature_annotation_consensus_status"],
        "feature_annotation_consensus_path": merge_provenance["feature_annotation_consensus_path"],
        "feature_annotation_consensus_reason": merge_provenance["feature_annotation_consensus_reason"],
        "feature_annotation_is_usable": merge_provenance["feature_annotation_is_usable"],
        "sample_metadata_id_column": sample_metadata_alignment["id_column"],
        "sample_metadata_alignment_status": sample_metadata_alignment["status"],
        "sample_metadata_row_count_input": sample_metadata_alignment["row_count_input"],
        "sample_metadata_row_count_aligned": sample_metadata_alignment["row_count_aligned"],
        "sample_metadata_extra_ids": sample_metadata_alignment["extra_metadata_ids"],
        "sample_metadata_recommended_columns": sample_metadata_alignment.get(
            "recommended_columns", {}
        ),
        "warning_count": len(warnings),
        "warnings": warnings,
    }

    if outdir is not None:
        preview_matrix_id = matrix_id or "PREVIEW_ANALYSIS_MATRIX"
        preview_run_id = run_id or (
            f"PREVIEW_BUILD_ANALYSIS_MATRIX_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
        preview["analysis_bundle"] = _build_analysis_bundle_preview(
            outdir=outdir,
            matrix_id=preview_matrix_id,
            run_id=preview_run_id,
            feature_annotation_path=merge_provenance["feature_annotation_consensus_path"],
        )

    return preview


def _build_analysis_merge_summary(
    *,
    matrix_id: str,
    matrix_path: Path,
    sample_metadata_path: Path,
    merge_provenance: dict[str, Any],
    sample_metadata_alignment: dict[str, Any],
    warnings: list[str],
    run_id: str,
    output_refs: list[str],
    input_refs: list[str],
    bundle_manifest_path: Path | None = None,
) -> dict[str, Any]:
    producer_version = _get_counter_producer_version()

    return {
        "schema_name": "AnalysisMergeSummary",
        "schema_version": "0.1.0",
        "status": "completed",
        "analysis_bundle_contract": {
            "contract_name": ANALYSIS_BUNDLE_CONTRACT_NAME,
            "contract_version": ANALYSIS_BUNDLE_CONTRACT_VERSION,
            "bundle_kind": ANALYSIS_BUNDLE_KIND,
            "producer": ANALYSIS_BUNDLE_PRODUCER,
            "producer_version": producer_version,
        },
        "matrix_id": matrix_id,
        "matrix_path": str(matrix_path.resolve()),
        "sample_metadata_path": str(sample_metadata_path.resolve()),
        "aligned_sample_metadata_path": sample_metadata_alignment.get("aligned_sample_metadata_path"),
        "run_id": run_id,
        "source_matrix_count": merge_provenance["source_matrix_count"],
        "source_quantifiers": merge_provenance["source_quantifiers"],
        "source_quantifier_versions": merge_provenance["source_quantifier_versions"],
        "source_aggregation_input_kinds": merge_provenance["source_aggregation_input_kinds"],
        "feature_annotation_consensus_status": merge_provenance["feature_annotation_consensus_status"],
        "feature_annotation_consensus_path": merge_provenance["feature_annotation_consensus_path"],
        "feature_annotation_consensus_reason": merge_provenance["feature_annotation_consensus_reason"],
        "feature_annotation_is_usable": merge_provenance["feature_annotation_is_usable"],
        "sample_metadata_alignment_status": sample_metadata_alignment["status"],
        "sample_metadata_id_column": sample_metadata_alignment["id_column"],
        "sample_metadata_row_count_input": sample_metadata_alignment["row_count_input"],
        "sample_metadata_row_count_aligned": sample_metadata_alignment["row_count_aligned"],
        "sample_metadata_extra_ids": sample_metadata_alignment["extra_metadata_ids"],
        "sample_metadata_recommended_columns": sample_metadata_alignment.get(
            "recommended_columns", {}
        ),
        "warning_count": len(warnings),
        "warnings": warnings,
        "analysis_bundle_manifest_path": (
            str(bundle_manifest_path.resolve()) if bundle_manifest_path else None
        ),
        "analysis_bundle_entrypoint_kind": "analysis_bundle_manifest",
        "analysis_bundle_artifact_paths": (
            analysis_bundle.get("artifact_paths", {}) if analysis_bundle else {}
        ),
        "input_refs": input_refs,
        "output_refs": output_refs,
    }


def _write_analysis_merge_summary(summary_path: Path, summary: dict[str, Any]) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)


def _write_analysis_merge_log(
    log_path: Path,
    *,
    matrix_id: str,
    source_matrix_count: int,
    source_quantifiers: list[str],
    feature_annotation_consensus_status: str,
    sample_metadata_alignment_status: str,
    warnings: list[str],
    analysis_bundle_manifest_path: str | None = None,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "[build_analysis_matrix]",
        "status=completed",
        f"matrix_id={matrix_id}",
        f"source_matrix_count={source_matrix_count}",
        f"source_quantifiers={source_quantifiers}",
        f"feature_annotation_consensus_status={feature_annotation_consensus_status}",
        f"sample_metadata_alignment_status={sample_metadata_alignment_status}",
        f"warning_count={len(warnings)}",
        f"analysis_bundle_contract_name={ANALYSIS_BUNDLE_CONTRACT_NAME}",
        f"analysis_bundle_contract_version={ANALYSIS_BUNDLE_CONTRACT_VERSION}",
        f"analysis_bundle_producer_version={_get_counter_producer_version()}",
    ]

    if analysis_bundle_manifest_path:
        lines.append(f"analysis_bundle_manifest_path={analysis_bundle_manifest_path}")
        lines.append("analysis_bundle_entrypoint_kind=analysis_bundle_manifest")

    if warnings:
        lines.append("warnings:")
        for w in warnings:
            lines.append(f"- {w}")
    else:
        lines.append("warnings: []")

    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _bundle_ref_path(path_value: str | Path | None, *, bundle_root: Path) -> str | None:
    """
    manifest 用の path 正規化。
    bundle_root 配下なら relative path、それ以外は absolute path で返す。
    """
    if path_value is None:
        return None

    path = path_value if isinstance(path_value, Path) else Path(str(path_value))
    if not path.is_absolute():
        path = bundle_root / path

    resolved = path.resolve()
    bundle_root_resolved = bundle_root.resolve()

    try:
        return str(resolved.relative_to(bundle_root_resolved))
    except ValueError:
        return str(resolved)


def _build_analysis_bundle_manifest(
    *,
    outdir: Path,
    matrix_id: str,
    run_id: str,
    matrix_spec_path: Path,
    execution_run_spec_path: Path,
    merged_matrix_path: Path,
    aligned_sample_metadata_path: str | None,
    analysis_merge_summary_path: Path,
    build_analysis_matrix_log_path: Path,
    feature_annotation_path: str | None,
) -> dict[str, Any]:
    """
    analysis handoff artifact 群を 1 枚の入口から辿れるようにする manifest を作る。
    v0.9.1-1:
      - build_analysis_matrix.py 側で manifest を出力
      - spec JSON 自体は現行どおり CLI 側で保存される前提のため、
        ここでは existence ではなく artifact ref registry として path を束ねる
    """
    bundle_root = outdir.resolve()

    artifacts = {
        "matrix_spec": {
            "path": _bundle_ref_path(matrix_spec_path, bundle_root=bundle_root),
            "kind": "spec",
            "schema_name": "MatrixSpec",
            "required": True,
        },
        "execution_run_spec": {
            "path": _bundle_ref_path(execution_run_spec_path, bundle_root=bundle_root),
            "kind": "spec",
            "schema_name": "ExecutionRunSpec",
            "required": True,
        },
        "merged_matrix": {
            "path": _bundle_ref_path(merged_matrix_path, bundle_root=bundle_root),
            "kind": "matrix",
            "required": True,
        },
        "aligned_sample_metadata": {
            "path": _bundle_ref_path(aligned_sample_metadata_path, bundle_root=bundle_root),
            "kind": "metadata",
            "required": True,
        },
        "analysis_merge_summary": {
            "path": _bundle_ref_path(analysis_merge_summary_path, bundle_root=bundle_root),
            "kind": "summary",
            "required": True,
        },
        "build_analysis_matrix_log": {
            "path": _bundle_ref_path(build_analysis_matrix_log_path, bundle_root=bundle_root),
            "kind": "log",
            "required": True,
        },
    }

    if feature_annotation_path:
        artifacts["feature_annotation"] = {
            "path": _bundle_ref_path(feature_annotation_path, bundle_root=bundle_root),
            "kind": "annotation",
            "required": False,
        }

    contract_block = _build_analysis_bundle_contract_block(
        producer_version=_get_counter_producer_version(),
    )

    return {
        **contract_block,
        "bundle_scope": "analysis",
        "matrix_id": matrix_id,
        "run_id": run_id,
        "bundle_root": str(bundle_root),
        "created_at": datetime.now(timezone.utc).astimezone().isoformat(),
        "artifacts": artifacts,
    }


def _write_analysis_bundle_manifest(manifest_path: Path, manifest: dict[str, Any]) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)


def _build_analysis_bundle_preview(
    *,
    outdir: Path,
    matrix_id: str,
    run_id: str,
    feature_annotation_path: str | None,
) -> dict[str, Any]:
    """
    dry-run / summary / log からも見やすいように、
    analysis bundle の entrypoint と artifact shape を組み立てる。
    実ファイルの existence check はせず、bundle contract の planned shape を返す。
    """
    outdir = outdir.resolve()
    counts_dir = outdir / "counts"
    logs_dir = outdir / "logs"
    specs_dir = outdir / "specs"
    metadata_dir = outdir / "metadata"
    results_dir = outdir / "results"

    manifest_path = results_dir / "analysis_bundle_manifest.json"

    manifest = _build_analysis_bundle_manifest(
        outdir=outdir,
        matrix_id=matrix_id,
        run_id=run_id,
        matrix_spec_path=specs_dir / "matrix.spec.json",
        execution_run_spec_path=specs_dir / "execution-run.spec.json",
        merged_matrix_path=counts_dir / "merged_gene_numreads.tsv",
        aligned_sample_metadata_path=metadata_dir / "aligned_sample_metadata.tsv",
        analysis_merge_summary_path=results_dir / "analysis_merge_summary.json",
        build_analysis_matrix_log_path=logs_dir / "build_analysis_matrix.log",
        feature_annotation_path=feature_annotation_path,
    )

    artifacts = manifest.get("artifacts", {}) or {}
    artifact_paths = {
        name: (entry.get("path") if isinstance(entry, dict) else entry)
        for name, entry in artifacts.items()
    }
    required_artifacts = [
        name
        for name, entry in artifacts.items()
        if isinstance(entry, dict) and entry.get("required") is True
    ]
    optional_artifacts = [
        name
        for name, entry in artifacts.items()
        if isinstance(entry, dict) and entry.get("required") is False
    ]

    producer_version = _get_counter_producer_version()

    return {
        "entrypoint_path": str(manifest_path.resolve()),
        "entrypoint_kind": "analysis_bundle_manifest",
        "bundle_root": str(outdir.resolve()),
        "contract": {
            "contract_name": ANALYSIS_BUNDLE_CONTRACT_NAME,
            "contract_version": ANALYSIS_BUNDLE_CONTRACT_VERSION,
            "bundle_kind": ANALYSIS_BUNDLE_KIND,
            "producer": ANALYSIS_BUNDLE_PRODUCER,
            "producer_version": producer_version,
        },
        "artifact_paths": artifact_paths,
        "required_artifacts": required_artifacts,
        "optional_artifacts": optional_artifacts,
        "manifest_preview": manifest,
    }


def build_analysis_matrix(
    matrix_specs: list[MatrixSpec],
    sample_metadata_path: Path,
    outdir: Path,
    matrix_id: str,
    run_id: str | None = None,
) -> tuple[MatrixSpec, ExecutionRunSpec]:
    if not matrix_specs:
        raise ValueError("matrix_specs must not be empty")

    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    counts_dir = outdir / "counts"
    logs_dir = outdir / "logs"
    specs_dir = outdir / "specs"
    metadata_dir = outdir / "metadata"
    results_dir = outdir / "results"
    counts_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    specs_dir.mkdir(exist_ok=True)
    metadata_dir.mkdir(exist_ok=True)
    results_dir.mkdir(exist_ok=True)

    manifest_path = results_dir / "analysis_bundle_manifest.json"

    started_at = datetime.now(timezone.utc).astimezone().isoformat()

    _validate_mergeable_matrix_specs(matrix_specs)

    normalized_tables: list[pd.DataFrame] = []
    source_assay_ids: list[str] = []
    source_specimen_ids: list[str] = []
    source_subject_ids: list[str] = []
    specimen_ids_in_order: list[str] = []

    for spec in matrix_specs:
        df, specimen_id = _load_matrix_for_merge(spec)
        normalized_tables.append(df)
        specimen_ids_in_order.append(specimen_id)
        source_assay_ids.extend(spec.source_assay_ids)
        source_specimen_ids.extend(spec.source_specimen_ids)
        source_subject_ids.extend(spec.source_subject_ids)

    merge_provenance = _build_merge_provenance(matrix_specs)
    sample_metadata_alignment = _validate_and_align_sample_metadata(
        sample_metadata_path=sample_metadata_path,
        specimen_ids_in_order=specimen_ids_in_order,
        metadata_dir=metadata_dir,
    )
    warnings = _collect_analysis_merge_warnings(
        merge_provenance=merge_provenance,
        sample_metadata_alignment=sample_metadata_alignment,
    )

    merged_df = pd.concat(normalized_tables, axis=1, join="outer").fillna(0)
    merged_df = merged_df.astype(int)

    matrix_path = counts_dir / "merged_gene_numreads.tsv"
    merged_df.to_csv(matrix_path, sep="\t")

    producer_version = _get_counter_producer_version()
    analysis_metadata = {
        "producer_app": "iwa_rnaseq_counter",
        "producer_version": producer_version,
        "merge_strategy": "column_bind_by_feature_id",
        "sample_metadata_path": str(sample_metadata_path.resolve()),
        "aligned_sample_metadata_path": sample_metadata_alignment["aligned_sample_metadata_path"],
        "sample_metadata_alignment_status": sample_metadata_alignment["status"],
        "sample_metadata_id_column": sample_metadata_alignment["id_column"],
        "sample_metadata_row_count_input": sample_metadata_alignment["row_count_input"],
        "sample_metadata_row_count_aligned": sample_metadata_alignment["row_count_aligned"],
        "sample_metadata_extra_ids": sample_metadata_alignment["extra_metadata_ids"],
        "sample_metadata_aligned_columns": sample_metadata_alignment["aligned_columns"],
        "sample_ids": source_specimen_ids,
        "source_matrix_count": merge_provenance["source_matrix_count"],
        "source_quantifiers": merge_provenance["source_quantifiers"],
        "source_quantifier_versions": merge_provenance["source_quantifier_versions"],
        "source_aggregation_input_kinds": merge_provenance["source_aggregation_input_kinds"],
        "feature_annotation_consensus_status": merge_provenance["feature_annotation_consensus_status"],
        "feature_annotation_consensus_path": merge_provenance["feature_annotation_consensus_path"],
        "feature_annotation_consensus_reason": merge_provenance["feature_annotation_consensus_reason"],
        "feature_annotation_is_usable": merge_provenance["feature_annotation_is_usable"],
        "source_feature_annotation_inspections": merge_provenance["source_feature_annotation_inspections"],
        "source_matrix_contexts": merge_provenance["source_matrix_contexts"],
        "analysis_bundle_manifest_path": str(manifest_path.resolve()),
        "analysis_bundle_entrypoint_kind": "analysis_bundle_manifest",
        "warnings": warnings,
    }

    analysis_spec = MatrixSpec(
        schema_name="MatrixSpec",
        schema_version="0.1.0",
        matrix_id=matrix_id,
        matrix_scope="analysis",
        matrix_kind="count_matrix",
        feature_type="gene",
        value_type="integer",
        normalization="raw",
        feature_id_system=matrix_specs[0].feature_id_system,
        sample_axis="specimen",
        matrix_path=str(matrix_path.resolve()),
        feature_annotation_path=(
            merge_provenance["feature_annotation_consensus_path"]
            if merge_provenance["feature_annotation_is_usable"]
            else None
        ),
        source_assay_ids=source_assay_ids,
        source_specimen_ids=source_specimen_ids,
        source_subject_ids=source_subject_ids,
        metadata=analysis_metadata,
        overlay={},
    )

    finished_at = datetime.now(timezone.utc).astimezone().isoformat()
    log_path = logs_dir / "build_analysis_matrix.log"

    run_spec = ExecutionRunSpec(
        schema_name="ExecutionRunSpec",
        schema_version="0.1.0",
        run_id=run_id or f"RUN_BUILD_ANALYSIS_MATRIX_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        app_name="iwa_rnaseq_counter",
        app_version=producer_version,
        started_at=started_at,
        input_refs=[spec.matrix_id for spec in matrix_specs],
        output_refs=[analysis_spec.matrix_id],
        parameters={
            "merge_strategy": "column_bind_by_feature_id",
            "sample_metadata_path": str(sample_metadata_path.resolve()),
            "aligned_sample_metadata_path": sample_metadata_alignment["aligned_sample_metadata_path"],
            "sample_metadata_alignment_status": sample_metadata_alignment["status"],
            "sample_metadata_id_column": sample_metadata_alignment["id_column"],
            "sample_metadata_row_count_input": sample_metadata_alignment["row_count_input"],
            "sample_metadata_row_count_aligned": sample_metadata_alignment["row_count_aligned"],
            "sample_metadata_extra_ids": sample_metadata_alignment["extra_metadata_ids"],
            "source_matrix_count": merge_provenance["source_matrix_count"],
            "source_quantifiers": merge_provenance["source_quantifiers"],
            "source_quantifier_versions": merge_provenance["source_quantifier_versions"],
            "source_aggregation_input_kinds": merge_provenance["source_aggregation_input_kinds"],
            "feature_annotation_consensus_status": merge_provenance["feature_annotation_consensus_status"],
            "feature_annotation_consensus_path": merge_provenance["feature_annotation_consensus_path"],
            "feature_annotation_consensus_reason": merge_provenance["feature_annotation_consensus_reason"],
            "feature_annotation_is_usable": merge_provenance["feature_annotation_is_usable"],
            "analysis_bundle_manifest_path": str(manifest_path.resolve()),
            "analysis_bundle_entrypoint_kind": "analysis_bundle_manifest",
            "warnings": warnings,
        },
        execution_backend="local",
        finished_at=finished_at,
        status="completed",
        log_path=str(log_path.resolve()),
    )

    analysis_bundle = _build_analysis_bundle_preview(
        outdir=outdir,
        matrix_id=matrix_id,
        run_id=run_spec.run_id,
        feature_annotation_path=merge_provenance["feature_annotation_consensus_path"],
    )

    summary = _build_analysis_merge_summary(
        matrix_id=matrix_id,
        matrix_path=matrix_path,
        sample_metadata_path=sample_metadata_path,
        merge_provenance=merge_provenance,
        sample_metadata_alignment=sample_metadata_alignment,
        warnings=warnings,
        run_id=run_spec.run_id,
        input_refs=run_spec.input_refs,
        output_refs=run_spec.output_refs,
        bundle_manifest_path=manifest_path,
        analysis_bundle=analysis_bundle,
    )

    summary_path = results_dir / "analysis_merge_summary.json"

    _write_analysis_merge_summary(summary_path, summary)
    _write_analysis_merge_log(
        log_path,
        matrix_id=matrix_id,
        source_matrix_count=merge_provenance["source_matrix_count"],
        source_quantifiers=merge_provenance["source_quantifiers"],
        feature_annotation_consensus_status=merge_provenance["feature_annotation_consensus_status"],
        sample_metadata_alignment_status=sample_metadata_alignment["status"],
        warnings=warnings,
        analysis_bundle_manifest_path=str(manifest_path.resolve()),
    )

    _write_analysis_bundle_manifest(
        manifest_path,
        analysis_bundle["manifest_preview"],
    )

    return analysis_spec, run_spec
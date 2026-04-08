import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from ..models.matrix import MatrixSpec
from ..models.execution_run import ExecutionRunSpec

logger = logging.getLogger(__name__)


def _canonical_source_context(spec: MatrixSpec) -> dict:
    """
    v0.9.0-1:
    source MatrixSpec から downstream が読みやすい canonical provenance を抽出する。
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


def _normalize_optional_path(value: str | None) -> str | None:
    if value is None:
        return None
    s = str(value).strip()
    return s if s else None


def _inspect_feature_annotation_file(path_value: str | None) -> dict[str, Any]:
    """
    v0.9.0-2:
    feature_annotation.tsv を analysis handoff 用 artifact として使ってよいかを軽く検査する。
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
    v0.9.0-1 最小版:
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


def _read_sample_metadata_table(sample_metadata_path: Path) -> pd.DataFrame:
    """
    v0.9.0-3:
    sample metadata を TSV / CSV のどちらでも読み込めるようにする。
    """
    if not sample_metadata_path.exists() or not sample_metadata_path.is_file():
        raise FileNotFoundError(f"sample metadata file not found: {sample_metadata_path}")

    suffix = sample_metadata_path.suffix.lower()

    # まず拡張子ベースで試す
    if suffix in {".tsv", ".txt"}:
        df = pd.read_csv(sample_metadata_path, sep="\t")
        if len(df.columns) == 1:
            # 誤って TSV で1列になったら CSV も試す
            df = pd.read_csv(sample_metadata_path)
        return df

    if suffix == ".csv":
        df = pd.read_csv(sample_metadata_path)
        if len(df.columns) == 1:
            # 誤って CSV で1列になったら TSV も試す
            df = pd.read_csv(sample_metadata_path, sep="\t")
        return df

    # suffix 不明なら CSV -> TSV の順で試す
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
    metadata_dir: Path,
) -> dict[str, Any]:
    """
    v0.9.0-3:
    merged matrix の specimen 列順に sample metadata を整列し、
    aligned artifact を保存する。

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

    metadata_dir.mkdir(parents=True, exist_ok=True)
    aligned_path = metadata_dir / "aligned_sample_metadata.tsv"
    aligned_df.to_csv(aligned_path, sep="\t", index=False)

    return {
        "status": "aligned",
        "id_column": id_col,
        "aligned_sample_metadata_path": str(aligned_path.resolve()),
        "row_count_input": int(len(df)),
        "row_count_aligned": int(len(aligned_df)),
        "missing_specimen_ids": missing_specimen_ids,
        "extra_metadata_ids": extra_metadata_ids,
        "duplicate_ids": duplicate_ids,
        "aligned_columns": [str(c) for c in aligned_df.columns],
    }


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
        # 複数 path だが列構造は一致
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


def build_analysis_matrix(
    matrix_specs: list[MatrixSpec],
    sample_metadata_path: Path,
    outdir: Path,
    matrix_id: str,
    run_id: str | None = None,
) -> tuple[MatrixSpec, ExecutionRunSpec]:
    if not matrix_specs:
        raise ValueError("matrix_specs must not be empty")

    _validate_mergeable_matrix_specs(matrix_specs)

    outdir.mkdir(parents=True, exist_ok=True)
    counts_dir = outdir / "counts"
    logs_dir = outdir / "logs"
    specs_dir = outdir / "specs"
    metadata_dir = outdir / "metadata"
    counts_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    specs_dir.mkdir(exist_ok=True)
    metadata_dir.mkdir(exist_ok=True)

    started_at = datetime.now(timezone.utc).astimezone().isoformat()

    merged_tables = []
    source_assay_ids: list[str] = []
    source_specimen_ids: list[str] = []
    source_subject_ids: list[str] = []

    for spec in matrix_specs:
        df = pd.read_csv(spec.matrix_path, sep="\t", index_col=0)

        specimen_id = spec.source_specimen_ids[0]

        # 単一 assay matrix を specimen_id 列名に正規化
        if df.shape[1] != 1:
            if specimen_id not in df.columns:
                raise ValueError(
                    f"MatrixSpec {spec.matrix_id} expected single-sample matrix or column named {specimen_id}, "
                    f"got columns={list(df.columns)}"
                )
            df = df[[specimen_id]]
        else:
            df.columns = [specimen_id]

        merged_tables.append(df)
        source_assay_ids.extend(spec.source_assay_ids)
        source_specimen_ids.extend(spec.source_specimen_ids)
        source_subject_ids.extend(spec.source_subject_ids)

    merged_df = pd.concat(merged_tables, axis=1, join="outer").fillna(0)
    merged_df = merged_df.astype(int)

    matrix_path = counts_dir / "merged_gene_numreads.tsv"
    merged_df.to_csv(matrix_path, sep="\t")

    merge_provenance = _build_merge_provenance(matrix_specs)

    sample_metadata_alignment = _validate_and_align_sample_metadata(
        sample_metadata_path=sample_metadata_path,
        specimen_ids_in_order=[str(c) for c in merged_df.columns.tolist()],
        metadata_dir=metadata_dir,
    )

    analysis_metadata = {
        "producer_app": "iwa_rnaseq_counter",
        "producer_version": "0.3.5",
        "merge_strategy": "column_bind_by_feature_id",
        "sample_metadata_path": str(sample_metadata_path),
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

    run_spec = ExecutionRunSpec(
        schema_name="ExecutionRunSpec",
        schema_version="0.1.0",
        run_id=run_id or f"RUN_BUILD_ANALYSIS_MATRIX_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        app_name="iwa_rnaseq_counter",
        app_version="0.3.5",
        started_at=started_at,
        input_refs=[spec.matrix_id for spec in matrix_specs],
        output_refs=[analysis_spec.matrix_id],
        parameters={
            "merge_strategy": "column_bind_by_feature_id",
            "sample_metadata_path": str(sample_metadata_path),
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
        },
        execution_backend="local",
        finished_at=finished_at,
        status="completed",
        log_path=str((logs_dir / "build_analysis_matrix.log").resolve()),
    )

    return analysis_spec, run_spec
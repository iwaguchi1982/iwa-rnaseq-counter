import logging
from datetime import datetime, timezone
from pathlib import Path

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
    feature_annotation_paths = _unique_preserve_order(
        [
            _normalize_optional_path(ctx.get("feature_annotation_path"))
            for ctx in source_contexts
            if _normalize_optional_path(ctx.get("feature_annotation_path"))
        ]
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

    if len(feature_annotation_paths) == 0:
        feature_annotation_consensus_status = "missing"
        feature_annotation_consensus_path = None
    elif len(feature_annotation_paths) == 1:
        feature_annotation_consensus_status = "consistent"
        feature_annotation_consensus_path = feature_annotation_paths[0]
    else:
        feature_annotation_consensus_status = "inconsistent"
        feature_annotation_consensus_path = feature_annotation_paths[0]

    return {
        "source_matrix_count": len(matrix_specs),
        "source_quantifiers": quantifiers,
        "source_quantifier_versions": quantifier_versions,
        "source_aggregation_input_kinds": aggregation_input_kinds,
        "source_feature_annotation_paths": feature_annotation_paths,
        "source_quantifier_index_paths": quantifier_index_paths,
        "source_tx2gene_paths": tx2gene_paths,
        "source_annotation_gtf_paths": annotation_gtf_paths,
        "feature_annotation_consensus_status": feature_annotation_consensus_status,
        "feature_annotation_consensus_path": feature_annotation_consensus_path,
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
    counts_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    specs_dir.mkdir(exist_ok=True)

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

    analysis_metadata = {
        "producer_app": "iwa_rnaseq_counter",
        "producer_version": "0.3.5",
        "merge_strategy": "column_bind_by_feature_id",
        "sample_metadata_path": str(sample_metadata_path),
        "sample_ids": source_specimen_ids,
        "source_matrix_count": merge_provenance["source_matrix_count"],
        "source_quantifiers": merge_provenance["source_quantifiers"],
        "source_quantifier_versions": merge_provenance["source_quantifier_versions"],
        "source_aggregation_input_kinds": merge_provenance["source_aggregation_input_kinds"],
        "feature_annotation_consensus_status": merge_provenance["feature_annotation_consensus_status"],
        "feature_annotation_consensus_path": merge_provenance["feature_annotation_consensus_path"],
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
        feature_annotation_path=merge_provenance["feature_annotation_consensus_path"],
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
            "source_matrix_count": merge_provenance["source_matrix_count"],
            "source_quantifiers": merge_provenance["source_quantifiers"],
            "source_quantifier_versions": merge_provenance["source_quantifier_versions"],
            "source_aggregation_input_kinds": merge_provenance["source_aggregation_input_kinds"],
            "feature_annotation_consensus_status": merge_provenance["feature_annotation_consensus_status"],
            "feature_annotation_consensus_path": merge_provenance["feature_annotation_consensus_path"],
        },
        execution_backend="local",
        finished_at=finished_at,
        status="completed",
        log_path=str((logs_dir / "build_analysis_matrix.log").resolve()),
    )

    return analysis_spec, run_spec
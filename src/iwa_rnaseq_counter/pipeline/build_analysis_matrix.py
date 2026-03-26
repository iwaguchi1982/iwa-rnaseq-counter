import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from ..models.matrix import MatrixSpec
from ..models.execution_run import ExecutionRunSpec

logger = logging.getLogger(__name__)


def build_analysis_matrix(
    matrix_specs: list[MatrixSpec],
    sample_metadata_path: Path,
    outdir: Path,
    matrix_id: str,
    run_id: str | None = None,
) -> tuple[MatrixSpec, ExecutionRunSpec]:
    if not matrix_specs:
        raise ValueError("matrix_specs must not be empty")

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

        if not spec.source_specimen_ids:
            raise ValueError(f"MatrixSpec {spec.matrix_id} has no source_specimen_ids")

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

    analysis_spec = MatrixSpec(
        schema_name="MatrixSpec",
        schema_version="0.1.0",
        matrix_id=matrix_id,
        matrix_scope="analysis",
        matrix_kind="count_matrix",
        feature_type="gene",
        value_type="integer",
        normalization="raw",
        feature_id_system="ensembl_gene_id",
        sample_axis="specimen",
        matrix_path=str(matrix_path.resolve()),
        feature_annotation_path=matrix_specs[0].feature_annotation_path,
        source_assay_ids=source_assay_ids,
        source_specimen_ids=source_specimen_ids,
        source_subject_ids=source_subject_ids,
        metadata={
            "producer_app": "iwa_rnaseq_counter",
            "producer_version": "0.3.0",
            "merge_strategy": "column_bind_by_feature_id",
            "sample_metadata_path": str(sample_metadata_path),
            "sample_ids": source_specimen_ids,
        },
        overlay={},
    )

    finished_at = datetime.now(timezone.utc).astimezone().isoformat()

    run_spec = ExecutionRunSpec(
        schema_name="ExecutionRunSpec",
        schema_version="0.1.0",
        run_id=run_id or f"RUN_BUILD_ANALYSIS_MATRIX_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        app_name="iwa_rnaseq_counter",
        app_version="0.3.0",
        started_at=started_at,
        input_refs=[spec.matrix_id for spec in matrix_specs],
        output_refs=[analysis_spec.matrix_id],
        parameters={
            "merge_strategy": "column_bind_by_feature_id",
            "sample_metadata_path": str(sample_metadata_path),
        },
        execution_backend="local",
        finished_at=finished_at,
        status="completed",
        log_path=str((logs_dir / "build_analysis_matrix.log").resolve()),
    )

    return analysis_spec, run_spec

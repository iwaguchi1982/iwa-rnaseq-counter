import logging
import pandas as pd
from datetime import datetime, timezone
from pathlib import Path

from iwa_rnaseq_counter.legacy.salmon_runner import run_salmon_quant
from iwa_rnaseq_counter.legacy.gene_aggregator import load_tx2gene_map, build_transcript_quant_table, aggregate_transcript_to_gene
from ..models.assay import AssaySpec
from ..models.matrix import MatrixSpec
from ..models.execution_run import ExecutionRunSpec

logger = logging.getLogger(__name__)

def run_counter_pipeline(
    assay_spec: AssaySpec,
    outdir: Path,
    threads: int = 4,
    run_id: str | None = None,
    profile: str = "local",
    quantifier: str = "salmon",
) -> tuple[MatrixSpec, ExecutionRunSpec]:
    started_at = datetime.now(timezone.utc).astimezone().isoformat()

    outdir.mkdir(parents=True, exist_ok=True)
    counts_dir = outdir / "counts"
    logs_dir = outdir / "logs"
    specs_dir = outdir / "specs"
    counts_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    specs_dir.mkdir(exist_ok=True)

    fastq_r1 = next((f.path for f in assay_spec.input_files if f.file_role.lower() == "fastq_r1"), None)
    fastq_r2 = next((f.path for f in assay_spec.input_files if f.file_role.lower() == "fastq_r2"), None)

    layout_final = "paired-end" if fastq_r2 else "single-end"
    r1_paths = [fastq_r1] if fastq_r1 else []
    r2_paths = [fastq_r2] if fastq_r2 else []
    all_paths = r1_paths if layout_final == "single-end" else []

    sample_df = pd.DataFrame([{
        "sample_id": assay_spec.specimen_id,
        "layout_final": layout_final,
        "r1_paths": r1_paths,
        "r2_paths": r2_paths,
        "all_paths": all_paths,
    }])

    salmon_index = None
    tx2gene = None
    if assay_spec.reference_resources:
        salmon_index = assay_spec.reference_resources.quantifier_index
        tx2gene = assay_spec.reference_resources.tx2gene_path

    if not salmon_index:
        raise ValueError("salmon_index is required in AssaySpec.reference_resources")

    if quantifier != "salmon":
        raise NotImplementedError(f"Only quantifier='salmon' is supported now, got: {quantifier!r}")

    run_result = run_salmon_quant(
        sample_df=sample_df,
        salmon_index_path=salmon_index,
        run_output_dir=str(outdir),
        strandedness_mode=assay_spec.strandedness or "Auto-detect",
        threads=threads,
    )

    outputs = run_result["outputs"]
    if not outputs or not outputs[0].get("is_success"):
        raise RuntimeError("Salmon quantification failed.")

    if not tx2gene:
        raise NotImplementedError("Transcript-level un-aggregated output is not currently handled without tx2gene.")

    tx2gene_df = load_tx2gene_map(tx2gene)
    t_nr_df = build_transcript_quant_table(outputs, value_type="NumReads")
    g_nr_df = aggregate_transcript_to_gene(t_nr_df, tx2gene_df)

    matrix_path = counts_dir / "gene_numreads.tsv"
    g_nr_df.to_csv(matrix_path, sep="\t")

    subject_id = assay_spec.metadata.get("subject_id")
    source_subject_ids = [subject_id] if subject_id else []

    matrix_spec = MatrixSpec(
        schema_name="MatrixSpec",
        schema_version="0.1.0",
        matrix_id=f"MAT_{assay_spec.assay_id}",
        matrix_scope="assay",
        matrix_kind="count_matrix",
        feature_type="gene",
        value_type="integer",
        normalization="raw",
        feature_id_system="ensembl_gene_id",
        sample_axis="specimen",
        matrix_path=str(matrix_path.resolve()),
        feature_annotation_path=str(tx2gene) if tx2gene else "",
        source_assay_ids=[assay_spec.assay_id],
        source_specimen_ids=[assay_spec.specimen_id],
        source_subject_ids=source_subject_ids,
        metadata={
            "producer_app": "iwa_rnaseq_counter",
            "producer_version": "0.3.5",
            "quantifier": quantifier,
            "salmon_index": str(salmon_index),
            "tx2gene_path": str(tx2gene),
            "sample_ids": [assay_spec.specimen_id],
            "feature_id_system_inferred": False,
            "feature_annotation_available": True,
        },
        overlay={},
    )

    finished_at = datetime.now(timezone.utc).astimezone().isoformat()

    run_spec = ExecutionRunSpec(
        schema_name="ExecutionRunSpec",
        schema_version="0.1.0",
        run_id=run_id or f"RUN_COUNTER_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        app_name="iwa_rnaseq_counter",
        app_version="0.3.5",
        input_refs=[assay_spec.assay_id],
        output_refs=[matrix_spec.matrix_id],
        parameters={
            "quantifier": quantifier,
            "threads": threads,
            "strandedness_mode": assay_spec.strandedness or "Auto-detect",
            "profile": profile,
            "salmon_index": str(salmon_index),
            "tx2gene_path": str(tx2gene),
        },
        execution_backend=profile,
        started_at=started_at,
        finished_at=finished_at,
        status="completed",
        log_path=str((logs_dir / "counter.log").resolve()),
    )
    
    return matrix_spec, run_spec

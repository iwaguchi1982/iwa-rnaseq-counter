import logging
import pandas as pd
from datetime import datetime, timezone
from pathlib import Path

from iwa_rnaseq_counter.pipeline.quantifiers.registry import get_quantifier
from iwa_rnaseq_counter.legacy.gene_aggregator import load_tx2gene_map, build_transcript_quant_table, aggregate_transcript_to_gene
from ..models.assay import AssaySpec
from ..models.matrix import MatrixSpec
from ..models.execution_run import ExecutionRunSpec

logger = logging.getLogger(__name__)

def _get_success_outputs(outputs: list[dict]) -> list[dict]:
    return [o for o in outputs if o.get("is_success")]


def _build_gene_numreads_matrix(run_result: dict, tx2gene_path: str | None) -> pd.DataFrame:
    aggregation_input_kind = run_result.get("aggregation_input_kind", "transcript_quant")
    success_outputs = _get_success_outputs(run_result.get("outputs", []))

    if not success_outputs:
        raise RuntimeError("No successful quantifier outputs found.")

    if aggregation_input_kind == "transcript_quant":
        if not tx2gene_path:
            raise NotImplementedError(
                "Transcript-level output requires tx2gene_path for gene-level aggregation."
            )
        tx2gene_df = load_tx2gene_map(tx2gene_path)
        t_nr_df = build_transcript_quant_table(success_outputs, value_type="NumReads")
        return aggregate_transcript_to_gene(t_nr_df, tx2gene_df)

    if aggregation_input_kind == "gene_counts":
        # Future: implement direct gene-level loading if backend provides it
        raise NotImplementedError("Direct gene_counts aggregation input is not yet implemented in runner.py")

    raise ValueError(f"Unknown aggregation_input_kind: {aggregation_input_kind!r}")


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

    quantifier_index = None
    tx2gene = None
    if assay_spec.reference_resources:
        quantifier_index = assay_spec.reference_resources.quantifier_index
        tx2gene = assay_spec.reference_resources.tx2gene_path

    if not quantifier_index:
        raise ValueError("quantifier_index is required in AssaySpec.reference_resources")

    # Quantifier execution via Registry (v0.7.0)
    quant = get_quantifier(quantifier)
    run_result = quant.run_quant(
        sample_df=sample_df,
        run_output_dir=outdir,
        threads=threads,
        strandedness_mode=assay_spec.strandedness or "Auto-detect",
        reference_config={
            "quantifier_index": quantifier_index,
            "tx2gene_path": tx2gene,
        },
    )

    # Aggregation entry point (v0.7.0)
    g_nr_df = _build_gene_numreads_matrix(run_result, tx2gene)

    matrix_path = counts_dir / "gene_numreads.tsv"
    g_nr_df.to_csv(matrix_path, sep="\t")
    
    # v0.5.1: Prepare standardized feature_annotation.tsv
    from iwa_rnaseq_counter.legacy.annotation_helper import prepare_feature_annotation, get_standard_annotation_path
    annotation_out = get_standard_annotation_path(outdir)
    # Ensure results dir exists
    (outdir / "results").mkdir(exist_ok=True)
    has_annotation = prepare_feature_annotation(tx2gene, annotation_out)

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
        # Use specimen as axis
        sample_axis="specimen",
        matrix_path=str(matrix_path.resolve()),
        feature_annotation_path=str(annotation_out.resolve()) if has_annotation else None,
        source_assay_ids=[assay_spec.assay_id],
        source_specimen_ids=[assay_spec.specimen_id],
        source_subject_ids=source_subject_ids,
        metadata={
            "producer_app": "iwa_rnaseq_counter",
            "producer_version": "0.3.5",
            "quantifier": run_result["quantifier"],
            "quantifier_version": run_result.get("quantifier_version"),
            "aggregation_input_kind": run_result.get(
                "aggregation_input_kind",
                "transcript_quant",
            ),
            "reference_context": run_result.get("reference_context", {}),
            "tx2gene_path": str(tx2gene),
            "sample_ids": [assay_spec.specimen_id],
            "feature_id_system_inferred": False,
            "feature_annotation_available": has_annotation,
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
            "quantifier": run_result["quantifier"],
            "quantifier_version": run_result.get("quantifier_version"),
            "aggregation_input_kind": run_result.get(
                "aggregation_input_kind",
                "transcript_quant",
            ),
            "threads": threads,
            "strandedness_mode": assay_spec.strandedness or "Auto-detect",
            "profile": profile,
            "reference_context": run_result.get("reference_context", {}),
            "tx2gene_path": str(tx2gene),
        },
        execution_backend=profile,
        started_at=started_at,
        finished_at=finished_at,
        status="completed",
        log_path=str((logs_dir / "run.log").resolve()),
    )
    
    return matrix_spec, run_spec

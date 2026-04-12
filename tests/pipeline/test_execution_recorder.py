import pytest
from pathlib import Path
from iwa_rnaseq_counter.pipeline.runner import run_counter_pipeline
from iwa_rnaseq_counter.models.assay import AssaySpec, InputFile, ReferenceResources
import json

def test_runner_records_failure_on_missing_index(tmp_path):
    # Setup outdir
    outdir = tmp_path / "run_fail"
    outdir.mkdir()
    
    # Create invalid AssaySpec (missing quantifier index)
    assay_spec = AssaySpec(
        schema_name="AssaySpec",
        schema_version="0.1.0",
        assay_id="TEST_ASSAY",
        specimen_id="TEST_SPECIMEN",
        assay_type="RNA-Seq",
        input_files=[
            InputFile(path="r1.fastq.gz", file_role="fastq_r1")
        ],
        reference_resources=ReferenceResources(
            quantifier_index=None # Missing!
        ),
        metadata={}
    )
    
    # Run and expect ValueError (re-thrown after recording)
    with pytest.raises(ValueError, match="requires quantifier_index"):
        run_counter_pipeline(
            assay_spec=assay_spec,
            outdir=outdir,
            quantifier="salmon"
        )
    
    # Verify ExecutionRunSpec was written
    spec_path = outdir / "specs" / "execution-run.spec.json"
    assert spec_path.exists()
    
    with open(spec_path, "r") as f:
        data = json.load(f)
        
    assert data["status"] == "failed"
    assert data["metadata"]["failure_stage"] == "quantifier"
    assert "requires quantifier_index" in data["metadata"]["failure_summary"]

def test_runner_records_success_structure(tmp_path, mocker):
    # Mock quantifier to avoid real execution
    mock_quant = mocker.patch("iwa_rnaseq_counter.pipeline.runner.get_quantifier")
    mock_inst = mock_quant.return_value
    mock_inst.name = "mock_salmon"
    
    from iwa_rnaseq_counter.pipeline.quantifiers.base import BackendCapabilities, BackendReferenceRequirements
    mock_inst.get_capabilities.return_value = BackendCapabilities(
        aggregation_input_kind="transcript_quant",
        has_transcript_quant=True,
        has_gene_counts=False,
        has_mapping_metrics=True,
        reference_requirements=BackendReferenceRequirements(
            quantifier_index="required",
            tx2gene="required",
            annotation_gtf="unused"
        )
    )
    
    # Mock run_quant
    mock_inst.run_quant.return_value = {
        "quantifier": "mock_salmon",
        "outputs": [{"sample_id": "TEST_SPECIMEN", "is_success": True, "gene_counts_path": "fake.tsv"}]
    }
    
    # Mock aggregation and annotation
    mocker.patch("iwa_rnaseq_counter.pipeline.runner._build_gene_numreads_matrix", return_value=mocker.MagicMock())
    mocker.patch("iwa_rnaseq_counter.legacy.annotation_helper.prepare_feature_annotation", return_value=True)
    mocker.patch("pandas.DataFrame.to_csv")
    
    outdir = tmp_path / "run_success"
    outdir.mkdir()
    
    assay_spec = AssaySpec(
        schema_name="AssaySpec",
        schema_version="0.1.0",
        assay_id="TEST_ASSAY",
        specimen_id="TEST_SPECIMEN",
        assay_type="RNA-Seq",
        input_files=[InputFile(path="r1.fastq.gz", file_role="fastq_r1")],
        reference_resources=ReferenceResources(quantifier_index="idx", tx2gene_path="t2g"),
        metadata={}
    )
    
    matrix_spec, run_spec = run_counter_pipeline(
        assay_spec=assay_spec,
        outdir=outdir,
        quantifier="salmon"
    )
    
    assert run_spec.status == "completed"
    assert run_spec.metadata["status_detail"] == "completed successfully"

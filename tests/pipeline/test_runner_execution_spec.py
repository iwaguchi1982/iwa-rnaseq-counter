import pytest
from pathlib import Path
from tempfile import TemporaryDirectory
import json

from iwa_rnaseq_counter.models.assay import AssaySpec, InputFile, ReferenceResources
from iwa_rnaseq_counter.pipeline.runner import run_counter_pipeline
from iwa_rnaseq_counter.models.execution_run import ExecutionRunSpec

def _fixture_assay_spec(fastq_r1: Path) -> AssaySpec:
    return AssaySpec(
        schema_name="AssaySpec",
        schema_version="0.1.0",
        assay_id="ASSAY_TEST_001",
        specimen_id="SPEC_TEST_001",
        assay_type="bulk_rnaseq",
        library_layout="single-end",
        input_files=[
            InputFile(file_role="fastq_r1", path=str(fastq_r1)),
        ],
        reference_resources=ReferenceResources(
            quantifier_index="/dummy/index",
            tx2gene_path="/dummy/tx2gene.tsv",
            annotation_gtf_path=None
        )
    )

def test_execution_run_spec_includes_preprocessing_steps(mocker):
    # Mock quantifier to avoid running salmon
    mock_quantifier = mocker.MagicMock()
    mock_quantifier.run_quant.return_value = {
        "quantifier": "mock_salmon",
        "quantifier_version": "1.0",
        "aggregation_input_kind": "transcript_quant",
        "reference_context": {},
        "outputs": [
            {
                "is_success": True,
                "sample_id": "SPEC_TEST_001",
                "quant_dir": "mock_path",
                "log_dir": "mock_path",
                "lib_type": "A"
            }
        ]
    }
    mocker.patch("iwa_rnaseq_counter.pipeline.runner.get_quantifier", return_value=mock_quantifier)
    
    # Mock transcript to gene aggregation output to avoid actually loading files
    mocker.patch("iwa_rnaseq_counter.pipeline.runner.load_tx2gene_map")
    mocker.patch("iwa_rnaseq_counter.pipeline.runner.build_transcript_quant_table")
    
    mock_agg = mocker.patch("iwa_rnaseq_counter.pipeline.runner.aggregate_transcript_to_gene")
    import pandas as pd
    mock_agg.return_value = pd.DataFrame({"SPEC_TEST_001": [1, 2]}, index=["gene1", "gene2"])
    
    mocker.patch("iwa_rnaseq_counter.legacy.annotation_helper.prepare_feature_annotation", return_value=True)
    mocker.patch("iwa_rnaseq_counter.legacy.annotation_helper.get_standard_annotation_path", return_value=Path("/dummy/annotation.tsv"))
    
    with TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        fastq = tmp / "test.fastq.gz"
        fastq.touch()
        
        assay = _fixture_assay_spec(fastq)
        
        matrix_spec, exec_spec = run_counter_pipeline(
            assay_spec=assay,
            outdir=tmp,
            threads=2,
            run_id="TEST_RUN_ID",
            quantifier="salmon"
        )
        
        assert isinstance(exec_spec, ExecutionRunSpec)
        
        # Verify preprocessing steps are present and correct shape
        assert exec_spec.preprocessing_steps is not None
        assert "qc" in exec_spec.preprocessing_steps
        assert "trimming" in exec_spec.preprocessing_steps
        
        step_qc = exec_spec.preprocessing_steps["qc"]
        assert step_qc.enabled is False
        assert step_qc.status == "not_run"
        
        # Verify matrix spec does NOT have qc info
        assert "preprocessing_steps" not in matrix_spec.metadata
        
        # Verify deep dumping
        dumped = exec_spec.to_dict()
        assert "preprocessing_steps" in dumped
        assert dumped["preprocessing_steps"]["qc"]["enabled"] is False
        assert dumped["preprocessing_steps"]["qc"]["status"] == "not_run"

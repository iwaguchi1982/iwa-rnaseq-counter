import pytest
from pathlib import Path
from iwa_rnaseq_counter.io.validate_counter_output import validate_counter_output

def test_validate_success_fixture(tmp_path):
    # Use existing fixture path if possible, or copy to tmp
    fixture_base = Path(__file__).parent.parent / "fixtures" / "contracts" / "counter" / "success_minimal"
    
    result = validate_counter_output(fixture_base)
    assert result.is_valid is True
    assert result.error_count == 0

def test_validate_failed_fixture(tmp_path):
    fixture_base = Path(__file__).parent.parent / "fixtures" / "contracts" / "counter" / "failed_minimal"
    
    result = validate_counter_output(fixture_base)
    # Failed run is "valid" if the record itself is complete
    assert result.is_valid is True
    assert result.error_count == 0

def test_validate_detects_missing_matrix(tmp_path):
    # Setup a "broken" success fixture in tmp
    outdir = tmp_path / "broken_success"
    outdir.mkdir()
    (outdir / "specs").mkdir()
    (outdir / "results").mkdir()
    (outdir / "logs").mkdir()
    
    # Write execution spec with status=completed but NO matrix file
    exec_spec = {
        "status": "completed",
        "schema_name": "ExecutionRunSpec",
        "schema_version": "0.1.0",
        "run_id": "TEST",
        "app_name": "app",
        "app_version": "1.0",
        "started_at": "now"
    }
    import json
    with open(outdir / "specs" / "execution-run.spec.json", "w") as f:
        json.dump(exec_spec, f)
    
    result = validate_counter_output(outdir)
    assert result.is_valid is False
    # Should report missing MatrixSpec and missing matrix file
    assert any(i.code == "missing_matrix_spec" for i in result.issues)
    assert any(i.code == "missing_matrix_file" for i in result.issues)

def test_validate_completed_with_errors_is_success_equivalent(tmp_path):
    # Setup a "completed_with_errors" fixture in tmp
    outdir = tmp_path / "completed_with_errors"
    outdir.mkdir()
    (outdir / "specs").mkdir()
    (outdir / "results").mkdir()
    (outdir / "logs").mkdir()
    
    import json
    # Write execution spec with status=completed_with_errors
    exec_spec = {
        "status": "completed_with_errors",
        "schema_name": "ExecutionRunSpec",
        "schema_version": "0.1.0",
        "run_id": "TEST_ERRORS",
        "app_name": "app",
        "app_version": "1.0",
        "started_at": "now",
        "output_refs": {"matrix": "MAT_001"}
    }
    with open(outdir / "specs" / "execution-run.spec.json", "w") as f:
        json.dump(exec_spec, f)
        
    # Write matrix spec
    matrix_spec = {
        "$schema_name": "MatrixSpec",
        "schema_version": "0.1.0",
        "matrix_id": "MAT_001",
        "matrix_kind": "count_matrix",
        "matrix_path": "results/gene_numreads.tsv",
        "metadata": {}
    }
    with open(outdir / "specs" / "matrix.spec.json", "w") as f:
        json.dump(matrix_spec, f)
        
    # Write dummy matrix file
    (outdir / "results" / "gene_numreads.tsv").write_text("dummy")
    # Write dummy log file
    (outdir / "logs" / "run.log").write_text("dummy")
    
    result = validate_counter_output(outdir)
    # Should be valid because it has the required matrix outcome
    assert result.is_valid is True
    assert result.error_count == 0

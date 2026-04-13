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

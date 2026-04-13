import pytest
from pathlib import Path
import json

def test_success_minimal_contract_check():
    """
    Ensure the success_minimal fixture satisfies the basic contract.
    If this fails, it means we either broke the validator or changed the contract intentionally.
    """
    fixture_dir = Path(__file__).parent.parent / "fixtures" / "contracts" / "counter" / "success_minimal"
    
    # 1. Check ExecutionRunSpec semantics
    with open(fixture_dir / "specs" / "execution-run.spec.json", "r") as f:
        data = json.load(f)
        assert data["status"] == "completed"
        assert len(data["output_refs"]) > 0
        assert data["metadata"]["output_generated"] is True
    
    # 2. Check MatrixSpec semantics
    with open(fixture_dir / "specs" / "matrix.spec.json", "r") as f:
        data = json.load(f)
        assert data["matrix_id"] == "MAT_001"
        assert data["matrix_kind"] == "count_matrix"

def test_failed_minimal_contract_check():
    """
    Ensure the failed_minimal fixture satisfies the basic contract.
    """
    fixture_dir = Path(__file__).parent.parent / "fixtures" / "contracts" / "counter" / "failed_minimal"
    
    with open(fixture_dir / "specs" / "execution-run.spec.json", "r") as f:
        data = json.load(f)
        assert data["status"] == "failed"
        assert data["metadata"]["failure_stage"] == "quantifier"
        assert len(data.get("output_refs", [])) == 0

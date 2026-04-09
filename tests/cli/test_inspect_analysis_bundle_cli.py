import subprocess
import sys
import json
from pathlib import Path

from tests.helpers.fixture_bundle import materialize_analysis_bundle_fixture

def _repo_root() -> Path:
    return Path(__file__).parent.parent.parent

def _fixture_root() -> Path:
    return Path(__file__).parent.parent / "fixtures" / "analysis_bundle"

def test_inspect_analysis_bundle_cli_human_readable_smoke(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    cmd = [
        sys.executable,
        str(_repo_root() / "cli.py"),
        "inspect-analysis-bundle",
        "--manifest",
        str(bundle_dir),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0
    assert "Analysis Bundle Inspect" in result.stdout
    assert "TEST_ANALYSIS_MATRIX" in result.stdout
    assert "TEST_RUN_001" in result.stdout
    assert "analysis_bundle" in result.stdout
    assert "rna_seq_analysis_bundle" in result.stdout

def test_inspect_analysis_bundle_cli_json_output(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    cmd = [
        sys.executable,
        str(_repo_root() / "cli.py"),
        "inspect-analysis-bundle",
        "--manifest",
        str(bundle_dir),
        "--json"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0
    
    # Verify valid JSON and contract keys
    data = json.loads(result.stdout)
    assert data["contract_name"] == "analysis_bundle"
    assert data["bundle_kind"] == "rna_seq_analysis_bundle"
    assert data["matrix_id"] == "TEST_ANALYSIS_MATRIX"
    assert data["run_id"] == "TEST_RUN_001"
    assert isinstance(data["matrix_shape"], dict)
    assert "analysis_bundle_manifest_path" in data

def test_inspect_analysis_bundle_cli_invalid_manifest_returns_nonzero():
    cmd = [
        sys.executable,
        str(_repo_root() / "cli.py"),
        "inspect-analysis-bundle",
        "--manifest",
        "non_existent_manifest.json",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 1
    assert "Error" in result.stderr

if __name__ == "__main__":
    from pathlib import Path
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        test_inspect_analysis_bundle_cli_human_readable_smoke(tmp_path)
        test_inspect_analysis_bundle_cli_json_output(tmp_path)
        test_inspect_analysis_bundle_cli_invalid_manifest_returns_nonzero()
    print("test_inspect_analysis_bundle_cli: ALL PASSED")

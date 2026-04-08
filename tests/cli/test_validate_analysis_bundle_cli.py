import sys
import subprocess
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from tests.helpers.fixture_bundle import materialize_analysis_bundle_fixture


def _fixture_root() -> Path:
    return Path(__file__).parent.parent / "fixtures" / "analysis_bundle"


def _repo_root() -> Path:
    return Path(__file__).parent.parent.parent


def test_cli_validate_analysis_bundle_valid(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    cmd = [
        sys.executable,
        str(_repo_root() / "cli.py"),
        "validate-analysis-bundle",
        "--manifest",
        str(bundle_dir),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0


def test_cli_validate_analysis_bundle_invalid(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "invalid_bundle")

    target = bundle_dir / "specs" / "matrix.spec.json"
    target.unlink()

    cmd = [
        sys.executable,
        str(_repo_root() / "cli.py"),
        "validate-analysis-bundle",
        "--manifest",
        str(bundle_dir),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 1

if __name__ == "__main__":
    from pathlib import Path
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        test_cli_validate_analysis_bundle_valid(tmp_path)
        test_cli_validate_analysis_bundle_invalid(tmp_path)
    print("test_validate_analysis_bundle_cli: ALL PASSED")

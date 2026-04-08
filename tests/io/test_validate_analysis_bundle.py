import sys
import json
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from iwa_rnaseq_counter.io.read_analysis_bundle import validate_analysis_bundle
from tests.helpers.fixture_bundle import materialize_analysis_bundle_fixture

def _fixture_root() -> Path:
    return Path(__file__).parent.parent / "fixtures" / "analysis_bundle"

def test_validate_analysis_bundle_returns_valid_for_valid_fixture(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    result = validate_analysis_bundle(bundle_dir)

    assert result.is_valid is True
    assert result.error_count == 0
    assert result.contract_info.is_supported is True

def test_validate_analysis_bundle_detects_missing_required_artifact(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "invalid_bundle")

    target = bundle_dir / "specs" / "matrix.spec.json"
    target.unlink()

    result = validate_analysis_bundle(bundle_dir)

    assert result.is_valid is False
    assert result.error_count >= 1
    assert any(
        issue.code in {"artifact_file_not_found", "missing_required_artifact"}
        and issue.artifact_name == "matrix_spec"
        for issue in result.issues
    )

def test_validate_analysis_bundle_detects_unsupported_contract(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "unsupported_contract_bundle")

    manifest_path = bundle_dir / "results" / "analysis_bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["contract_version"] = "2.0.0"  # Supported is 1.x
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    result = validate_analysis_bundle(bundle_dir)

    assert result.is_valid is False
    assert any(issue.code == "unsupported_contract_major" for issue in result.issues)

if __name__ == "__main__":
    from pathlib import Path
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        test_validate_analysis_bundle_returns_valid_for_valid_fixture(tmp_path)
        test_validate_analysis_bundle_detects_missing_required_artifact(tmp_path)
        test_validate_analysis_bundle_detects_unsupported_contract(tmp_path)
    print("test_validate_analysis_bundle: ALL PASSED")

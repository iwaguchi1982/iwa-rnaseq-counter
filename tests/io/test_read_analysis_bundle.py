import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from iwa_rnaseq_counter.io.read_analysis_bundle import read_analysis_bundle
from tests.helpers.fixture_bundle import materialize_analysis_bundle_fixture

def _fixture_root() -> Path:
    return Path(__file__).parent.parent / "fixtures" / "analysis_bundle"

def test_read_analysis_bundle_from_valid_fixture(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    bundle = read_analysis_bundle(bundle_dir)

    assert bundle.contract_info.is_supported is True
    assert bundle.matrix_spec.matrix_id == "TEST_ANALYSIS_MATRIX"
    assert bundle.execution_run_spec.run_id == "TEST_RUN_001"
    assert bundle.analysis_merge_summary["schema_name"] == "AnalysisMergeSummary"

    # Check aligned metadata
    assert list(bundle.aligned_sample_metadata["specimen_id"]) == ["SP001", "SP002"]
    
    # By default merged matrix is not loaded
    assert bundle.merged_matrix is None

def test_read_analysis_bundle_from_manifest_path(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    manifest_path = bundle_dir / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)

    assert bundle.contract_info.is_supported is True
    assert bundle.matrix_spec.matrix_id == "TEST_ANALYSIS_MATRIX"
    assert bundle.execution_run_spec.run_id == "TEST_RUN_001"

def test_read_analysis_bundle_loads_merged_matrix_when_requested(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    bundle = read_analysis_bundle(bundle_dir, load_merged_matrix=True)

    assert bundle.merged_matrix is not None
    assert list(bundle.merged_matrix.columns) == ["SP001", "SP002"]
    assert bundle.merged_matrix.shape == (3, 2)

if __name__ == "__main__":
    from pathlib import Path
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        test_read_analysis_bundle_from_valid_fixture(tmp_path)
        test_read_analysis_bundle_from_manifest_path(tmp_path)
        test_read_analysis_bundle_loads_merged_matrix_when_requested(tmp_path)
    print("test_read_analysis_bundle: ALL PASSED")

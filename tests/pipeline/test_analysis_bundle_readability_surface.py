import sys
import json
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from tests.helpers.fixture_bundle import materialize_analysis_bundle_fixture

def _fixture_root() -> Path:
    return Path(__file__).parent.parent / "fixtures" / "analysis_bundle"

def test_analysis_merge_summary_contains_consumer_readability_fields(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    summary_path = bundle_dir / "results" / "analysis_merge_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    # v0.10.3 Enhanced Fields
    assert "matrix_shape" in summary
    assert "sample_axis" in summary
    assert "feature_id_system" in summary
    assert "column_order_specimen_ids" in summary
    assert "source_quantifier_summary" in summary
    assert "feature_annotation_status" in summary
    assert "warning_summary" in summary

def test_matrixspec_metadata_contains_consumer_handoff_fields(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    spec_path = bundle_dir / "specs" / "matrix.spec.json"
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    metadata = spec["metadata"]

    # v0.10.3 Enhanced Fields
    assert "column_order_specimen_ids" in metadata
    assert "sample_count" in metadata
    assert "feature_count" in metadata
    assert "sample_axis_kind" in metadata
    assert "analysis_bundle_manifest_path" in metadata
    assert "analysis_bundle_entrypoint_kind" in metadata
    assert "warning_summary" in metadata

def test_executionrunspec_parameters_contains_consumer_summary_fields(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    run_path = bundle_dir / "specs" / "execution-run.spec.json"
    run_spec = json.loads(run_path.read_text(encoding="utf-8"))
    params = run_spec["parameters"]

    # v0.10.3 Enhanced Fields
    assert "matrix_shape" in params
    assert "sample_count" in params
    assert "feature_count" in params
    assert "sample_axis" in params
    assert "warning_count" in params
    assert "warning_summary" in params
    assert "source_quantifier_summary" in params
    assert "analysis_bundle_manifest_path" in params

if __name__ == "__main__":
    from pathlib import Path
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        test_analysis_merge_summary_contains_consumer_readability_fields(tmp_path)
        test_matrixspec_metadata_contains_consumer_handoff_fields(tmp_path)
        test_executionrunspec_parameters_contains_consumer_summary_fields(tmp_path)
    print("test_analysis_bundle_readability_surface: ALL PASSED")

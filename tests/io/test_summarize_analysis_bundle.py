import sys
import json
from pathlib import Path
import tempfile

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from iwa_rnaseq_counter.io.read_analysis_bundle import (
    read_analysis_bundle,
    summarize_analysis_bundle_for_consumer,
)
from tests.helpers.fixture_bundle import materialize_analysis_bundle_fixture

def _fixture_root() -> Path:
    return Path(__file__).parent.parent / "fixtures" / "analysis_bundle"

def test_summarize_analysis_bundle_from_valid_fixture(tmp_path):
    src = _fixture_root() / "valid_minimal_bundle"
    bundle_dir = materialize_analysis_bundle_fixture(src, tmp_path / "valid_bundle")

    # Read the bundle
    bundle = read_analysis_bundle(bundle_dir)
    
    # Summarize
    summary = summarize_analysis_bundle_for_consumer(bundle)

    # Check basic fields
    assert summary["contract_name"] == "analysis_bundle"
    assert summary["matrix_id"] == "TEST_ANALYSIS_MATRIX"
    assert summary["run_id"] == "TEST_RUN_001"
    assert summary["matrix_shape"]["feature_count"] == 3
    assert summary["matrix_shape"]["sample_count"] == 2
    assert summary["sample_axis"] == "specimen"
    assert summary["feature_id_system"] == "ensembl_gene_id"
    assert summary["column_order_specimen_ids"] == ["SP001", "SP002"]
    
    # Check normalization
    assert "analysis_bundle_manifest_path" in summary
    assert str(summary["analysis_bundle_manifest_path"]).endswith("analysis_bundle_manifest.json")

def test_summarize_analysis_bundle_priority_fallback(tmp_path):
    # Dummy bundle as dict
    bundle_dict = {
        "manifest": {
            "contract_name": "analysis_bundle",
            "run_id": "MANIFEST_RUN"
        },
        "analysis_merge_summary": {
            "matrix_id": "SUMMARY_MATRIX",
            "matrix_shape": {"feature_count": 100, "sample_count": 10}
        },
        "matrix_spec": None, # Missing matrix spec
        "execution_run_spec": {
            "run_id": "EXEC_RUN",
            "parameters": {
                "warning_count": 5
            }
        }
    }
    
    summary = summarize_analysis_bundle_for_consumer(bundle_dict)
    
    # matrix_id should be from summary
    assert summary["matrix_id"] == "SUMMARY_MATRIX"
    
    # run_id should be from spec preferentially
    assert summary["run_id"] == "EXEC_RUN"
    
    # warning_summary should fall back to execution_parameters.warning_count
    assert summary["warning_summary"]["warning_count"] == 5
    assert summary["warning_summary"]["has_warnings"] is True

if __name__ == "__main__":
    from pathlib import Path
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        test_summarize_analysis_bundle_from_valid_fixture(tmp_path)
        test_summarize_analysis_bundle_priority_fallback(tmp_path)
    print("test_summarize_analysis_bundle: ALL PASSED")

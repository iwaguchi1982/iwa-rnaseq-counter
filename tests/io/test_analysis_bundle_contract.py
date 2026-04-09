import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from iwa_rnaseq_counter.io.read_analysis_bundle import (
    read_analysis_bundle,
    summarize_analysis_bundle_for_consumer,
    get_bundle_matrix_shape,
    get_bundle_warning_summary,
    get_bundle_source_quantifier_summary,
    get_bundle_feature_annotation_status,
    get_bundle_column_order_specimen_ids,
    get_bundle_manifest_path,
    get_bundle_matrix_id,
    get_bundle_run_id,
    get_bundle_sample_axis,
    get_bundle_feature_id_system,
    get_bundle_sample_metadata_alignment_status
)

def _fixture_root() -> Path:
    return Path(__file__).parent.parent / "fixtures" / "analysis_bundle"

def _get_valid_bundle():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    return read_analysis_bundle(manifest_path)

def test_summarize_analysis_bundle_contract_has_stable_top_level_keys():
    bundle = _get_valid_bundle()
    summary = summarize_analysis_bundle_for_consumer(bundle)
    
    expected_keys = {
        "contract_name", "contract_version", "bundle_kind", "producer", "producer_version",
        "matrix_id", "run_id", "matrix_shape", "sample_axis", "feature_id_system",
        "column_order_specimen_ids", "source_quantifier_summary", "feature_annotation_status",
        "sample_metadata_alignment_status", "warning_summary", "analysis_bundle_manifest_path"
    }
    
    actual_keys = set(summary.keys())
    # Ensure all expected keys are present
    missing = expected_keys - actual_keys
    assert not missing, f"Missing contract keys: {missing}"
    
    # Assert types
    assert isinstance(summary["matrix_shape"], dict)
    assert isinstance(summary["column_order_specimen_ids"], list)
    assert isinstance(summary["source_quantifier_summary"], dict)
    assert isinstance(summary["analysis_bundle_manifest_path"], str)

def test_summarize_analysis_bundle_warning_summary_fallback_from_warning_count():
    # Contract: Even if warning_summary is missing in raw JSON, it should be built from warning_count
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {}},
        "execution_run_spec": {"parameters": {"warning_count": 3}},
        "analysis_merge_summary": {}
    }
    summary = summarize_analysis_bundle_for_consumer(bundle_dict)
    ws = summary["warning_summary"]
    assert isinstance(ws, dict)
    assert ws["has_warnings"] is True
    assert ws["warning_count"] == 3
    assert isinstance(ws["messages"], list)

def test_accessors_contract_from_valid_fixture():
    bundle = _get_valid_bundle()
    
    assert isinstance(get_bundle_matrix_shape(bundle), dict)
    assert isinstance(get_bundle_column_order_specimen_ids(bundle), list)
    assert isinstance(get_bundle_source_quantifier_summary(bundle), dict)
    assert isinstance(get_bundle_manifest_path(bundle), str)
    
    # Specific values for valid_minimal_bundle
    assert get_bundle_matrix_id(bundle) == "TEST_ANALYSIS_MATRIX"
    assert get_bundle_run_id(bundle) == "TEST_RUN_001"

def test_accessors_contract_normalizes_empty_collections():
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {}},
        "execution_run_spec": {"parameters": {}},
        "analysis_merge_summary": {}
    }
    assert get_bundle_column_order_specimen_ids(bundle_dict) == []
    assert get_bundle_source_quantifier_summary(bundle_dict) == {}
    assert get_bundle_warning_summary(bundle_dict) is None

def test_summary_helper_and_accessors_are_contract_consistent():
    bundle = _get_valid_bundle()
    summary = summarize_analysis_bundle_for_consumer(bundle)
    
    # Contract: summary helper MUST be consistent with individual accessors
    assert summary["matrix_id"] == get_bundle_matrix_id(bundle)
    assert summary["run_id"] == get_bundle_run_id(bundle)
    assert summary["matrix_shape"] == get_bundle_matrix_shape(bundle)
    assert summary["sample_axis"] == get_bundle_sample_axis(bundle)
    assert summary["feature_id_system"] == get_bundle_feature_id_system(bundle)
    assert summary["column_order_specimen_ids"] == get_bundle_column_order_specimen_ids(bundle)
    assert summary["source_quantifier_summary"] == get_bundle_source_quantifier_summary(bundle)
    assert summary["feature_annotation_status"] == get_bundle_feature_annotation_status(bundle)
    assert summary["sample_metadata_alignment_status"] == get_bundle_sample_metadata_alignment_status(bundle)
    assert summary["warning_summary"] == get_bundle_warning_summary(bundle)
    assert summary["analysis_bundle_manifest_path"] == get_bundle_manifest_path(bundle)

if __name__ == "__main__":
    test_summarize_analysis_bundle_contract_has_stable_top_level_keys()
    test_summarize_analysis_bundle_warning_summary_fallback_from_warning_count()
    test_accessors_contract_from_valid_fixture()
    test_accessors_contract_normalizes_empty_collections()
    test_summary_helper_and_accessors_are_contract_consistent()
    print("test_analysis_bundle_contract: ALL PASSED")

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

def test_get_bundle_matrix_id():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)
    # Contract check: analysis_summary priority
    assert get_bundle_matrix_id(bundle) == "TEST_ANALYSIS_MATRIX"

def test_get_bundle_run_id():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)
    assert get_bundle_run_id(bundle) == "TEST_RUN_001"

def test_get_bundle_matrix_shape_contract():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)
    shape = get_bundle_matrix_shape(bundle)
    assert isinstance(shape, dict)
    assert shape["feature_count"] == 3
    assert shape["sample_count"] == 2

def test_get_bundle_matrix_shape_falls_back_to_metadata_and_parameters():
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {"feature_count": 100, "sample_count": 50}},
        "execution_run_spec": {"parameters": {}},
        "analysis_merge_summary": {}
    }
    shape = get_bundle_matrix_shape(bundle_dict)
    assert shape["feature_count"] == 100
    assert shape["sample_count"] == 50

def test_get_bundle_warning_summary_builds_from_warning_count():
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {}},
        "execution_run_spec": {"parameters": {"warning_count": 5}},
        "analysis_merge_summary": {}
    }
    summary = get_bundle_warning_summary(bundle_dict)
    assert summary["has_warnings"] is True
    assert summary["warning_count"] == 5
    assert isinstance(summary["messages"], list)

def test_get_bundle_source_quantifier_summary_contract():
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {}},
        "execution_run_spec": {"parameters": {}},
        "analysis_merge_summary": {}
    }
    summary = get_bundle_source_quantifier_summary(bundle_dict)
    assert isinstance(summary, dict)
    assert summary == {}

def test_get_bundle_feature_annotation_status_returns_mapping_or_none():
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {"feature_annotation_status": {"is_usable": True}}},
        "execution_run_spec": {"parameters": {}},
        "analysis_merge_summary": {}
    }
    status = get_bundle_feature_annotation_status(bundle_dict)
    assert status["is_usable"] is True
    
    bundle_dict["matrix_spec"]["metadata"] = {}
    assert get_bundle_feature_annotation_status(bundle_dict) is None


def test_get_bundle_column_order_specimen_ids_contract():
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {"column_order_specimen_ids": ("S1", "S2")}},
        "execution_run_spec": {"parameters": {}},
        "analysis_merge_summary": {}
    }
    ids = get_bundle_column_order_specimen_ids(bundle_dict)
    assert isinstance(ids, list)
    assert ids == ["S1", "S2"]

def test_get_bundle_manifest_path_contract():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)
    path = get_bundle_manifest_path(bundle)
    assert isinstance(path, str)
    assert path.endswith("analysis_bundle_manifest.json")

def test_summary_helper_and_accessors_consistency_check():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)
    
    summary = summarize_analysis_bundle_for_consumer(bundle)
    
    # Contract: Consistency check between helper and individual accessors
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
    test_get_bundle_matrix_id()
    test_get_bundle_run_id()
    test_get_bundle_matrix_shape_contract()
    test_get_bundle_matrix_shape_falls_back_to_metadata_and_parameters()
    test_get_bundle_warning_summary_builds_from_warning_count()
    test_get_bundle_source_quantifier_summary_contract()
    test_get_bundle_feature_annotation_status_returns_mapping_or_none()
    test_get_bundle_column_order_specimen_ids_contract()
    test_get_bundle_manifest_path_contract()
    test_summary_helper_and_accessors_consistency_check()
    print("test_analysis_bundle_accessors: ALL PASSED")

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from iwa_rnaseq_counter.io.read_analysis_bundle import (
    read_analysis_bundle,
    get_bundle_matrix_shape,
    get_bundle_warning_summary,
    get_bundle_source_quantifier_summary,
    get_bundle_feature_annotation_status,
    get_bundle_column_order_specimen_ids,
    get_bundle_manifest_path,
    get_bundle_matrix_id,
    get_bundle_run_id
)

def _fixture_root() -> Path:
    return Path(__file__).parent.parent / "fixtures" / "analysis_bundle"

def test_get_bundle_matrix_id():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)
    assert get_bundle_matrix_id(bundle) == "TEST_ANALYSIS_MATRIX"

def test_get_bundle_run_id():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)
    assert get_bundle_run_id(bundle) == "TEST_RUN_001"

def test_get_bundle_matrix_shape_prefers_analysis_summary():
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
    assert summary["messages"] == []

def test_get_bundle_source_quantifier_summary_normalizes_empty_dict():
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {}},
        "execution_run_spec": {"parameters": {}},
        "analysis_merge_summary": {}
    }
    assert get_bundle_source_quantifier_summary(bundle_dict) == {}

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

def test_get_bundle_column_order_specimen_ids_normalizes_to_list():
    bundle_dict = {
        "manifest": {},
        "matrix_spec": {"metadata": {"column_order_specimen_ids": ("S1", "S2")}},
        "execution_run_spec": {"parameters": {}},
        "analysis_merge_summary": {}
    }
    ids = get_bundle_column_order_specimen_ids(bundle_dict)
    assert isinstance(ids, list)
    assert ids == ["S1", "S2"]

def test_get_bundle_manifest_path_prefers_paths_manifest_path():
    manifest_path = _fixture_root() / "valid_minimal_bundle" / "results" / "analysis_bundle_manifest.json"
    bundle = read_analysis_bundle(manifest_path)
    path = get_bundle_manifest_path(bundle)
    assert isinstance(path, str)
    assert path.endswith("analysis_bundle_manifest.json")

if __name__ == "__main__":
    test_get_bundle_matrix_id()
    test_get_bundle_run_id()
    test_get_bundle_matrix_shape_prefers_analysis_summary()
    test_get_bundle_matrix_shape_falls_back_to_metadata_and_parameters()
    test_get_bundle_warning_summary_builds_from_warning_count()
    test_get_bundle_source_quantifier_summary_normalizes_empty_dict()
    test_get_bundle_feature_annotation_status_returns_mapping_or_none()
    test_get_bundle_column_order_specimen_ids_normalizes_to_list()
    test_get_bundle_manifest_path_prefers_paths_manifest_path()
    print("test_analysis_bundle_accessors: ALL PASSED")

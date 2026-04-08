import sys
import json
from pathlib import Path
import tempfile

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from iwa_rnaseq_counter.io.read_analysis_bundle import validate_analysis_bundle

def _create_mock_bundle(base_dir: Path):
    bundle_root = base_dir / "bundle"
    results_dir = bundle_root / "results"
    specs_dir = bundle_root / "specs"
    counts_dir = bundle_root / "counts"
    metadata_dir = bundle_root / "metadata"
    logs_dir = bundle_root / "logs"

    for d in [results_dir, specs_dir, counts_dir, metadata_dir, logs_dir]:
        d.mkdir(parents=True)

    # Dummy files
    (specs_dir / "matrix.spec.json").write_text(json.dumps({"schema_name": "MatrixSpec"}))
    (specs_dir / "execution-run.spec.json").write_text(json.dumps({"schema_name": "ExecutionRunSpec"}))
    (counts_dir / "merged.tsv").write_text("gene_id\ts1\nG1\t10\n")
    (metadata_dir / "metadata.tsv").write_text("sample_id\tcondition\ns1\tA\n")
    (results_dir / "summary.json").write_text(json.dumps({"schema_name": "AnalysisMergeSummary"}))
    (logs_dir / "run.log").write_text("hello log")

    manifest = {
        "schema_name": "AnalysisBundleManifest",
        "schema_version": "0.2.0",
        "contract_name": "analysis_bundle",
        "contract_version": "1.0.0",
        "bundle_kind": "rna_seq_analysis_bundle",
        "bundle_root": str(bundle_root),
        "artifacts": {
            "matrix_spec": {"path": "specs/matrix.spec.json", "required": True},
            "execution_run_spec": {"path": "specs/execution-run.spec.json", "required": True},
            "merged_matrix": {"path": "counts/merged.tsv", "required": True},
            "aligned_sample_metadata": {"path": "metadata/metadata.tsv", "required": True},
            "analysis_merge_summary": {"path": "results/summary.json", "required": True},
            "build_analysis_matrix_log": {"path": "logs/run.log", "required": True},
        }
    }
    manifest_path = results_dir / "analysis_bundle_manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    return manifest_path

def test_validate_healthy_bundle():
    with tempfile.TemporaryDirectory() as tmp:
        manifest_path = _create_mock_bundle(Path(tmp))
        result = validate_analysis_bundle(manifest_path)
        assert result.is_valid is True
        assert result.error_count == 0

def test_validate_missing_artifact():
    with tempfile.TemporaryDirectory() as tmp:
        manifest_path = _create_mock_bundle(Path(tmp))
        # Corrupt manifest
        data = json.loads(manifest_path.read_text())
        del data["artifacts"]["matrix_spec"]
        manifest_path.write_text(json.dumps(data))

        result = validate_analysis_bundle(manifest_path)
        assert result.is_valid is False
        assert result.error_count == 1
        assert any(i.code == "missing_required_artifact" for i in result.issues)

def test_validate_missing_file():
    with tempfile.TemporaryDirectory() as tmp:
        manifest_path = _create_mock_bundle(Path(tmp))
        # Delete file
        (manifest_path.parent.parent / "specs/matrix.spec.json").unlink()

        result = validate_analysis_bundle(manifest_path)
        assert result.is_valid is False
        assert result.error_count == 1
        assert any(i.code == "artifact_file_not_found" for i in result.issues)

def test_validate_unsupported_contract():
    with tempfile.TemporaryDirectory() as tmp:
        manifest_path = _create_mock_bundle(Path(tmp))
        # Change version
        data = json.loads(manifest_path.read_text())
        data["contract_version"] = "2.0.0"
        manifest_path.write_text(json.dumps(data))

        result = validate_analysis_bundle(manifest_path)
        assert result.is_valid is False
        assert result.error_count == 1
        assert any(i.code == "unsupported_contract_major" for i in result.issues)

if __name__ == "__main__":
    test_validate_healthy_bundle()
    test_validate_missing_artifact()
    test_validate_missing_file()
    test_validate_unsupported_contract()
    print("test_validate_analysis_bundle: ALL PASSED")

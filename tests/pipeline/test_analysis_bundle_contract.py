import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from iwa_rnaseq_counter.pipeline.build_analysis_matrix import (
    _build_analysis_bundle_manifest,
)

def test_build_analysis_bundle_manifest_includes_contract_fields(tmp_path):
    manifest = _build_analysis_bundle_manifest(
        outdir=tmp_path,
        matrix_id="TEST_MAT",
        run_id="TEST_RUN",
        matrix_spec_path=tmp_path / "specs/matrix.spec.json",
        execution_run_spec_path=tmp_path / "specs/execution-run.spec.json",
        merged_matrix_path=tmp_path / "counts/merged_gene_numreads.tsv",
        aligned_sample_metadata_path=None,
        analysis_merge_summary_path=tmp_path / "results/analysis_merge_summary.json",
        build_analysis_matrix_log_path=tmp_path / "logs/build_analysis_matrix.log",
        feature_annotation_path=None,
    )

    assert manifest["contract_name"] == "analysis_bundle"
    assert manifest["contract_version"] == "1.0.0"
    assert manifest["bundle_kind"] == "rna_seq_analysis_bundle"
    assert manifest["producer"] == "iwa_rnaseq_counter"
    assert "producer_version" in manifest

if __name__ == "__main__":
    from pathlib import Path
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        test_build_analysis_bundle_manifest_includes_contract_fields(Path(tmp_dir))
    print("test_build_analysis_bundle_manifest_includes_contract_fields: PASSED")

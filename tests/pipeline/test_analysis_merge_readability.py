import sys
from pathlib import Path
import json
import tempfile
import pandas as pd

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from iwa_rnaseq_counter.models.matrix import MatrixSpec
from iwa_rnaseq_counter.pipeline.build_analysis_matrix import build_analysis_matrix

def test_build_analysis_matrix_readability_fields():
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        
        # Create mock input specs
        spec1_path = tmp_dir / "s1/specs/matrix.spec.json"
        spec1_path.parent.mkdir(parents=True)
        m1_path = tmp_dir / "s1/counts/gene_numreads.tsv"
        m1_path.parent.mkdir(parents=True)
        m1_path.write_text("gene_id\ts1\nG1\t10\nG2\t20\n")
        
        # Feature Annotation
        anno_path = tmp_dir / "feature_annotation.tsv"
        anno_path.write_text("feature_id\tgene_symbol\nG1\tGene1\nG2\tGene2\n")

        spec1 = MatrixSpec(
            schema_name="MatrixSpec", schema_version="0.1.0",
            matrix_id="SPEC1", matrix_scope="assay", matrix_kind="count_matrix",
            feature_type="gene", value_type="integer", normalization="raw",
            feature_id_system="ensembl_gene_id", sample_axis="specimen",
            matrix_path=str(m1_path), source_assay_ids=["A1"],
            source_specimen_ids=["s1"], source_subject_ids=["SUBJ1"],
            feature_annotation_path=str(anno_path),
            metadata={"producer_app": "test", "quantifier": "salmon"}
        )
        
        spec2_path = tmp_dir / "s2/specs/matrix.spec.json"
        spec2_path.parent.mkdir(parents=True)
        m2_path = tmp_dir / "s2/counts/gene_numreads.tsv"
        m2_path.parent.mkdir(parents=True)
        m2_path.write_text("gene_id\ts2\nG1\t15\nG2\t25\n")
        
        spec2 = MatrixSpec(
            schema_name="MatrixSpec", schema_version="0.1.0",
            matrix_id="SPEC2", matrix_scope="assay", matrix_kind="count_matrix",
            feature_type="gene", value_type="integer", normalization="raw",
            feature_id_system="ensembl_gene_id", sample_axis="specimen",
            matrix_path=str(m2_path), source_assay_ids=["A2"],
            source_specimen_ids=["s2"], source_subject_ids=["SUBJ2"],
            feature_annotation_path=str(anno_path),
            metadata={"producer_app": "test", "quantifier": "salmon"}
        )
        
        # Metadata
        metadata_path = tmp_dir / "metadata.tsv"
        metadata_path.write_text("specimen_id\tsubject_id\tcondition\tgroup\tbatch\ns1\tS1\tA\tG1\tB1\ns2\tS2\tB\tG1\tB1\n")
        
        outdir = tmp_dir / "out"
        
        # Execution
        matrix_spec, exec_run_spec = build_analysis_matrix(
            matrix_specs=[spec1, spec2],
            sample_metadata_path=metadata_path,
            outdir=outdir,
            matrix_id="ANALYSIS_MAT",
            run_id="ANALYSIS_RUN"
        )
        
        # 1. Check MatrixSpec Metadata
        meta = matrix_spec.metadata
        assert meta["sample_count"] == 2
        assert meta["feature_count"] == 2
        assert meta["column_order_specimen_ids"] == ["s1", "s2"]
        assert meta["sample_axis_kind"] == "specimen"
        
        # 2. Check ExecutionRunSpec Parameters
        params = exec_run_spec.parameters
        assert params["warning_count"] == 0
        assert params["source_quantifier_summary"]["values"] == ["salmon"]
        
        # 3. Check AnalysisMergeSummary
        summary_path = outdir / "results" / "analysis_merge_summary.json"
        assert summary_path.exists()
        summary = json.loads(summary_path.read_text())
        
        assert summary["matrix_shape"] == {"feature_count": 2, "sample_count": 2}
        assert summary["sample_axis"] == "specimen"
        assert summary["feature_id_system"] == "ensembl_gene_id"
        assert summary["column_order_specimen_ids"] == ["s1", "s2"]
        assert summary["source_quantifier_summary"]["values"] == ["salmon"]
        assert summary["warning_summary"]["count"] == 0
        
        print("test_build_analysis_matrix_readability_fields: ALL PASSED")

if __name__ == "__main__":
    test_build_analysis_matrix_readability_fields()

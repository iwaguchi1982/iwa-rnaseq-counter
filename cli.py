import argparse
import json
import logging
import sys
from pathlib import Path

# Add local src and root src to sys.path
sys.path.insert(0, str(Path(__file__).parent.resolve() / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent.resolve() / "src"))

from iwa_rnaseq_counter.io.read_assay_spec import read_assay_spec
from iwa_rnaseq_counter.models.assay import AssaySpec, ReferenceResources
from iwa_rnaseq_counter.io.read_matrix_spec import read_matrix_spec
from iwa_rnaseq_counter.io.write_matrix_spec import write_matrix_spec
from iwa_rnaseq_counter.io.write_execution_run_spec import write_execution_run_spec
from iwa_rnaseq_counter.pipeline.runner import run_counter_pipeline
from iwa_rnaseq_counter.pipeline.build_analysis_matrix import (
    build_analysis_matrix,
    preview_build_analysis_matrix,
)
from iwa_rnaseq_counter.io.read_analysis_bundle import (
    validate_analysis_bundle,
    read_analysis_bundle,
    summarize_analysis_bundle_for_consumer,
)


def setup_logging(level: str, logfile: Path | None = None) -> None:
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {level}")

    handlers = [logging.StreamHandler(sys.stdout)]
    if logfile is not None:
        logfile.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(logfile, encoding="utf-8"))

    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
        force=True,
    )


def _normalize_quantifier_name(quantifier: str | None) -> str:
    return str(quantifier or "salmon").strip().lower()


def _validate_reference_requirements_for_quantifier(
    *,
    quantifier: str,
    quantifier_index: str | None = None,
    tx2gene_path: str | None = None,
    annotation_gtf_path: str | None = None,
    context: str = "run",
) -> None:
    """
    v0.8.4:
    assay / batch / GUI で揃えたい required reference ルールを
    CLI 側でも同じ形で使う。
    """
    q = _normalize_quantifier_name(quantifier)

    if not quantifier_index:
        raise ValueError(
            f"{context}: quantifier_index is required. "
            "Pass --quantifier-index (or legacy --salmon-index), or define it in the spec."
        )

    if q in {"salmon", "kallisto"} and not tx2gene_path:
        raise ValueError(
            f"{context}: tx2gene_path is required for quantifier '{q}'. "
            "Pass --tx2gene or define it in the spec/sample sheet context."
        )

    if q == "hisat2" and not annotation_gtf_path:
        raise ValueError(
            f"{context}: annotation_gtf_path is required for quantifier 'hisat2'. "
            "Pass --annotation-gtf or define it in the spec/sample sheet context."
        )


def _apply_run_assay_reference_overrides(
    assay_spec: AssaySpec,
    *,
    quantifier: str,
    quantifier_index: str | None = None,
    tx2gene_path: str | None = None,
    annotation_gtf_path: str | None = None,
) -> AssaySpec:
    """
    run-assay 用:
    AssaySpec.reference_resources に対して CLI override を適用し、
    その後 quantifier 別 required reference を検証する。
    """
    ref = assay_spec.reference_resources
    if ref is None:
        ref = ReferenceResources()

    if quantifier_index:
        ref.quantifier_index = str(quantifier_index)
    if tx2gene_path:
        ref.tx2gene_path = str(tx2gene_path)
    if annotation_gtf_path:
        ref.annotation_gtf_path = str(annotation_gtf_path)

    assay_spec.reference_resources = ref

    _validate_reference_requirements_for_quantifier(
        quantifier=quantifier,
        quantifier_index=ref.quantifier_index,
        tx2gene_path=ref.tx2gene_path,
        annotation_gtf_path=ref.annotation_gtf_path,
        context="run-assay",
    )

    return assay_spec


def main():
    parser = argparse.ArgumentParser(description="iwa_rnaseq_counter CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_run = subparsers.add_parser("run-assay", help="Run single AssaySpec and emit assay-level MatrixSpec")
    p_run.add_argument("--assay-spec", required=True, type=Path)
    p_run.add_argument("--outdir", required=True, type=Path)
    p_run.add_argument("--run-id", type=str)
    p_run.add_argument("--profile", type=str, default="local")
    p_run.add_argument(
        "--quantifier",
        type=str,
        default="salmon",
        help="Quantifier backend: salmon, star, hisat2, kallisto",
    )
    p_run.add_argument("--quantifier-index", type=Path, help="Override quantifier index path")
    p_run.add_argument("--salmon-index", type=Path, help="Legacy alias for --quantifier-index")
    p_run.add_argument("--tx2gene", type=Path, help="Override tx2gene path")
    p_run.add_argument("--annotation-gtf", type=Path, help="Override annotation GTF path (required for hisat2)")
    p_run.add_argument("--threads", type=int, default=4)
    p_run.add_argument("--dry-run", action="store_true")
    p_run.add_argument("--log-level", type=str, default="INFO")

    p_batch = subparsers.add_parser("run-batch", help="Run multiple AssaySpecs from a sample sheet")
    p_batch.add_argument("--sample-sheet", required=True, type=Path)
    p_batch.add_argument("--quantifier-index", type=Path, help="Generic quantifier index path (preferred)")
    p_batch.add_argument("--salmon-index", type=Path, help="Legacy alias for --quantifier-index")
    p_batch.add_argument("--tx2gene", type=Path, help="tx2gene path (required for salmon / kallisto)")
    p_batch.add_argument("--annotation-gtf", type=Path, help="Annotation GTF path (required for HISAT2)")
    p_batch.add_argument("--strandedness", type=str, default="Auto-detect")
    p_batch.add_argument("--outdir", required=True, type=Path)
    p_batch.add_argument("--run-id", type=str)
    p_batch.add_argument("--profile", type=str, default="local")
    p_batch.add_argument("--quantifier", type=str, default="salmon")
    p_batch.add_argument("--threads", type=int, default=4)
    p_batch.add_argument("--dry-run", action="store_true")
    p_batch.add_argument("--log-level", type=str, default="INFO")

    p_merge = subparsers.add_parser("build-analysis-matrix", help="Merge assay-level MatrixSpecs into reporter-ready analysis MatrixSpec")
    p_merge.add_argument("--matrix-spec", required=True, nargs="+", type=Path)
    p_merge.add_argument("--sample-metadata", required=True, type=Path)
    p_merge.add_argument("--outdir", required=True, type=Path)
    p_merge.add_argument("--matrix-id", required=True, type=str)
    p_merge.add_argument("--run-id", type=str)
    p_merge.add_argument("--dry-run", action="store_true")
    p_merge.add_argument("--log-level", type=str, default="INFO")

    p_gui = subparsers.add_parser("run-gui-backend", help="Execute the monolithic GUI pipeline via CLI for job running")
    p_gui.add_argument("--config", required=True, type=Path)
    p_gui.add_argument("--sample-sheet", required=True, type=Path)
    p_gui.add_argument("--outdir", required=True, type=Path)
    p_gui.add_argument("--started-at", type=str, required=True)
    p_gui.add_argument("--log-level", type=str, default="INFO")

    p_validate = subparsers.add_parser("validate-analysis-bundle", help="Validate an existing analysis bundle for integrity and contract compatibility")
    p_validate.add_argument("--manifest", required=True, type=Path, help="Path to manifest or bundle directory")
    p_validate.add_argument("--json", action="store_true", help="Output validation results as JSON")
    p_validate.add_argument("--log-level", type=str, default="INFO")

    p_inspect = subparsers.add_parser("inspect-analysis-bundle", help="Inspect an analysis bundle and print its summary as JSON")
    p_inspect.add_argument("--manifest", required=True, type=Path, help="Path to manifest or bundle directory")
    p_inspect.add_argument("--json", action="store_true", help="Output bundle summary as machine-readable JSON")
    p_inspect.add_argument("--log-level", type=str, default="WARNING")

    args = parser.parse_args()

    if args.command == "run-assay":
        setup_logging(args.log_level, args.outdir / "logs" / "counter.log")
        logger = logging.getLogger(__name__)

        assay_spec = read_assay_spec(args.assay_spec)

        quantifier_index_arg = args.quantifier_index or args.salmon_index

        try:
            assay_spec = _apply_run_assay_reference_overrides(
                assay_spec,
                quantifier=args.quantifier,
                quantifier_index=str(quantifier_index_arg) if quantifier_index_arg else None,
                tx2gene_path=str(args.tx2gene) if args.tx2gene else None,
                annotation_gtf_path=str(args.annotation_gtf) if args.annotation_gtf else None,
            )
        except Exception as e:
            logger.error(f"Invalid run-assay reference configuration: {e}")
            sys.exit(1)

        if args.dry_run:
            logger.info("Dry run complete. AssaySpec validated successfully.")
            return

        matrix_spec, exec_spec = run_counter_pipeline(
            assay_spec=assay_spec,
            outdir=args.outdir,
            threads=args.threads,
            run_id=args.run_id,
            profile=args.profile,
            quantifier=args.quantifier,
        )

        write_matrix_spec(matrix_spec, args.outdir / "specs" / "matrix.spec.json")
        write_execution_run_spec(exec_spec, args.outdir / "specs" / "execution-run.spec.json")
        logger.info("run-assay completed")

    elif args.command == "run-batch":
        setup_logging(args.log_level, args.outdir / "logs" / "batch_counter.log")
        logger = logging.getLogger(__name__)

        from iwa_rnaseq_counter.io.read_sample_sheet import read_sample_sheet

        quantifier_index_arg = args.quantifier_index or args.salmon_index

        try:
            _validate_reference_requirements_for_quantifier(
                quantifier=args.quantifier,
                quantifier_index=str(quantifier_index_arg) if quantifier_index_arg else None,
                tx2gene_path=str(args.tx2gene) if args.tx2gene else None,
                annotation_gtf_path=str(args.annotation_gtf) if args.annotation_gtf else None,
                context="run-batch",
            )
        except Exception as e:
            parser.error(str(e))

        try:
            assay_specs = read_sample_sheet(
                sample_sheet_path=args.sample_sheet,
                quantifier_index_path=str(quantifier_index_arg) if quantifier_index_arg else None,
                tx2gene_path=str(args.tx2gene) if args.tx2gene else None,
                strandedness=args.strandedness,
                annotation_gtf_path=str(args.annotation_gtf) if args.annotation_gtf else None,
            )
        except Exception as e:
            logger.error(f"Failed to read sample sheet: {e}")
            sys.exit(1)

        if args.dry_run:
            logger.info(f"Dry run complete. Validated {len(assay_specs)} AssaySpecs.")
            return

        out_matrix_specs = []
        for spec in assay_specs:
            logger.info(f"Processing assay: {spec.assay_id}")
            assay_outdir = args.outdir / spec.specimen_id
            try:
                matrix_spec, exec_spec = run_counter_pipeline(
                    assay_spec=spec,
                    outdir=assay_outdir,
                    threads=args.threads,
                    run_id=args.run_id,
                    profile=args.profile,
                    quantifier=args.quantifier,
                )
                write_matrix_spec(matrix_spec, assay_outdir / "specs" / "matrix.spec.json")
                write_execution_run_spec(exec_spec, assay_outdir / "specs" / "execution-run.spec.json")
                out_matrix_specs.append(matrix_spec)
            except Exception as e:
                logger.error(f"Pipeline failed for {spec.assay_id}: {e}")
                
        logger.info(f"run-batch completed. Processed {len(out_matrix_specs)} assays successfully.")

    elif args.command == "build-analysis-matrix":
        setup_logging(args.log_level, args.outdir / "logs" / "build_analysis_matrix.log")
        logger = logging.getLogger(__name__)

        matrix_specs = [read_matrix_spec(p) for p in args.matrix_spec]
        manifest_path = args.outdir / "results" / "analysis_bundle_manifest.json"

        if args.dry_run:
            preview = preview_build_analysis_matrix(
                matrix_specs=matrix_specs,
                sample_metadata_path=args.sample_metadata,
                outdir=args.outdir,
                matrix_id=args.matrix_id,
                run_id=args.run_id,
            )
            bundle_contract = preview.get("analysis_bundle", {}).get("contract", {})
            if bundle_contract:
                logger.info("analysis bundle contract preview:")
                logger.info(f"  contract_name: {bundle_contract.get('contract_name')}")
                logger.info(f"  contract_version: {bundle_contract.get('contract_version')}")
                logger.info(f"  bundle_kind: {bundle_contract.get('bundle_kind')}")
                logger.info(f"  producer: {bundle_contract.get('producer')}")
                logger.info(f"  producer_version: {bundle_contract.get('producer_version')}")

            logger.info(json.dumps(preview, indent=2, ensure_ascii=False))
            return

        matrix_spec, exec_spec = build_analysis_matrix(
            matrix_specs=matrix_specs,
            sample_metadata_path=args.sample_metadata,
            outdir=args.outdir,
            matrix_id=args.matrix_id,
            run_id=args.run_id,
        )

        write_matrix_spec(matrix_spec, args.outdir / "specs" / "matrix.spec.json")
        write_execution_run_spec(exec_spec, args.outdir / "specs" / "execution-run.spec.json")

        logger.info("build-analysis-matrix completed")
        logger.info("analysis bundle entrypoint: %s", str(manifest_path.resolve()))
        logger.info("analysis bundle contract: analysis_bundle@1.x (rna_seq_analysis_bundle)")
        logger.info("matrix spec: %s", str((args.outdir / "specs" / "matrix.spec.json").resolve()))
        logger.info(
            "execution run spec: %s",
            str((args.outdir / "specs" / "execution-run.spec.json").resolve()),
        )


    elif args.command == "run-gui-backend":
        setup_logging(args.log_level, args.outdir / "logs" / "run.log")
        logger = logging.getLogger(__name__)
        
        import pandas as pd
        from iwa_rnaseq_counter.pipeline.gui_backend import run_gui_backend_pipeline
        
        try:
            with open(args.config, "r") as f:
                config_data = json.load(f)
            import ast
            def parse_list(x):
                try:
                    if isinstance(x, str) and x.startswith("[") and x.endswith("]"):
                        return ast.literal_eval(x)
                    return x
                except:
                    return x

            sample_df = pd.read_csv(args.sample_sheet)
            for col in ["r1_paths", "r2_paths", "all_paths"]:
                if col in sample_df.columns:
                    sample_df[col] = sample_df[col].apply(parse_list)

            run_gui_backend_pipeline(args.outdir, config_data, sample_df, args.started_at)
        except Exception as e:
            logger.error(f"GUI Backend pipeline failed: {e}")
            sys.exit(1)

    elif args.command == "validate-analysis-bundle":
        setup_logging(args.log_level)
        logger = logging.getLogger(__name__)

        result = validate_analysis_bundle(args.manifest)

        if args.json:
            # We need a custom serializer for Path objects
            def default_serializer(obj):
                if isinstance(obj, Path):
                    return str(obj)
                if hasattr(obj, "__dict__"):
                    return obj.__dict__
                return str(obj)
            
            print(json.dumps(result, indent=2, ensure_ascii=False, default=default_serializer))
        else:
            status_str = "VALID" if result.is_valid else "INVALID"
            logger.info("=" * 60)
            logger.info(f"ANALYSIS BUNDLE VALIDATION: {status_str}")
            logger.info("=" * 60)
            logger.info(f"manifest: {result.manifest_path}")
            logger.info(f"bundle_root: {result.bundle_root}")
            contract = result.contract_info
            logger.info(f"contract: {contract.contract_name} {contract.contract_version} ({contract.bundle_kind})")
            logger.info(f"errors: {result.error_count}")
            logger.info(f"warnings: {result.warning_count}")
            
            if result.issues:
                logger.info("-" * 60)
                for issue in result.issues:
                    level_label = f"[{issue.level.upper()}]"
                    artifact_label = f" ({issue.artifact_name})" if issue.artifact_name else ""
                    logger.info(f"{level_label} {issue.code}{artifact_label}: {issue.message}")
            
            logger.info("=" * 60)

        if not result.is_valid:
            sys.exit(1)

    elif args.command == "inspect-analysis-bundle":
        setup_logging(args.log_level)
        
        try:
            bundle = read_analysis_bundle(args.manifest)
            summary = summarize_analysis_bundle_for_consumer(bundle)
            
            if args.json:
                def default_serializer(obj):
                    if isinstance(obj, Path):
                        return str(obj)
                    return str(obj)
                
                print(json.dumps(summary, indent=2, ensure_ascii=False, default=default_serializer))
            else:
                print("=" * 60)
                print("Analysis Bundle Inspect")
                print("=" * 60)
                print(f"Manifest:          {summary.get('analysis_bundle_manifest_path')}")
                print(f"Contract:          {summary.get('contract_name')}@{summary.get('contract_version')}")
                print(f"Bundle Kind:       {summary.get('bundle_kind')}")
                print(f"Producer:          {summary.get('producer')}@{summary.get('producer_version')}")
                print(f"Run ID:            {summary.get('run_id')}")
                print(f"Matrix ID:         {summary.get('matrix_id')}")
                
                shape = summary.get("matrix_shape", {})
                print(f"Matrix Shape:      features={shape.get('feature_count')}, samples={shape.get('sample_count')}")
                
                print(f"Sample Axis:       {summary.get('sample_axis')}")
                print(f"Feature ID System: {summary.get('feature_id_system')}")
                
                def _fmt(val):
                    if val is None: return "None"
                    if isinstance(val, (dict, list)):
                        return json.dumps(val, ensure_ascii=False)
                    return str(val)

                print(f"Column Order Specimen IDs:        {_fmt(summary.get('column_order_specimen_ids'))}")
                print(f"Source Quantifier Summary:        {_fmt(summary.get('source_quantifier_summary'))}")
                print(f"Feature Annotation Status:        {_fmt(summary.get('feature_annotation_status'))}")
                print(f"Sample Metadata Alignment Status: {_fmt(summary.get('sample_metadata_alignment_status'))}")
                print(f"Warning Summary:                  {_fmt(summary.get('warning_summary'))}")
                print("=" * 60)
                
        except Exception as e:
            print(f"Error: Failed to inspect bundle: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()

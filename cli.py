import argparse
import logging
import sys
from pathlib import Path

# Add src to sys.path to allow importing from iwa_rnaseq_counter package
sys.path.insert(0, str(Path(__file__).parent.resolve() / "src"))

from iwa_rnaseq_counter.io.read_assay_spec import read_assay_spec
from iwa_rnaseq_counter.io.read_matrix_spec import read_matrix_spec
from iwa_rnaseq_counter.io.write_matrix_spec import write_matrix_spec
from iwa_rnaseq_counter.io.write_execution_run_spec import write_execution_run_spec
from iwa_rnaseq_counter.pipeline.runner import run_counter_pipeline
from iwa_rnaseq_counter.pipeline.build_analysis_matrix import build_analysis_matrix


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


def main():
    parser = argparse.ArgumentParser(description="iwa_rnaseq_counter CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_run = subparsers.add_parser("run-assay", help="Run single AssaySpec and emit assay-level MatrixSpec")
    p_run.add_argument("--assay-spec", required=True, type=Path)
    p_run.add_argument("--outdir", required=True, type=Path)
    p_run.add_argument("--run-id", type=str)
    p_run.add_argument("--profile", type=str, default="local")
    p_run.add_argument("--quantifier", type=str, default="salmon")
    p_run.add_argument("--threads", type=int, default=4)
    p_run.add_argument("--dry-run", action="store_true")
    p_run.add_argument("--log-level", type=str, default="INFO")

    p_batch = subparsers.add_parser("run-batch", help="Run multiple AssaySpecs from a sample sheet")
    p_batch.add_argument("--sample-sheet", required=True, type=Path)
    p_batch.add_argument("--salmon-index", required=True, type=Path)
    p_batch.add_argument("--tx2gene", required=True, type=Path)
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
    p_merge.add_argument("--log-level", type=str, default="INFO")

    args = parser.parse_args()

    if args.command == "run-assay":
        setup_logging(args.log_level, args.outdir / "logs" / "counter.log")
        logger = logging.getLogger(__name__)
        assay_spec = read_assay_spec(args.assay_spec)

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
        
        try:
            assay_specs = read_sample_sheet(
                sample_sheet_path=args.sample_sheet,
                salmon_index_path=str(args.salmon_index),
                tx2gene_path=str(args.tx2gene),
                strandedness=args.strandedness,
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


if __name__ == "__main__":
    main()

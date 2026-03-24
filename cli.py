import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.resolve()))

from src.iwa_rnaseq_counter.io.read_assay_spec import read_assay_spec
from src.iwa_rnaseq_counter.io.write_matrix_spec import write_matrix_spec
from src.iwa_rnaseq_counter.io.write_execution_run_spec import write_execution_run_spec
from src.iwa_rnaseq_counter.pipeline.runner import run_counter_pipeline

def main():
    parser = argparse.ArgumentParser(description="iwa_rnaseq_counter CLI")
    parser.add_argument("--assay-spec", required=True, type=Path, help="AssaySpec JSON path")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--run-id", type=str, help="Run ID")
    parser.add_argument("--profile", type=str, default="default", help="Execution profile")
    parser.add_argument("--quantifier", type=str, default="salmon", help="Quantifier")
    parser.add_argument("--threads", type=int, default=4, help="Threads")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--write-run-spec-only", action="store_true")
    parser.add_argument("--log-level", type=str, default="INFO")

    args = parser.parse_args()

    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {args.log_level}")
    logging.basicConfig(level=numeric_level, format="%(asctime)s [%(levelname)s] %(message)s")
    logger = logging.getLogger(__name__)

    logger.info(f"Reading AssaySpec from {args.assay_spec}")
    try:
        assay_spec = read_assay_spec(args.assay_spec)
    except Exception as e:
        logger.error(f"Failed to read AssaySpec: {e}")
        sys.exit(1)

    if args.dry_run:
        logger.info("Dry run complete. Validated AssaySpec successfully.")
        return

    logger.info("Starting pipeline execution")
    try:
        matrix_spec, exec_spec = run_counter_pipeline(assay_spec, args.outdir, args.threads)
        
        write_matrix_spec(matrix_spec, args.outdir / "specs" / "matrix.spec.json")
        write_execution_run_spec(exec_spec, args.outdir / "specs" / "execution-run.spec.json")
        
        logger.info(f"Pipeline completed perfectly. MatrixSpec written to {args.outdir / 'specs' / 'matrix.spec.json'}")
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

import argparse
import logging
import sys
from pathlib import Path

# Add local src and root src to sys.path
sys.path.insert(0, str(Path(__file__).parent.resolve() / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent.resolve() / "src"))

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
   
    # [v0.6.0 C-02]
    # CLI の表面で quantifier の既定値が salmon に固定されている。
    # v0.6.0 では salmon を唯一実装のままでも、
    # 入口としては「最初の backend 実装例」に寄せたい。
    p_run.add_argument("--quantifier", type=str, default="salmon")
    p_run.add_argument("--threads", type=int, default=4)
    p_run.add_argument("--dry-run", action="store_true")
    p_run.add_argument("--log-level", type=str, default="INFO")

    p_batch = subparsers.add_parser("run-batch", help="Run multiple AssaySpecs from a sample sheet")
    p_batch.add_argument("--sample-sheet", required=True, type=Path)
    
    # [v0.6.0 C-03 / C-08]
    # CLI 引数名そのものが --salmon-index であり、
    # reference/index 入力契約が Salmon 固有語彙のまま露出している。
    # v0.6.0 では表面名を backend 非依存に寄せるか、
    # 少なくとも CLI facade で抽象名 -> backend 固有名へ変換したい。
    p_batch.add_argument("--salmon-index", required=True, type=Path)
    p_batch.add_argument("--tx2gene", required=True, type=Path)
    p_batch.add_argument("--strandedness", type=str, default="Auto-detect")
    p_batch.add_argument("--outdir", required=True, type=Path)
    p_batch.add_argument("--run-id", type=str)
    p_batch.add_argument("--profile", type=str, default="local")
   
    # [v0.6.0 C-02]
    # run-batch 側も quantifier 既定値が salmon 固定。
    # 実装の現状としては妥当でも、CLI 契約としては v0.6.0 で整理対象。
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
   
    # [v0.6.0 C-03 / C-09]
    # run-gui-backend 自体は抽象的な名前だが、
    # 実際には --config の中身が Salmon 固有契約である点が本質。
    # CLI 名よりも config schema の抽象化が優先。
    p_gui = subparsers.add_parser("run-gui-backend", help="Execute the monolithic GUI pipeline via CLI for job running")
    p_gui.add_argument("--config", required=True, type=Path)
    p_gui.add_argument("--sample-sheet", required=True, type=Path)
    p_gui.add_argument("--outdir", required=True, type=Path)
    p_gui.add_argument("--started-at", type=str, required=True)
    p_gui.add_argument("--log-level", type=str, default="INFO")

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
            # [v0.6.0 C-09]
            # CLI から runner へ quantifier 名をそのまま流している。
            # 将来的にはここで registry / resolver 解決を挟むか、
            # 少なくとも runner 側の責務と CLI 側の責務を分けたい。
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
            # [v0.6.0 C-03 / C-08]
            # read_sample_sheet の入力契約も salmon_index_path を直接要求している。
            # ここは CLI 表面だけでなく内部 I/O 契約にも Salmon 語彙が残っている証拠。
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


    elif args.command == "run-gui-backend":
        setup_logging(args.log_level, args.outdir / "logs" / "run.log")
        logger = logging.getLogger(__name__)
        
        import json
        import pandas as pd
        from iwa_rnaseq_counter.pipeline.gui_backend import run_gui_backend_pipeline
        
        try:
            with open(args.config, "r") as f:
                config_data = json.load(f)
            # [v0.6.0 C-03 / C-09]
            # run-gui-backend は config JSON を受け取るだけに見えるが、
            # その config_data の中身が現状 salmon_index_path / tx2gene_path 前提。
            # つまり CLI の抽象度より config schema の抽象度の方が低い状態。
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


if __name__ == "__main__":
    main()

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
    # --- CLI entrypoint ---
    # cli.py は counter の各実行モードに対する統一入口である。
    # assay 単位実行、batch 実行、analysis matrix 構築、GUI backend 実行を
    # サブコマンドとして束ねる。
    parser = argparse.ArgumentParser(description="iwa_rnaseq_counter CLI")
    # --- Subcommand parser construction ---
    # 利用シーンごとにサブコマンドを分け、
    # それぞれに必要な入力契約を argparse 上で定義する。
    subparsers = parser.add_subparsers(dest="command", required=True)
    # --- run-assay parser ---
    # 単一 AssaySpec を受け取り、assay 単位の matrix/spec を出力する入口。
    # 主に spec 駆動の単体実行や検証用途を想定している。
    p_run = subparsers.add_parser("run-assay", help="Run single AssaySpec and emit assay-level MatrixSpec")
    p_run.add_argument("--assay-spec", required=True, type=Path)
    p_run.add_argument("--outdir", required=True, type=Path)
    p_run.add_argument("--run-id", type=str)
    p_run.add_argument("--profile", type=str, default="local")
   
    # [v0.6.0 C-02]
    # CLI は quantifier を受け取れる形になっているが、
    # 現在の既定値は salmon であり、実装もまずは Salmon backend を基準にしている。
    # 入口としては将来の複数 backend を見越しつつ、現段階では既定値で運用する。
    p_run.add_argument("--quantifier", type=str, default="salmon")
    p_run.add_argument("--threads", type=int, default=4)
    p_run.add_argument("--dry-run", action="store_true")
    p_run.add_argument("--log-level", type=str, default="INFO")

    # --- run-batch parser ---
    # sample sheet から複数 AssaySpec 相当を解釈し、
    # specimen 単位で連続実行する batch 入口である。
    p_batch = subparsers.add_parser("run-batch", help="Run multiple AssaySpecs from a sample sheet")
    p_batch.add_argument("--sample-sheet", required=True, type=Path)
    
    # [v0.6.0 C-03 / C-08]
    # quantifier 自体は抽象化の入口を持ち始めたが、
    # reference/index 側の CLI 契約には依然として --salmon-index という
    # Salmon 固有語彙が残っている。
    # 今後は CLI 表面名と内部 I/O 契約の両方を段階的に整理したい。
    p_batch.add_argument("--salmon-index", required=True, type=Path)
    p_batch.add_argument("--tx2gene", required=True, type=Path)
    p_batch.add_argument("--strandedness", type=str, default="Auto-detect")
    p_batch.add_argument("--outdir", required=True, type=Path)
    p_batch.add_argument("--run-id", type=str)
    p_batch.add_argument("--profile", type=str, default="local")
   
    # [v0.6.0 C-02]
    # run-batch でも quantifier を受け取れるようにしている。
    # ただし既定値は引き続き salmon であり、
    # まずは既存実装を壊さずに複数 backend の差し込み口だけを確保している。
    p_batch.add_argument("--quantifier", type=str, default="salmon")
    p_batch.add_argument("--threads", type=int, default=4)
    p_batch.add_argument("--dry-run", action="store_true")
    p_batch.add_argument("--log-level", type=str, default="INFO")

    # --- build-analysis-matrix parser ---
    # assay 単位の MatrixSpec 群を統合し、
    # reporter 側で扱いやすい analysis 単位の matrix/spec を構築する入口
    p_merge = subparsers.add_parser("build-analysis-matrix", help="Merge assay-level MatrixSpecs into reporter-ready analysis MatrixSpec")
    p_merge.add_argument("--matrix-spec", required=True, nargs="+", type=Path)
    p_merge.add_argument("--sample-metadata", required=True, type=Path)
    p_merge.add_argument("--outdir", required=True, type=Path)
    p_merge.add_argument("--matrix-id", required=True, type=str)
    p_merge.add_argument("--run-id", type=str)
    p_merge.add_argument("--log-level", type=str, default="INFO")
   
    # [v0.6.0 C-03 / C-09]
    # run-gui-backend 自体は抽象的な名前だが、
    # 実際に受け取る config schema には reference 側の Salmon 語彙がまだ残っている。
    # 一方で quantifier / quantifier_version もこの経路で受け渡すようになっており、
    # 今後は config schema 全体の backend 非依存化を進めたい。
    # --- run-gui-backend parser ---
    # GUI から生成された run artifact を受け取り、
    # background job として GUI backend pipeline を再実行する入口。
    # app.py が作成した config/sample sheet を CLI 経由で backend へ橋渡しする。
    p_gui = subparsers.add_parser("run-gui-backend", help="Execute the monolithic GUI pipeline via CLI for job running")
    p_gui.add_argument("--config", required=True, type=Path)
    p_gui.add_argument("--sample-sheet", required=True, type=Path)
    p_gui.add_argument("--outdir", required=True, type=Path)
    p_gui.add_argument("--started-at", type=str, required=True)
    p_gui.add_argument("--log-level", type=str, default="INFO")

    args = parser.parse_args()

    # --- run-assay execution path ---
    # 単一 AssaySpec を読み込み、runner に処理を委譲して
    # assay-level artifact を出力する実行経路。
    if args.command == "run-assay":
        # --- Logging setup ---
        # 各サブコマンドで共通利用する logging 初期化関数。
        # 標準出力と必要に応じた logfile の両方へ同一フォーマットで出力する。
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
            # CLI は受け取った quantifier 名をそのまま runner に渡す。
            # backend 解決そのものは runner 側の registry に委譲し、
            # CLI は入力契約の受け口に留める方針である。
            quantifier=args.quantifier,
        )

        write_matrix_spec(matrix_spec, args.outdir / "specs" / "matrix.spec.json")
        write_execution_run_spec(exec_spec, args.outdir / "specs" / "execution-run.spec.json")
        logger.info("run-assay completed")

    elif args.command == "run-batch":
        # --- run-batch execution path ---
        # sample sheet を AssaySpec 群へ変換し、
        # 各 assay を順に runner へ渡して batch 実行する経路。
        setup_logging(args.log_level, args.outdir / "logs" / "batch_counter.log")
        logger = logging.getLogger(__name__)

        from iwa_rnaseq_counter.io.read_sample_sheet import read_sample_sheet
        
        try:
            # [v0.6.0 C-03 / C-08]
            # read_sample_sheet の入力契約も salmon_index_path を直接要求している。
            # ここは CLI 表面だけでなく内部 I/O 契約にも Salmon 語彙が残っている証拠。
            # [v0.6.0 C-03 / C-08]
            # sample sheet 読み込みの内部契約にも salmon_index_path が残っている。
            # つまり quantifier 実行の抽象化とは別に、
            # reference/resource 側の I/O 契約整理がまだ必要である。
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

        # --- Per-assay batch execution ---
        # 読み込まれた各 assay を specimen 単位の出力ディレクトリへ実行し、
        # 成功したものだけを結果一覧へ積み上げる。
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
        # --- build-analysis-matrix execution path ---
        # 複数の assay-level MatrixSpec を読み込み、
        # analysis 単位の matrix と execution spec を生成する経路。
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
        # --- GUI backend execution path ---
        # GUI が事前に保存した config と sample sheet を読み込み、
        # GUI backend pipeline にそのまま受け渡して background 実行する。
        setup_logging(args.log_level, args.outdir / "logs" / "run.log")
        logger = logging.getLogger(__name__)
        
        import json
        import pandas as pd
        from iwa_rnaseq_counter.pipeline.gui_backend import run_gui_backend_pipeline
        
        try:
            # --- GUI config loading ---
            # app.py が保存した run_config.json を読み込み、
            # backend 実行に必要な設定を復元する。
            with open(args.config, "r") as f:
                config_data = json.load(f)
            # [v0.6.0 C-03 / C-09]
            # run-gui-backend は config JSON を受け取り、そのまま gui_backend へ渡す。
            # 現在の config schema には salmon_index_path / tx2gene_path などの
            # Salmon 語彙が残る一方、quantifier / quantifier_version もこの経路で通る。
            # 今後は config schema 自体の backend 非依存化を進めたい。
            # --- Sample sheet deserialization ---
            # CSV へ保存された path list 列を Python の list へ戻し、
            # backend pipeline が扱える形へ整形する。
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

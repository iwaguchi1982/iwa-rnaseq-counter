#
# runner.py
# 全体の役割：アッセイ（単一サンプル）単位のパイプライン実行オーケストレーター。AssaySpec を受け取り、定量計算、遺伝子集計、および標準規格の Spec 出力を担当。
# 抽象化の維持：具体的な計算ロジックは quantifier registry 経由で各 backend に委譲し、この層では共通の処理フローを維持する。
# ---
# 主要な処理ステップ
# 1. 出力ディレクトリの初期化:
# アッセイごとの実行結果を格納するための counts, logs, specs ディレクトリを自動作成し、成果物の責務を分離。
# 
# 2. 入力データの正規化:
# AssaySpec から FASTQ ファイルのパスを抽出し、Single-end/Paired-end の判定と backend が解釈可能な DataFrame 形式への変換を行う。
# 
# 3. リソースの解決と検証:
# 実行に必要な quantifier index や tx2gene マップを AssaySpec から取り出す。実行前にインデックスの存在を最低限チェックし、エラーを未然に防ぐ。
# 
# 4. 定量計算の実行:
# 指定された Quantifier（Salmon 等）を呼び出し、共通インターフェースを通じて定量処理を実行。runner 自身は各ツールのコマンド詳細を直接は関知しない。
# 
# 5. 遺伝子レベルへの集約と保存:
# 転写産物レベルの定量結果を遺伝子レベルのカウント行列（raw count）へ集約。結果を TSV 形式で counts ディレクトリに保存。
# 
# 6. 標準 Spec (Contract) の生成:
# 解析結果のメタデータを MatrixSpec として、実行記録を ExecutionRunSpec として構築。
# これにより、後続の解析ステップや統合処理がツールの種類に依存せず、型安全に結果を利用できるようにする。
# 

import logging
import pandas as pd
from datetime import datetime, timezone
from pathlib import Path
from iwa_rnaseq_counter.pipeline.quantifiers.registry import get_quantifier
from iwa_rnaseq_counter.legacy.gene_aggregator import load_tx2gene_map, build_transcript_quant_table, aggregate_transcript_to_gene
from ..models.assay import AssaySpec
from ..models.matrix import MatrixSpec
from ..models.execution_run import ExecutionRunSpec

# --- Assay-level pipeline orchestration ---
# runner.py は AssaySpec を受け取り、
# assay 単位の quantification 実行、gene-level 集約、spec 出力情報の組み立てを担当する。
# backend 実装の解決は quantifier registry に委譲し、この層では共通フローを維持する。
logger = logging.getLogger(__name__)

# --- Assay pipeline entrypoint ---
# 1 assay 分の入力を受け取り、counts / logs / specs の出力先を初期化し、
# downstream で利用する MatrixSpec / ExecutionRunSpec を返す。
def run_counter_pipeline(
    assay_spec: AssaySpec,
    outdir: Path,
    threads: int = 4,
    run_id: str | None = None,
    profile: str = "local",
    quantifier: str = "salmon",
) -> tuple[MatrixSpec, ExecutionRunSpec]:
    started_at = datetime.now(timezone.utc).astimezone().isoformat()

    # --- Output directory preparation ---
    # assay 実行に必要な出力ディレクトリを先に確保する。
    # counts, logs, specs をここで分けておくことで、artifact の責務を明確にする。
    outdir.mkdir(parents=True, exist_ok=True)
    counts_dir = outdir / "counts"
    logs_dir = outdir / "logs"
    specs_dir = outdir / "specs"
    counts_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    specs_dir.mkdir(exist_ok=True)

    # --- Input file extraction ---
    # AssaySpec から FASTQ 入力を取り出し、
    # quantifier 実行に必要な最小構造へ正規化する。
    fastq_r1 = next((f.path for f in assay_spec.input_files if f.file_role.lower() == "fastq_r1"), None)
    fastq_r2 = next((f.path for f in assay_spec.input_files if f.file_role.lower() == "fastq_r2"), None)

    layout_final = "paired-end" if fastq_r2 else "single-end"
    r1_paths = [fastq_r1] if fastq_r1 else []
    r2_paths = [fastq_r2] if fastq_r2 else []
    all_paths = r1_paths if layout_final == "single-end" else []

    # --- Minimal sample dataframe assembly ---
    # legacy quantification 実装が受け取れる形式に合わせて、
    # 1 assay / 1 specimen 分の sample_df を最小構成で組み立てる。
    sample_df = pd.DataFrame([{
        "sample_id": assay_spec.specimen_id,
        "layout_final": layout_final,
        "r1_paths": r1_paths,
        "r2_paths": r2_paths,
        "all_paths": all_paths,
    }])

    # --- Reference resource extraction ---
    # AssaySpec に含まれる reference resource から、
    # quantifier 実行と gene-level 集約に必要な参照情報を取り出す。
    salmon_index = None
    tx2gene = None
    if assay_spec.reference_resources:
        # [v0.6.0 C-08]
        # reference_resources から quantifier 実行用 index と tx2gene を取り出している。
        # quantifier 実行の抽象化は進んだが、reference 側の命名と契約には
        # まだ Salmon 由来の概念が残っている。
        # 今後は backend 非依存な reference/index/resource 契約へ段階的に寄せたい。
        salmon_index = assay_spec.reference_resources.quantifier_index
        tx2gene = assay_spec.reference_resources.tx2gene_path

    # --- Required reference validation ---
    # 現行 backend 実行では quantifier index が必須であるため、
    # 実行前に最低限の前提条件をここで確認する。
    if not salmon_index:
        raise ValueError("salmon_index is required in AssaySpec.reference_resources")

    # --- Quantifier execution ---
    # 指定された quantifier 名から backend 実装を解決し、
    # 共通インターフェース経由で quantification を実行する。
    # runner 自体は backend 実装の詳細を直接知らない。
    quant = get_quantifier(quantifier)

    run_result = quant.run_quant(
        sample_df=sample_df,
        run_output_dir=outdir,
        threads=threads,
        strandedness_mode=assay_spec.strandedness or "Auto-detect",
        reference_config={
            "quantifier_index": salmon_index,
            "tx2gene_path": tx2gene,
        },
    )

    # --- Quantifier result validation ---
    # backend 実行結果から per-sample 出力を取り出し、
    # downstream 集約へ進める最低条件を確認する
    outputs = run_result["outputs"]
    # TODO(v0.6.x): 失敗メッセージも backend 非依存な表現へ寄せる余地がある。
    if not outputs or not outputs[0].get("is_success"):
        raise RuntimeError("Salmon quantification failed.")

    if not tx2gene:
        raise NotImplementedError("Transcript-level un-aggregated output is not currently handled without tx2gene.")

    # --- Gene-level aggregation ---
    # quantifier 出力から transcript-level の NumReads を取得し、
    # tx2gene を使って gene-level count matrix へ集約する。
    tx2gene_df = load_tx2gene_map(tx2gene)
    t_nr_df = build_transcript_quant_table(outputs, value_type="NumReads")
    g_nr_df = aggregate_transcript_to_gene(t_nr_df, tx2gene_df)

    # --- Count matrix export ---
    # reporter や downstream spec で参照できるように、
    # gene-level raw count matrix を assay 出力として保存する。
    matrix_path = counts_dir / "gene_numreads.tsv"
    g_nr_df.to_csv(matrix_path, sep="\t")
    
    # --- Feature annotation preparation ---
    # reporter 側で使う feature_annotation.tsv を準備する。
    # annotation が無い場合でも matrix 自体は保持し、表示時に feature_id fallback を許容する。
    from iwa_rnaseq_counter.legacy.annotation_helper import prepare_feature_annotation, get_standard_annotation_path
    annotation_out = get_standard_annotation_path(outdir)
    # Ensure results dir exists
    (outdir / "results").mkdir(exist_ok=True)
    has_annotation = prepare_feature_annotation(tx2gene, annotation_out)

    subject_id = assay_spec.metadata.get("subject_id")
    source_subject_ids = [subject_id] if subject_id else []

    # --- MatrixSpec assembly ---
    # assay 実行で得られた gene-level matrix を、
    # downstream が読める標準 spec として表現する。
    matrix_spec = MatrixSpec(
        schema_name="MatrixSpec",
        schema_version="0.1.0",
        matrix_id=f"MAT_{assay_spec.assay_id}",
        matrix_scope="assay",
        matrix_kind="count_matrix",
        feature_type="gene",
        value_type="integer",
        normalization="raw",
        feature_id_system="ensembl_gene_id",
        # Use specimen as axis
        sample_axis="specimen",
        matrix_path=str(matrix_path.resolve()),
        feature_annotation_path=str(annotation_out.resolve()) if has_annotation else None,
        source_assay_ids=[assay_spec.assay_id],
        source_specimen_ids=[assay_spec.specimen_id],
        source_subject_ids=source_subject_ids,
        metadata={
            "producer_app": "iwa_rnaseq_counter",
            "producer_version": "0.3.5",
            "quantifier": quantifier,
            # [v0.6.0 C-09]
            # backend 情報が metadata に散発的に入っている。
            # quantifier 名や index 参照は正式な backend metadata 領域へ整理したい。
            # 特に "salmon_index" というキー名は backend 固有なので、
            # v0.6.0 で汎化候補。
            "salmon_index": str(salmon_index),
            "tx2gene_path": str(tx2gene),
            "sample_ids": [assay_spec.specimen_id],
            "feature_id_system_inferred": False,
            "feature_annotation_available": has_annotation,
        },
        overlay={},
    )

    finished_at = datetime.now(timezone.utc).astimezone().isoformat()

    run_spec = ExecutionRunSpec(
        schema_name="ExecutionRunSpec",
        schema_version="0.1.0",
        run_id=run_id or f"RUN_COUNTER_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        app_name="iwa_rnaseq_counter",
        app_version="0.3.5",
        input_refs=[assay_spec.assay_id],
        output_refs=[matrix_spec.matrix_id],
        parameters={
            "quantifier": quantifier,
            "threads": threads,
            "strandedness_mode": assay_spec.strandedness or "Auto-detect",
            "profile": profile,
            # [v0.6.0 C-09]
            # parameters には run-time option と execution context の一部が入る。
            # 現在は quantifier 名や reference 参照も含めているが、
            # 将来的には run option と backend identity / resource refs の境界を整理したい。
            "salmon_index": str(salmon_index),
            "tx2gene_path": str(tx2gene),
        },
        execution_backend=profile,
        started_at=started_at,
        finished_at=finished_at,
        status="completed",
        log_path=str((logs_dir / "run.log").resolve()),
    )
    # --- Pipeline return values ---
    # assay-level の matrix spec と execution run spec を返し、
    # 呼び出し側で永続化または downstream 処理へ接続できるようにする。
    return matrix_spec, run_spec

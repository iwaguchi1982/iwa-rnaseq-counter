#
# gui_artifact_export.py
# 全体の役割：GUI 実行で得られた成果物（Run artifact）を、標準的な Spec 形式（MatrixSpec, ExecutionRunSpec）へ変換・排出するアダプター層。
# 責務：GUI 固有の結果構造を、後続の統計解析ツール（Reporter 等）が共通契約として解釈できる形式へ正規化する。
# ---
# 主要な機能ブロック
# 1. MatrixSpec アダプター:
# GUI 実行結果から解析単位の MatrixSpec を構築。GUI 側で保持していないアッセイ ID をサンプル ID から合成し、行列のメタデータを補完。
# アノテーションの可用性や、生成元のアプリ情報などのトレーサビリティ情報を付与。
# 
# 2. ExecutionRunSpec アダプター:
# GUI 実行全体を一つの実行記録として表現。実行時の入力モード、スレッド数、使用した Quantifier などのコンテキストを記録。
# 成功・失敗サンプル数に基づいた実行ステータスの判定と、ログパスの紐付けを担当。
# 
# 3. Spec 書き出しエントリーポイント (write_gui_supporting_inputs):
# 実行ディレクトリ内の specs/ フォルダを確保し、MatrixSpec および ExecutionRunSpec を JSON 形式で永続化。
# GUI 処理の最終ステップとして、解析結果が標準エコシステムに接続可能な状態であることを保証する。
# 

import logging
from datetime import datetime, timezone
from pathlib import Path

from iwa_rnaseq_counter.models.matrix import MatrixSpec
from iwa_rnaseq_counter.models.execution_run import ExecutionRunSpec
from iwa_rnaseq_counter.io.write_matrix_spec import write_matrix_spec
from iwa_rnaseq_counter.io.write_execution_run_spec import write_execution_run_spec

# --- GUI artifact adapter ---
# gui_artifact_export.py は、GUI 実行で得られた run artifact を
# 標準 Spec 系へ橋渡しする adapter 層である。
# core 実行ロジックではなく、GUI 結果を downstream が読める契約へ正規化する責務を持つ。
logger = logging.getLogger(__name__)

# --- MatrixSpec adapter from GUI result ---
# GUI 実行結果から、analysis 単位の MatrixSpec を組み立てる。
# GUI 側には assay ID が明示的に存在しないため、必要な識別子はここで補う。
def build_matrix_spec_from_gui_result(
    run_dir: Path,
    run_summary: dict,
    matrix_rel_path: str = "results/gene_numreads.csv",
    feature_annotation_path: str | None = None,
    feature_annotation_available: bool = False
) -> MatrixSpec:
    # GUI 出力ディレクトリ内の gene-level matrix への実体パスを解決する。
    matrix_abs_path = (run_dir / matrix_rel_path).resolve()
    
    # GUI 実行では assay ID が独立していないため、
    # successful sample_id から downstream 用の assay ID を合成する。
    sample_ids = run_summary.get("sample_ids_success", [])
    assay_ids = [f"ASSAY_{sid}" for sid in sample_ids]
    
    # GUI run 名をもとに、analysis 単位の matrix_id を安定して組み立てる。
    run_name = run_summary.get("run_name", "UNKNOWN_RUN")
    matrix_id = f"MAT_{run_name}"
    
    # --- Matrix metadata assembly ---
    # GUI 実行由来であること、入力モード、sample_ids、annotation 可用性など、
    # downstream 表示や追跡に必要な最小 metadata をここで付与する。
    metadata = {
        "producer_app": "iwa_rnaseq_counter",
        "producer_version": "0.3.5",
        "run_origin": "gui",
        "gui_mode": run_summary.get("discovery_mode", "unknown"),
        "sample_ids": sample_ids,
        "input_dir": run_summary.get("input_dir", ""),
        "feature_id_system_inferred": False,
        "feature_annotation_available": feature_annotation_available,
        
        # [v0.6.0 C-09]
        # MatrixSpec.metadata は比較的 backend 非依存に保ちたい。
        # 現時点では GUI 由来の実行文脈だけを最小限に載せ、
        # backend / quantifier 情報をここで過剰に再解釈しない方が安全である。
        # どこまでを MatrixSpec.metadata に持たせるかは今後も整理対象
    }
    
    # GUI 結果を downstream 共通契約の MatrixSpec として返す。
    return MatrixSpec(
        schema_name="MatrixSpec",
        schema_version="0.1.0",
        matrix_id=matrix_id,
        matrix_scope="analysis",
        matrix_kind="count_matrix",
        feature_type="gene",
        value_type="integer",
        normalization="raw",
        feature_id_system="ensembl_gene_id",
        sample_axis="specimen",
        matrix_path=str(matrix_abs_path),
        feature_annotation_path=feature_annotation_path,
        source_assay_ids=assay_ids,
        source_specimen_ids=sample_ids,
        source_subject_ids=sample_ids,
        metadata=metadata,
        overlay={}
    )

# --- ExecutionRunSpec adapter from GUI result ---
# GUI 実行全体を 1 つの実行記録として表現する ExecutionRunSpec を組み立てる。
# downstream から見たときに、GUI 実行も通常 run と同じ形式で追跡できるようにする。
def build_execution_run_spec_from_gui_result(
    run_dir: Path,
    run_summary: dict,
    matrix_spec: MatrixSpec,
    started_at: str,
    log_rel_path: str = "logs/run.log"
) -> ExecutionRunSpec:
    # GUI 実行ログの実体パスを解決する。
    log_abs_path = (run_dir / log_rel_path).resolve()

    # GUI run 名を execution run ID の基礎として利用する。
    # 名前が無い場合は GUI 実行用の fallback ID を生成する。
    run_name = run_summary.get("run_name", f"RUN_GUI_{datetime.now().strftime('%Y%m%d_%H%M%S')}")

    # --- Execution parameter assembly ---
    # GUI 実行時の入力モード、sample 数、strandedness、threads など、
    # 実行時コンテキストとして追跡したい項目を parameters にまとめる。
    parameters = {
        "run_origin": "gui",
        "input_mode": run_summary.get("discovery_mode", "unknown"),
        "input_dir": str(run_summary.get("input_dir", "")),
        "sample_count": run_summary.get("sample_count", 0),
        "success_count": run_summary.get("success_count", 0),
        "strandedness_mode": run_summary.get("strandedness", {}).get("mode", "Auto-detect")
        if isinstance(run_summary.get("strandedness"), dict) else "Auto-detect",
        "threads": run_summary.get("threads", 4),

        # [v0.6.0 C-05 / C-09]
        # GUI adapter は ExecutionRunSpec.parameters に quantifier 情報を記録している。
        # 現在は gui_backend が構築した run_summary 由来の実行文脈をそのまま受け取り、
        # adapter 層で backend 名を固定値として再定義しないようにしている。
        # 将来的には execution context の注入境界をさらに整理したい。
        "quantifier": run_summary.get("quantifier", "salmon"),
        "quantifier_version": run_summary.get("quantifier_version", "unknown"),
    }
    
    # MatrixSpec 側と対応が取れるように、successful sample_id から input_refs 用の assay ID を合成する。
    sample_ids = run_summary.get("sample_ids_success", [])
    assay_ids = [f"ASSAY_{sid}" for sid in sample_ids]
    
    # GUI 実行を標準的な実行記録として返す。
    # これにより、GUI 実行結果も spec ベースの downstream 処理に接続できる。
    return ExecutionRunSpec(
        schema_name="ExecutionRunSpec",
        schema_version="0.1.0",
        run_id=run_name,
        app_name="iwa_rnaseq_counter",
        app_version="0.3.5",
        input_refs=assay_ids,
        output_refs=[matrix_spec.matrix_id],
        parameters=parameters,

        # [v0.6.0 C-09]
        # execution_backend="local-gui" は実行形態を表す値として使っている。
        # ただし quantifier/backend 種別とは別概念なので、
        # execution_backend の意味定義は今後さらに明確にしたい。
        execution_backend="local-gui",
        started_at=started_at,
        finished_at=datetime.now(timezone.utc).astimezone().isoformat(),
        status="completed" if run_summary.get("failure_count", 1) == 0 else "completed_with_errors",
        log_path=str(log_abs_path)
    )

# --- GUI spec export entrypoint ---
# GUI 実行結果から MatrixSpec / ExecutionRunSpec を生成し、
# specs/ 配下へ保存するためのメイン入口である。
def write_gui_supporting_inputs(
    run_dir: Path,
    run_summary: dict,
    matrix_rel_path: str = "results/gene_numreads.csv",
    log_rel_path: str = "logs/run.log",
    started_at: str | None = None,
    feature_annotation_path: str | None = None,
    feature_annotation_available: bool = False
) -> None:
    # started_at が明示されない場合は、spec 生成時刻を fallback として補う。
    if started_at is None:
        started_at = datetime.now(timezone.utc).astimezone().isoformat()
    
    # --- Spec output directory preparation ---
    # GUI adapter が生成する spec 類の保存先を確保する。
    specs_dir = run_dir / "specs"
    specs_dir.mkdir(parents=True, exist_ok=True)

    # [v0.6.0 C-03 / C-09]
    # この関数は GUI 実行結果を標準 Spec 群へ橋渡しする adapter の本体である。
    # backend 実行そのものではなく、run_summary と annotation 情報を
    # downstream 共通契約へ写像する責務を持つ。
    # backend 由来 metadata の注入境界は今後さらに整理したい。
    
    # --- MatrixSpec generation and persistence ---
    # GUI 結果から MatrixSpec を生成し、specs/ 配下へ保存する。
    matrix_spec = build_matrix_spec_from_gui_result(
        run_dir, 
        run_summary, 
        matrix_rel_path,
        feature_annotation_path=feature_annotation_path,
        feature_annotation_available=feature_annotation_available
    )
    write_matrix_spec(matrix_spec, specs_dir / "matrix.spec.json")
    
    # --- ExecutionRunSpec generation and persistence ---
    # GUI 実行全体を表す ExecutionRunSpec を生成し、specs/ 配下へ保存する。
    exec_run_spec = build_execution_run_spec_from_gui_result(run_dir, run_summary, matrix_spec, started_at, log_rel_path)
    write_execution_run_spec(exec_run_spec, specs_dir / "execution-run.spec.json")
    
    # GUI adapter による spec 生成が正常完了したことを記録する。
    logger.info("Successfully generated GUI adapter specs in specs/ folder")

#
# execution_run.py
# 全体の役割：1 回の実行（Run）を追跡するための標準契約「ExecutionRunSpec」のデータモデル定義。
# 役割の分離：実行内容そのものではなく、「何が、いつ、どのような設定で実行され、結果はどうだったか」という実行記録としての責務を持つ。
# ---
# 主要な構成要素
# 1. ExecutionRunSpec モデル (dataclass):
# アッセイ実行、バッチ実行、GUI 実行を問わず、1 回の実行単位を共通の形式で表現。アプリ名、バージョン、一意な実行 ID（run_id）を保持。
# 
# 2. 実行の起源と成果物の紐付け (Provenance):
# 開始・終了時刻と合わせて、入力（input_refs）および出力（output_refs）の参照リストを保持。解析の連鎖（リネージ）を後追い可能にする。
# 
# 3. 実行時パラメータ (Parameters):
# 実行時に与えられたオプションや、使用した Quantifier、リソース参照（スレッド数やインデックスパス）などを保持。再現性のための実行文脈を記録。
# 
# 4. 実行状態とモニタリング:
# 実行ステータス（pending, completed, failed 等）とログファイルへのパス（log_path）を保持。実行の監視と、障害発生時の追跡を容易にする。
# 
# 5. 実行環境の識別 (Execution Backend):
# local, local-gui など、どのような環境・形態で実行されたかを識別するラベルを保持。
# 
# 6. シリアライズ (to_dict):
# 保存・永続化に適した dict 表現への変換を担当。$schema_name や $schema_version といった共通メタデータ規格に従って整形。
# 

from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any

# --- Execution run specification model ---
# ExecutionRunSpec は 1 回の実行を追跡するための標準契約である。
# 入力参照、出力参照、実行時パラメータ、状態、ログ位置などを保持し、
# downstream から「何がどう実行されたか」を後追いできるようにする。
@dataclass

# --- Run record contract ---
# このモデルは assay 実行、batch 実行、GUI 実行などを問わず、
# 1 run を共通の shape で表現するための dataclass である。
# 実行内容そのものではなく、実行記録としての責務を持つ。
class ExecutionRunSpec:
    schema_name: str
    schema_version: str
    run_id: str
    app_name: str
    app_version: str

    # --- Run timing and provenance ---
    # この run がいつ開始され、何を入力として受け取り、
    # 何を出力したかを追跡するための基本情報を保持する。
    started_at: str
    input_refs: List[str] = field(default_factory=list)
    output_refs: List[str] = field(default_factory=list)

    # [v0.6.0 C-09]
    # parameters は run-time option の受け皿だが、
    # 現状は quantifier / salmon_index / tx2gene_path なども流入している。
    # 今後は「実行オプション」と「backend identity / resource refs」の
    # 置き場所をより明確に整理したい。

    # --- Runtime parameter slot ---
    # 実行時に与えた option や実行文脈の一部を保持する拡張スロット。
    # ただし backend identity や resource reference まで無秩序に混ぜないよう注意が必要。
    parameters: Dict[str, Any] = field(default_factory=dict)

    # [v0.6.0 C-09]
    # execution_backend は "local" / "local-gui" のように
    # 実行形態を表す値として使っている。
    # quantifier 情報とは別軸なので、今後も意味を明確に分けて扱いたい。

    # --- Execution environment label ---
    # この run がどの実行形態で走ったかを表すラベルを保持する。
    # ただし quantifier/backend 種別とは別概念なので、意味を混ぜないようにしたい。
    execution_backend: Optional[str] = None

    # --- Run status and logging ---
    # 実行完了時刻、実行状態、ログ参照先を保持し、
    # 実行監視や障害時の追跡に使えるようにする。
    finished_at: Optional[str] = None

    status: str = "pending"
    log_path: str = ""

    # [v0.6.0 C-09]
    # metadata は将来の拡張先として有力である。
    # ただし parameters の肥大化対策として使う場合でも、
    # backend 固有情報をどの規約で載せるかを決めてから扱うべきである。

    # --- Extensible run metadata ---
    # 共通フィールドでは表しきれない補足的な実行情報を保持する拡張スロット。
    # ただし parameters の代替ゴミ箱にはせず、意味づけを保ったまま使いたい。
    metadata: Dict[str, Any] = field(default_factory=dict)

    # --- Overlay slot ---
    # 派生的・補助的な上書き情報や追加情報を保持するための拡張スロット。
    # core run contract を崩さずに局所的な情報を持たせたい場合に使う。
    overlay: Dict[str, Any] = field(default_factory=dict)

    # --- Serialization to stored spec ---
    # ExecutionRunSpec を永続化しやすい dict 表現へ変換する。
    # schema key の整形と None 除去のみを行い、意味解釈はここで増やさない。
    def to_dict(self) -> Dict[str, Any]:
        # [v0.6.0 C-09]
        # serializer 自体は汎用のままでよい。
        # backend 差分の吸収や特殊処理をここで始めるのではなく、
        # writer / adapter 側で整理する方針を維持したい。
        d = asdict(self)
        d["$schema_name"] = d.pop("schema_name")
        d["$schema_version"] = d.pop("schema_version")
        return {k: v for k, v in d.items() if v is not None}

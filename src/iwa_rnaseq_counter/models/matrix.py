#
# matrix.py
# 全体の役割：行列データ（Matrix artifact）を扱うための標準契約「MatrixSpec」のデータモデル定義。
# 役割の分離：計算コードから「データがどこにあり、どのような構造か」という情報を分離し、Reporter などの下流ツールが型安全にデータを扱えるようにする。
# ---
# 主要な構成要素
# 1. MatrixSpec モデル (dataclass):
# アッセイ単位および解析単位の行列成果物を共通の形式で表現。実データのパス（matrix_path）に加えて、軸（sample_axis）や特徴量（feature_type）の定義を保持。
# 
# 2. レポーター向けアノテーション参照:
# feature_annotation.tsv（UI表示用アノテーション）への参照を保持。tx2gene などの一次リソースとは明確に区別し、ダウンストリーム側の表示責務を支える。
# 
# 3. トレーサビリティ (Provenance):
# 行列の由来となる Assay ID, Specimen ID, Subject ID のリストを保持。最終的な解析結果から元のサンプルや実験への追跡を可能にする。
# 
# 4. 拡張スロット (metadata / overlay):
# 基本契約（Contract）に含まれない補足情報や、一時的な上書き情報を保持。特定のバックエンド依存情報の流入を制御しつつ、柔軟なメタデータ拡張を許容。
# 
# 5. シリアライズ・復元 (to_dict / from_dict):
# 共通契約（$schema_name, $schema_version）に基づいた永続化と復元を担当。保存された JSON ファイルとメモリ内のデータモデルを安定して橋渡しする。
# 

from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any


# --- Matrix specification model ---
# MatrixSpec は downstream が matrix artifact を読むための標準契約である。
# 実データ本体の path、軸情報、feature 情報、由来情報を保持し、
# 実行 backend の詳細そのものは持ち込みすぎない方針で設計する。

# --- Matrix artifact contract ---
# このモデルは assay-level / analysis-level を問わず、
# matrix 系成果物を共通の shape で表現するための dataclass である。
# runner / gui adapter / reporter 側の受け渡しの土台になる。
@dataclass
class MatrixSpec:
    schema_name: str
    schema_version: str
    matrix_id: str
    matrix_scope: Optional[str]
    matrix_kind: str
    feature_type: str
    value_type: str
    normalization: str
    feature_id_system: str
    sample_axis: str
    matrix_path: str

    # [v0.6.0 C-08 / C-09]
    # このフィールド自体は backend 非依存で維持したい。
    # 呼び出し側が tx2gene_path や backend 固有 annotation 資源を
    # 直接ここへ流し込まないことを明確に守る必要がある。
    # feature_annotation_path は reporter 表示用 annotation の参照に限定する。

    # --- Reporter-facing annotation reference ---
    # reporter 表示で使う feature annotation への参照を保持する。
    # tx2gene そのものではなく、UI / report 向けに整形済みの annotation 成果物を指す。
    feature_annotation_path: Optional[str] = None

    # --- Provenance references ---
    # この matrix がどの assay / specimen / subject に由来するかを保持する。
    # downstream で traceability を保つための基本参照情報である。
    source_assay_ids: List[str] = field(default_factory=list)
    source_specimen_ids: List[str] = field(default_factory=list)
    source_subject_ids: List[str] = field(default_factory=list)

    # [v0.6.0 C-09]
    # metadata は拡張余地として必要だが、
    # salmon_index / tx2gene_path / backend 固有 key を無秩序に積む場所にはしたくない。
    # 今後は「共通 metadata」と「backend 由来 metadata」の整理方針を明確にしたい。

    # --- Extensible matrix metadata ---
    # 共通契約では表現しきれない補足情報を保持する拡張スロット。
    # ただし backend 固有情報の無秩序な流入先にはしない。
    metadata: Dict[str, Any] = field(default_factory=dict)

    # --- Overlay slot ---
    # 派生的・補助的な上書き情報や付加情報を載せるための拡張スロット。
    # core contract を汚さずに局所的な補助情報を持たせたい場合に使う。
    overlay: Dict[str, Any] = field(default_factory=dict)

    # --- Deserialization from stored spec ---
    # 保存済みの dict 表現から MatrixSpec を復元する。
    # reader 側ではこの入口を通して共通契約へ戻す。
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MatrixSpec":

        # [v0.6.0 C-09]
        # ここは受け皿としては汎用のままでよい。
        # 重要なのは upstream writer が backend 固有 key をどう入れるかであり、
        # from_dict 側で backend 差分吸収を始めない方が安全である。
        return cls(
            schema_name=data.get("$schema_name", ""),
            schema_version=data.get("$schema_version", ""),
            matrix_id=data.get("matrix_id", ""),
            matrix_scope=data.get("matrix_scope"),
            matrix_kind=data.get("matrix_kind", ""),
            feature_type=data.get("feature_type", ""),
            value_type=data.get("value_type", ""),
            normalization=data.get("normalization", ""),
            feature_id_system=data.get("feature_id_system", ""),
            sample_axis=data.get("sample_axis", ""),
            matrix_path=data.get("matrix_path", ""),
            feature_annotation_path=data.get("feature_annotation_path") or None,
            source_assay_ids=data.get("source_assay_ids", []),
            source_specimen_ids=data.get("source_specimen_ids", []),
            source_subject_ids=data.get("source_subject_ids", []),
            metadata=data.get("metadata", {}),
            overlay=data.get("overlay", {}),
        )

    # --- Serialization to stored spec ---
    # MatrixSpec を永続化しやすい dict 表現へ変換する。
    # schema key の整形と None 除去のみを行い、意味解釈はここで増やさない。
    def to_dict(self) -> Dict[str, Any]:
        # [v0.6.0 C-09]
        # serializer 自体は汎用のままでよい。
        # ここで特定 backend の都合を吸収し始めるとモデル層が汚れるため、
        # backend 依存の整形は writer / adapter 側で扱う方針を維持したい。
        d = asdict(self)
        d["$schema_name"] = d.pop("schema_name")
        d["$schema_version"] = d.pop("schema_version")
        # Remove None values
        return {k: v for k, v in d.items() if v is not None}

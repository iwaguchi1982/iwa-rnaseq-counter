from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any

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
    # ただし呼び出し側で tx2gene_path や backend 固有 annotation 資源を
    # 代入しないことを今後も明確に守る必要がある。
    # feature_annotation_path は reporter 表示用 annotation の参照に限定する。
    feature_annotation_path: Optional[str] = None

    """
    レポーター対応アノテーションファイル（v0.5.0契約）への参照。
    - 有効なアノテーションが存在する場合：feature_annotation.tsvへの実際の絶対パスまたは相対パス。
    - アノテーションが存在しない場合：null（None）または空文字列 ""。
    - [重要]：このフィールドの代わりにtx2geneパスを使用しない。
    """
    source_assay_ids: List[str] = field(default_factory=list)
    source_specimen_ids: List[str] = field(default_factory=list)
    source_subject_ids: List[str] = field(default_factory=list)

    # [v0.6.0 C-09]
    # metadata は拡張余地として必要だが、
    # salmon_index / tx2gene_path / backend 固有 key を無秩序に積む場所にはしたくない。
    # v0.6.0 では「共通 metadata」と「backend 由来 metadata」の整理方針を決めたい。
    metadata: Dict[str, Any] = field(default_factory=dict)
    overlay: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MatrixSpec":

        # [v0.6.0 C-09]
        # ここは受け皿としては汎用のままでよい。
        # 重要なのは upstream writer が backend 固有 key をどう入れるか。
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

    def to_dict(self) -> Dict[str, Any]:
        # [v0.6.0 C-09]
        # serializer 自体は汎用のままでよい。
        # ここで特定 backend の都合を吸収し始めるとモデル層が汚れるので注意。
        d = asdict(self)
        d["$schema_name"] = d.pop("schema_name")
        d["$schema_version"] = d.pop("schema_version")
        # Remove None values
        return {k: v for k, v in d.items() if v is not None}

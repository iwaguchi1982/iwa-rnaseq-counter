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
    feature_annotation_path: str
    source_assay_ids: List[str] = field(default_factory=list)
    source_specimen_ids: List[str] = field(default_factory=list)
    source_subject_ids: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    overlay: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MatrixSpec":
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
            feature_annotation_path=data.get("feature_annotation_path", ""),
            source_assay_ids=data.get("source_assay_ids", []),
            source_specimen_ids=data.get("source_specimen_ids", []),
            source_subject_ids=data.get("source_subject_ids", []),
            metadata=data.get("metadata", {}),
            overlay=data.get("overlay", {}),
        )

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["$schema_name"] = d.pop("schema_name")
        d["$schema_version"] = d.pop("schema_version")
        # Remove None values
        return {k: v for k, v in d.items() if v is not None}

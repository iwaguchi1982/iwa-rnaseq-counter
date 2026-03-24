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

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["$schema_name"] = d.pop("schema_name")
        d["$schema_version"] = d.pop("schema_version")
        # Remove None values
        return {k: v for k, v in d.items() if v is not None}

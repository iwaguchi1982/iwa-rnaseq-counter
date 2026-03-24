from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any

@dataclass
class ExecutionRunSpec:
    schema_name: str
    schema_version: str
    run_id: str
    app_name: str
    app_version: str
    started_at: str
    input_refs: List[str] = field(default_factory=list)
    output_refs: List[str] = field(default_factory=list)
    parameters: Dict[str, Any] = field(default_factory=dict)
    execution_backend: Optional[str] = None
    finished_at: Optional[str] = None

    status: str = "pending"
    log_path: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)
    overlay: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["$schema_name"] = d.pop("schema_name")
        d["$schema_version"] = d.pop("schema_version")
        return {k: v for k, v in d.items() if v is not None}

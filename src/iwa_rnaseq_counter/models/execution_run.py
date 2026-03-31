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

    # [v0.6.0 C-09]
    # metadata は将来の拡張先として有力。
    # ただし parameters の肥大化対策として使う場合でも、
    # backend 固有情報の整理ルールを決めてから移すべき。
    metadata: Dict[str, Any] = field(default_factory=dict)
    overlay: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        # [v0.6.0 C-09]
        # serializer は単純維持でよい。
        # backend 差分吸収をここで始めるのではなく、writer / builder 側で整理する
        d = asdict(self)
        d["$schema_name"] = d.pop("schema_name")
        d["$schema_version"] = d.pop("schema_version")
        return {k: v for k, v in d.items() if v is not None}

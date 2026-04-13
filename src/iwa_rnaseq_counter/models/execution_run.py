from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any

from .execution_step import ExecutionStepRecord

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

    preprocessing_steps: Optional[Dict[str, ExecutionStepRecord]] = None

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["$schema_name"] = d.pop("schema_name")
        d["$schema_version"] = d.pop("schema_version")
        
        # Ensure deep dict conversion for preprocessing_steps if provided
        if "preprocessing_steps" in d and d["preprocessing_steps"]:
            # Need to re-build as the default asdict behavior isn't enough if we filter Nones
            # We explicitly use the custom to_dict of the ExecutionStepRecord
            steps_dict = {}
            for k, v in self.preprocessing_steps.items():
                steps_dict[k] = v.to_dict()
            d["preprocessing_steps"] = steps_dict
            
        return {k: v for k, v in d.items() if v is not None}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ExecutionRunSpec":
        steps_raw = data.get("preprocessing_steps")
        preprocessing_steps = None
        if steps_raw and isinstance(steps_raw, dict):
            preprocessing_steps = {
                k: ExecutionStepRecord.from_dict(v) if isinstance(v, dict) else v
                for k, v in steps_raw.items()
            }
            
        return cls(
            schema_name=data.get("$schema_name", data.get("schema_name", "")),
            schema_version=data.get("$schema_version", data.get("schema_version", "")),
            run_id=data.get("run_id", ""),
            app_name=data.get("app_name", ""),
            app_version=data.get("app_version", ""),
            started_at=data.get("started_at", ""),
            input_refs=data.get("input_refs", []) or [],
            output_refs=data.get("output_refs", []) or [],
            parameters=data.get("parameters", {}) or {},
            execution_backend=data.get("execution_backend"),
            finished_at=data.get("finished_at"),
            status=data.get("status", "pending"),
            log_path=data.get("log_path", ""),
            metadata=data.get("metadata", {}) or {},
            overlay=data.get("overlay", {}) or {},
            preprocessing_steps=preprocessing_steps
        )

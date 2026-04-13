from dataclasses import dataclass, field, asdict
from typing import Optional, Dict, Any

@dataclass
class ExecutionStepRecord:
    enabled: bool
    status: str
    tool_name: Optional[str] = None
    warning_count: Optional[int] = None
    error_summary: Optional[str] = None
    log_ref: Optional[str] = None
    report_ref: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert the record to a dictionary, filtering out None values."""
        d = asdict(self)
        return {k: v for k, v in d.items() if v is not None}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ExecutionStepRecord":
        return cls(
            enabled=data.get("enabled", False),
            status=data.get("status", "not_run"),
            tool_name=data.get("tool_name"),
            warning_count=data.get("warning_count"),
            error_summary=data.get("error_summary"),
            log_ref=data.get("log_ref"),
            report_ref=data.get("report_ref"),
            metadata=data.get("metadata", {}) or {}
        )

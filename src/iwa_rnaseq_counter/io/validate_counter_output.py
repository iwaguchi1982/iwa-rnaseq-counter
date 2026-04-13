import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Optional
import pandas as pd

from ..io.read_matrix_spec import read_matrix_spec
from ..models.execution_run import ExecutionRunSpec
from ..models.matrix import MatrixSpec

logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class CounterContractValidationIssue:
    level: str  # "error" | "warning"
    code: str
    message: str
    artifact_name: Optional[str] = None
    path: Optional[str] = None

@dataclass
class CounterContractValidationResult:
    outdir: Path
    is_valid: bool
    error_count: int
    warning_count: int
    issues: List[CounterContractValidationIssue] = field(default_factory=list)
    
    def summary(self) -> str:
        status = "VALID" if self.is_valid else "INVALID"
        return f"Validation {status}: {self.error_count} errors, {self.warning_count} warnings."

def _append_issue(
    issues: List[CounterContractValidationIssue],
    level: str,
    code: str,
    message: str,
    artifact_name: Optional[str] = None,
    path: Optional[str] = None
):
    issues.append(CounterContractValidationIssue(level, code, message, artifact_name, path))

def _read_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def validate_counter_output(outdir: str | Path) -> CounterContractValidationResult:
    outdir = Path(outdir).resolve()
    issues: List[CounterContractValidationIssue] = []
    
    # 1. Basic Structure Check
    specs_dir = outdir / "specs"
    results_dir = outdir / "results"
    logs_dir = outdir / "logs"
    
    for d in [specs_dir, results_dir, logs_dir]:
        if not d.exists() or not d.is_dir():
            _append_issue(issues, "error", "missing_directory", f"Required directory missing: {d.name}", path=str(d))

    # 2. Key Artifact Existence
    exec_spec_path = specs_dir / "execution-run.spec.json"
    matrix_spec_path = specs_dir / "matrix.spec.json"
    run_log_path = logs_dir / "run.log"
    
    # ExecutionRunSpec is mandatory in both success/failed
    if not exec_spec_path.exists():
        _append_issue(issues, "error", "missing_execution_spec", "execution-run.spec.json is missing", artifact_name="execution_run_spec")
    
    if not run_log_path.exists():
        _append_issue(issues, "warning", "missing_run_log", "run.log is missing", artifact_name="run_log")

    # 3. Spec Content and Consistency
    exec_run_spec = None
    if exec_spec_path.exists():
        try:
            data = _read_json(exec_spec_path)
            exec_run_spec = ExecutionRunSpec.from_dict(data) # Assuming we added from_dict to model or helper
        except Exception as e:
            # Fallback if from_dict doesn't exist yet or fails
            try:
                data = _read_json(exec_spec_path)
                # Minimal mock check
                if data.get("schema_name") != "ExecutionRunSpec":
                    _append_issue(issues, "error", "invalid_schema", f"Expected ExecutionRunSpec, got {data.get('schema_name')}", artifact_name="execution_run_spec")
            except Exception as e2:
                _append_issue(issues, "error", "invalid_json", f"Failed to load ExecutionRunSpec: {e2}", artifact_name="execution_run_spec")

    status = getattr(exec_run_spec, "status", "unknown") if exec_run_spec else "unknown"
    
    if status in ("completed", "completed_with_errors"):
        # Success expectations
        if not matrix_spec_path.exists():
            _append_issue(issues, "error", "missing_matrix_spec", "MatrixSpec is missing despite 'completed' status", artifact_name="matrix_spec")
        
        # Check matrix file existence (could be .tsv or .csv)
        matrix_file = results_dir / "gene_numreads.tsv"
        if not matrix_file.exists():
            matrix_file = results_dir / "gene_numreads.csv"
            
        if not matrix_file.exists():
            _append_issue(issues, "error", "missing_matrix_file", "Count matrix file is missing in results/", artifact_name="matrix_file")
        else:
            # Consistency check between ExecutionRunSpec and MatrixSpec
            if matrix_spec_path.exists():
                try:
                    ms = read_matrix_spec(matrix_spec_path)
                    if exec_run_spec and ms.matrix_id not in exec_run_spec.output_refs:
                        _append_issue(issues, "warning", "id_mismatch", f"Matrix ID {ms.matrix_id} not found in ExecutionRunSpec.output_refs", artifact_name="matrix_spec")
                except Exception as e:
                    _append_issue(issues, "error", "invalid_matrix_spec", f"Failed to load MatrixSpec: {e}", artifact_name="matrix_spec")

    elif status == "failed":
        # Failure expectations
        meta = getattr(exec_run_spec, "metadata", {})
        if not meta.get("failure_summary"):
            _append_issue(issues, "warning", "missing_failure_summary", "Failure summary missing in metadata", artifact_name="execution_run_spec")
        if not meta.get("failure_stage"):
            _append_issue(issues, "warning", "missing_failure_stage", "Failure stage missing in metadata", artifact_name="execution_run_spec")

    error_count = sum(1 for i in issues if i.level == "error")
    warning_count = sum(1 for i in issues if i.level == "warning")
    
    return CounterContractValidationResult(
        outdir=outdir,
        is_valid=(error_count == 0),
        error_count=error_count,
        warning_count=warning_count,
        issues=issues
    )

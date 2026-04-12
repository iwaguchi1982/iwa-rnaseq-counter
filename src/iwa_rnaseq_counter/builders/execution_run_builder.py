from datetime import datetime, timezone
from typing import Any, Dict, List, Optional
from iwa_rnaseq_counter.models.execution_run import ExecutionRunSpec
from iwa_rnaseq_counter.models.execution_step import ExecutionStepRecord

def build_execution_run_spec_for_success(
    run_id: str,
    app_version: str,
    started_at: str,
    input_refs: List[str],
    output_refs: List[str],
    parameters: Dict[str, Any],
    execution_backend: str,
    log_path: str,
    status: str = "completed",
    preprocessing_steps: Optional[Dict[str, ExecutionStepRecord]] = None,
    metadata: Optional[Dict[str, Any]] = None
) -> ExecutionRunSpec:
    """
    Builds a standard ExecutionRunSpec for a successful run.
    """
    final_metadata = metadata.copy() if metadata else {}
    final_metadata.setdefault("status_detail", "completed successfully")
    final_metadata.setdefault("output_generated", True)
    
    return ExecutionRunSpec(
        schema_name="ExecutionRunSpec",
        schema_version="0.1.0",
        run_id=run_id,
        app_name="iwa_rnaseq_counter",
        app_version=app_version,
        started_at=started_at,
        input_refs=input_refs,
        output_refs=output_refs,
        parameters=parameters,
        execution_backend=execution_backend,
        finished_at=datetime.now(timezone.utc).astimezone().isoformat(),
        status=status,
        log_path=log_path,
        preprocessing_steps=preprocessing_steps,
        metadata=final_metadata
    )

def build_execution_run_spec_for_failure(
    run_id: str,
    app_version: str,
    started_at: str,
    input_refs: List[str],
    parameters: Dict[str, Any],
    execution_backend: str,
    log_path: str,
    failure_stage: str,
    failure_summary: str,
    error_messages: Optional[List[str]] = None,
    preprocessing_steps: Optional[Dict[str, ExecutionStepRecord]] = None,
    metadata: Optional[Dict[str, Any]] = None
) -> ExecutionRunSpec:
    """
    Builds a structured ExecutionRunSpec for a failed run.
    """
    final_metadata = metadata.copy() if metadata else {}
    final_metadata.setdefault("failure_stage", failure_stage)
    final_metadata.setdefault("failure_summary", failure_summary)
    final_metadata.setdefault("error_messages", error_messages or [failure_summary])
    final_metadata.setdefault("output_generated", False)
    
    return ExecutionRunSpec(
        schema_name="ExecutionRunSpec",
        schema_version="0.1.0",
        run_id=run_id,
        app_name="iwa_rnaseq_counter",
        app_version=app_version,
        started_at=started_at,
        input_refs=input_refs,
        output_refs=[],  # No outputs on failure
        parameters=parameters,
        execution_backend=execution_backend,
        finished_at=datetime.now(timezone.utc).astimezone().isoformat(),
        status="failed",
        log_path=log_path,
        preprocessing_steps=preprocessing_steps,
        metadata=final_metadata
    )

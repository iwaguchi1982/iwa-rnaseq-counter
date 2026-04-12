import json
from pathlib import Path
from datetime import datetime, timezone
from iwa_rnaseq_counter.builders.execution_run_builder import build_execution_run_spec_for_success, build_execution_run_spec_for_failure
from iwa_rnaseq_counter.models.execution_step import ExecutionStepRecord

def generate_samples():
    dev_docs = Path("/home/manager/iwa_bio_analysis_orchestra/dev_docs")
    dev_docs.mkdir(parents=True, exist_ok=True)
    
    # Success Sample
    success_spec = build_execution_run_spec_for_success(
        run_id="RUN_20240412_SUCCESS",
        app_version="0.3.5",
        started_at="2024-04-12T10:00:00+09:00",
        input_refs=["ASSAY_001"],
        output_refs=["MAT_001"],
        parameters={
            "quantifier": "salmon",
            "threads": 4,
            "profile": "local"
        },
        execution_backend="local",
        log_path="/path/to/logs/run.log",
        preprocessing_steps={
            "qc": ExecutionStepRecord(enabled=True, status="completed", warning_count=0),
            "trimming": ExecutionStepRecord(enabled=False, status="not_run")
        }
    )
    
    with open(dev_docs / "execution_run_success_sample.json", "w") as f:
        json.dump(success_spec.to_dict(), f, indent=2)
        
    # Failure Sample (Quantifier stage)
    failed_spec = build_execution_run_spec_for_failure(
        run_id="RUN_20240412_FAILED",
        app_version="0.3.5",
        started_at="2024-04-12T11:00:00+09:00",
        input_refs=["ASSAY_002"],
        parameters={
            "quantifier": "star",
            "threads": 8,
            "profile": "hpc"
        },
        execution_backend="hpc",
        log_path="/path/to/logs/run.log",
        failure_stage="quantifier",
        failure_summary="STAR binary not found in PATH",
        error_messages=["Executable 'STAR' not found", "Search path: /usr/bin:/bin"],
        preprocessing_steps={
            "qc": ExecutionStepRecord(enabled=True, status="completed"),
            "trimming": ExecutionStepRecord(enabled=True, status="failed", error_summary="Disk full")
        }
    )
    
    with open(dev_docs / "execution_run_failed_sample.json", "w") as f:
        json.dump(failed_spec.to_dict(), f, indent=2)

if __name__ == "__main__":
    generate_samples()

import pytest
from iwa_rnaseq_counter.models.execution_step import ExecutionStepRecord

def test_execution_step_record_warning_to_dict():
    record = ExecutionStepRecord(
        enabled=True,
        status="warning",
        tool_name="fastqc",
        warning_count=3,
        error_summary=None,
        log_ref="file:///path/to/log.txt",
        metadata={"version": "0.12.1"}
    )
    
    d = record.to_dict()
    assert d["enabled"] is True
    assert d["status"] == "warning"
    assert d["tool_name"] == "fastqc"
    assert d["warning_count"] == 3
    assert d["log_ref"] == "file:///path/to/log.txt"
    assert d["metadata"] == {"version": "0.12.1"}
    
    # Check that Nones are omitted
    assert "error_summary" not in d
    assert "report_ref" not in d

def test_execution_step_record_failed_to_dict():
    record = ExecutionStepRecord(
        enabled=True,
        status="failed",
        tool_name="trimmomatic",
        error_summary="Trimming failed due to invalid adapter sequences",
        log_ref="file:///path/to/trimming_error.log"
    )
    
    d = record.to_dict()
    assert d["enabled"] is True
    assert d["status"] == "failed"
    assert d["tool_name"] == "trimmomatic"
    assert d["error_summary"] == "Trimming failed due to invalid adapter sequences"
    assert d["log_ref"] == "file:///path/to/trimming_error.log"
    
    assert "warning_count" not in d
    assert "report_ref" not in d

def test_execution_step_record_default_factory():
    record = ExecutionStepRecord(
        enabled=False,
        status="not_run"
    )
    assert record.enabled is False
    assert record.status == "not_run"
    assert record.metadata == {}
    
    d = record.to_dict()
    assert "tool_name" not in d
    assert "warning_count" not in d

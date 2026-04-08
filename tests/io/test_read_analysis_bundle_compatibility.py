import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from iwa_rnaseq_counter.io.read_analysis_bundle import (
    evaluate_analysis_bundle_contract,
)

def test_evaluate_analysis_bundle_contract_accepts_major_1():
    manifest = {
        "contract_name": "analysis_bundle",
        "contract_version": "1.2.0",
        "bundle_kind": "rna_seq_analysis_bundle",
        "producer": "iwa_rnaseq_counter",
        "producer_version": "0.10.1",
    }

    info = evaluate_analysis_bundle_contract(manifest)

    assert info.is_supported is True
    assert info.compatibility_status == "supported"

def test_evaluate_analysis_bundle_contract_rejects_missing_version():
    manifest = {
        "contract_name": "analysis_bundle",
        "bundle_kind": "rna_seq_analysis_bundle",
    }

    info = evaluate_analysis_bundle_contract(manifest)

    assert info.is_supported is False
    assert info.compatibility_status == "missing_contract_version"

def test_evaluate_analysis_bundle_contract_rejects_bundle_kind():
    manifest = {
        "contract_name": "analysis_bundle",
        "contract_version": "1.0.0",
        "bundle_kind": "other_bundle",
    }

    info = evaluate_analysis_bundle_contract(manifest)

    assert info.is_supported is False
    assert info.compatibility_status == "unsupported_bundle_kind"

def test_evaluate_analysis_bundle_contract_rejects_major_2():
    manifest = {
        "contract_name": "analysis_bundle",
        "contract_version": "2.0.0",
        "bundle_kind": "rna_seq_analysis_bundle",
    }

    info = evaluate_analysis_bundle_contract(manifest)

    assert info.is_supported is False
    assert info.compatibility_status == "unsupported_contract_major"

if __name__ == "__main__":
    test_evaluate_analysis_bundle_contract_accepts_major_1()
    test_evaluate_analysis_bundle_contract_rejects_missing_version()
    test_evaluate_analysis_bundle_contract_rejects_bundle_kind()
    test_evaluate_analysis_bundle_contract_rejects_major_2()
    print("test_read_analysis_bundle_compatibility: ALL PASSED")

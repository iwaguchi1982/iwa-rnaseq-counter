from __future__ import annotations
from typing import Any, TypedDict


class QCResult(TypedDict, total=False):
    sample_id: str
    status: str
    icon: str
    alerts: list[str]
    metric_mode: str
    mapping_rate: float | None
    decoy_rate: float | None
    filter_rate: float | None
    processed_reads: int | None
    mapped_reads: int | None


def evaluate_sample_qc(
    sample_id: str,
    num_processed: int,
    num_mapped: int,
    num_decoy: int,
    num_filter: int
) -> QCResult:
    """
    Salmon の実行統計に基づき、生物学的な品質を判定する。
    """
    total = max(num_processed, 1)
    mapping_rate = (num_mapped / total) * 100
    decoy_rate = (num_decoy / total) * 100
    filter_rate = (num_filter / total) * 100

    alerts = []
    status = "OK"
    icon = "✅"

    if mapping_rate < 1.0:
        status = "Critical"
        icon = "🔴"
        alerts.append(f"Low mapping rate ({mapping_rate:.2f}%)")
    elif mapping_rate < 10.0:
        status = "Warning"
        icon = "🟡"
        alerts.append(f"Sub-optimal mapping rate ({mapping_rate:.2f}%)")

    if decoy_rate > 20.0:
        if status == "OK":
            status = "Warning"
            icon = "🟡"
        alerts.append(f"High decoy discard ({decoy_rate:.1f}%)")

    if filter_rate > 30.0:
        if status == "OK":
            status = "Warning"
            icon = "🟡"
        alerts.append(f"High low-score filter ({filter_rate:.1f}%)")

    return {
        "sample_id": sample_id,
        "status": status,
        "icon": icon,
        "alerts": alerts,
        "metric_mode": "mapping_stats",
        "mapping_rate": mapping_rate,
        "decoy_rate": decoy_rate,
        "filter_rate": filter_rate,
        "processed_reads": num_processed,
        "mapped_reads": num_mapped,
    }


def output_has_mapping_metrics(output: dict[str, Any]) -> bool:
    return any(
        output.get(k) is not None
        for k in ["num_processed", "num_mapped", "mapping_rate"]
    )


def summarize_output_qc(output: dict[str, Any]) -> QCResult:
    """
    backend 非依存の UI 用 summary。
    mapping 統計がある backend は evaluate_sample_qc を使い、
    ない backend は artifact/status ベースの軽量表示に落とす。
    """
    sample_id = str(output.get("sample_id", "unknown"))

    if output_has_mapping_metrics(output):
        return evaluate_sample_qc(
            sample_id=sample_id,
            num_processed=int(output.get("num_processed") or 0),
            num_mapped=int(output.get("num_mapped") or 0),
            num_decoy=int(output.get("num_decoy") or 0),
            num_filter=int(output.get("num_filter") or 0),
        )

    if output.get("is_success"):
        alerts: list[str] = []
        if output.get("gene_counts_path"):
            alerts.append("Gene-level output available")
        if output.get("transcript_quant_path") or output.get("quant_path"):
            alerts.append("Transcript-level output available")

        return {
            "sample_id": sample_id,
            "status": "OK",
            "icon": "✅",
            "alerts": alerts,
            "metric_mode": "artifact_only",
            "mapping_rate": None,
            "decoy_rate": None,
            "filter_rate": None,
            "processed_reads": None,
            "mapped_reads": None,
        }

    return {
        "sample_id": sample_id,
        "status": "Critical",
        "icon": "❌",
        "alerts": [output.get("error_reason", "Execution failed")],
        "metric_mode": "error",
        "mapping_rate": None,
        "decoy_rate": None,
        "filter_rate": None,
        "processed_reads": None,
        "mapped_reads": None,
    }

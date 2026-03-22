from __future__ import annotations
from typing import TypedDict


class QCResult(TypedDict):
    sample_id: str
    status: str  # "OK", "Warning", "Critical"
    icon: str    # "✅", "🟡", "🔴"
    alerts: list[str]
    mapping_rate: float
    decoy_rate: float
    filter_rate: float


def evaluate_sample_qc(
    sample_id: str,
    num_processed: int,
    num_mapped: int,
    num_decoy: int,
    num_filter: int
) -> QCResult:
    """
    Salmon の実行統計に基づき、生物学的な品質を判定する。
    
    判定基準:
    - mapping_rate < 1% -> critical
    - mapping_rate < 10% -> warning
    - decoy_rate > 20% -> warning
    - filter_rate > 30% -> warning
    """
    total = max(num_processed, 1)
    mapping_rate = (num_mapped / total) * 100
    decoy_rate = (num_decoy / total) * 100
    filter_rate = (num_filter / total) * 100
    
    alerts = []
    status = "OK"
    icon = "✅"
    
    # Mapping Rate Check
    if mapping_rate < 1.0:
        status = "Critical"
        icon = "🔴"
        alerts.append(f"Low mapping rate ({mapping_rate:.2f}%)")
    elif mapping_rate < 10.0:
        if status != "Critical":
            status = "Warning"
            icon = "🟡"
        alerts.append(f"Sub-optimal mapping rate ({mapping_rate:.2f}%)")
        
    # Decoy Rate Check
    if decoy_rate > 20.0:
        if status == "OK":
            status = "Warning"
            icon = "🟡"
        alerts.append(f"High decoy discard ({decoy_rate:.1f}%)")
        
    # Filter Rate Check
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
        "mapping_rate": mapping_rate,
        "decoy_rate": decoy_rate,
        "filter_rate": filter_rate
    }

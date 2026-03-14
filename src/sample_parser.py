from __future__ import annotations

from pathlib import Path

import pandas as pd


def group_fastq_by_sample(fastq_files: list[str]) -> dict:
    groups: dict[str, dict] = {}
    for path_str in fastq_files:
        path = Path(path_str)
        sample_id = _extract_sample_id(path.name)
        lane = _extract_lane(path.name) or "single_lane"
        read = _extract_read(path.name)

        groups.setdefault(sample_id, {"files": [], "read_pairs": {}})
        groups[sample_id]["files"].append(path_str)
        groups[sample_id]["read_pairs"].setdefault(lane, {})
        if read:
            groups[sample_id]["read_pairs"][lane][read] = path_str
        else:
            groups[sample_id]["read_pairs"][lane]["SE"] = path_str
    return groups


def infer_sample_layout(sample_groups: dict) -> dict:
    layout_map: dict[str, str] = {}
    for sample_id, info in sample_groups.items():
        lane_maps = info.get("read_pairs", {})
        is_paired = any("R1" in reads or "R2" in reads for reads in lane_maps.values())
        layout_map[sample_id] = "paired-end" if is_paired else "single-end"
    return layout_map


def detect_lane_groups(sample_groups: dict) -> dict:
    lane_map: dict[str, dict] = {}
    for sample_id, info in sample_groups.items():
        lanes = sorted(info.get("read_pairs", {}).keys())
        lane_map[sample_id] = {
            "lane_count": len(lanes),
            "lanes": lanes,
        }
    return lane_map


def build_sample_table(sample_groups: dict, layout_map: dict, lane_map: dict) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for sample_id, info in sample_groups.items():
        layout = layout_map.get(sample_id, "unknown")
        lane_info = lane_map.get(sample_id, {"lane_count": 0, "lanes": []})
        files = sorted(info.get("files", []))
        r1_paths = [v.get("R1") for v in info.get("read_pairs", {}).values() if v.get("R1")]
        r2_paths = [v.get("R2") for v in info.get("read_pairs", {}).values() if v.get("R2")]
        status = "ok"
        notes = ""
        if layout == "paired-end" and len(r1_paths) != len(r2_paths):
            status = "error"
            notes = "R1 / R2 count mismatch"
        rows.append(
            {
                "sample_id": sample_id,
                "layout_predicted": layout,
                "layout_final": layout,
                "lane_count": lane_info["lane_count"],
                "file_count": len(files),
                "status": status,
                "r1_paths": r1_paths,
                "r2_paths": r2_paths,
                "all_paths": files,
                "notes": notes,
            }
        )
    return pd.DataFrame(rows)


def apply_sample_table_edits(sample_df: pd.DataFrame, edited_sample_df: pd.DataFrame | None) -> pd.DataFrame:
    if edited_sample_df is None:
        return sample_df
    return edited_sample_df.copy()


def _extract_sample_id(filename: str) -> str:
    name = filename
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        name = name.replace(suffix, "")
    for token in ("_R1", "_R2", "_1", "_2"):
        name = name.replace(token, "")
    for i in range(1, 9):
        name = name.replace(f"_L00{i}", "")
    return name.strip("_-")


def _extract_lane(filename: str) -> str | None:
    for i in range(1, 9):
        token = f"L00{i}"
        if token in filename:
            return token
    return None


def _extract_read(filename: str) -> str | None:
    if "_R1" in filename or "_1" in filename:
        return "R1"
    if "_R2" in filename or "_2" in filename:
        return "R2"
    return None

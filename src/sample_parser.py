import re
import pandas as pd
from pathlib import Path

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
        notes_list = []
        if layout == "paired-end" and len(r1_paths) != len(r2_paths):
            status = "error"
            notes_list.append("R1 / R2 count mismatch")
        
        if lane_info["lane_count"] > 1:
            lanes_str = ", ".join(lane_info["lanes"])
            notes_list.append(f"{lane_info['lane_count']} lanes merged ({lanes_str})")
            
        rows.append(
            {
                "sample_id": sample_id,
                "group": "",
                "layout_predicted": layout,
                "layout_final": layout,
                "lane_count": lane_info["lane_count"],
                "file_count": len(files),
                "status": status,
                "r1_paths": r1_paths,
                "r2_paths": r2_paths,
                "all_paths": files,
                "notes": "; ".join(notes_list),
            }
        )
    return pd.DataFrame(rows)


def parse_sample_sheet(csv_path: str, input_dir: str) -> pd.DataFrame:
    """
    CSV からサンプル表を作成する。
    仕様: sample_id, r1_path, r2_path, layout, group
    r1_path, r2_path は input_dir からの相対パス、または絶対パス。
    """
    df_csv = pd.read_csv(csv_path)
    # 必要列の正規化
    col_map = {c.lower(): c for c in df_csv.columns}
    
    rows = []
    for _, row in df_csv.iterrows():
        sid = str(row.get(col_map.get("sample_id"), "")).strip()
        r1 = str(row.get(col_map.get("r1_path"), "")).strip()
        r2 = str(row.get(col_map.get("r2_path"), "")).strip()
        layout = str(row.get(col_map.get("layout"), "paired-end")).strip().lower()
        group = str(row.get(col_map.get("group"), "")).strip()
        
        # パス解決
        def resolve(p):
            if not p or p == "nan": return None
            pth = Path(p)
            if not pth.is_absolute():
                pth = Path(input_dir) / pth
            return str(pth)

        r1_resolved = resolve(r1)
        r2_resolved = resolve(r2)
        
        all_p = [p for p in [r1_resolved, r2_resolved] if p]
        
        # layout 文字列の統一
        if layout in ("pe", "paired"): layout = "paired-end"
        if layout in ("se", "single"): layout = "single-end"
        
        rows.append({
            "sample_id": sid,
            "group": group,
            "layout_predicted": layout,
            "layout_final": layout,
            "lane_count": 1,
            "file_count": len(all_p),
            "status": "ok",
            "r1_paths": [r1_resolved] if r1_resolved else [],
            "r2_paths": [r2_resolved] if r2_resolved else [],
            "all_paths": all_p,
            "notes": "Loaded from CSV",
        })
    
    df = pd.DataFrame(rows)
    # apply_sample_table_edits と同様のバリデーションを適用
    df = df.apply(_recalculate_sample_row_status, axis=1)
    df = _mark_duplicate_sample_ids(df)
    return df


def apply_sample_table_edits(sample_df: pd.DataFrame, edited_sample_df: pd.DataFrame | None) -> pd.DataFrame:
    if edited_sample_df is None:
        return sample_df

    df = edited_sample_df.copy()
    
    # 1. 基礎的なクリーンアップ
    df["sample_id"] = df["sample_id"].astype(str).str.strip()
    if "group" in df.columns:
        df["group"] = df["group"].astype(str).str.strip()
    df["layout_final"] = df["layout_final"].astype(str).str.strip().str.lower()
    
    # 2. 状態の再計算
    df = df.apply(_recalculate_sample_row_status, axis=1)
    
    # 3. 重複チェック
    df = _mark_duplicate_sample_ids(df)
    
    return df


def _recalculate_sample_row_status(row: pd.Series) -> pd.Series:
    sid = str(row["sample_id"])
    layout = row["layout_final"]
    r1 = row["r1_paths"]
    r2 = row["r2_paths"]
    all_p = row["all_paths"]
    
    notes = []
    status = "ok"
    
    if not sid or sid == "nan" or sid == "None":
        status = "error"
        notes.append("Empty sample ID")
    
    if layout not in ("paired-end", "single-end"):
        status = "error"
        notes.append(f"Invalid layout: {layout}")
    
    # パス存在チェック
    missing = []
    for p in all_p:
        if not p or not Path(p).exists():
            missing.append(Path(p).name if p else "None")
    
    if missing:
        status = "error"
        notes.append(f"File not found: {', '.join(missing)}")

    if layout == "paired-end":
        if not r1 or not r2:
            status = "error"
            notes.append("R1 or R2 missing for PE")
        elif len(r1) != len(r2):
            status = "error"
            notes.append("R1/R2 count mismatch")
    elif layout == "single-end":
        if not all_p and not r1:
            status = "error"
            notes.append("No FASTQ for SE")
            
    # 元々の notes も一部引き継ぐ (lane 統合情報など)
    orig_notes = str(row.get("notes", ""))
    if "merged" in orig_notes:
        notes.append(orig_notes.split(";")[0]) # 最初の一つを保持
        
    row["status"] = status
    row["notes"] = "; ".join(sorted(set(notes)))
    return row


def _mark_duplicate_sample_ids(df: pd.DataFrame) -> pd.DataFrame:
    ids = df["sample_id"]
    is_dup = ids.duplicated(keep=False) & (ids != "") & (ids != "nan")
    for idx in df.index[is_dup]:
        df.at[idx, "status"] = "error"
        notes = str(df.at[idx, "notes"])
        if "Duplicate" not in notes:
            df.at[idx, "notes"] = (notes + "; Duplicate ID").strip("; ")
    return df


def _extract_sample_id(filename: str) -> str:
    # 拡張子を除去
    name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", filename, flags=re.IGNORECASE)
    # リード識別子を除去 (_R1, _R2, _1, _2)
    name = re.sub(r"(_R[12]|_[12])$", "", name)
    # レーン識別子を除去 (_L001 ~ _L008)
    name = re.sub(r"_L00[1-8]", "", name)
    return name.strip("_-")


def _extract_lane(filename: str) -> str | None:
    match = re.search(r"L00[1-8]", filename)
    return match.group(0) if match else None


def _extract_read(filename: str) -> str | None:
    if re.search(r"(_R1|[^0-9]1)\.(fastq|fq)", filename, re.IGNORECASE):
        return "R1"
    if re.search(r"(_R2|[^0-9]2)\.(fastq|fq)", filename, re.IGNORECASE):
        return "R2"
    return None

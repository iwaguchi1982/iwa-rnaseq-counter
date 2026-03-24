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


METADATA_COLUMNS = [
    "group",
    "condition",
    "replicate",
    "batch",
    "pair_id",
    "note",
    "display_name",
    "color",
    "exclude"
]

STANDARD_SAMPLE_COLUMNS = [
    "sample_id",
    "group",
    "condition",
    "replicate",
    "batch",
    "pair_id",
    "note",
    "display_name",
    "color",
    "exclude",
    "layout_predicted",
    "layout_final",
    "lane_count",
    "file_count",
    "status",
    "notes",
    "input_source",
    "r1_paths",
    "r2_paths",
    "all_paths",
]


def normalize_sample_sheet_columns(df: pd.DataFrame) -> pd.DataFrame:
    """列名を標準名に寄せる。"""
    df = df.copy()
    df.columns = [str(c).lower().strip() for c in df.columns]
    return df


def add_missing_metadata_columns(df: pd.DataFrame) -> pd.DataFrame:
    """任意 metadata 列が無ければ空列で補完する。"""
    df = df.copy()
    for col in set(METADATA_COLUMNS) - set(df.columns):
        df[col] = ""
    return df


def normalize_layout_value(layout: str) -> str:
    """layout を single-end / paired-end に正規化する。"""
    l = str(layout).strip().lower()
    if l in ("pe", "paired"): return "paired-end"
    if l in ("se", "single"): return "single-end"
    return l


def normalize_exclude_value(value) -> bool:
    """exclude を bool に正規化する。"""
    val = str(value).strip().lower()
    return val == "true" if val not in ("", "none", "nan") else False


def resolve_sample_sheet_path(path_value: str, input_dir: str | None = None, csv_dir: str | Path | None = None) -> str | None:
    """相対 / 絶対 path を解決する。優先順位: 絶対パス > input_dir相対 > csv_dir相対 > CWD相対"""
    pv = str(path_value).strip()
    if not pv or pv == "nan" or pv == "None":
        return None
    pth = Path(pv)
    if pth.is_absolute():
        return str(pth)
    
    # 1. input_dir があればそこから探す
    if input_dir:
        candidate = Path(input_dir) / pth
        if candidate.exists():
            return str(candidate)
            
    # 2. csv_dir があればそこから探す
    if csv_dir:
        candidate = Path(csv_dir) / pth
        if candidate.exists():
            return str(candidate)
            
    # 3. CWD からの相対パスで探す
    if pth.exists():
        return str(pth)
        
    # 4. どちらでも見つからない場合のフォールバック
    if input_dir:
        # 入力ディレクトリが指定されている場合は、そこからの相対パスであると期待するのが基本
        return str(Path(input_dir) / pth)
    return str(pth)


def standardize_sample_df_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    sample_df の列を STANDARD_SAMPLE_COLUMNS に揃える。
    """
    # 欠落している標準列があれば空で埋める
    for col in STANDARD_SAMPLE_COLUMNS:
        if col not in df.columns:
            if col in ("r1_paths", "r2_paths", "all_paths"):
                df[col] = [[] for _ in range(len(df))]
            elif col == "exclude":
                df[col] = False
            else:
                df[col] = ""
            
    # 存在する列を STANDARD_SAMPLE_COLUMNS の順序で抽出
    final_cols = [c for c in STANDARD_SAMPLE_COLUMNS if c in df.columns]
    return df[final_cols]


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
            
        row = {
            "sample_id": sample_id,
            "layout_predicted": layout,
            "layout_final": layout,
            "lane_count": lane_info["lane_count"],
            "file_count": len(files),
            "status": status,
            "r1_paths": r1_paths,
            "r2_paths": r2_paths,
            "all_paths": files,
            "notes": "; ".join(notes_list),
            "input_source": "auto_detect",
        }
        for col in METADATA_COLUMNS:
            row[col] = ""
            
        rows.append(row)
    
    df = pd.DataFrame(rows)
    return standardize_sample_df_columns(df)


def parse_sample_sheet(csv_path: str, input_dir: str | None = None) -> pd.DataFrame:
    """
    CSV からサンプル表を作成する。
    仕様: sample_id, r1_path, r2_path, layout, group, condition, replicate, batch, pair_id, note
    """
    df_csv = pd.read_csv(csv_path)
    df_csv = normalize_sample_sheet_columns(df_csv)
    df_csv = add_missing_metadata_columns(df_csv)
    
    csv_dir = Path(csv_path).parent
    
    rows = []
    for _, row in df_csv.iterrows():
        sid = str(row.get("sample_id", "")).strip()
        layout = normalize_layout_value(row.get("layout", "paired-end"))
        
        # フォールバック処理
        cond = str(row.get("condition", "")).strip()
        grp = str(row.get("group", "")).strip()
        if not cond or cond == "nan" or cond == "None":
            cond = grp
            
        dname = str(row.get("display_name", "")).strip()
        if not dname or dname == "nan" or dname == "None":
            dname = sid
            
        exclude = normalize_exclude_value(row.get("exclude", ""))
        
        r1_resolved = resolve_sample_sheet_path(row.get("r1_path", ""), input_dir, csv_dir)
        r2_resolved = resolve_sample_sheet_path(row.get("r2_path", ""), input_dir, csv_dir)
        all_p = [p for p in [r1_resolved, r2_resolved] if p is not None]
        
        row_data = {
            "sample_id": sid,
            "layout_predicted": layout,
            "layout_final": layout,
            "lane_count": 1,
            "file_count": len(all_p),
            "status": "ok",
            "r1_paths": [r1_resolved] if r1_resolved else [],
            "r2_paths": [r2_resolved] if r2_resolved else [],
            "all_paths": all_p,
            "notes": "Loaded from CSV",
            "input_source": "sample_sheet",
            "condition": cond,
            "display_name": dname,
            "exclude": exclude,
        }
        
        # その他のメタデータを取得
        for col in METADATA_COLUMNS:
            if col not in ("condition", "display_name", "exclude"):
                val = row.get(col, "")
                row_data[col] = str(val).strip() if pd.notna(val) else ""
            
        rows.append(row_data)
    
    df = pd.DataFrame(rows)
    df = df.apply(_recalculate_sample_row_status, axis=1)
    df = _mark_duplicate_sample_ids(df)
    return standardize_sample_df_columns(df)


def apply_sample_table_edits(sample_df: pd.DataFrame, edited_sample_df: pd.DataFrame | None) -> pd.DataFrame:
    if edited_sample_df is None:
        return sample_df

    df = edited_sample_df.copy()
    
    # 1. 基礎的なクリーンアップ
    df["sample_id"] = df["sample_id"].astype(str).str.strip()
    df["layout_final"] = df["layout_final"].astype(str).str.strip().str.lower()
    
    # メタデータ列のクリーンアップ
    for col in METADATA_COLUMNS:
        if col in df.columns and col != "exclude":
            df[col] = df[col].astype(str).str.strip()
            
    # 値の正規化・フォールバック
    if "exclude" in df.columns:
        df["exclude"] = df["exclude"].apply(normalize_exclude_value)
        
    if "condition" in df.columns and "group" in df.columns:
        df["condition"] = df.apply(
            lambda r: r["group"] if not str(r["condition"]).strip() or str(r["condition"]).strip() == "nan" else str(r["condition"]).strip(), 
            axis=1
        )
        
    if "display_name" in df.columns and "sample_id" in df.columns:
        df["display_name"] = df.apply(
            lambda r: r["sample_id"] if not str(r["display_name"]).strip() or str(r["display_name"]).strip() == "nan" else str(r["display_name"]).strip(), 
            axis=1
        )
    
    # 2. 状態の再計算
    df = df.apply(_recalculate_sample_row_status, axis=1)
    
    # 3. 重複チェック
    df = _mark_duplicate_sample_ids(df)
    
    return standardize_sample_df_columns(df)


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

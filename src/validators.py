from __future__ import annotations

from pathlib import Path
import os
import re
import shutil

import pandas as pd

from src.strandedness import validate_strandedness_selection


def validate_input_directory(input_dir: str) -> dict:
    if not input_dir:
        return _invalid("入力ディレクトリを指定してください。")
    path = Path(input_dir)
    if not path.exists():
        return _invalid("入力ディレクトリが存在しません。")
    if not path.is_dir():
        return _invalid("入力パスがディレクトリではありません。")
    
    # 読み取り可能確認 (Section 19.2 (2))
    if not os.access(path, os.R_OK):
        return _invalid("入力ディレクトリへの読み取り権限がありません。")
        
    return _valid()


def validate_analysis_name(name: str) -> dict:
    if not name or not name.strip():
        return _invalid("解析名を入力してください。")
    
    # ファイル名として使える文字 (英数字, _, -) かつ 空白を含まない
    if not re.match(r"^[a-zA-Z0-9_\-]+$", name):
        return _invalid("解析名には英数字、アンダースコア(_)、ハイフン(-)のみ使用してください。")
    
    return _valid()


def validate_output_directory(output_dir: str) -> dict:
    if not output_dir:
        return _invalid("出力ディレクトリを指定してください。")
    path = Path(output_dir)
    
    # 親ディレクトリの存在確認
    parent = path if path.exists() else path.parent
    if not parent.exists():
        return _invalid("出力先の親ディレクトリが存在しません。")
    
    # 書き込みテスト (書き込み権限の確認)
    # 既存ディレクトリ、または親が存在すれば作成を試みる想定
    try:
        path.mkdir(parents=True, exist_ok=True)
        test_file = path / ".write_test"
        test_file.touch()
        test_file.unlink()
    except Exception as e:
        return _invalid(f"出力ディレクトリへの書き込み権限がありません: {e}")
        
    return _valid()


def validate_salmon_index(salmon_index_path: str) -> dict:
    if not salmon_index_path:
        return _invalid("Salmon index が指定されていません。")
    path = Path(salmon_index_path)
    if not path.exists() or not path.is_dir():
        return _invalid("Salmon index パスが有効ではありません。")
    
    # Salmon インデックスディレクトリには通常 version.json または versionInfo.json が含まれる
    version_files = ["version.json", "versionInfo.json", "ref_sigs.json"]
    if not any((path / f).exists() for f in version_files):
        return _invalid(f"Salmon インデックスとして正しくないか、不完全なディレクトリです: {path}")
    
    # 構造上の必須ファイル (最低限)
    if not (path / "pos.bin").exists() and not (path / "seq.bin").exists():
         return _invalid("Salmon インデックス内の必須ファイル (pos.bin 等) が見つかりません。")
    
    return _valid()


def validate_tx2gene_file(tx2gene_path: str) -> dict:
    if not tx2gene_path:
        return _invalid("tx2gene ファイルが指定されていません。")
    path = Path(tx2gene_path)
    if not path.exists() or not path.is_file():
        return _invalid("tx2gene ファイルが存在しません。")
    
    # 内容の簡易検証 (読めるか、2列以上あるか)
    try:
        sep = "," if tx2gene_path.endswith(".csv") else "\t"
        df_stub = pd.read_csv(tx2gene_path, sep=sep, nrows=2)
        if len(df_stub.columns) < 2:
            return _invalid("tx2gene ファイルには少なくとも2列 (transcript_id, gene_id) 必要です。")
    except Exception as e:
        return _invalid(f"tx2gene ファイルを読み込めません: {e}")

    return _valid()


def validate_sample_structure(sample_df: pd.DataFrame | None) -> dict:
    if sample_df is None or sample_df.empty:
        return _invalid("サンプル情報がありません。")
    
    errors = []
    # 重複チェック
    if sample_df["sample_id"].duplicated().any():
        dups = sample_df[sample_df["sample_id"].duplicated()]["sample_id"].unique()
        errors.append(f"サンプル名が重複しています: {', '.join(dups)}")
    
    # 空文字チェック
    if (sample_df["sample_id"].str.strip() == "").any():
        errors.append("サンプル名が空の行があります。")
        
    # エラー行チェック
    if "status" in sample_df.columns and (sample_df["status"] == "error").any():
        errors.append("サンプル構造に不整合またはエラーが含まれています。一覧を確認してください。")
    
    if errors:
        return {"is_valid": False, "errors": errors, "warnings": []}
    return _valid()


def validate_sample_sheet_schema(sample_df: pd.DataFrame) -> dict:
    if sample_df is None or sample_df.empty:
        return _invalid("サンプルシートが空です。")
    
    errors = []
    # sample_id 欠損
    if (sample_df["sample_id"].str.strip() == "").any() or sample_df["sample_id"].isna().any():
        errors.append("sample_id が欠損している行があります。")
        
    # sample_id 重複
    ids = sample_df["sample_id"].astype(str).str.strip()
    valid_ids = ids[(ids != "") & (ids != "nan")]
    if valid_ids.duplicated().any():
        dups = valid_ids[valid_ids.duplicated()].unique()
        errors.append(f"sample_id が重複しています: {', '.join(dups)}")
        
    # layout_final が single-end / paired-end 以外
    if "layout_final" in sample_df.columns:
        invalid_layout = sample_df[~sample_df["layout_final"].isin(["paired-end", "single-end"])]
        if not invalid_layout.empty:
            samples = invalid_layout['sample_id'].tolist()
            errors.append(f"無効な layout が指定されています ({', '.join(samples)})")
            
    # paired-end で r1_paths / r2_paths 空
    if "layout_final" in sample_df.columns and "r1_paths" in sample_df.columns and "r2_paths" in sample_df.columns:
        pe_df = sample_df[sample_df["layout_final"] == "paired-end"]
        if not pe_df.empty:
            if pe_df["r1_paths"].apply(lambda x: isinstance(x, list) and len(x) == 0).any():
                errors.append("paired-end サンプルに r1_path が欠損しています。")
            if pe_df["r2_paths"].apply(lambda x: isinstance(x, list) and len(x) == 0).any():
                errors.append("paired-end サンプルに r2_path が欠損しています。")
                
    # single-end で r1_paths 空
    if "layout_final" in sample_df.columns and "r1_paths" in sample_df.columns:
        se_df = sample_df[sample_df["layout_final"] == "single-end"]
        if not se_df.empty:
            if se_df["r1_paths"].apply(lambda x: isinstance(x, list) and len(x) == 0).any():
                errors.append("single-end サンプルに r1_path が欠損しています。")

    if errors:
        return {"is_valid": False, "errors": errors, "warnings": []}
    return _valid()


def validate_sample_paths_from_sheet(sample_df: pd.DataFrame) -> dict:
    if sample_df is None or sample_df.empty:
        return _valid()
        
    errors = []
    for _, row in sample_df.iterrows():
        sid = row.get("sample_id", "Unknown")
        all_p = row.get("all_paths", [])
        if not isinstance(all_p, list):
            continue
        missing = [p for p in all_p if not p or not Path(p).exists()]
        if missing:
            errors.append(f"[{sid}] FASTQファイルが見つかりません: {', '.join([Path(p).name if p else 'None' for p in missing])}")
            
    if errors:
        return {"is_valid": False, "errors": errors, "warnings": []}
    return _valid()


def validate_sample_metadata_completeness(sample_df: pd.DataFrame) -> dict:
    if sample_df is None or sample_df.empty:
        return _valid()
        
    warnings = []
    check_cols = ["group", "condition", "replicate", "batch", "pair_id", "display_name"]
    
    for col in check_cols:
        if col in sample_df.columns:
            empty_mask = (sample_df[col].astype(str).str.strip() == "") | (sample_df[col].astype(str).str.strip() == "nan") | sample_df[col].isna()
            if empty_mask.any():
                samples = sample_df[empty_mask]["sample_id"].tolist()
                sample_str = ", ".join(samples[:3]) + ("..." if len(samples) > 3 else "")
                warnings.append(f"メタデータ '{col}' が空のサンプルがあります ({sample_str})。将来の連携時に警告される可能性があります。")
                
    if warnings:
        return {"is_valid": True, "errors": [], "warnings": warnings}
    return _valid()


def validate_fastq_detected(sample_df: pd.DataFrame | None) -> dict:
    if sample_df is None or sample_df.empty:
        return _invalid("FASTQ ファイルが検出されていません。")
    return _valid()


def validate_run_conditions(
    input_dir: str,
    output_dir: str,
    sample_df: pd.DataFrame | None,
    salmon_index_path: str,
    tx2gene_path: str,
    strandedness_mode: str,
    strandedness_result: dict | None,
) -> dict:
    checks = {
        "fastq_detected": validate_fastq_detected(sample_df),
        "sample_structure": validate_sample_structure(sample_df),
        "salmon_index": validate_salmon_index(salmon_index_path),
        "tx2gene": validate_tx2gene_file(tx2gene_path),
        "strandedness": validate_strandedness_selection(strandedness_mode, strandedness_result),
        "output_dir": validate_output_directory(output_dir),
        "salmon_binary": validate_salmon_binary(),
    }
    
    if sample_df is not None and not sample_df.empty and "input_source" in sample_df.columns:
        if (sample_df["input_source"] == "sample_sheet").any():
            checks["schema"] = validate_sample_sheet_schema(sample_df)
            checks["paths"] = validate_sample_paths_from_sheet(sample_df)
            checks["metadata_warnings"] = validate_sample_metadata_completeness(sample_df)

    errors = []
    warnings = []
    for result in checks.values():
        errors.extend(result.get("errors", []))
        warnings.extend(result.get("warnings", []))

    return {
        "is_valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "checks": {name: result.get("is_valid", False) for name, result in checks.items()},
    }


def validate_salmon_binary() -> dict:
    if shutil.which("salmon") is None:
        return _invalid("salmon コマンドが見つかりません。PATH を確認してください。")
    return _valid()


def _valid() -> dict:
    return {"is_valid": True, "errors": [], "warnings": []}


def _invalid(message: str) -> dict:
    return {"is_valid": False, "errors": [message], "warnings": []}

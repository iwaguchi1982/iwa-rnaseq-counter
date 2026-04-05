from __future__ import annotations

from pathlib import Path
import os
import re
import shutil

import pandas as pd

from .strandedness import validate_strandedness_selection

# --- Validation utilities ---
# validators.py は GUI / CLI から受け取った入力を検証し、
# 実行前に満たすべき前提条件を共通ルールとしてまとめる層である。
# 単体の値検証と、run 可否の総合判定の両方を担当する。


# --- Basic input/output validation ---
# 入力ディレクトリ、解析名、出力ディレクトリなど、
# 実行前に必ず成立していてほしい基本条件を検証する。

# 入力ディレクトリが存在し、ディレクトリであり、読み取り可能であることを確認する。
# FASTQ 探索や sample sheet 解決の前提条件をここで担保する。
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

# 解析名が空でなく、ファイル名や run 名として安全に使える文字だけで構成されることを確認する。
def validate_analysis_name(name: str) -> dict:
    if not name or not name.strip():
        return _invalid("解析名を入力してください。")
    
    # ファイル名として使える文字 (英数字, _, -) かつ 空白を含まない
    if not re.match(r"^[a-zA-Z0-9_\-]+$", name):
        return _invalid("解析名には英数字、アンダースコア(_)、ハイフン(-)のみ使用してください。")
    
    return _valid()

# 出力先が作成可能であり、最低限の書き込みテストを通ることを確認する。
# 実行途中で permission error を起こしにくくするための事前検証である。
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

# [v0.6.0 C-04 / C-08]
# validator 名と入力契約にはまだ Salmon 固有語彙が残っている。
# さらに判定ロジック自体も Salmon index の directory structure を直接知っている。
# quantifier 実行の抽象化とは別に、reference 側の validator 契約も
# 今後段階的に facade 化・backend 分離していく必要がある。
# - app.py から見える表面名は backend 非依存へ寄せる
# - 実体の判定ロジックは backend adapter / reference validator 側へ隔離
# を検討したい。
# --- Reference resource validation ---
# quantifier 実行と gene-level 集約に必要な reference/resource を検証する。
# 現時点では Salmon 前提の命名と directory structure が残ってい
def validate_salmon_index(salmon_index_path: str) -> dict:
    if not salmon_index_path:
        return _invalid("Salmon index が指定されていません。")
    path = Path(salmon_index_path)
    if not path.exists() or not path.is_dir():
        return _invalid("Salmon index パスが有効ではありません。")
    
    # [v0.6.0 C-08]
    # ここは単なる名前の問題ではなく、Salmon index 特有の directory structure を前提にしている。
    # version.json / versionInfo.json / ref_sigs.json を見ているため、f
    # reference/index 抽象化の際はロジックの移設先を先に決める必要がある。
    # Salmon index として最低限期待するメタデータ系ファイルの存在を確認する。
    # ここは「名前」ではなく、Salmon index 構造そのものに依存した検証である。
    version_files = ["version.json", "versionInfo.json", "ref_sigs.json"]
    if not any((path / f).exists() for f in version_files):
        return _invalid(f"Salmon インデックスとして正しくないか、不完全なディレクトリです: {path}")
    
    # [v0.6.0 C-08]
    # pos.bin / seq.bin の存在確認も Salmon 固有。
    # validator facade の名前変更だけで済ませず、
    # backend ごとの index completeness check に落とす必要あり。
    # 実体ファイルの存在を確認し、不完全な index directory を弾く。
    # backend 非依存 validator に寄せる場合も、この completeness check はどこかで必要になる。
    if not (path / "pos.bin").exists() and not (path / "seq.bin").exists():
         return _invalid("Salmon インデックス内の必須ファイル (pos.bin 等) が見つかりません。")
    
    return _valid()

# tx2gene mapping file が存在し、少なくとも transcript_id / gene_id 相当の
# 2 列以上を持つことを確認する。
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

# --- Sample structure validation ---
# sample_df の内容が downstream 実行に耐えられるかを確認する。
# 重複 sample_id、空文字、構造エラー行などをここで弾く。
# GUI / CSV 由来の sample_df が実行可能な最小整合性を持つかを確認する。
def validate_sample_structure(sample_df: pd.DataFrame | None) -> dict:
    if sample_df is None or sample_df.empty:
        return _invalid("サンプル情報がありません。")
    
    # --- Sample sheet and path validation ---
    # sample sheet 由来の入力について、
    # schema、path 解決、metadata の最低限の整合性を補助的に検証する。
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

# --- Aggregated run readiness validation ---
# 個別 validator の結果をまとめて、
# 実行前に run を開始してよいかどうかを総合判定する。
# UI 側ではこの結果をそのまま run 可否判定と詳細表示に使う。
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
            missing_paths = [str(Path(p).absolute()) if p else 'None' for p in missing]
            errors.append(f"[{sid}] FASTQファイルが見つかりません: {', '.join(missing_paths)}")
            
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


# [v0.6.0 C-03 / C-04 / C-08]
# ここは app.py / CLI と validator 群の境界にある総合判定関数である。
# quantifier 実行の抽象化は進んだ一方で、
# reference 側の引数名・check key・validator 呼び出し先にはまだ Salmon 語彙が残っている。
# 今後は facade 化と backend 非依存 naming を段階的に進めたい。

# v0.6.0 では
# - 引数名の抽象化
# - checks 辞書キーの backend 非依存化
# - backend 固有判定の委譲先整理
# を段階的に進める必要がある。
def validate_run_conditions(
    input_dir: str,
    output_dir: str,
    sample_df: pd.DataFrame | None,
    salmon_index_path: str,
    tx2gene_path: str,
    strandedness_mode: str,
    strandedness_result: dict | None,
) -> dict:
    # 各観点の validator をまとめて実行し、
    # エラーと警告を集約可能な形で保持する。
    checks = {
        "fastq_detected": validate_fastq_detected(sample_df),
        "sample_structure": validate_sample_structure(sample_df),
        # [v0.6.0 C-04 / C-08]
        # check key も validator 呼び出し先も Salmon 固有。
        # UI 表示や validation details の表示項目に影響する可能性があるため、
        # 置換時は app.py 側の参照も合わせて確認する。

        # reference 側の検証は現時点では Salmon index validator を直接利用している。
        # 今後ここを generic reference validator へ差し替えられる形にしたい。
        "salmon_index": validate_salmon_index(salmon_index_path),
        "tx2gene": validate_tx2gene_file(tx2gene_path),
        "strandedness": validate_strandedness_selection(strandedness_mode, strandedness_result),
        "output_dir": validate_output_directory(output_dir),
        # [v0.6.0 C-04]
        # 実行前チェックに salmon binary の存在確認が直結している。
        # backend 抽象化後は validate_quantifier_binary() 相当の facade か、
        # backend 実装側の preflight check へ寄せたい。

        # 実行環境チェックとして salmon binary の存在もここで確認している。
        # 将来的には backend ごとの preflight check へ分離する余地がある。
        "salmon_binary": validate_salmon_binary(),
    }
    
    if sample_df is not None and not sample_df.empty and "input_source" in sample_df.columns:
        if (sample_df["input_source"] == "sample_sheet").any():
            checks["schema"] = validate_sample_sheet_schema(sample_df)
            checks["paths"] = validate_sample_paths_from_sheet(sample_df)
            checks["metadata_warnings"] = validate_sample_metadata_completeness(sample_df)

    # 個別 check の結果を一つの run 判定へ集約し、
    # UI が扱いやすい errors / warnings / checks 形式へ正規化する。
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

# [v0.6.0 C-04]
# ここは validator 層に置かれているが、性質としては backend 実行環境の preflight check に近い。
# quantifier が複数になった場合は、backend 実装または facade 側へ責務を寄せる候補である。

# --- Backend binary validation ---
# 現行 backend 実行に必要な CLI バイナリが利用可能かを確認する。
# 今は Salmon 固定だが、将来的には backend ごとの preflight check に整理したい。
def validate_salmon_binary() -> dict:
    if shutil.which("salmon") is None:
        return _invalid("salmon コマンドが見つかりません。PATH を確認してください。")
    return _valid()


def _valid() -> dict:
    return {"is_valid": True, "errors": [], "warnings": []}


def _invalid(message: str) -> dict:
    return {"is_valid": False, "errors": [message], "warnings": []}

#
# read_sample_sheet.py
# 全体の役割：ユーザー提供のサンプルシート（CSV）を読み込み、内部モデルである AssaySpec のリスト形式へ変換するパーサー。
# 責務：表形式の入力データを解析し、各行を独立したアッセイ定義へと写像することで、バッチ処理の基盤を提供。
# ---
# 主要な機能ブロック
# 1. 共通リソースの注入:
# 全サンプルに共通して適用される定量用インデックスや tx2gene マップのパスを受け取り、各 AssaySpec の reference_resources 領域へ埋め込む。
# 
# 2. CSV 解析とバリデーション:
# csv.DictReader を用いて各行を走査。必須列（sample_id, r1_path）の存在確認や、exclude 列による実行対象からの除外判定を行う。
# 
# 3. アッセイモデルの構築:
# 各行のデータから AssaySpec オブジェクトを生成。ライブラリレイアウト（Single/Paired）の明示的な指定または R2 パスの有無による自動推定を担当。
# 
# 4. メタデータの自動収集:
# 解析に直接使われない任意の列（Group, Condition, Replicate 等）を抽出し、metadata 辞書に格納。実験デザイン情報を後続の工程まで引き継ぐ。
# 
# 5. パスと構造の整合性チェック:
# 必須ファイルの欠損や無効なレイアウト指定を検出し、パイプライン実行前に適切なエラーを発生させることでデータ整合性を担保。
# 

import csv
from pathlib import Path

from ..models.assay import AssaySpec, InputFile, ReferenceResources

def read_sample_sheet(
    sample_sheet_path: Path,
    quantifier_index_path: str | None = None,
    tx2gene_path: str | None = None,
    strandedness: str = "Auto-detect",
    salmon_index_path: str | None = None,
) -> list[AssaySpec]:
    resolved_quantifier_index = quantifier_index_path or salmon_index_path

    if not sample_sheet_path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {sample_sheet_path}")

    ref_res = None
    if resolved_quantifier_index or tx2gene_path:
        ref_res = ReferenceResources(
            quantifier_index=resolved_quantifier_index,
            tx2gene_path=tx2gene_path,
    )

    if not sample_sheet_path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {sample_sheet_path}")

    ref_res = None
    if salmon_index_path or tx2gene_path:
        ref_res = ReferenceResources(
            quantifier_index=salmon_index_path,
            tx2gene_path=tx2gene_path,
        )

    assays: list[AssaySpec] = []
    
    with sample_sheet_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        
        # Header validation
        header = reader.fieldnames if reader.fieldnames else []
        required_cols = {"sample_id", "r1_path"}
        missing_cols = required_cols - set(header)
        if missing_cols:
            raise ValueError(f"Sample sheet is missing required columns: {', '.join(missing_cols)}")

        for row in reader:
            # Handle exclude (support true/false, 0/1, yes/no)
            raw_exclude = str(row.get("exclude", "")).strip().lower()
            exclude = raw_exclude in ("true", "1", "yes", "y", "on")
            if exclude:
                continue

            sample_id = str(row.get("sample_id", "")).strip()
            if not sample_id:
                continue # Skip empty rows if sample_id is missing

            r1_path = str(row.get("r1_path", "")).strip()
            r2_path = str(row.get("r2_path", "")).strip()
            
            if not r1_path:
                raise ValueError(f"r1_path is required for sample_id: {sample_id}")

            # Layout inference
            raw_layout = str(row.get("layout", "")).strip().lower()
            if raw_layout in ("pe", "paired", "paired-end"):
                layout = "paired-end"
            elif raw_layout in ("se", "single", "single-end"):
                layout = "single-end"
            else:
                # Infer from r2_path if layout is empty or unknown
                layout = "paired-end" if r2_path else "single-end"
                
            if layout == "paired-end" and not r2_path:
                raise ValueError(f"r2_path is required for paired-end sample_id: {sample_id}")

            input_files = [InputFile(file_role="fastq_r1", path=r1_path)]
            if layout == "paired-end" and r2_path:
                input_files.append(InputFile(file_role="fastq_r2", path=r2_path))

            # Move all other columns to metadata
            metadata = {k: v.strip() for k, v in row.items() if k and k not in ("sample_id", "r1_path", "r2_path", "layout", "exclude") and v is not None}
            metadata["subject_id"] = metadata.get("subject_id", sample_id) # Set subject_id fallback if not provided

            assay = AssaySpec(
                schema_name="AssaySpec",
                schema_version="0.1.0",
                assay_id=f"ASSAY_{sample_id}",
                specimen_id=sample_id,
                assay_type="bulk_rnaseq",
                library_layout=layout,
                strandedness=strandedness,
                reference_resources=ref_res,
                input_files=input_files,
                metadata=metadata,
                overlay={}
            )
            assays.append(assay)

    return assays

from __future__ import annotations

import pandas as pd
import json
import subprocess
import tempfile
import shutil
from pathlib import Path


def infer_strandedness(sample_df: pd.DataFrame, input_dir: str, salmon_index_path: str) -> dict:
    """
    Salmon をライブラリタイプ A (Auto) で小規模実行し、strandedness を推定する。
    """
    if sample_df is None or sample_df.empty:
        return {
            "mode": "unknown",
            "confidence": "low",
            "reason": "No sample information available",
            "is_valid": False,
            "warnings": ["サンプル情報がないため strandedness を推定できません。"],
        }

    # 代表サンプルの選定
    sample_row = _pick_representative_sample(sample_df)
    if not sample_row:
        return {
            "mode": "unknown",
            "confidence": "low",
            "reason": "Compatible sample not found for probing",
            "is_valid": False,
            "warnings": ["推定に使用できる有効なサンプルが見つかりません。"],
        }

    temp_dir = Path(tempfile.mkdtemp(prefix="salmon_probe_"))
    try:
        # Salmon バイナリの存在確認
        if shutil.which("salmon") is None:
            return {
                "mode": "unknown",
                "confidence": "low",
                "reason": "Salmon binary not found",
                "is_valid": False,
                "warnings": ["salmon コマンドが見つかりません。PATH を確認してください。"],
                "probe_cmd": [],
                "stderr": "",
            }

        cmd = _build_probe_command(sample_row, salmon_index_path, temp_dir)
        # 実行 (小規模・短時間で済ませるために --validateMappings のみ利用)
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        # 結果の解析
        lib_counts_json = temp_dir / "lib_format_counts.json"
        
        # 実行失敗または結果ファイルがない場合
        if res.returncode != 0 or not lib_counts_json.exists():
            return {
                "mode": "unknown",
                "confidence": "low",
                "reason": f"Salmon probe failed (rc={res.returncode})",
                "is_valid": False,
                "warnings": ["strandedness 自動推定に失敗しました。"],
                "probe_cmd": cmd,
                "stderr": res.stderr,
            }

        result = _parse_salmon_libtype(lib_counts_json, res.stdout, res.stderr, sample_row["sample_id"])
        result["probe_cmd"] = cmd
        result["stdout"] = res.stdout
        result["stderr"] = res.stderr
        return result

    except Exception as e:
        return {
            "mode": "unknown",
            "confidence": "low",
            "reason": f"Error during probing: {str(e)}",
            "is_valid": False,
            "warnings": [f"推定実行中にエラーが発生しました: {e}"],
            "probe_cmd": [],
            "stderr": str(e),
        }
    finally:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)


def validate_strandedness_selection(strandedness_mode: str, strandedness_result: dict | None) -> dict:
    if strandedness_mode == "Auto-detect":
        if not strandedness_result:
            return {"is_valid": False, "errors": ["Auto-detect が選択されていますが、推定結果がありません。"], "warnings": []}
        if not strandedness_result.get("is_valid", False):
            return {"is_valid": False, "errors": [f"Auto-detect に失敗しました: {strandedness_result.get('reason')}"], "warnings": []}
        if strandedness_result.get("mode") == "unknown":
            return {"is_valid": False, "errors": ["strandedness を特定できませんでした。手動で指定してください。"], "warnings": []}
        
        confidence = strandedness_result.get("confidence", "low")
        warnings = []
        if confidence == "low":
            warnings.append("strandedness 推定の信頼度が低いです。結果を確認してください。")
        return {"is_valid": True, "errors": [], "warnings": warnings}
    
    return {"is_valid": True, "errors": [], "warnings": []}


def _pick_representative_sample(sample_df: pd.DataFrame) -> dict | None:
    # status が ok なものを優先
    ok_samples = sample_df[sample_df["status"] == "ok"]
    if not ok_samples.empty:
        return ok_samples.iloc[0].to_dict()
    if not sample_df.empty:
        return sample_df.iloc[0].to_dict()
    return None


def _build_probe_command(sample_row: dict, salmon_index_path: str, temp_dir: Path) -> list[str]:
    # 最小限の実行コマンド
    # -l A で自動判定。進捗は不要なので --quiet
    cmd = [
        "salmon", "quant",
        "-i", str(salmon_index_path),
        "-l", "A",
        "-o", str(temp_dir),
        "-p", "2", # 2 threads for probe
        "--validateMappings",
        "--quiet"
    ]
    
    layout = sample_row.get("layout_final", sample_row.get("layout_predicted", ""))
    r1 = ",".join(sample_row.get("r1_paths", []))
    r2 = ",".join(sample_row.get("r2_paths", []))
    all_p = ",".join(sample_row.get("all_paths", []))

    if layout == "paired-end":
        cmd.extend(["-1", r1, "-2", r2])
    else:
        cmd.extend(["-r", all_p or r1])
    
    return cmd


def _parse_salmon_libtype(json_path: Path, stdout: str, stderr: str, sample_id: str) -> dict:
    if not json_path.exists():
        # JSONがない場合は stderr 等から A -> U/SR/SF を探す
        # v0.1.0 では JSON 必須とするが、保険で unknown 扱い
        return {
            "mode": "unknown",
            "confidence": "low",
            "reason": "lib_format_counts.json not found",
            "is_valid": False,
            "warnings": ["Salmon のライブラリ判定結果ファイルが見つかりませんでした。"]
        }

    try:
        with open(json_path) as f:
            data = json.load(f)
        
        # Salmon inferred type (e.g., "IU", "ISR", "ISF", "U", "SR", "SF")
        inferred = data.get("expected_format", "unknown")
        
        mode = _normalize_salmon_inferred_type(inferred)
        
        # 信頼度の計算 (一例として 80% 以上のリードが一致しているか)
        total = data.get("num_mappings", 0)
        compatible = data.get("compatible_fragment_ratio", 0.0)
        
        confidence = "high" if compatible > 0.8 else ("medium" if compatible > 0.5 else "low")
        
        return {
            "mode": mode,
            "confidence": confidence,
            "reason": f"Inferred via Salmon probe on {sample_id} ({inferred}, frag_ratio={compatible:.2f})",
            "is_valid": True,
            "warnings": []
        }
    except Exception as e:
        return {
            "mode": "unknown",
            "confidence": "low",
            "reason": f"Parse error: {e}",
            "is_valid": False,
            "warnings": [f"判定結果の解析に失敗しました: {e}"]
        }


def _normalize_salmon_inferred_type(raw_type: str) -> str:
    # Salmon code -> forward/reverse/unstranded
    if "SR" in raw_type or "SF" in raw_type:
        # paired or single containing SR or SF
        if "SR" in raw_type: return "reverse"
        if "SF" in raw_type: return "forward"
    if "U" in raw_type:
        return "unstranded"
    return "unknown"

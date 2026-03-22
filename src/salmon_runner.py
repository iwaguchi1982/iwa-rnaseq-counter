from __future__ import annotations

import subprocess
import shutil
import json
from pathlib import Path
from typing import Any

import pandas as pd


def _normalize_strandedness_mode(strandedness_mode: str, layout: str) -> str:
    """
    UI側の strandedness 指定を Salmon 用の libType に寄せる。

    PE (paired-end) の場合:
    - Auto-detect -> A
    - unstranded -> IU
    - forward -> ISF
    - reverse -> ISR

    SE (single-end) の場合:
    - Auto-detect -> A
    - unstranded -> U
    - forward -> SF
    - reverse -> SR
    """
    is_pe = (layout == "paired-end")
    mode_map = {
        "Auto-detect": "A",
        "unstranded": "IU" if is_pe else "U",
        "forward": "ISF" if is_pe else "SF",
        "reverse": "ISR" if is_pe else "SR",
    }
    return mode_map.get(strandedness_mode, "A")


def _ensure_directory(path: str | Path) -> Path:
    """ディレクトリを作成して Path で返す。"""
    path_obj = Path(path)
    path_obj.mkdir(parents=True, exist_ok=True)
    return path_obj


def _to_path_list(value: Any) -> list[str]:
    """
    sample_df 内の r1_paths / r2_paths を list[str] に正規化する。

    想定入力:
    - list[str]
    - tuple[str, ...]
    - None
    """
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value if str(v).strip()]
    if isinstance(value, tuple):
        return [str(v) for v in value if str(v).strip()]
    # 文字列1本だけ入ってしまった場合にも最低限対応
    if isinstance(value, str) and value.strip():
        return [value]
    return []


def build_salmon_command(
    sample_row: dict,
    salmon_index_path: str,
    output_dir: str,
    strandedness_mode: str,
    threads: int = 4,
) -> list[str]:
    """
    Salmon 実行コマンドを構築する。

    Parameters
    ----------
    sample_row : dict
        少なくとも以下を含む想定:
        - sample_id
        - layout_final
        - r1_paths
        - r2_paths
        - all_paths
    salmon_index_path : str
        既存 Salmon index のパス
    output_dir : str
        サンプル別 Salmon 出力ディレクトリ
    strandedness_mode : str
        UI の strandedness 指定
    threads : int
        使用スレッド数

    Returns
    -------
    list[str]
        subprocess.run() にそのまま渡せるコマンド
    """
    sample_id = str(sample_row.get("sample_id", "unknown_sample"))
    # layout_final がなければ layout_predicted を使う
    layout_final = str(sample_row.get("layout_final", sample_row.get("layout_predicted", ""))).strip()

    r1_paths = _to_path_list(sample_row.get("r1_paths"))
    r2_paths = _to_path_list(sample_row.get("r2_paths"))
    all_paths = _to_path_list(sample_row.get("all_paths"))

    cmd = [
        "salmon",
        "quant",
        "-i", str(salmon_index_path),
        "-l", "A", # Placeholder, will be updated
        "-p", str(threads),
        "--validateMappings",
        "-o", str(output_dir),
    ]

    lib_type = _normalize_strandedness_mode(strandedness_mode, layout_final)
    # lib_type をコマンドに挿入
    cmd[cmd.index("-l") + 1] = lib_type

    if layout_final == "paired-end":
        if not r1_paths or not r2_paths:
            raise ValueError(
                f"Sample '{sample_id}' is marked as paired-end, but r1_paths or r2_paths is missing."
            )
        cmd.extend(["-1", ",".join(r1_paths)])
        cmd.extend(["-2", ",".join(r2_paths)])
    elif layout_final == "single-end":
        if all_paths:
            cmd.extend(["-r", ",".join(all_paths)])
        elif r1_paths:
            # layout が single-end でも r1_paths に入っている場合がある
            cmd.extend(["-r", ",".join(r1_paths)])
        else:
            raise ValueError(
                f"Sample '{sample_id}' is marked as single-end, but no readable FASTQ paths were found."
            )
    else:
        raise ValueError(
            f"Sample '{sample_id}' has unknown layout_final='{layout_final}'. "
            "Expected 'single-end' or 'paired-end'."
        )

    return cmd


def _write_command_log(log_path: Path, command: list[str]) -> None:
    """実行コマンドをログに追記する。"""
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write("[COMMAND]\n")
        fh.write(" ".join(command))
        fh.write("\n\n")


def _write_process_log(log_path: Path, stdout: str, stderr: str, returncode: int) -> None:
    """stdout / stderr / returncode をログ保存する。"""
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(f"[RETURN CODE]\n{returncode}\n\n")
        fh.write("[STDOUT]\n")
        fh.write(stdout if stdout else "(empty)\n")
        fh.write("\n[STDERR]\n")
        fh.write(stderr if stderr else "(empty)\n")
        fh.write("\n")


def collect_salmon_quant_paths(run_output_dir: str) -> list[dict]:
    """
    run_output_dir 配下から quant.sf を収集する。

    想定構造:
    run_output_dir/
      salmon/
        SampleA/
          quant.sf
        SampleB/
          quant.sf
    """
    root = Path(run_output_dir)
    results: list[dict] = []

    # run_output_dir/salmon/SAMPLE/quant.sf を探す
    # ただし user は output_dir に直接 salmon/ を掘る場合もあるので
    # iwa-rnaseq-counter.py の構造に合わせる
    salmon_dir = root / "salmon"
    if not salmon_dir.exists():
        return []

    for quant_path in salmon_dir.glob("*/quant.sf"):
        sample_id = quant_path.parent.name
        aux_info_dir = quant_path.parent / "aux_info"
        meta_info_json = aux_info_dir / "meta_info.json"

        results.append(
            {
                "sample_id": sample_id,
                "quant_path": str(quant_path),
                "aux_info_dir": str(aux_info_dir),
                "meta_info_json": str(meta_info_json) if meta_info_json.exists() else None,
                "is_success": quant_path.exists(),
            }
        )

    return sorted(results, key=lambda x: x["sample_id"])


def run_salmon_quant(
    sample_df: pd.DataFrame,
    salmon_index_path: str,
    run_output_dir: str,
    strandedness_mode: str,
    threads: int = 4,
) -> dict:
    """
    sample_df の各サンプルについて Salmon quant を順次実行する。

    Parameters
    ----------
    sample_df : pandas.DataFrame
        少なくとも以下の列を持つ想定:
        - sample_id
        - layout_final (なければ layout_predicted)
        - r1_paths
        - r2_paths
        - all_paths
    salmon_index_path : str
        既存 Salmon index のパス
    run_output_dir : str
        run の出力ディレクトリ
    strandedness_mode : str
        UI で選ばれた strandedness
    threads : int
        Salmon に渡すスレッド数

    Returns
    -------
    dict
        {
            "is_success": bool,
            "errors": list[str],
            "warnings": list[str],
            "log_summary": str,
            "outputs": list[dict],
        }
    """
    if sample_df is None or len(sample_df) == 0:
        return {
            "is_success": False,
            "errors": ["sample_df is empty."],
            "warnings": [],
            "log_summary": "No samples available for Salmon execution.",
            "outputs": [],
        }

    if not isinstance(sample_df, pd.DataFrame):
        return {
            "is_success": False,
            "errors": ["sample_df must be a pandas.DataFrame."],
            "warnings": [],
            "log_summary": "Invalid sample table type.",
            "outputs": [],
        }

    run_root = Path(run_output_dir)
    salmon_root = _ensure_directory(run_root / "salmon")
    logs_root = _ensure_directory(run_root / "logs")
    master_log_path = logs_root / "run.log"

    errors: list[str] = []
    warnings: list[str] = []
    outputs: list[dict] = []
    summary_lines: list[str] = []

    # マスターログ初期化
    with master_log_path.open("w", encoding="utf-8") as f:
        f.write(f"Salmon Run Started: {pd.Timestamp.now()}\n")
        f.write(f"Index: {salmon_index_path}\n\n")

    for _, row in sample_df.iterrows():
        sample_row = row.to_dict()
        sample_id = str(sample_row.get("sample_id", "unknown_sample")).strip() or "unknown_sample"

        sample_outdir = _ensure_directory(salmon_root / sample_id)
        sample_log_path = logs_root / f"{sample_id}_salmon.log"

        try:
            command = build_salmon_command(
                sample_row=sample_row,
                salmon_index_path=salmon_index_path,
                output_dir=str(sample_outdir),
                strandedness_mode=strandedness_mode,
                threads=threads,
            )
            _write_command_log(sample_log_path, command)

            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
            )
            completed_stdout = completed.stdout
            completed_stderr = completed.stderr
            completed_returncode = completed.returncode
            
            quant_path = sample_outdir / "quant.sf"
            # 実在大・非空であることを確認
            is_success = (
                completed.returncode == 0 
                and quant_path.exists() 
                and quant_path.stat().st_size > 0
            )

            _write_process_log(
                log_path=sample_log_path,
                stdout=completed_stdout,
                stderr=completed_stderr,
                returncode=completed_returncode,
            )

            if is_success:
                # meta_info.json から統計情報を取得試行
                mapping_rate = 0.0
                num_mapped = 0
                num_processed = 0
                num_decoy = 0
                num_filter = 0
                meta_path = sample_outdir / "aux_info" / "meta_info.json"
                if meta_path.exists():
                    try:
                        with open(meta_path) as f:
                            meta_data = json.load(f)
                            mapping_rate = meta_data.get("percent_mapped", 0.0)
                            num_mapped = meta_data.get("num_mapped", 0)
                            num_processed = meta_data.get("num_processed", 0)
                            num_decoy = meta_data.get("num_decoy_fragments", 0)
                            num_filter = meta_data.get("num_fragments_filtered_vm", 0)
                    except Exception:
                        pass

                outputs.append(
                    {
                        "sample_id": sample_id,
                        "quant_path": str(quant_path),
                        "aux_info_dir": str(sample_outdir / "aux_info"),
                        "log_path": str(sample_log_path),
                        "is_success": True,
                        "mapping_rate": mapping_rate,
                        "num_mapped": num_mapped,
                        "num_processed": num_processed,
                        "num_decoy": num_decoy,
                        "num_filter": num_filter,
                    }
                )
                summary_lines.append(f"[OK] {sample_id}: {num_mapped:,} reads mapped ({mapping_rate:.2f}%)")
            else:
                error_msg = f"Salmon failed for sample '{sample_id}'. See log: {sample_log_path}"
                errors.append(error_msg)
                outputs.append(
                    {
                        "sample_id": sample_id,
                        "quant_path": str(quant_path) if quant_path.exists() else None,
                        "aux_info_dir": str(sample_outdir / "aux_info"),
                        "log_path": str(sample_log_path),
                        "is_success": False,
                    }
                )
                summary_lines.append(f"[FAILED] {sample_id}")

        except Exception as exc:
            error_msg = f"Exception while running Salmon for sample '{sample_id}': {exc}"
            errors.append(error_msg)
            summary_lines.append(f"[ERROR] {sample_id}: {exc}")

            with sample_log_path.open("a", encoding="utf-8") as fh:
                fh.write("[EXCEPTION]\n")
                fh.write(str(exc))
                fh.write("\n")

            outputs.append(
                {
                    "sample_id": sample_id,
                    "quant_path": None,
                    "aux_info_dir": str(sample_outdir / "aux_info"),
                    "log_path": str(sample_log_path),
                    "is_success": False,
                }
            )

    log_summary = "\n".join(summary_lines) if summary_lines else "No Salmon execution summary available."

    return {
        "is_success": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "log_summary": log_summary,
        "outputs": outputs,
        "master_log_path": str(master_log_path),
    }



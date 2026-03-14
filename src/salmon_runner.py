import subprocess
from pathlib import Path

import pandas as pd


def build_salmon_command(
    sample_row: dict,
    salmon_index_path: str,
    output_dir: str,
    strandedness_mode: str,
) -> list[str]:
    sample_id = sample_row["sample_id"]
    sample_out = Path(output_dir) / sample_id
    cmd = ["salmon", "quant", "-i", salmon_index_path, "-o", str(sample_out), "-l", _to_salmon_libtype(strandedness_mode)]

    r1_paths = sample_row.get("r1_paths", [])
    r2_paths = sample_row.get("r2_paths", [])
    all_paths = sample_row.get("all_paths", [])

    if sample_row.get("layout_final") == "paired-end":
        cmd.extend(["-1", ",".join(r1_paths), "-2", ",".join(r2_paths)])
    else:
        cmd.extend(["-r", ",".join(all_paths)])
    return cmd



def run_salmon_quant(
    sample_df: pd.DataFrame,
    salmon_index_path: str,
    run_output_dir: str,
    strandedness_mode: str,
) -> dict:
    sample_outputs = []
    errors = []
    
    # Create required directories
    run_out_path = Path(run_output_dir)
    salmon_root = run_out_path / "salmon"
    log_root = run_out_path / "logs"
    salmon_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    for _, row in sample_df.iterrows():
        sample_id = row["sample_id"]
        sample_out = salmon_root / sample_id
        log_path = log_root / f"{sample_id}_salmon.log"
        
        cmd = build_salmon_command(
            sample_row=row.to_dict(),
            salmon_index_path=salmon_index_path,
            output_dir=str(salmon_root),
            strandedness_mode=strandedness_mode
        )
        
        try:
            with open(log_path, "w") as log_file:
                process = subprocess.Popen(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    text=True
                )
                process.wait()
                
            is_success = process.returncode == 0
            if not is_success:
                errors.append(f"Salmon failed for sample {sample_id} (exit code: {process.returncode})")
            
            sample_outputs.append(
                {
                    "sample_id": sample_id,
                    "quant_path": str(sample_out / "quant.sf"),
                    "aux_info_dir": str(sample_out / "aux_info"),
                    "is_success": is_success,
                    "log_path": str(log_path),
                }
            )
        except Exception as e:
            errors.append(f"Unexpected error running Salmon for {sample_id}: {str(e)}")
            sample_outputs.append(
                {
                    "sample_id": sample_id,
                    "quant_path": "",
                    "aux_info_dir": "",
                    "is_success": False,
                    "log_path": str(log_path) if 'log_path' in locals() else "",
                }
            )

    is_overall_success = len(errors) == 0 and len(sample_outputs) > 0
    return {
        "is_success": is_overall_success,
        "errors": errors,
        "warnings": [],
        "log_summary": "Salmon quantification completed." if is_overall_success else "Salmon quantification failed for some samples.",
        "outputs": sample_outputs,
    }


def collect_salmon_quant_paths(run_output_dir: str) -> list[dict]:
    quant_paths = []
    salmon_root = Path(run_output_dir) / "salmon"
    if not salmon_root.exists():
        return quant_paths
    for quant_file in salmon_root.rglob("quant.sf"):
        quant_paths.append(
            {
                "sample_id": quant_file.parent.name,
                "quant_path": str(quant_file),
            }
        )
    return quant_paths


def _to_salmon_libtype(strandedness_mode: str) -> str:
    mapping = {
        "Auto-detect": "A",
        "unstranded": "U",
        "forward": "SF",
        "reverse": "SR",
    }
    return mapping.get(strandedness_mode, "A")

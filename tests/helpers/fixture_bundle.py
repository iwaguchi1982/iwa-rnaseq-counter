import json
import shutil
from pathlib import Path

PLACEHOLDER = "__BUNDLE_ROOT__"

def materialize_analysis_bundle_fixture(src_dir: Path, dst_dir: Path) -> Path:
    """
    Copies a static bundle fixture to a temporary directory and replaces 
    all instances of the __BUNDLE_ROOT__ placeholder with the actual destination path.
    """
    if dst_dir.exists():
        shutil.rmtree(dst_dir)
    
    shutil.copytree(src_dir, dst_dir)
    bundle_root = dst_dir.resolve()

    for path in bundle_root.rglob("*"):
        if not path.is_file():
            continue

        # Process text files that might contain the placeholder
        if path.suffix.lower() in {".json", ".tsv", ".log"}:
            try:
                text = path.read_text(encoding="utf-8")
                if PLACEHOLDER in text:
                    path.write_text(
                        text.replace(PLACEHOLDER, str(bundle_root)),
                        encoding="utf-8",
                    )
            except Exception:
                # Skip files that aren't valid text (e.g. accidentally binary or encoding issues)
                pass

    return bundle_root

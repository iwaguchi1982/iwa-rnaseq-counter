import json
from pathlib import Path
from ..models.assay import AssaySpec

def read_assay_spec(file_path: str | Path) -> AssaySpec:
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    schema_name = data.get("$schema_name")
    if schema_name != "AssaySpec":
        raise ValueError(f"Expected $schema_name 'AssaySpec', got '{schema_name}'")
        
    return AssaySpec.from_dict(data)

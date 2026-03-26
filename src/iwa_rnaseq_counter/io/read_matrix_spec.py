import json
from pathlib import Path

from ..models.matrix import MatrixSpec


def read_matrix_spec(file_path: str | Path) -> MatrixSpec:
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    schema_name = data.get("$schema_name")
    if schema_name != "MatrixSpec":
        raise ValueError(f"Expected $schema_name 'MatrixSpec', got {schema_name!r}")

    return MatrixSpec.from_dict(data)

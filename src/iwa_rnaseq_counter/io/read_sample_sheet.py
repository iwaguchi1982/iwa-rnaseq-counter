import csv
from pathlib import Path

from ..models.assay import AssaySpec, InputFile, ReferenceResources

def read_sample_sheet(
    sample_sheet_path: Path,
    salmon_index_path: str | None = None,
    tx2gene_path: str | None = None,
    strandedness: str = "Auto-detect",
) -> list[AssaySpec]:
    """
    Reads a sample sheet CSV and yields AssaySpec objects.
    """
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
        for row in reader:
            exclude = str(row.get("exclude", "")).strip().lower() == "true"
            if exclude:
                continue

            sample_id = str(row.get("sample_id", "")).strip()
            if not sample_id:
                raise ValueError("sample_id is required in sample sheet")

            layout = str(row.get("layout", "")).strip().lower()
            if layout in ("pe", "paired"):
                layout = "paired-end"
            elif layout in ("se", "single"):
                layout = "single-end"
            else:
                layout = "paired-end" # default fallback
                
            r1_path = str(row.get("r1_path", "")).strip()
            r2_path = str(row.get("r2_path", "")).strip()

            if not r1_path:
                raise ValueError(f"r1_path is required for sample_id: {sample_id}")
            if layout == "paired-end" and not r2_path:
                raise ValueError(f"r2_path is required for paired-end sample_id: {sample_id}")

            input_files = [InputFile(file_role="fastq_r1", path=r1_path)]
            if layout == "paired-end" and r2_path:
                input_files.append(InputFile(file_role="fastq_r2", path=r2_path))

            # Move all other columns to metadata
            metadata = {k: v.strip() for k, v in row.items() if k not in ("sample_id", "r1_path", "r2_path", "layout", "exclude") and v.strip()}
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

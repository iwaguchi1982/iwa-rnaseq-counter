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
    
    Required columns: sample_id, r1_path
    Optional columns: r2_path, layout, exclude, etc.
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

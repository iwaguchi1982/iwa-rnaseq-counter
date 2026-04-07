from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any

@dataclass
class InputFile:
    file_role: str
    path: str

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "InputFile":
        return cls(
            file_role=data.get("file_role", ""),
            path=data.get("path", "")
        )

@dataclass
class ReferenceResources:
    genome_build: Optional[str] = None
    annotation_release: Optional[str] = None
    quantifier_index: Optional[str] = None
    tx2gene_path: Optional[str] = None
    annotation_gtf_path: Optional[str] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ReferenceResources":
        return cls(
            genome_build=data.get("genome_build"),
            annotation_release=data.get("annotation_release"),
            quantifier_index=data.get("quantifier_index"),
            tx2gene_path=data.get("tx2gene_path"),
            annotation_gtf_path=data.get("annotation_gtf_path")
        )

@dataclass
class AssaySpec:
    schema_name: str
    schema_version: str
    assay_id: str
    specimen_id: str
    assay_type: str
    platform: Optional[str] = None
    library_strategy: Optional[str] = None
    library_layout: Optional[str] = None
    strandedness: Optional[str] = None
    measurement_target: Optional[str] = None
    processing_level: Optional[str] = None
    reference_resources: Optional[ReferenceResources] = None
    input_files: List[InputFile] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    overlay: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AssaySpec":
        ref_res = data.get("reference_resources")
        reference_resources = ReferenceResources.from_dict(ref_res) if ref_res else None
        
        input_files = [InputFile.from_dict(f) for f in data.get("input_files", [])]
        
        return cls(
            schema_name=data.get("$schema_name", ""),
            schema_version=data.get("$schema_version", ""),
            assay_id=data.get("assay_id", ""),
            specimen_id=data.get("specimen_id", ""),
            assay_type=data.get("assay_type", ""),
            platform=data.get("platform"),
            library_strategy=data.get("library_strategy"),
            library_layout=data.get("library_layout"),
            strandedness=data.get("strandedness"),
            measurement_target=data.get("measurement_target"),
            processing_level=data.get("processing_level"),
            reference_resources=reference_resources,
            input_files=input_files,
            metadata=data.get("metadata", {}),
            overlay=data.get("overlay", {})
        )

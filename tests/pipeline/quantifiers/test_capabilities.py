import pytest
from iwa_rnaseq_counter.pipeline.quantifiers.registry import get_quantifier
from iwa_rnaseq_counter.pipeline.quantifiers.base import BackendCapabilities

def test_salmon_capabilities():
    quant = get_quantifier("salmon")
    cap = quant.get_capabilities()
    
    assert isinstance(cap, BackendCapabilities)
    assert cap.aggregation_input_kind == "transcript_quant"
    assert cap.has_transcript_quant is True
    assert cap.has_gene_counts is False
    assert cap.has_mapping_metrics is True
    
    req = cap.reference_requirements
    assert req.quantifier_index == "required"
    assert req.tx2gene == "required"
    assert req.annotation_gtf == "unused"

def test_star_capabilities():
    quant = get_quantifier("star")
    cap = quant.get_capabilities()
    
    assert isinstance(cap, BackendCapabilities)
    assert cap.aggregation_input_kind == "gene_counts"
    assert cap.has_transcript_quant is False
    assert cap.has_gene_counts is True
    assert cap.has_mapping_metrics is False
    
    req = cap.reference_requirements
    assert req.quantifier_index == "required"
    assert req.tx2gene == "unused"
    assert req.annotation_gtf == "unused"

def test_hisat2_capabilities():
    quant = get_quantifier("hisat2")
    cap = quant.get_capabilities()
    
    assert isinstance(cap, BackendCapabilities)
    assert cap.aggregation_input_kind == "gene_counts"
    assert cap.has_transcript_quant is False
    assert cap.has_gene_counts is True
    assert cap.has_mapping_metrics is False
    
    req = cap.reference_requirements
    assert req.quantifier_index == "required"
    assert req.tx2gene == "unused"
    assert req.annotation_gtf == "required"

def test_unsupported_quantifier():
    with pytest.raises(NotImplementedError, match="Unsupported quantifier"):
        get_quantifier("unsupported_backend")

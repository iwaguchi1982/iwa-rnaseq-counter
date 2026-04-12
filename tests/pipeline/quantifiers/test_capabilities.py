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
    assert cap.requires_tx2gene is True
    assert cap.requires_annotation_gtf is False

def test_star_capabilities():
    quant = get_quantifier("star")
    cap = quant.get_capabilities()
    
    assert isinstance(cap, BackendCapabilities)
    assert cap.aggregation_input_kind == "gene_counts"
    assert cap.has_transcript_quant is False
    assert cap.has_gene_counts is True
    assert cap.has_mapping_metrics is False
    assert cap.requires_tx2gene is False
    assert cap.requires_annotation_gtf is False

def test_hisat2_capabilities():
    quant = get_quantifier("hisat2")
    cap = quant.get_capabilities()
    
    assert isinstance(cap, BackendCapabilities)
    assert cap.aggregation_input_kind == "gene_counts"
    assert cap.has_transcript_quant is False
    assert cap.has_gene_counts is True
    assert cap.has_mapping_metrics is False
    assert cap.requires_tx2gene is False
    assert cap.requires_annotation_gtf is True

def test_unsupported_quantifier():
    with pytest.raises(NotImplementedError, match="Unsupported quantifier"):
        get_quantifier("unsupported_backend")

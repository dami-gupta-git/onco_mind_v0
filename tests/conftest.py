"""Pytest configuration and fixtures."""

import pytest


@pytest.fixture
def sample_variant_input():
    """Sample variant input for testing."""
    from oncomind.models.variant import VariantInput

    return VariantInput(
        gene="BRAF",
        variant="V600E",
        tumor_type="Melanoma",
    )


@pytest.fixture
def sample_evidence():
    """Sample evidence for testing."""
    from oncomind.models.evidence import CIViCEvidence, Evidence

    civic_ev = CIViCEvidence(
        evidence_type="Predictive",
        evidence_level="A",
        evidence_direction="Supports",
        clinical_significance="Sensitivity/Response",
        disease="Melanoma",
        drugs=["Vemurafenib", "Dabrafenib"],
        description="BRAF V600E mutation confers sensitivity to BRAF inhibitors in melanoma.",
        source="PubMed",
        rating=5,
    )

    return Evidence(
        variant_id="BRAF:V600E",
        gene="BRAF",
        variant="V600E",
        civic=[civic_ev],
        clinvar=[],
        cosmic=[],
    )


@pytest.fixture
def mock_llm_response():
    """Mock LLM response for testing."""
    return """{
        "summary": "BRAF V600E is a well-established actionable mutation in melanoma with FDA-approved targeted therapies.",
        "rationale": "Multiple FDA-approved BRAF inhibitors (vemurafenib, dabrafenib, encorafenib) exist for this mutation in melanoma. Strong clinical evidence from multiple phase III trials.",
        "recommended_therapies": [
            {
                "drug_name": "Vemurafenib",
                "evidence_level": "FDA-approved",
                "approval_status": "Approved",
                "clinical_context": "First-line therapy"
            },
            {
                "drug_name": "Dabrafenib + Trametinib",
                "evidence_level": "FDA-approved",
                "approval_status": "Approved",
                "clinical_context": "First-line therapy"
            }
        ],
        "references": ["Chapman PB et al. NEJM 2011", "Hauschild A et al. Lancet 2012"]
    }"""

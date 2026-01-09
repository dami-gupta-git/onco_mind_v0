"""Integration tests for FDA label data flow.

Verifies that all fields (especially initial_approval_date) are correctly
passed through the entire data pipeline:
1. OpenFDA API fetch -> FDALabelEvidence (directly)
2. FDALabelEvidence -> backend serialization (for Streamlit UI)

Run with: pytest tests/integration/api/test_fda_label_data_flow.py -v
"""

import pytest
from unittest.mock import patch, MagicMock

from oncomind.api.fda_label_service import (
    fetch_drug_label_from_openfda,
    get_fda_labels_for_drugs,
)
from oncomind.models.evidence.fda import FDALabelEvidence, ClinicalStudyEvidence


class TestFDALabelDataFlow:
    """Test that all fields flow through the entire pipeline."""

    def test_fetch_returns_initial_approval_date(self):
        """Verify fetch_drug_label_from_openfda extracts initial_approval_date."""
        # This test hits the real API - skip in CI if needed
        result = fetch_drug_label_from_openfda("capivasertib")

        assert result is not None, "Should return data for capivasertib"
        assert "initial_approval_date" in result, "Should have initial_approval_date key"
        assert result["initial_approval_date"] is not None, "initial_approval_date should not be None"
        # Format should be YYYY-MM-DD
        assert len(result["initial_approval_date"]) == 10, "Date should be in YYYY-MM-DD format"
        assert result["initial_approval_date"][4] == "-", "Date should have dash at position 4"

    def test_fda_label_evidence_has_initial_approval_date(self):
        """Verify FDALabelEvidence model has initial_approval_date field."""
        evidence = FDALabelEvidence(
            drug="test_drug",
            gene="TEST",
            initial_approval_date="2023-11-14",
        )

        assert evidence.initial_approval_date == "2023-11-14"

    def test_fda_label_evidence_model_dump_includes_initial_approval_date(self):
        """Verify FDALabelEvidence.model_dump() includes initial_approval_date."""
        evidence = FDALabelEvidence(
            drug="test_drug",
            gene="TEST",
            initial_approval_date="2023-11-14",
        )

        result = evidence.model_dump()

        assert "initial_approval_date" in result
        assert result["initial_approval_date"] == "2023-11-14"

    def test_backend_serialization_includes_initial_approval_date(self):
        """Verify backend serialization includes initial_approval_date.

        This simulates what backend.py does when serializing fda_labels.
        """
        # Create FDALabelEvidence with all fields
        evidence = FDALabelEvidence(
            drug="capivasertib",
            gene="AKT1",
            brand_name="TRUQAP",
            initial_approval_date="2023-11-14",
            clinical_studies=ClinicalStudyEvidence(
                trial_name="CAPItello-291",
                nct_id="NCT04305496",
            ),
        )

        # Simulate backend serialization (from backend.py lines 252-288)
        serialized = {
            "drug": evidence.drug,
            "gene": evidence.gene,
            "brand_name": evidence.brand_name,
            "generic_name": evidence.generic_name,
            "manufacturer": evidence.manufacturer,
            "indications_and_usage": evidence.indications_and_usage,
            "initial_approval_date": evidence.initial_approval_date,
            "clinical_studies": {
                "trial_name": evidence.clinical_studies.trial_name,
                "nct_id": evidence.clinical_studies.nct_id,
            } if evidence.clinical_studies else None,
            "last_label_update": evidence.last_label_update,
            "update_reason": evidence.update_reason,
        }

        assert serialized["initial_approval_date"] == "2023-11-14"
        assert serialized["clinical_studies"]["trial_name"] == "CAPItello-291"


class TestFDALabelEvidenceFields:
    """Ensure FDALabelEvidence has all required fields."""

    def test_required_fields_present(self):
        """Verify key fields exist in FDALabelEvidence."""
        required_fields = [
            "drug",
            "gene",
            "brand_name",
            "generic_name",
            "manufacturer",
            "indications_and_usage",
            "initial_approval_date",
            "approved_indications",
            "clinical_studies",
            "mechanism_of_action",
            "adverse_reactions",
            "locus_variant_match",
        ]

        fda_evidence = FDALabelEvidence(drug="test", gene="TEST")

        for field in required_fields:
            assert hasattr(fda_evidence, field), \
                f"FDALabelEvidence missing field: {field}"


class TestLocusVariantMatchPopulation:
    """Test populate_locus_variant_match for computing locus_variant_match."""

    def test_populate_locus_variant_match_exact(self):
        """Test locus_variant_match for exact variant match."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match

        # Sotorasib is approved for KRAS G12C
        evidence = FDALabelEvidence(
            drug="sotorasib",
            gene="KRAS",
            indications_and_usage="Treatment of KRAS G12C-mutated NSCLC",
        )

        # Populate with exact matching variant
        populate_locus_variant_match([evidence], query_variant="G12C")

        # Should match at variant level (exact match)
        assert evidence.locus_variant_match is not None
        assert evidence.locus_variant_match.level == "variant"
        assert evidence.locus_variant_match.scope == "specific"

    def test_populate_locus_variant_match_codon_level(self):
        """Test locus_variant_match when query variant is at same codon but different AA."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match

        # Sotorasib is approved for KRAS G12C
        evidence = FDALabelEvidence(
            drug="sotorasib",
            gene="KRAS",
            indications_and_usage="Treatment of KRAS G12C-mutated NSCLC",
        )

        # Query with G12D (same codon, different variant)
        populate_locus_variant_match([evidence], query_variant="G12D")

        # Should match at codon level (not variant level)
        assert evidence.locus_variant_match is not None
        assert evidence.locus_variant_match.level == "codon"
        assert evidence.locus_variant_match.scope == "specific"

    def test_populate_locus_variant_match_gene_level(self):
        """Test locus_variant_match for gene-level approvals."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match

        # Capivasertib is approved for AKT1 alterations (gene-level)
        evidence = FDALabelEvidence(
            drug="capivasertib",
            gene="AKT1",
            indications_and_usage="Treatment of AKT1 alteration-positive breast cancer",
        )

        # Query with E17K
        populate_locus_variant_match([evidence], query_variant="E17K")

        # Should match at gene level
        assert evidence.locus_variant_match is not None
        assert evidence.locus_variant_match.level == "gene"
        assert evidence.locus_variant_match.scope == "specific"

    def test_populate_locus_variant_match_no_query_variant(self):
        """Test locus_variant_match when no query variant is provided."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match

        # Gene-level approval
        evidence = FDALabelEvidence(
            drug="capivasertib",
            gene="AKT1",
            indications_and_usage="Treatment of AKT1 alteration-positive breast cancer",
        )

        # No query variant
        populate_locus_variant_match([evidence], query_variant=None)

        # Should have gene-level match with unspecified scope
        assert evidence.locus_variant_match is not None
        assert evidence.locus_variant_match.level == "gene"
        assert evidence.locus_variant_match.scope == "unspecified"

    def test_populate_locus_variant_match_no_match(self):
        """Test locus_variant_match when query doesn't match the approval."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match

        # KRAS G12C approval
        evidence = FDALabelEvidence(
            drug="sotorasib",
            gene="KRAS",
            indications_and_usage="Treatment of KRAS G12C-mutated NSCLC",
        )

        # Query with completely different variant (different codon)
        populate_locus_variant_match([evidence], query_variant="Q61H")

        # No match - Q61H is not covered by G12C approval
        assert evidence.locus_variant_match is None


class TestGetFDALabelsForDrugsReturnsEvidence:
    """Test that get_fda_labels_for_drugs returns FDALabelEvidence directly."""

    @patch('oncomind.api.fda_label_service.load_fda_labels_json')
    def test_returns_fda_label_evidence(self, mock_load):
        """Verify get_fda_labels_for_drugs returns FDALabelEvidence objects."""
        # Mock the cache to return data
        mock_load.return_value = {
            "sotorasib": {
                "brand_name": ["LUMAKRAS"],
                "generic_name": ["sotorasib"],
                "manufacturer_name": ["Amgen Inc."],
                "indications_and_usage": ["Treatment of KRAS G12C-mutated NSCLC"],
                "initial_approval_date": "2021-05-28",
                "clinical_studies": {
                    "trial_name": "CodeBreaK 100",
                    "nct_id": "NCT03600883",
                },
            }
        }

        results = get_fda_labels_for_drugs({"sotorasib"}, gene="KRAS", fetch_missing=False)

        assert len(results) == 1
        assert isinstance(results[0], FDALabelEvidence)
        assert results[0].drug == "sotorasib"
        assert results[0].gene == "KRAS"
        assert results[0].brand_name == "LUMAKRAS"
        assert results[0].initial_approval_date == "2021-05-28"
        assert results[0].clinical_studies is not None
        assert results[0].clinical_studies.trial_name == "CodeBreaK 100"

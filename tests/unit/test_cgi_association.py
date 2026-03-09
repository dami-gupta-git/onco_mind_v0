"""Tests for CGI association field handling.

CGI biomarkers have an 'association' field with values:
- "Responsive" - drug sensitivity
- "Resistant" - drug resistance
- "No Responsive" - NOT sensitivity (should not be treated as sensitivity)
- "Increased Toxicity" variants - toxicity signals (not sensitivity/resistance)

These tests verify that only explicit "Responsive" is treated as sensitivity.
"""

import pytest
from oncomind.models.evidence.evidence import Evidence, VariantIdentifiers
from oncomind.models.evidence.cgi import CGIBiomarkerEvidence


def make_identifiers() -> VariantIdentifiers:
    """Create minimal variant identifiers for testing."""
    return VariantIdentifiers(
        variant_id="EGFR:L858R",
        gene="EGFR",
        variant="L858R",
    )


def make_biomarker(association: str, drug: str = "TestDrug") -> CGIBiomarkerEvidence:
    """Create a CGI biomarker with the given association."""
    return CGIBiomarkerEvidence(
        gene="EGFR",
        alteration="EGFR:L858R",
        drug=drug,
        drug_status="Approved",
        association=association,
        evidence_level="FDA guidelines",
        source="CGI",
        tumor_type="NSCLC",
        tumor_type_full="Non-Small Cell Lung Cancer",
        fda_approved=True,
        locus_match="variant",
    )


class TestCGIAssociationHandling:
    """Tests for CGI association field interpretation."""

    def test_responsive_is_sensitivity(self):
        """'Responsive' association should be treated as sensitivity."""
        biomarker = make_biomarker("Responsive", "Erlotinib")
        evidence = Evidence(identifiers=make_identifiers(), cgi_biomarkers=[biomarker])

        # Check that therapeutic data has response_type="Sensitivity"
        therapeutic = evidence.get_therapeutic_evidence()
        erlotinib_entries = [t for t in therapeutic if t.drug_name == "Erlotinib"]
        assert len(erlotinib_entries) == 1
        assert erlotinib_entries[0].response_type == "Sensitivity"

    def test_resistant_is_resistance(self):
        """'Resistant' association should be treated as resistance."""
        biomarker = make_biomarker("Resistant", "Erlotinib")
        evidence = Evidence(identifiers=make_identifiers(), cgi_biomarkers=[biomarker])

        therapeutic = evidence.get_therapeutic_evidence()
        erlotinib_entries = [t for t in therapeutic if t.drug_name == "Erlotinib"]
        assert len(erlotinib_entries) == 1
        assert erlotinib_entries[0].response_type == "Resistance"

    def test_no_responsive_is_not_sensitivity(self):
        """'No Responsive' should NOT be treated as sensitivity.

        This is a critical test - the old code would mark this as sensitivity
        because it didn't contain 'RESIST'. The fix requires explicit 'RESPONSIVE'.
        """
        biomarker = make_biomarker("No Responsive", "Gefitinib")
        evidence = Evidence(identifiers=make_identifiers(), cgi_biomarkers=[biomarker])

        therapeutic = evidence.get_therapeutic_evidence()
        gefitinib_entries = [t for t in therapeutic if t.drug_name == "Gefitinib"]
        assert len(gefitinib_entries) == 1
        # Should be None, NOT "Sensitivity"
        assert gefitinib_entries[0].response_type is None

    def test_increased_toxicity_is_not_sensitivity(self):
        """'Increased Toxicity' should NOT be treated as sensitivity.

        Toxicity signals are about adverse effects, not drug efficacy.
        """
        biomarker = make_biomarker("Increased Toxicity", "Irinotecan")
        evidence = Evidence(identifiers=make_identifiers(), cgi_biomarkers=[biomarker])

        therapeutic = evidence.get_therapeutic_evidence()
        irinotecan_entries = [t for t in therapeutic if t.drug_name == "Irinotecan"]
        assert len(irinotecan_entries) == 1
        # Should be None, NOT "Sensitivity"
        assert irinotecan_entries[0].response_type is None

    def test_increased_toxicity_myelosuppression(self):
        """'Increased Toxicity (Myelosupression)' should NOT be sensitivity."""
        biomarker = make_biomarker("Increased Toxicity (Myelosupression)", "Cisplatin")
        evidence = Evidence(identifiers=make_identifiers(), cgi_biomarkers=[biomarker])

        therapeutic = evidence.get_therapeutic_evidence()
        cisplatin_entries = [t for t in therapeutic if t.drug_name == "Cisplatin"]
        assert len(cisplatin_entries) == 1
        assert cisplatin_entries[0].response_type is None

    def test_drug_has_sensitivity_responsive_only(self):
        """_drug_has_sensitivity_for_gene should only match explicit 'Responsive'."""
        # Create evidence with "No Responsive" - should NOT count as sensitivity
        no_responsive = make_biomarker("No Responsive", "Gefitinib")
        evidence = Evidence(
            identifiers=make_identifiers(), cgi_biomarkers=[no_responsive]
        )

        # Should return False - "No Responsive" is not sensitivity
        assert evidence._drug_has_sensitivity_for_gene("gefitinib") is False

        # Now with actual "Responsive"
        responsive = make_biomarker("Responsive", "Erlotinib")
        evidence2 = Evidence(
            identifiers=make_identifiers(), cgi_biomarkers=[responsive]
        )

        # Should return True
        assert evidence2._drug_has_sensitivity_for_gene("erlotinib") is True

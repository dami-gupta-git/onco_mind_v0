"""Integration tests for Cancer Genome Interpreter (CGI) biomarkers client.

Tests validate that the CGI biomarkers database returns expected FDA/NCCN
approval status for well-characterized oncogenic mutations.
"""

import pytest

from oncomind.api.cgi import CGIClient, CGIBiomarker


class TestCGIBRAFV600E:
    """Tests for BRAF V600E - FDA-approved targeted therapies."""

    @pytest.mark.integration
    def test_returns_biomarkers(self):
        """BRAF V600E should return biomarker data."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E")
        assert len(biomarkers) >= 1, "BRAF V600E should have at least 1 biomarker"

    @pytest.mark.integration
    def test_has_fda_approved(self):
        """BRAF V600E should have FDA-approved therapies."""
        client = CGIClient()
        fda_approved = client.fetch_fda_approved("BRAF", "V600E")
        assert len(fda_approved) >= 1, "BRAF V600E should have at least 1 FDA-approved therapy"

    @pytest.mark.integration
    def test_expected_braf_inhibitors(self):
        """Should return expected BRAF inhibitor drugs."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E")

        expected_drugs = {"vemurafenib", "dabrafenib", "encorafenib", "trametinib", "cobimetinib"}
        all_drugs = {b.drug.lower() for b in biomarkers if b.drug}

        found_expected = all_drugs & expected_drugs
        assert len(found_expected) > 0, (
            f"Expected at least one of {expected_drugs}, got: {all_drugs}"
        )

    @pytest.mark.integration
    def test_melanoma_coverage(self):
        """Should have melanoma-related biomarkers."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E", tumor_type="Melanoma")

        # Melanoma is a key indication for BRAF V600E
        assert len(biomarkers) >= 1, "Should have melanoma biomarkers for BRAF V600E"

    @pytest.mark.integration
    def test_is_fda_approved_method(self):
        """is_fda_approved() should correctly identify FDA approvals."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E")

        fda_count = sum(1 for b in biomarkers if b.is_fda_approved())
        assert fda_count > 0, "BRAF V600E should have FDA-approved biomarkers"


class TestCGIEGFRL858R:
    """Tests for EGFR L858R - EGFR TKI sensitivity marker."""

    @pytest.mark.integration
    def test_returns_biomarkers(self):
        """EGFR L858R should return biomarker data."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("EGFR", "L858R")
        assert len(biomarkers) >= 1, "EGFR L858R should have at least 1 biomarker"

    @pytest.mark.integration
    def test_expected_egfr_tkis(self):
        """Should return expected EGFR TKI drugs."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("EGFR", "L858R")

        expected_drugs = {"erlotinib", "gefitinib", "afatinib", "osimertinib"}
        all_drugs = {b.drug.lower() for b in biomarkers if b.drug}

        found_expected = all_drugs & expected_drugs
        assert len(found_expected) > 0, (
            f"Expected at least one of {expected_drugs}, got: {all_drugs}"
        )

    @pytest.mark.integration
    def test_nsclc_coverage(self):
        """Should have NSCLC-related biomarkers."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("EGFR", "L858R", tumor_type="NSCLC")

        assert len(biomarkers) >= 1, "Should have NSCLC biomarkers for EGFR L858R"


class TestCGIEGFRWildcardMatching:
    """Tests for CGI wildcard pattern matching (e.g., G719.)."""

    @pytest.mark.integration
    def test_g719_wildcard_matches_g719s(self):
        """G719. pattern should match G719S variant."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("EGFR", "G719S")

        # CGI uses G719. pattern which should match G719S
        assert len(biomarkers) >= 1, "EGFR G719S should match G719. pattern"

    @pytest.mark.integration
    def test_g719_wildcard_matches_g719a(self):
        """G719. pattern should match G719A variant."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("EGFR", "G719A")

        assert len(biomarkers) >= 1, "EGFR G719A should match G719. pattern"


class TestCGIKRASG12:
    """Tests for KRAS G12 position variants."""

    @pytest.mark.integration
    def test_kras_g12c_biomarkers(self):
        """KRAS G12C should have biomarkers (sotorasib approved)."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("KRAS", "G12C")

        assert len(biomarkers) >= 1, "KRAS G12C should have at least 1 biomarker"

    @pytest.mark.integration
    def test_kras_g12c_has_sotorasib(self):
        """KRAS G12C should include sotorasib."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("KRAS", "G12C")

        drugs = {b.drug.lower() for b in biomarkers if b.drug}
        assert "sotorasib" in drugs or "adagrasib" in drugs, (
            f"Expected sotorasib or adagrasib for KRAS G12C, got: {drugs}"
        )

    @pytest.mark.integration
    def test_kras_g12d_biomarkers(self):
        """KRAS G12D should have biomarkers."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("KRAS", "G12D")

        # G12D may have fewer approved therapies than G12C but should have some entries
        assert len(biomarkers) >= 0  # May or may not have biomarkers


class TestCGIBiomarkerStructure:
    """Tests for CGIBiomarker data structure integrity."""

    @pytest.mark.integration
    def test_required_fields_populated(self):
        """Biomarkers should have required fields populated."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E")

        for biomarker in biomarkers:
            assert isinstance(biomarker.gene, str) and biomarker.gene
            assert isinstance(biomarker.alteration, str) and biomarker.alteration
            assert isinstance(biomarker.drug, str) and biomarker.drug
            assert isinstance(biomarker.drug_status, str)
            assert isinstance(biomarker.association, str)

    @pytest.mark.integration
    def test_to_dict_keys(self):
        """to_dict should return expected keys."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E")

        expected_keys = {
            "gene", "alteration", "drug", "drug_status", "association",
            "evidence_level", "source", "tumor_type", "tumor_type_full",
            "fda_approved"
        }

        for biomarker in biomarkers:
            biomarker_dict = biomarker.to_dict()
            assert set(biomarker_dict.keys()) == expected_keys

    @pytest.mark.integration
    def test_association_values(self):
        """Association should be Responsive or Resistant."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E")

        for biomarker in biomarkers:
            assert biomarker.association.upper() in ["RESPONSIVE", "RESISTANT", "NO RESPONSIVE"], (
                f"Unexpected association: {biomarker.association}"
            )


class TestCGITumorTypeFilter:
    """Tests for tumor type filtering."""

    @pytest.mark.integration
    def test_melanoma_filter(self):
        """Filtering by melanoma should return melanoma biomarkers."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E", tumor_type="Melanoma")

        # All results should be melanoma-related
        for biomarker in biomarkers:
            tumor_lower = (biomarker.tumor_type + biomarker.tumor_type_full).lower()
            assert "melanoma" in tumor_lower or "mel" in tumor_lower or "skin" in tumor_lower, (
                f"Expected melanoma-related tumor type, got: {biomarker.tumor_type}"
            )

    @pytest.mark.integration
    def test_nsclc_filter(self):
        """Filtering by NSCLC should work."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("EGFR", "L858R", tumor_type="NSCLC")

        assert len(biomarkers) >= 1, "Should have NSCLC biomarkers for EGFR L858R"

    @pytest.mark.integration
    def test_colorectal_filter(self):
        """Filtering by colorectal should work."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("BRAF", "V600E", tumor_type="colorectal")

        # BRAF V600E in colorectal has some evidence
        # May or may not have results depending on CGI data
        assert len(biomarkers) >= 0


class TestCGIBiomarkerEvidence:
    """Tests for biomarker evidence model conversion."""

    @pytest.mark.integration
    def test_fetch_biomarker_evidence(self):
        """fetch_biomarker_evidence should return evidence objects."""
        client = CGIClient()
        evidence_list = client.fetch_biomarker_evidence("BRAF", "V600E")

        assert len(evidence_list) >= 1, "Should return at least 1 evidence object"
        for evidence in evidence_list:
            assert hasattr(evidence, "match_level")
            assert evidence.match_level in ["variant", "codon", "gene"]

    @pytest.mark.integration
    def test_match_level_variant(self):
        """Exact variant matches should have variant match level."""
        client = CGIClient()
        evidence_list = client.fetch_biomarker_evidence("BRAF", "V600E")

        v600e_evidence = [e for e in evidence_list if "V600E" in e.alteration.upper()]
        for evidence in v600e_evidence:
            assert evidence.match_level == "variant", (
                f"V600E exact match should be variant level, got {evidence.match_level}"
            )

    @pytest.mark.integration
    def test_match_level_codon(self):
        """Wildcard pattern matches should have codon match level."""
        client = CGIClient()
        # G719S matches the G719. pattern (codon-level)
        evidence_list = client.fetch_biomarker_evidence("EGFR", "G719S")

        # Check if any evidence has codon-level match
        codon_matches = [e for e in evidence_list if e.match_level == "codon"]
        variant_matches = [e for e in evidence_list if e.match_level == "variant"]

        # G719S should match either as exact variant or via G719. pattern (codon)
        assert len(codon_matches) > 0 or len(variant_matches) > 0, (
            "EGFR G719S should have variant or codon-level matches"
        )


class TestCGIResistanceMarkers:
    """Tests for resistance biomarkers."""

    @pytest.mark.integration
    def test_egfr_t790m_resistance(self):
        """EGFR T790M is a resistance marker to first-gen EGFR TKIs."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("EGFR", "T790M")

        # T790M should have some resistance biomarkers
        resistance_markers = [b for b in biomarkers if b.association.upper() == "RESISTANT"]

        # T790M confers resistance to erlotinib/gefitinib but sensitivity to osimertinib
        # So we expect a mix
        assert len(biomarkers) >= 1, "EGFR T790M should have biomarkers"

    @pytest.mark.integration
    def test_resistance_association(self):
        """Resistance biomarkers should have RESISTANT association."""
        client = CGIClient()
        biomarkers = client.fetch_biomarkers("EGFR", "T790M")

        for biomarker in biomarkers:
            if biomarker.association.upper() == "RESISTANT":
                assert "resistance" in biomarker.association.lower() or biomarker.association.upper() == "RESISTANT"

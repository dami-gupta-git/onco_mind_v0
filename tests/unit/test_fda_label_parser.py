"""Unit tests for FDA label parser.

Tests the FDALabelParser class and match_variant_to_indications function.
"""

import pytest

from oncomind.api.fda_label_parser import FDALabelParser, match_variant_to_indications
from oncomind.models.evidence.fda_biomarker import (
    BiomarkerRequirement,
    SpecificityLevel,
)


class TestFDALabelParserIMJUDO:
    """Tests for IMJUDO label parsing - the key REQUIRED_NEGATIVE case."""

    @pytest.fixture
    def parser(self):
        return FDALabelParser()

    @pytest.fixture
    def imjudo_label(self):
        """IMJUDO label with EGFR/ALK REQUIRED_NEGATIVE indication."""
        imjudo_indication = """
        1 INDICATIONS AND USAGE IMJUDO is a cytotoxic T-lymphocyte-associated antigen 4 (CTLA-4)
        blocking antibody indicated: in combination with durvalumab, for the treatment of adult
        patients with unresectable hepatocellular carcinoma (uHCC). in combination with
        durvalumab and platinum-based chemotherapy for the treatment of adult patients with
        metastatic non-small cell lung cancer (NSCLC) with no sensitizing epidermal growth factor
        receptor (EGFR) mutation or anaplastic lymphoma kinase (ALK) genomic tumor aberrations.
        1.1 Hepatocellular Carcinoma IMJUDO, in combination with durvalumab, is indicated
        for the treatment of adult patients with unresectable hepatocellular carcinoma (uHCC).
        1.2 Non-Small Cell Lung Cancer (NSCLC) IMJUDO, in combination with durvalumab and
        platinum-based chemotherapy, is indicated for the treatment of adult patients with
        metastatic NSCLC with no sensitizing epidermal growth factor receptor (EGFR) mutations
        or anaplastic lymphoma kinase (ALK) genomic tumor aberrations.
        """
        return {
            "openfda": {
                "generic_name": ["tremelimumab"],
                "brand_name": ["IMJUDO"],
            },
            "indications_and_usage": [imjudo_indication],
        }

    def test_parses_egfr_required_negative(self, parser, imjudo_label):
        """IMJUDO should have EGFR with REQUIRED_NEGATIVE requirement."""
        indications = parser.parse_label(imjudo_label)

        egfr_indications = [i for i in indications if i.gene == "EGFR"]
        assert len(egfr_indications) >= 1

        egfr = egfr_indications[0]
        assert egfr.requirement == BiomarkerRequirement.REQUIRED_NEGATIVE

    def test_parses_alk_required_negative(self, parser, imjudo_label):
        """IMJUDO should have ALK with REQUIRED_NEGATIVE requirement."""
        indications = parser.parse_label(imjudo_label)

        alk_indications = [i for i in indications if i.gene == "ALK"]
        assert len(alk_indications) >= 1

        alk = alk_indications[0]
        assert alk.requirement == BiomarkerRequirement.REQUIRED_NEGATIVE

    def test_egfr_t790m_excluded_from_imjudo(self, parser, imjudo_label):
        """EGFR T790M in NSCLC should be EXCLUDED from IMJUDO."""
        indications = parser.parse_label(imjudo_label)

        results = match_variant_to_indications(indications, "EGFR", "T790M", "NSCLC")

        # Should find IMJUDO but marked as excluded
        imjudo_results = [
            r for r in results if r["drug"] and "tremelimumab" in r["drug"].lower()
        ]
        assert len(imjudo_results) >= 1

        imjudo = imjudo_results[0]
        assert imjudo["matches"] is False
        assert imjudo["match_type"] == "excluded"

    def test_egfr_l858r_excluded_from_imjudo(self, parser, imjudo_label):
        """EGFR L858R in NSCLC should be EXCLUDED from IMJUDO."""
        indications = parser.parse_label(imjudo_label)

        results = match_variant_to_indications(indications, "EGFR", "L858R", "NSCLC")

        imjudo_results = [
            r for r in results if r["drug"] and "tremelimumab" in r["drug"].lower()
        ]
        assert len(imjudo_results) >= 1

        imjudo = imjudo_results[0]
        assert imjudo["matches"] is False
        assert imjudo["match_type"] == "excluded"


class TestFDALabelParserOsimertinib:
    """Tests for Osimertinib/TAGRISSO label parsing."""

    @pytest.fixture
    def parser(self):
        return FDALabelParser()

    @pytest.fixture
    def osimertinib_label(self):
        """Osimertinib label with T790M and L858R indications."""
        osimertinib_indication = """
        TAGRISSO is a kinase inhibitor indicated for:
        The first-line treatment of adult patients with metastatic NSCLC whose tumors have
        EGFR exon 19 deletions or exon 21 L858R mutations, as detected by an FDA-approved test.
        The treatment of adult patients with metastatic EGFR T790M mutation-positive NSCLC,
        as detected by an FDA-approved test, whose disease has progressed on or after EGFR TKI therapy.
        """
        return {
            "openfda": {
                "generic_name": ["osimertinib"],
                "brand_name": ["TAGRISSO"],
            },
            "indications_and_usage": [osimertinib_indication],
        }

    def test_parses_t790m_indication(self, parser, osimertinib_label):
        """Osimertinib should have T790M variant-level indication."""
        indications = parser.parse_label(osimertinib_label)

        egfr_indications = [i for i in indications if i.gene == "EGFR"]
        assert len(egfr_indications) >= 1

        # Find the T790M indication
        t790m_indications = [
            i
            for i in egfr_indications
            if i.specified_variants and "T790M" in i.specified_variants
        ]
        assert len(t790m_indications) >= 1

        t790m = t790m_indications[0]
        assert t790m.requirement == BiomarkerRequirement.REQUIRED_POSITIVE
        assert t790m.specificity == SpecificityLevel.VARIANT

    def test_parses_l858r_indication(self, parser, osimertinib_label):
        """Osimertinib should have L858R variant-level indication."""
        indications = parser.parse_label(osimertinib_label)

        egfr_indications = [i for i in indications if i.gene == "EGFR"]

        # Find the L858R indication
        l858r_indications = [
            i
            for i in egfr_indications
            if i.specified_variants and "L858R" in i.specified_variants
        ]
        assert len(l858r_indications) >= 1

        l858r = l858r_indications[0]
        assert l858r.requirement == BiomarkerRequirement.REQUIRED_POSITIVE

    def test_t790m_matches_osimertinib(self, parser, osimertinib_label):
        """EGFR T790M should match osimertinib with exact match."""
        indications = parser.parse_label(osimertinib_label)

        results = match_variant_to_indications(indications, "EGFR", "T790M", "NSCLC")

        # Should find osimertinib as a match
        osi_results = [
            r for r in results if r["drug"] and "osimertinib" in r["drug"].lower()
        ]
        assert len(osi_results) >= 1

        osi = osi_results[0]
        assert osi["matches"] is True
        assert osi["match_type"] == "exact"

    def test_l858r_matches_osimertinib(self, parser, osimertinib_label):
        """EGFR L858R should match osimertinib."""
        indications = parser.parse_label(osimertinib_label)

        results = match_variant_to_indications(indications, "EGFR", "L858R", "NSCLC")

        osi_results = [
            r for r in results if r["drug"] and "osimertinib" in r["drug"].lower()
        ]
        assert len(osi_results) >= 1

        osi = osi_results[0]
        assert osi["matches"] is True


class TestMatchVariantToIndications:
    """Tests for match_variant_to_indications function."""

    @pytest.fixture
    def parser(self):
        return FDALabelParser()

    def test_returns_empty_for_different_gene(self, parser):
        """Should not return results for a different gene."""
        label = {
            "openfda": {
                "generic_name": ["braftarget"],
                "brand_name": ["BRAFBRAND"],
            },
            "indications_and_usage": ["For BRAF V600E melanoma."],
        }

        indications = parser.parse_label(label)

        # Query for EGFR (not BRAF)
        results = match_variant_to_indications(indications, "EGFR", "L858R")

        # Should be empty - different gene
        assert len(results) == 0

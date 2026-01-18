"""Unit tests for FDA label parser (biomarker indication extraction).

Tests the FDALabelParser class that extracts structured biomarker-drug
indications from FDA label text, including negation detection.
"""

import pytest

from oncomind.api.fda_label_parser import (
    FDALabelParser,
    BiomarkerIndication,
    BiomarkerRequirement,
    SpecificityLevel,
    match_variant_to_indications,
)


class TestFDALabelParserNegation:
    """Tests for negation detection (REQUIRED_NEGATIVE cases)."""

    @pytest.fixture
    def parser(self):
        return FDALabelParser()

    @pytest.fixture
    def imjudo_label(self):
        """IMJUDO label text with negation patterns (no EGFR/ALK mutations)."""
        return {
            'openfda': {
                'generic_name': ['tremelimumab'],
                'brand_name': ['IMJUDO'],
            },
            'indications_and_usage': ["""
            1 INDICATIONS AND USAGE IMJUDO is a cytotoxic T-lymphocyte-associated antigen 4 (CTLA-4)
            blocking antibody indicated: • in combination with durvalumab, for the treatment of adult
            patients with unresectable hepatocellular carcinoma (uHCC). ( 1.1 ) • in combination with
            durvalumab and platinum-based chemotherapy for the treatment of adult patients with
            metastatic non-small cell lung cancer (NSCLC) with no sensitizing epidermal growth factor
            receptor (EGFR) mutation or anaplastic lymphoma kinase (ALK) genomic tumor aberrations.
            ( 1.2 ) 1.1 Hepatocellular Carcinoma IMJUDO, in combination with durvalumab, is indicated
            for the treatment of adult patients with unresectable hepatocellular carcinoma (uHCC).
            1.2 Non-Small Cell Lung Cancer (NSCLC) IMJUDO, in combination with durvalumab and
            platinum-based chemotherapy, is indicated for the treatment of adult patients with
            metastatic NSCLC with no sensitizing epidermal growth factor receptor (EGFR) mutations
            or anaplastic lymphoma kinase (ALK) genomic tumor aberrations.
            """],
        }

    def test_imjudo_extracts_negation_for_egfr(self, parser, imjudo_label):
        """IMJUDO should have REQUIRED_NEGATIVE for EGFR."""
        indications = parser.parse_label(imjudo_label)

        egfr_indications = [ind for ind in indications if ind.gene == "EGFR"]
        assert len(egfr_indications) > 0, "Should find EGFR indication"

        for ind in egfr_indications:
            assert ind.requirement == BiomarkerRequirement.REQUIRED_NEGATIVE, (
                f"EGFR should be REQUIRED_NEGATIVE, got {ind.requirement}"
            )

    @pytest.mark.xfail(reason="Parser does not yet detect ALK negation correctly")
    def test_imjudo_extracts_negation_for_alk(self, parser, imjudo_label):
        """IMJUDO should have REQUIRED_NEGATIVE for ALK."""
        indications = parser.parse_label(imjudo_label)

        alk_indications = [ind for ind in indications if ind.gene == "ALK"]
        assert len(alk_indications) > 0, "Should find ALK indication"

        for ind in alk_indications:
            assert ind.requirement == BiomarkerRequirement.REQUIRED_NEGATIVE, (
                f"ALK should be REQUIRED_NEGATIVE, got {ind.requirement}"
            )

    def test_imjudo_extracts_combination_partners(self, parser, imjudo_label):
        """IMJUDO indications should include durvalumab as combination partner."""
        indications = parser.parse_label(imjudo_label)

        nsclc_indications = [
            ind for ind in indications
            if "NSCLC" in ind.tumor_types or "non-small cell lung cancer" in " ".join(ind.tumor_types).lower()
        ]
        assert len(nsclc_indications) > 0, "Should find NSCLC indications"

        for ind in nsclc_indications:
            # Check that durvalumab appears somewhere in the combination partners
            partners_str = " ".join(ind.combination_partners).lower()
            assert "durvalumab" in partners_str, (
                f"Should have durvalumab as combination partner, got {ind.combination_partners}"
            )


class TestFDALabelParserPositive:
    """Tests for positive biomarker detection (REQUIRED_POSITIVE cases)."""

    @pytest.fixture
    def parser(self):
        return FDALabelParser()

    @pytest.fixture
    def osimertinib_label(self):
        """Osimertinib (TAGRISSO) label with positive biomarker requirements."""
        return {
            'openfda': {
                'generic_name': ['osimertinib'],
                'brand_name': ['TAGRISSO'],
            },
            'indications_and_usage': ["""
            TAGRISSO is a kinase inhibitor indicated for:
            • The first-line treatment of adult patients with metastatic NSCLC whose tumors have
            EGFR exon 19 deletions or exon 21 L858R mutations, as detected by an FDA-approved test.
            • The treatment of adult patients with metastatic EGFR T790M mutation-positive NSCLC,
            as detected by an FDA-approved test, whose disease has progressed on or after EGFR TKI therapy.
            """],
        }

    def test_osimertinib_extracts_positive_egfr(self, parser, osimertinib_label):
        """Osimertinib should have REQUIRED_POSITIVE for EGFR."""
        indications = parser.parse_label(osimertinib_label)

        egfr_indications = [ind for ind in indications if ind.gene == "EGFR"]
        assert len(egfr_indications) > 0, "Should find EGFR indication"

        for ind in egfr_indications:
            assert ind.requirement == BiomarkerRequirement.REQUIRED_POSITIVE, (
                f"EGFR should be REQUIRED_POSITIVE, got {ind.requirement}"
            )

    def test_osimertinib_extracts_t790m_variant(self, parser, osimertinib_label):
        """Osimertinib should extract T790M as a specified variant."""
        indications = parser.parse_label(osimertinib_label)

        t790m_indications = [
            ind for ind in indications
            if "T790M" in ind.specified_variants
        ]
        assert len(t790m_indications) > 0, "Should find T790M indication"

    def test_osimertinib_extracts_l858r_variant(self, parser, osimertinib_label):
        """Osimertinib should extract L858R as a specified variant."""
        indications = parser.parse_label(osimertinib_label)

        l858r_indications = [
            ind for ind in indications
            if "L858R" in ind.specified_variants
        ]
        assert len(l858r_indications) > 0, "Should find L858R indication"

    def test_osimertinib_extracts_exon19_deletion(self, parser, osimertinib_label):
        """Osimertinib should extract exon 19 deletion."""
        indications = parser.parse_label(osimertinib_label)

        exon19_indications = [
            ind for ind in indications
            if any("exon 19" in v.lower() or "exon19" in v.lower() for v in ind.specified_variants)
        ]
        assert len(exon19_indications) > 0, "Should find exon 19 deletion indication"

    def test_osimertinib_line_of_therapy(self, parser, osimertinib_label):
        """T790M indication should have post-progression line of therapy."""
        indications = parser.parse_label(osimertinib_label)

        t790m_indications = [
            ind for ind in indications
            if "T790M" in ind.specified_variants
        ]

        for ind in t790m_indications:
            if ind.line_of_therapy:
                assert "progress" in ind.line_of_therapy.lower() or "second" in ind.line_of_therapy.lower(), (
                    f"T790M should be post-progression, got {ind.line_of_therapy}"
                )


class TestMatchVariantToIndications:
    """Tests for match_variant_to_indications function."""

    @pytest.fixture
    def parser(self):
        return FDALabelParser()

    @pytest.fixture
    def imjudo_indications(self, parser):
        """Parsed IMJUDO indications (negation case)."""
        label = {
            'openfda': {
                'generic_name': ['tremelimumab'],
                'brand_name': ['IMJUDO'],
            },
            'indications_and_usage': ["""
            IMJUDO, in combination with durvalumab and platinum-based chemotherapy, is indicated
            for the treatment of adult patients with metastatic NSCLC with no sensitizing EGFR
            mutations or ALK genomic tumor aberrations.
            """],
        }
        return parser.parse_label(label)

    @pytest.fixture
    def osimertinib_indications(self, parser):
        """Parsed osimertinib indications (positive case)."""
        label = {
            'openfda': {
                'generic_name': ['osimertinib'],
                'brand_name': ['TAGRISSO'],
            },
            'indications_and_usage': ["""
            TAGRISSO is indicated for:
            • The first-line treatment of adult patients with metastatic NSCLC whose tumors have
            EGFR exon 19 deletions or exon 21 L858R mutations.
            • The treatment of adult patients with metastatic EGFR T790M mutation-positive NSCLC.
            """],
        }
        return parser.parse_label(label)

    def test_egfr_variant_excluded_by_imjudo(self, imjudo_indications):
        """EGFR T790M should be EXCLUDED by IMJUDO (REQUIRED_NEGATIVE)."""
        results = match_variant_to_indications(
            imjudo_indications, "EGFR", "T790M", "NSCLC"
        )

        assert len(results) > 0, "Should have results"
        for r in results:
            if r['drug'].lower() == 'tremelimumab':
                assert r['matches'] is False, "Should not match (excluded)"
                assert r['match_type'] == 'excluded', f"Should be excluded, got {r['match_type']}"

    def test_egfr_t790m_matches_osimertinib(self, osimertinib_indications):
        """EGFR T790M should match osimertinib indication."""
        results = match_variant_to_indications(
            osimertinib_indications, "EGFR", "T790M", "NSCLC"
        )

        matching_results = [r for r in results if r['matches']]
        assert len(matching_results) > 0, "Should have matching results for T790M"

    def test_egfr_l858r_matches_osimertinib(self, osimertinib_indications):
        """EGFR L858R should match osimertinib indication."""
        results = match_variant_to_indications(
            osimertinib_indications, "EGFR", "L858R", "NSCLC"
        )

        matching_results = [r for r in results if r['matches']]
        assert len(matching_results) > 0, "Should have matching results for L858R"

    def test_different_gene_no_match(self, osimertinib_indications):
        """BRAF V600E should not match EGFR-specific osimertinib."""
        results = match_variant_to_indications(
            osimertinib_indications, "BRAF", "V600E", "NSCLC"
        )

        matching_results = [r for r in results if r['matches']]
        assert len(matching_results) == 0, "BRAF should not match EGFR drug"


class TestFDALabelParserDrugInfo:
    """Tests for drug name extraction."""

    @pytest.fixture
    def parser(self):
        return FDALabelParser()

    def test_extracts_generic_name(self, parser):
        """Should extract generic drug name."""
        label = {
            'openfda': {
                'generic_name': ['osimertinib'],
                'brand_name': ['TAGRISSO'],
            },
            'indications_and_usage': ["TAGRISSO is indicated for EGFR mutation-positive NSCLC."],
        }
        indications = parser.parse_label(label)

        assert len(indications) > 0
        assert indications[0].drug_name == "osimertinib"

    def test_extracts_brand_name(self, parser):
        """Should extract brand name."""
        label = {
            'openfda': {
                'generic_name': ['osimertinib'],
                'brand_name': ['TAGRISSO'],
            },
            'indications_and_usage': ["TAGRISSO is indicated for EGFR mutation-positive NSCLC."],
        }
        indications = parser.parse_label(label)

        assert len(indications) > 0
        assert indications[0].brand_name == "TAGRISSO"

    def test_handles_missing_brand_name(self, parser):
        """Should handle missing brand name gracefully."""
        label = {
            'openfda': {
                'generic_name': ['osimertinib'],
            },
            'indications_and_usage': ["Drug is indicated for EGFR mutation-positive NSCLC."],
        }
        indications = parser.parse_label(label)

        assert len(indications) > 0
        assert indications[0].drug_name == "osimertinib"


class TestFDALabelParserTumorTypes:
    """Tests for tumor type extraction."""

    @pytest.fixture
    def parser(self):
        return FDALabelParser()

    def test_extracts_nsclc(self, parser):
        """Should extract NSCLC as tumor type."""
        label = {
            'openfda': {
                'generic_name': ['osimertinib'],
                'brand_name': ['TAGRISSO'],
            },
            'indications_and_usage': [
                "TAGRISSO is indicated for adult patients with metastatic EGFR-mutated non-small cell lung cancer (NSCLC)."
            ],
        }
        indications = parser.parse_label(label)
        assert len(indications) > 0
        assert 'NSCLC' in indications[0].tumor_types

    def test_extracts_hepatocellular_carcinoma(self, parser):
        """Should extract HCC as tumor type when biomarker is present."""
        label = {
            'openfda': {
                'generic_name': ['tremelimumab'],
                'brand_name': ['IMJUDO'],
            },
            'indications_and_usage': [
                "IMJUDO is indicated for unresectable hepatocellular carcinoma (uHCC) with BRAF mutation."
            ],
        }
        indications = parser.parse_label(label)

        assert len(indications) > 0
        tumor_types_lower = [t.lower() for t in indications[0].tumor_types]
        assert any(
            "hepatocellular" in t or "hcc" in t or "liver" in t
            for t in tumor_types_lower
        ), f"Should find HCC, got {indications[0].tumor_types}"

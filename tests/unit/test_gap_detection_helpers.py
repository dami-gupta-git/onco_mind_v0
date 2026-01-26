"""Unit tests for gap_detection/helpers.py.

Tests utility functions like MatchCounts, count_with_levels, normalize_source,
compute_overall_quality, compute_research_priority, etc.
"""

import pytest
from unittest.mock import MagicMock

from oncomind.insight_builder.gap_detection.helpers import (
    MatchCounts,
    count_with_levels,
    count_fda_match_levels,
    has_pathogenic_signal,
    normalize_source,
    sort_characterized_by_category,
    compute_overall_quality,
    compute_research_priority,
    get_primary_drug,
    get_top_sensitive_drugs,
    get_top_cooccurring_gene,
    has_strong_cooccurrence,
    get_top_codependencies,
)
from oncomind.models.extracted.evidence_gaps import (
    EvidenceGap, GapCategory, GapSeverity, CharacterizedAspect
)


# =============================================================================
# TEST MatchCounts
# =============================================================================

class TestMatchCounts:
    """Tests for the MatchCounts dataclass."""

    def test_default_values(self):
        """Default values should be zero or empty."""
        counts = MatchCounts()
        assert counts.total == 0
        assert counts.variant == 0
        assert counts.codon == 0
        assert counts.gene == 0
        assert counts.tumor == 0
        assert counts.other_cancers == set()

    def test_matches_on_str_empty(self):
        """Empty counts should return '0 matches'."""
        counts = MatchCounts()
        assert counts.matches_on_str == "0 matches"

    def test_matches_on_str_single_level(self):
        """Single level counts should show correctly."""
        counts = MatchCounts(variant=3)
        assert counts.matches_on_str == "3 variant"

    def test_matches_on_str_multiple_levels(self):
        """Multiple level counts should be comma-separated."""
        counts = MatchCounts(variant=2, codon=1, gene=3)
        assert "2 variant" in counts.matches_on_str
        assert "1 codon" in counts.matches_on_str
        assert "3 gene" in counts.matches_on_str

    def test_tumor_breakdown_str(self):
        """Tumor breakdown should show tumor vs other counts."""
        counts = MatchCounts(total=5, tumor=3)
        breakdown = counts.tumor_breakdown_str
        assert "3 tumor" in breakdown
        assert "2 other" in breakdown

    def test_tumor_breakdown_str_no_other(self):
        """Tumor breakdown with no other cancers."""
        counts = MatchCounts(total=3, tumor=3)
        breakdown = counts.tumor_breakdown_str
        assert breakdown == "3 tumor"

    def test_add_method(self):
        """Add method should combine counts correctly."""
        counts1 = MatchCounts(total=3, variant=2, gene=1, tumor=2)
        counts2 = MatchCounts(total=2, variant=1, codon=1, tumor=1)
        counts2.other_cancers.add("Melanoma")

        result = counts1.add(counts2)

        assert result.total == 5
        assert result.variant == 3
        assert result.codon == 1
        assert result.gene == 1
        assert result.tumor == 3
        assert "Melanoma" in result.other_cancers
        # add() should return self for chaining
        assert result is counts1


# =============================================================================
# TEST count_with_levels
# =============================================================================

class TestCountWithLevels:
    """Tests for count_with_levels function."""

    def test_empty_list(self):
        """Empty list should return empty counts."""
        counts = count_with_levels([])
        assert counts.total == 0

    def test_counts_by_locus_match(self):
        """Should count items by their locus_match property."""
        item_variant = MagicMock()
        item_variant.locus_match = "variant"
        item_variant.tumor_type = None

        item_codon = MagicMock()
        item_codon.locus_match = "codon"
        item_codon.tumor_type = None

        item_gene = MagicMock()
        item_gene.locus_match = "gene"
        item_gene.tumor_type = None

        counts = count_with_levels([item_variant, item_codon, item_gene, item_gene])

        assert counts.total == 4
        assert counts.variant == 1
        assert counts.codon == 1
        assert counts.gene == 2

    def test_counts_tumor_matches_via_property(self):
        """Should use tumor_match property if available."""
        item1 = MagicMock()
        item1.locus_match = "variant"
        item1.tumor_type = "NSCLC"
        item1.tumor_match = True  # Computed property

        item2 = MagicMock()
        item2.locus_match = "gene"
        item2.tumor_type = "Melanoma"
        item2.tumor_match = False  # Different tumor

        counts = count_with_levels([item1, item2], tumor_type="NSCLC")

        assert counts.tumor == 1
        assert "Melanoma" in counts.other_cancers

    def test_tracks_other_cancers(self):
        """Should track non-matching cancer types."""
        item1 = MagicMock()
        item1.locus_match = "variant"
        item1.tumor_type = "Melanoma"
        item1.tumor_match = False

        item2 = MagicMock()
        item2.locus_match = "gene"
        item2.disease = "Colon Cancer"
        item2.tumor_type = None
        item2.tumor_match = False

        counts = count_with_levels([item1, item2], tumor_type="NSCLC")

        assert "Melanoma" in counts.other_cancers


# =============================================================================
# TEST normalize_source
# =============================================================================

class TestNormalizeSource:
    """Tests for normalize_source function."""

    def test_civic_normalization(self):
        """CIViC variants should normalize."""
        assert normalize_source("VICC/civic") == "CIViC"
        assert normalize_source("vicc/civic") == "CIViC"
        assert normalize_source("civic") == "CIViC"

    def test_cgi_normalization(self):
        """CGI variants should normalize."""
        assert normalize_source("cgi") == "CGI"
        assert normalize_source("VICC/cgi") == "CGI"

    def test_oncokb_normalization(self):
        """OncoKB variants should normalize."""
        assert normalize_source("VICC/oncokb") == "OncoKB"
        assert normalize_source("oncokb") == "OncoKB"

    def test_vicc_prefix_title_cases(self):
        """VICC/ prefix sources should be title-cased."""
        # Unknown VICC sources get title-cased
        assert normalize_source("VICC/moa") == "Moa"
        assert normalize_source("VICC/pmkb") == "Pmkb"

    def test_unknown_source_passthrough(self):
        """Unknown sources should pass through unchanged."""
        assert normalize_source("SomeNewSource") == "SomeNewSource"


# =============================================================================
# TEST sort_characterized_by_category
# =============================================================================

class TestSortCharacterizedByCategory:
    """Tests for sort_characterized_by_category function."""

    def test_sorts_by_category_order(self):
        """Should sort by clinical importance (CLINICAL first)."""
        items = [
            CharacterizedAspect(aspect="Functional", basis="test", category=GapCategory.FUNCTIONAL),
            CharacterizedAspect(aspect="Clinical", basis="test", category=GapCategory.CLINICAL),
            CharacterizedAspect(aspect="Drug", basis="test", category=GapCategory.DRUG_RESPONSE),
        ]

        sorted_items = sort_characterized_by_category(items)

        assert sorted_items[0].category == GapCategory.CLINICAL
        assert sorted_items[1].category == GapCategory.DRUG_RESPONSE
        assert sorted_items[2].category == GapCategory.FUNCTIONAL

    def test_none_category_goes_last(self):
        """Items with no category should sort last."""
        items = [
            CharacterizedAspect(aspect="No Cat", basis="test", category=None),
            CharacterizedAspect(aspect="Clinical", basis="test", category=GapCategory.CLINICAL),
        ]

        sorted_items = sort_characterized_by_category(items)

        assert sorted_items[0].category == GapCategory.CLINICAL
        assert sorted_items[1].category is None


# =============================================================================
# TEST compute_overall_quality
# =============================================================================

class TestComputeOverallQuality:
    """Tests for compute_overall_quality function."""

    def test_no_gaps_comprehensive(self):
        """No gaps should return comprehensive."""
        result = compute_overall_quality([], 5)
        assert result == "comprehensive"

    def test_critical_gap_caps_at_limited(self):
        """CRITICAL gap should cap quality at 'limited'."""
        critical_gap = EvidenceGap(
            category=GapCategory.VALIDATION,
            severity=GapSeverity.CRITICAL,
            description="Critical gap",
        )
        result = compute_overall_quality([critical_gap], 20)
        assert result == "limited"

    def test_high_gap_caps_at_moderate(self):
        """HIGH gap should cap quality at 'moderate'."""
        high_gap = EvidenceGap(
            category=GapCategory.DISCORDANT,
            severity=GapSeverity.HIGH,
            description="High gap",
        )
        result = compute_overall_quality([high_gap], 20)
        assert result == "moderate"

    def test_significant_gap_caps_at_moderate(self):
        """SIGNIFICANT gap should cap quality at 'moderate'."""
        significant_gap = EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.SIGNIFICANT,
            description="Significant gap",
        )
        result = compute_overall_quality([significant_gap], 20)
        assert result == "moderate"

    def test_many_well_characterized_improves_quality(self):
        """Many well-characterized aspects should improve quality."""
        minor_gap = EvidenceGap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.MINOR,
            description="Minor gap",
        )
        # With 0 well-characterized
        result_none = compute_overall_quality([minor_gap], 0)
        # With 10 well-characterized
        result_many = compute_overall_quality([minor_gap], 10)

        quality_order = ["comprehensive", "moderate", "limited", "minimal"]
        assert quality_order.index(result_many) <= quality_order.index(result_none)


# =============================================================================
# TEST compute_research_priority
# =============================================================================

class TestComputeResearchPriority:
    """Tests for compute_research_priority function."""

    def test_pathogenic_cancer_gene_very_high_priority(self):
        """Pathogenic signal in cancer gene with critical gap should be very_high priority."""
        critical_gap = EvidenceGap(
            category=GapCategory.VALIDATION,
            severity=GapSeverity.CRITICAL,
            description="Critical gap",
        )
        mock_evidence = MagicMock()
        mock_evidence.identifiers.gene = "BRAF"
        mock_evidence.identifiers.variant = "V600E"  # A hotspot variant

        result = compute_research_priority(
            evidence=mock_evidence,
            gaps=[critical_gap],
            overall_quality="limited",
            is_cancer_gene=True,
            has_pathogenic=True,
        )

        # BRAF V600E is a hotspot, so with critical gap it's very_high priority
        assert result == "very_high"

    def test_minimal_quality_with_gaps_high_priority(self):
        """Minimal quality with critical gaps should be high priority."""
        critical_gap = EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.CRITICAL,
            description="Critical gap",
        )
        mock_evidence = MagicMock()
        mock_evidence.identifiers.gene = "UNKNOWNGENE"
        mock_evidence.identifiers.variant = "X999Y"  # Not a hotspot

        result = compute_research_priority(
            evidence=mock_evidence,
            gaps=[critical_gap],
            overall_quality="minimal",
            is_cancer_gene=True,
            has_pathogenic=False,
        )

        # With critical gaps, should be at least high priority
        assert result in ("high", "very_high")

    def test_comprehensive_quality_low_priority(self):
        """Comprehensive quality should be low priority."""
        mock_evidence = MagicMock()
        mock_evidence.identifiers.gene = "BRAF"
        mock_evidence.identifiers.variant = "V600E"  # Must be a string for regex

        result = compute_research_priority(
            evidence=mock_evidence,
            gaps=[],
            overall_quality="comprehensive",
            is_cancer_gene=True,
            has_pathogenic=True,
        )

        assert result == "low"


# =============================================================================
# TEST has_pathogenic_signal
# =============================================================================

class TestHasPathogenicSignal:
    """Tests for has_pathogenic_signal function."""

    @pytest.fixture
    def mock_evidence(self):
        """Create mock evidence with default values."""
        evidence = MagicMock()
        evidence.functional.alphamissense_prediction = None
        evidence.functional.cadd_score = None
        evidence.functional.polyphen2_prediction = None
        evidence.functional.sift_prediction = None
        evidence.functional.snpeff_effect = None
        evidence.civic_assertions = []
        evidence.fda_biomarker_evidence = []
        evidence.clinvar_entries = []
        evidence.clinvar_significance = None
        return evidence

    def test_alphamissense_pathogenic(self, mock_evidence):
        """AlphaMissense pathogenic should return True."""
        mock_evidence.functional.alphamissense_prediction = "pathogenic"
        assert has_pathogenic_signal(mock_evidence) is True

    def test_alphamissense_likely_pathogenic(self, mock_evidence):
        """AlphaMissense likely_pathogenic should return True."""
        mock_evidence.functional.alphamissense_prediction = "likely_pathogenic"
        assert has_pathogenic_signal(mock_evidence) is True

    def test_high_cadd_score(self, mock_evidence):
        """High CADD score (>=20) should return True."""
        mock_evidence.functional.cadd_score = 25.0
        assert has_pathogenic_signal(mock_evidence) is True

    def test_low_cadd_score(self, mock_evidence):
        """Low CADD score (<20) should return False."""
        mock_evidence.functional.cadd_score = 15.0
        assert has_pathogenic_signal(mock_evidence) is False

    def test_polyphen_probably_damaging(self, mock_evidence):
        """PolyPhen2 probably_damaging should return True."""
        mock_evidence.functional.polyphen2_prediction = "probably_damaging"
        assert has_pathogenic_signal(mock_evidence) is True

    def test_sift_deleterious(self, mock_evidence):
        """SIFT deleterious should return True."""
        mock_evidence.functional.sift_prediction = "deleterious"
        assert has_pathogenic_signal(mock_evidence) is True

    def test_stop_gained(self, mock_evidence):
        """Stop gained variant should return True."""
        mock_evidence.functional.snpeff_effect = "stop_gained"
        assert has_pathogenic_signal(mock_evidence) is True

    def test_frameshift(self, mock_evidence):
        """Frameshift variant should return True."""
        mock_evidence.functional.snpeff_effect = "frameshift_variant"
        assert has_pathogenic_signal(mock_evidence) is True

    def test_civic_assertions(self, mock_evidence):
        """CIViC assertions should return True."""
        mock_evidence.civic_assertions = [MagicMock()]
        assert has_pathogenic_signal(mock_evidence) is True

    def test_fda_biomarker_evidence(self, mock_evidence):
        """FDA biomarker evidence should return True."""
        mock_evidence.fda_biomarker_evidence = [MagicMock()]
        assert has_pathogenic_signal(mock_evidence) is True

    def test_clinvar_pathogenic(self, mock_evidence):
        """ClinVar pathogenic entry should return True."""
        entry = MagicMock()
        entry.clinical_significance = "Pathogenic"
        mock_evidence.clinvar_entries = [entry]
        assert has_pathogenic_signal(mock_evidence) is True

    def test_no_signals(self, mock_evidence):
        """No pathogenic signals should return False."""
        assert has_pathogenic_signal(mock_evidence) is False


# =============================================================================
# TEST get_primary_drug
# =============================================================================

class TestGetPrimaryDrug:
    """Tests for get_primary_drug function."""

    def test_returns_fda_drug_first(self):
        """Should return FDA-approved drug first."""
        from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement

        fda_ev = MagicMock()
        fda_ev.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        fda_ev.drug_name = "Vemurafenib"
        fda_ev.brand_name = "Zelboraf"

        evidence = MagicMock()
        evidence.fda_biomarker_evidence = [fda_ev]
        evidence.civic_assertions = []

        result = get_primary_drug(evidence)
        assert result == "Zelboraf"

    def test_falls_back_to_civic_therapy(self):
        """Should fall back to CIViC therapy if no FDA."""
        from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement

        civic_assertion = MagicMock()
        civic_assertion.therapies = ["Trametinib", "Dabrafenib"]

        evidence = MagicMock()
        evidence.fda_biomarker_evidence = []
        evidence.civic_assertions = [civic_assertion]

        result = get_primary_drug(evidence)
        assert result == "Trametinib"

    def test_returns_none_if_no_drugs(self):
        """Should return None if no drug data."""
        evidence = MagicMock()
        evidence.fda_biomarker_evidence = []
        evidence.civic_assertions = []

        result = get_primary_drug(evidence)
        assert result is None


# =============================================================================
# TEST get_top_sensitive_drugs
# =============================================================================

class TestGetTopSensitiveDrugs:
    """Tests for get_top_sensitive_drugs function."""

    def test_returns_empty_if_no_depmap(self):
        """Should return empty list if no DepMap evidence."""
        evidence = MagicMock()
        evidence.depmap_evidence = None

        result = get_top_sensitive_drugs(evidence)
        assert result == []

    def test_returns_empty_if_no_tumor_matched_models(self):
        """Should return empty if no tumor-matched cell line models."""
        cell_line = MagicMock()
        cell_line.has_mutation = True
        cell_line.primary_disease = "Melanoma"

        drug_sens = MagicMock()
        drug_sens.drug_name = "Trametinib"

        depmap = MagicMock()
        depmap.drug_sensitivities = [drug_sens]
        depmap.cell_line_models = [cell_line]
        depmap.get_top_sensitive_drugs.return_value = [drug_sens]

        evidence = MagicMock()
        evidence.depmap_evidence = depmap

        # Tumor type doesn't match
        result = get_top_sensitive_drugs(evidence, tumor_type="NSCLC")
        assert result == []


# =============================================================================
# TEST get_top_cooccurring_gene
# =============================================================================

class TestGetTopCooccurringGene:
    """Tests for get_top_cooccurring_gene function."""

    def test_returns_top_cooccurring_gene(self):
        """Should return top co-occurring gene from cBioPortal."""
        cooc = MagicMock()
        cooc.gene = "TP53"

        cbioportal = MagicMock()
        cbioportal.co_occurring = [cooc]

        evidence = MagicMock()
        evidence.cbioportal_evidence = cbioportal

        result = get_top_cooccurring_gene(evidence)
        assert result == "TP53"

    def test_returns_none_if_no_cbioportal(self):
        """Should return None if no cBioPortal evidence."""
        evidence = MagicMock()
        evidence.cbioportal_evidence = None

        result = get_top_cooccurring_gene(evidence)
        assert result is None


# =============================================================================
# TEST has_strong_cooccurrence
# =============================================================================

class TestHasStrongCooccurrence:
    """Tests for has_strong_cooccurrence function."""

    def test_returns_true_above_threshold(self):
        """Should return True if co-occurrence >= threshold."""
        cooc = MagicMock()
        cooc.pct = 40.0

        cbioportal = MagicMock()
        cbioportal.co_occurring = [cooc]

        evidence = MagicMock()
        evidence.cbioportal_evidence = cbioportal

        result = has_strong_cooccurrence(evidence, threshold_pct=30.0)
        assert result is True

    def test_returns_false_below_threshold(self):
        """Should return False if co-occurrence < threshold."""
        cooc = MagicMock()
        cooc.pct = 20.0

        cbioportal = MagicMock()
        cbioportal.co_occurring = [cooc]

        evidence = MagicMock()
        evidence.cbioportal_evidence = cbioportal

        result = has_strong_cooccurrence(evidence, threshold_pct=30.0)
        assert result is False


# =============================================================================
# TEST get_top_codependencies
# =============================================================================

class TestGetTopCodependencies:
    """Tests for get_top_codependencies function."""

    def test_returns_top_codependent_genes(self):
        """Should return top 3 co-dependent genes from DepMap."""
        codep1 = MagicMock()
        codep1.gene = "MDM2"
        codep2 = MagicMock()
        codep2.gene = "BCL2"
        codep3 = MagicMock()
        codep3.gene = "MCL1"
        codep4 = MagicMock()
        codep4.gene = "XIAP"

        depmap = MagicMock()
        depmap.co_dependencies = [codep1, codep2, codep3, codep4]

        evidence = MagicMock()
        evidence.depmap_evidence = depmap

        result = get_top_codependencies(evidence)
        assert len(result) == 3
        assert result == ["MDM2", "BCL2", "MCL1"]

    def test_returns_empty_if_no_depmap(self):
        """Should return empty list if no DepMap evidence."""
        evidence = MagicMock()
        evidence.depmap_evidence = None

        result = get_top_codependencies(evidence)
        assert result == []

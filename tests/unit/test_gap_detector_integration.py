"""Integration tests for detect_evidence_gaps and helper functions.

Tests: has_pathogenic_signal, compute_overall_quality,
and the detect_evidence_gaps integration function.
"""

from unittest.mock import MagicMock

from oncomind.insight_builder.gap_detector import (
    GapDetectionContext,
    detect_evidence_gaps,
)
from oncomind.insight_builder.gap_detection.helpers import (
    has_pathogenic_signal as _has_pathogenic_signal,
    compute_overall_quality as _compute_overall_quality,
)
from oncomind.models.extracted.evidence_gaps import (
    GapCategory,
    GapSeverity,
    EvidenceGap,
)

# =============================================================================
# TEST _has_pathogenic_signal
# =============================================================================


class TestHasPathogenicSignal:
    """Tests for _has_pathogenic_signal helper function."""

    def test_alphamissense_pathogenic(self, mock_evidence):
        """AlphaMissense pathogenic prediction should return True."""
        mock_evidence.functional.alphamissense_prediction = "pathogenic"

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_cadd_high_score(self, mock_evidence):
        """High CADD score should return True."""
        mock_evidence.functional.cadd_score = 25.0

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_cadd_low_score(self, mock_evidence):
        """Low CADD score should return False."""
        mock_evidence.functional.alphamissense_prediction = None
        mock_evidence.functional.polyphen2_prediction = None
        mock_evidence.functional.sift_prediction = None
        mock_evidence.functional.snpeff_effect = None
        mock_evidence.functional.cadd_score = 15.0
        mock_evidence.civic_assertions = []
        mock_evidence.fda_biomarker_evidence = []
        mock_evidence.clinvar_entries = []
        mock_evidence.clinvar_significance = None

        assert _has_pathogenic_signal(mock_evidence) is False

    def test_polyphen2_damaging(self, mock_evidence):
        """PolyPhen2 damaging should return True."""
        mock_evidence.functional.polyphen2_prediction = "probably_damaging"

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_civic_assertions(self, mock_evidence):
        """CIViC assertions should return True."""
        mock_evidence.civic_assertions = [MagicMock()]

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_fda_biomarker_evidence(self, mock_evidence):
        """FDA biomarker evidence should return True."""
        mock_evidence.fda_biomarker_evidence = [MagicMock()]

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_clinvar_pathogenic(self, mock_evidence):
        """ClinVar pathogenic entry should return True."""
        entry = MagicMock()
        entry.clinical_significance = "Pathogenic"
        mock_evidence.clinvar_entries = [entry]

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_truncating_variant(self, mock_evidence):
        """Truncating variant should return True."""
        mock_evidence.functional.snpeff_effect = "stop_gained"

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_frameshift(self, mock_evidence):
        """Frameshift variant should return True."""
        mock_evidence.functional.snpeff_effect = "frameshift_variant"

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_no_signal(self, mock_evidence):
        """No pathogenic signals should return False."""
        mock_evidence.functional.alphamissense_prediction = None
        mock_evidence.functional.alphamissense_score = None
        mock_evidence.functional.cadd_score = None
        mock_evidence.functional.polyphen2_prediction = None
        mock_evidence.functional.sift_prediction = None
        mock_evidence.functional.snpeff_effect = None
        mock_evidence.civic_assertions = []
        mock_evidence.fda_biomarker_evidence = []
        mock_evidence.clinvar_entries = []
        mock_evidence.clinvar_significance = None

        assert _has_pathogenic_signal(mock_evidence) is False


# =============================================================================
# TEST detect_evidence_gaps (integration of all checks)
# =============================================================================


class TestDetectEvidenceGapsIntegration:
    """Integration tests for the main detect_evidence_gaps function."""

    def test_well_studied_variant_few_gaps(self, mock_evidence):
        """Well-studied variant should have few gaps."""
        mock_evidence.identifiers.gene = "BRAF"
        mock_evidence.identifiers.variant = "V600E"
        mock_evidence.context.tumor_type = "Melanoma"
        mock_evidence.context.gene_role = "oncogene"
        mock_evidence.functional.alphamissense_score = 0.98
        mock_evidence.functional.alphamissense_prediction = "pathogenic"
        mock_evidence.functional.cadd_score = 30.0

        from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement

        civic_assertion = MagicMock()
        civic_assertion.disease = "Melanoma"
        civic_assertion.locus_variant_match = None
        mock_evidence.civic_assertions = [civic_assertion]

        fda_ev = MagicMock()
        fda_ev.tumor_types = ["Melanoma"]
        fda_ev.locus_variant_match = None
        fda_ev.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        mock_evidence.fda_biomarker_evidence = [fda_ev]
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[fda_ev])

        cgi_biomarker = MagicMock()
        cgi_biomarker.tumor_type = "Melanoma"
        cgi_biomarker.locus_variant_match = None
        mock_evidence.cgi_biomarkers = [cgi_biomarker]

        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = [MagicMock() for _ in range(20)]

        gaps = detect_evidence_gaps(mock_evidence)

        assert gaps.overall_evidence_quality in ("comprehensive", "moderate")
        assert len(gaps.well_characterized) > len(gaps.gaps)

    def test_unknown_variant_many_gaps(self, mock_evidence):
        """Unknown variant should have many gaps."""
        mock_evidence.identifiers.gene = "UNKNOWNGENE"
        mock_evidence.identifiers.variant = "X999Y"
        mock_evidence.literature_searched = True
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        mock_evidence.functional.alphamissense_score = None
        mock_evidence.functional.alphamissense_prediction = None
        mock_evidence.functional.cadd_score = None
        mock_evidence.functional.polyphen2_prediction = None
        mock_evidence.functional.sift_prediction = None
        mock_evidence.context.gene_role = None
        mock_evidence.context.pathway = None
        mock_evidence.clinvar_significance = None

        gaps = detect_evidence_gaps(mock_evidence)

        assert gaps.overall_evidence_quality in ("limited", "minimal")
        assert len(gaps.gaps) > 0

    def test_returns_evidence_gaps_object(self, mock_evidence):
        """Should return EvidenceGaps object with all fields."""
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])
        gaps = detect_evidence_gaps(mock_evidence)

        assert hasattr(gaps, "gaps")
        assert hasattr(gaps, "overall_evidence_quality")
        assert hasattr(gaps, "well_characterized")
        assert hasattr(gaps, "well_characterized_detailed")
        assert hasattr(gaps, "poorly_characterized")
        assert hasattr(gaps, "research_priority")


# =============================================================================
# TEST _compute_overall_quality
# =============================================================================


class TestComputeOverallQuality:
    """Tests for _compute_overall_quality function."""

    def test_no_gaps_comprehensive(self):
        """No gaps should return comprehensive."""
        result = _compute_overall_quality([], 5)
        assert result == "comprehensive"

    def test_well_characterized_offsets_gaps(self):
        """Well-characterized aspects should offset gap penalties."""
        minor_gap = EvidenceGap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.MINOR,
            description="Test gap",
        )
        result_no_credit = _compute_overall_quality([minor_gap], 0)
        result_with_credit = _compute_overall_quality([minor_gap], 5)

        quality_order = ["comprehensive", "moderate", "limited", "minimal"]
        assert quality_order.index(result_with_credit) <= quality_order.index(
            result_no_credit
        )

    def test_many_well_characterized_with_minor_gap_yields_comprehensive(self):
        """Many well-characterized aspects should yield comprehensive with only MINOR gaps."""
        minor_gap = EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.MINOR,
            description="Test gap",
        )
        result = _compute_overall_quality([minor_gap], 10)
        assert result == "comprehensive"

    def test_significant_gap_caps_quality_at_moderate(self):
        """SIGNIFICANT gaps should cap quality at 'moderate' regardless of well-characterized count."""
        significant_gap = EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.SIGNIFICANT,
            description="Test gap",
        )
        result = _compute_overall_quality([significant_gap], 10)
        assert result == "moderate"

    def test_critical_gaps_cap_quality_at_limited(self):
        """CRITICAL gaps should cap quality at 'limited' regardless of well-characterized count."""
        critical_gap = EvidenceGap(
            category=GapCategory.VALIDATION,
            severity=GapSeverity.CRITICAL,
            description="Critical validation gap",
        )
        result = _compute_overall_quality([critical_gap], 20)
        assert result == "limited"

    def test_high_gaps_cap_quality_at_moderate(self):
        """HIGH gaps should cap quality at 'moderate' regardless of well-characterized count."""
        high_gap = EvidenceGap(
            category=GapCategory.DISCORDANT,
            severity=GapSeverity.HIGH,
            description="High severity gap",
        )
        result = _compute_overall_quality([high_gap], 20)
        assert result == "moderate"

    def test_many_gaps_yields_poor_quality(self):
        """Many gaps should yield poor quality even with some well-characterized."""
        gaps = [
            EvidenceGap(
                category=GapCategory.CLINICAL,
                severity=GapSeverity.CRITICAL,
                description="Gap 1",
            ),
            EvidenceGap(
                category=GapCategory.FUNCTIONAL,
                severity=GapSeverity.SIGNIFICANT,
                description="Gap 2",
            ),
            EvidenceGap(
                category=GapCategory.PRECLINICAL,
                severity=GapSeverity.SIGNIFICANT,
                description="Gap 3",
            ),
            EvidenceGap(
                category=GapCategory.VALIDATION,
                severity=GapSeverity.CRITICAL,
                description="Gap 4",
            ),
        ]
        result = _compute_overall_quality(gaps, 3)
        assert result in ("limited", "minimal")

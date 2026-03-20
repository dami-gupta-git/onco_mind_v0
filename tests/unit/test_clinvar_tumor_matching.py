"""Unit tests for ClinVar tumor type matching.

Tests the annotation of cancer_type_match on ClinVarEvidence entries,
and the tumor_match filter in discordance detection.
"""

import pytest
from unittest.mock import MagicMock

from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.base import EvidenceLevel, compute_cancer_specificity
from oncomind.insight_builder.gap_detection.discordance import (
    detect_discordant_evidence_internal,
)


# =============================================================================
# Helpers
# =============================================================================


def make_clinvar_entry(
    clinical_significance: str,
    locus_match: str = "variant",
    conditions: list[str] | None = None,
    tumor_match: bool | None = None,
) -> MagicMock:
    """Create a mock ClinVarEvidence entry."""
    entry = MagicMock()
    entry.clinical_significance = clinical_significance
    entry.locus_match = locus_match
    entry.conditions = conditions or []
    entry.tumor_match = tumor_match
    return entry


def annotate_clinvar_entries(
    entries: list[ClinVarEvidence], resolved_tumor: str | None
) -> None:
    """Apply cancer_type_match annotation — mirrors the aggregator loop."""
    for entry in entries:
        if entry.conditions:
            best: str | None = None
            for condition in entry.conditions:
                specificity = compute_cancer_specificity(condition, resolved_tumor)
                if specificity == "cancer_specific":
                    best = "cancer_specific"
                    break
                if specificity == "pan_cancer" and best != "cancer_specific":
                    best = "pan_cancer"
                elif best is None:
                    best = specificity
            if best is not None:
                entry.cancer_type_match = EvidenceLevel(level=best)


# =============================================================================
# Tests: cancer_type_match annotation
# =============================================================================


class TestClinVarTumorAnnotation:
    """Tests for the ClinVar cancer_type_match annotation logic."""

    def test_exact_match_sets_cancer_specific(self):
        """Condition matching queried tumor sets cancer_specific."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Lung adenocarcinoma"],
        )
        annotate_clinvar_entries([entry], "NSCLC")
        assert entry.tumor_match is True

    def test_no_match_sets_other_cancer(self):
        """Condition not matching queried tumor sets the condition name, not cancer_specific."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Hereditary breast ovarian cancer syndrome"],
        )
        annotate_clinvar_entries([entry], "NSCLC")
        assert entry.tumor_match is False

    def test_pan_cancer_condition(self):
        """Pan-cancer condition (e.g. 'Cancer') sets pan_cancer level, tumor_match is False."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Cancer"],
        )
        annotate_clinvar_entries([entry], "NSCLC")
        assert entry.cancer_type_match is not None
        assert entry.cancer_type_match.level == "pan_cancer"
        # pan_cancer is not cancer_specific, so tumor_match is False (not None)
        assert entry.tumor_match is False

    def test_multiple_conditions_one_matches(self):
        """If any condition matches, entry is cancer_specific."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Hereditary breast ovarian cancer syndrome", "Lung adenocarcinoma"],
        )
        annotate_clinvar_entries([entry], "NSCLC")
        assert entry.tumor_match is True

    def test_multiple_conditions_none_match(self):
        """If no condition matches, entry is not cancer_specific."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Hereditary breast ovarian cancer syndrome", "Li-Fraumeni syndrome"],
        )
        annotate_clinvar_entries([entry], "NSCLC")
        assert entry.tumor_match is False

    def test_cancer_specific_takes_priority_over_pan_cancer(self):
        """cancer_specific takes priority over pan_cancer when both are present."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Cancer", "Lung adenocarcinoma"],
        )
        annotate_clinvar_entries([entry], "NSCLC")
        assert entry.tumor_match is True
        assert entry.cancer_type_match.level == "cancer_specific"

    def test_no_conditions_leaves_cancer_type_match_unset(self):
        """Entry with no conditions leaves cancer_type_match as None."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=[],
        )
        annotate_clinvar_entries([entry], "NSCLC")
        assert entry.cancer_type_match is None
        assert entry.tumor_match is None

    def test_no_tumor_type_leaves_cancer_type_match_unset(self):
        """No resolved tumor type: cancer_type_match reflects condition without a match."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Lung adenocarcinoma"],
        )
        annotate_clinvar_entries([entry], None)
        # compute_cancer_specificity with no queried_tumor returns the source disease name
        assert entry.cancer_type_match is not None
        assert entry.cancer_type_match.level == "Lung adenocarcinoma"
        assert entry.tumor_match is False

    def test_tumor_type_mappings_alias_match(self):
        """TUMOR_TYPE_MAPPINGS aliases resolve correctly (e.g. NSCLC -> lung adenocarcinoma)."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Non-small cell lung carcinoma"],
        )
        annotate_clinvar_entries([entry], "NSCLC")
        assert entry.tumor_match is True


# =============================================================================
# Tests: discordance filter respects tumor_match
# =============================================================================


class TestDiscordanceClinVarTumorFilter:
    """Tests that discordance detection skips ClinVar entries irrelevant to the tumor."""

    def _make_evidence(self, entries: list) -> MagicMock:
        evidence = MagicMock()
        evidence.identifiers.gene = "CHEK2"
        evidence.identifiers.variant = "I157T"
        evidence.fda_biomarker_evidence = []
        evidence.cgi_biomarkers = []
        evidence.preclinical_biomarkers = []
        evidence.early_phase_biomarkers = []
        evidence.vicc_evidence = []
        evidence.civic_assertions = []
        evidence.civic_evidence = []
        evidence.clinvar_entries = entries
        return evidence

    def test_tumor_mismatched_entries_do_not_cause_conflict(self):
        """Pathogenic/benign conflict from other-tumor entries should not be flagged."""
        pathogenic = make_clinvar_entry(
            "Pathogenic", locus_match="variant", tumor_match=False
        )
        benign = make_clinvar_entry(
            "Benign", locus_match="variant", tumor_match=False
        )
        evidence = self._make_evidence([pathogenic, benign])

        high, significant = detect_discordant_evidence_internal(evidence)

        assert len(high) == 0
        assert len(significant) == 0

    def test_tumor_matched_entries_still_flag_conflict(self):
        """Pathogenic/benign conflict from matching-tumor entries is still flagged."""
        pathogenic = make_clinvar_entry(
            "Pathogenic", locus_match="variant", tumor_match=True
        )
        benign = make_clinvar_entry(
            "Benign", locus_match="variant", tumor_match=True
        )
        evidence = self._make_evidence([pathogenic, benign])

        high, significant = detect_discordant_evidence_internal(evidence)

        assert len(high) == 1
        assert "pathogenic" in high[0].lower()
        assert "benign" in high[0].lower()

    def test_none_tumor_match_entries_still_flag_conflict(self):
        """Entries with tumor_match=None (no conditions) are not filtered out."""
        pathogenic = make_clinvar_entry(
            "Pathogenic", locus_match="variant", tumor_match=None
        )
        benign = make_clinvar_entry(
            "Benign", locus_match="variant", tumor_match=None
        )
        evidence = self._make_evidence([pathogenic, benign])

        high, significant = detect_discordant_evidence_internal(evidence)

        assert len(high) == 1

    def test_mixed_tumor_match_partial_filter(self):
        """Only False tumor_match entries are skipped; True and None are included."""
        # pathogenic entry matches tumor
        pathogenic = make_clinvar_entry(
            "Pathogenic", locus_match="variant", tumor_match=True
        )
        # benign entry is from a different tumor type
        benign = make_clinvar_entry(
            "Benign", locus_match="variant", tumor_match=False
        )
        evidence = self._make_evidence([pathogenic, benign])

        high, significant = detect_discordant_evidence_internal(evidence)

        # benign is filtered out, so no conflict
        assert len(high) == 0

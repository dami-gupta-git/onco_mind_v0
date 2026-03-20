"""Integration tests for ClinVar tumor type matching.

Tests the full pipeline: ClinVarEvidence objects annotated with cancer_type_match,
their tumor_match computed property, and the effect on discordance detection.
Uses real model objects (no mocks).
"""

import pytest

from oncomind.models.evidence import Evidence
from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.evidence import VariantIdentifiers
from oncomind.models.evidence.base import EvidenceLevel, compute_cancer_specificity
from oncomind.insight_builder.gap_detection.discordance import (
    detect_discordant_evidence_internal,
)


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


class TestClinVarTumorMatchingIntegration:
    """Integration tests using real Evidence and ClinVarEvidence objects."""

    def test_chek2_breast_cancer_entry_not_matched_for_nsclc(self):
        """CHEK2 breast cancer ClinVar entry should not match when querying NSCLC."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic/Likely pathogenic",
            review_status="criteria provided, multiple submitters, no conflicts",
            conditions=["Hereditary breast ovarian cancer syndrome"],
            variation_id="141834",
        )
        annotate_clinvar_entries([entry], "NSCLC")

        assert entry.tumor_match is False

    def test_egfr_lung_entry_matched_for_nsclc(self):
        """Lung cancer ClinVar entry should match when querying NSCLC."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            review_status="criteria provided, single submitter",
            conditions=["Non-small cell lung carcinoma"],
            variation_id="12345",
        )
        annotate_clinvar_entries([entry], "NSCLC")

        assert entry.tumor_match is True

    def test_pan_cancer_entry_tumor_match_is_false(self):
        """Pan-cancer ClinVar condition sets pan_cancer level; tumor_match is False (not cancer_specific)."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Cancer"],
            variation_id="99999",
        )
        annotate_clinvar_entries([entry], "NSCLC")

        assert entry.cancer_type_match.level == "pan_cancer"
        assert entry.tumor_match is False

    def test_no_conditions_leaves_tumor_match_none(self):
        """ClinVar entry without conditions has tumor_match=None."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=[],
            variation_id="00001",
        )
        annotate_clinvar_entries([entry], "NSCLC")

        assert entry.cancer_type_match is None
        assert entry.tumor_match is None

    def test_discordance_not_triggered_for_other_tumor_entries(self):
        """Pathogenic/benign conflict from non-matching tumor should not be flagged."""
        pathogenic = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Hereditary breast ovarian cancer syndrome"],
            variation_id="111",
            locus_variant_match=EvidenceLevel(level="variant"),
        )
        benign = ClinVarEvidence(
            clinical_significance="Benign",
            conditions=["Hereditary breast ovarian cancer syndrome"],
            variation_id="111",
            locus_variant_match=EvidenceLevel(level="variant"),
        )
        annotate_clinvar_entries([pathogenic, benign], "Colorectal Cancer")

        assert pathogenic.tumor_match is False
        assert benign.tumor_match is False

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="CHEK2:I157T", gene="CHEK2", variant="I157T"
            ),
            clinvar_entries=[pathogenic, benign],
        )

        high, significant = detect_discordant_evidence_internal(evidence)
        assert len(high) == 0

    def test_discordance_triggered_for_matched_tumor_entries(self):
        """Pathogenic/benign conflict for matching tumor is still flagged."""
        pathogenic = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Colorectal cancer"],
            variation_id="222",
            locus_variant_match=EvidenceLevel(level="variant"),
        )
        benign = ClinVarEvidence(
            clinical_significance="Benign",
            conditions=["Colorectal cancer"],
            variation_id="222",
            locus_variant_match=EvidenceLevel(level="variant"),
        )
        annotate_clinvar_entries([pathogenic, benign], "Colorectal Cancer")

        assert pathogenic.tumor_match is True
        assert benign.tumor_match is True

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="CHEK2:I157T", gene="CHEK2", variant="I157T"
            ),
            clinvar_entries=[pathogenic, benign],
        )

        high, significant = detect_discordant_evidence_internal(evidence)
        assert len(high) == 1
        assert "pathogenic" in high[0].lower()
        assert "benign" in high[0].lower()

    def test_tumor_match_serializes_into_model_dump(self):
        """tumor_match is a computed_field and must appear in model_dump output."""
        entry = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Lung adenocarcinoma"],
            variation_id="333",
        )
        annotate_clinvar_entries([entry], "NSCLC")

        dumped = entry.model_dump()
        assert "tumor_match" in dumped
        assert dumped["tumor_match"] is True

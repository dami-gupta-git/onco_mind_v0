"""Unit tests for gap_detection/discordance.py.

Tests discordant evidence detection between different sources.
"""

import pytest
from unittest.mock import MagicMock

from oncomind.insight_builder.gap_detection.discordance import (
    check_discordant_evidence,
    detect_discordant_evidence_internal,
)
from oncomind.insight_builder.gap_detection.context import GapDetectionContext
from oncomind.models.extracted.evidence_gaps import GapCategory, GapSeverity
from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement


@pytest.fixture
def mock_evidence():
    """Create a mock Evidence object with empty evidence lists."""
    evidence = MagicMock()
    evidence.identifiers.gene = "BRAF"
    evidence.identifiers.variant = "V600E"
    evidence.fda_biomarker_evidence = []
    evidence.cgi_biomarkers = []
    evidence.preclinical_biomarkers = []
    evidence.early_phase_biomarkers = []
    evidence.vicc_evidence = []
    evidence.civic_assertions = []
    evidence.civic_evidence = []
    evidence.clinvar_entries = []
    return evidence


@pytest.fixture
def base_context():
    """Create a base GapDetectionContext for testing."""
    return GapDetectionContext(
        gene="BRAF",
        variant="V600E",
        tumor_type="Melanoma",
        is_cancer_gene=True,
        has_pathogenic_signal=True,
    )


# =============================================================================
# TEST detect_discordant_evidence_internal
# =============================================================================


class TestDetectDiscordantEvidenceInternal:
    """Tests for detect_discordant_evidence_internal function."""

    def test_no_conflicts_empty_evidence(self, mock_evidence):
        """Empty evidence should have no conflicts."""
        high, significant = detect_discordant_evidence_internal(mock_evidence)
        assert len(high) == 0
        assert len(significant) == 0

    def test_same_source_no_conflict(self, mock_evidence):
        """Same source with both sensitive and resistant should NOT flag conflict."""
        # CIViC evidence with sensitive
        civic_sens = MagicMock()
        civic_sens.drugs = ["dabrafenib"]
        civic_sens.is_sensitivity = True
        civic_sens.is_resistance = False
        civic_sens.locus_match = "variant"

        # CIViC evidence with resistance (different context)
        civic_resist = MagicMock()
        civic_resist.drugs = ["dabrafenib"]
        civic_resist.is_sensitivity = False
        civic_resist.is_resistance = True
        civic_resist.locus_match = "variant"

        mock_evidence.civic_evidence = [civic_sens, civic_resist]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        # Same source conflicts are NOT flagged (context-dependent responses)
        assert len(high) == 0
        assert len(significant) == 0

    def test_cross_source_conflict_cgi_vs_civic(self, mock_evidence):
        """CGI sensitive vs CIViC resistant should flag HIGH conflict."""
        # CGI says sensitive at variant level
        cgi = MagicMock()
        cgi.drug = "vemurafenib"
        cgi.association = "Sensitivity"
        cgi.locus_match = "variant"
        mock_evidence.cgi_biomarkers = [cgi]

        # CIViC says resistant at variant level
        civic = MagicMock()
        civic.drugs = ["vemurafenib"]
        civic.is_resistance = True
        civic.is_sensitivity = False
        civic.locus_match = "variant"
        mock_evidence.civic_evidence = [civic]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        assert len(high) == 1
        assert "vemurafenib" in high[0].lower()
        assert "sensitive" in high[0].lower()
        assert "resistant" in high[0].lower()

    def test_fda_vs_vicc_conflict_flagged_as_significant(self, mock_evidence):
        """FDA sensitive vs VICC resistant should flag SIGNIFICANT (not HIGH)."""
        # Use a variant/drug combo NOT in sensitizing variants list
        mock_evidence.identifiers.gene = "PIK3CA"
        mock_evidence.identifiers.variant = "H1047R"

        # FDA says sensitive at variant level (hypothetical drug not in sensitizing list)
        fda = MagicMock()
        fda.drug_name = "testdrug123"
        fda.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        fda.locus_match = "variant"
        mock_evidence.fda_biomarker_evidence = [fda]

        # VICC says resistant at variant level
        vicc = MagicMock()
        vicc.drugs = ["testdrug123"]
        vicc.response_type = "RESISTANCE"
        vicc.source = None
        vicc.locus_match = "variant"
        mock_evidence.vicc_evidence = [vicc]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        # FDA vs VICC is SIGNIFICANT, not HIGH (VICC less reliable)
        assert len(high) == 0
        assert len(significant) == 1
        assert "fda approves" in significant[0].lower()

    def test_fda_vs_vicc_skipped_for_known_sensitizing_variant(self, mock_evidence):
        """FDA sensitive vs VICC resistant should be SKIPPED for known sensitizing variants."""
        # Set up BRAF V600E which is a known sensitizing variant for dabrafenib
        mock_evidence.identifiers.gene = "BRAF"
        mock_evidence.identifiers.variant = "V600E"

        # FDA says sensitive at variant level
        fda = MagicMock()
        fda.drug_name = "dabrafenib"
        fda.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        fda.locus_match = "variant"
        mock_evidence.fda_biomarker_evidence = [fda]

        # VICC says resistant at variant level (acquired resistance is expected)
        vicc = MagicMock()
        vicc.drugs = ["dabrafenib"]
        vicc.response_type = "RESISTANCE"
        vicc.source = None
        vicc.locus_match = "variant"
        mock_evidence.vicc_evidence = [vicc]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        # Should be skipped because acquired resistance is expected
        assert len(high) == 0
        assert len(significant) == 0

    def test_clinvar_conflict_variant_level(self, mock_evidence):
        """ClinVar with both pathogenic and benign at variant level should flag conflict."""
        entry_pathogenic = MagicMock()
        entry_pathogenic.clinical_significance = "Pathogenic"
        entry_pathogenic.locus_match = "variant"

        entry_benign = MagicMock()
        entry_benign.clinical_significance = "Benign"
        entry_benign.locus_match = "variant"

        mock_evidence.clinvar_entries = [entry_pathogenic, entry_benign]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        assert len(high) == 1
        assert "clinvar" in high[0].lower()
        assert "pathogenic" in high[0].lower() and "benign" in high[0].lower()

    def test_clinvar_conflict_codon_level(self, mock_evidence):
        """ClinVar with both pathogenic and benign at codon level should flag conflict."""
        entry_pathogenic = MagicMock()
        entry_pathogenic.clinical_significance = "Pathogenic"
        entry_pathogenic.locus_match = "codon"

        entry_benign = MagicMock()
        entry_benign.clinical_significance = "Benign"
        entry_benign.locus_match = "codon"

        mock_evidence.clinvar_entries = [entry_pathogenic, entry_benign]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        assert len(high) == 1
        assert "codon level" in high[0].lower()

    def test_clinvar_gene_level_no_conflict(self, mock_evidence):
        """ClinVar with pathogenic and benign at gene level should NOT flag conflict."""
        entry_pathogenic = MagicMock()
        entry_pathogenic.clinical_significance = "Pathogenic"
        entry_pathogenic.locus_match = "gene"

        entry_benign = MagicMock()
        entry_benign.clinical_significance = "Benign"
        entry_benign.locus_match = "gene"

        mock_evidence.clinvar_entries = [entry_pathogenic, entry_benign]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        # Gene-level is not flagged (different variants legitimately differ)
        assert len(high) == 0
        assert len(significant) == 0

    def test_clinvar_benign_vs_civic_oncogenic_conflict(self, mock_evidence):
        """ClinVar benign at variant level vs CIViC oncogenic should flag conflict."""
        entry_benign = MagicMock()
        entry_benign.clinical_significance = "Benign"
        entry_benign.locus_match = "variant"
        mock_evidence.clinvar_entries = [entry_benign]

        # CIViC oncogenic assertion at variant level
        civic_assertion = MagicMock()
        civic_assertion.assertion_type = "ONCOGENIC"
        civic_assertion.significance = None
        civic_assertion.therapies = None
        civic_assertion.locus_match = "variant"
        mock_evidence.civic_assertions = [civic_assertion]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        assert len(high) == 1
        assert "clinvar" in high[0].lower()
        assert "benign" in high[0].lower()
        assert "civic" in high[0].lower()

    def test_clinvar_benign_vs_civic_predictive_conflict(self, mock_evidence):
        """ClinVar benign at variant level vs CIViC predictive with therapy should flag conflict."""
        entry_benign = MagicMock()
        entry_benign.clinical_significance = "Benign"
        entry_benign.locus_match = "variant"
        mock_evidence.clinvar_entries = [entry_benign]

        # CIViC predictive assertion with therapy at variant level
        civic_assertion = MagicMock()
        civic_assertion.assertion_type = "PREDICTIVE"
        civic_assertion.significance = None
        civic_assertion.therapies = ["Vemurafenib"]
        civic_assertion.locus_match = "variant"
        mock_evidence.civic_assertions = [civic_assertion]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        assert len(high) == 1
        assert "clinvar" in high[0].lower()
        assert "benign" in high[0].lower()

    def test_no_conflict_gene_level_only(self, mock_evidence):
        """Gene-level evidence conflicts should NOT be flagged."""
        # CGI says sensitive at gene level
        cgi = MagicMock()
        cgi.drug = "vemurafenib"
        cgi.association = "Sensitivity"
        cgi.locus_match = "gene"
        mock_evidence.cgi_biomarkers = [cgi]

        # CIViC says resistant at gene level
        civic = MagicMock()
        civic.drugs = ["vemurafenib"]
        civic.is_resistance = True
        civic.is_sensitivity = False
        civic.locus_match = "gene"
        mock_evidence.civic_evidence = [civic]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        # Gene-level conflicts are NOT flagged
        assert len(high) == 0
        assert len(significant) == 0

    def test_excludes_required_negative_fda(self, mock_evidence):
        """FDA evidence with REQUIRED_NEGATIVE should be excluded from sensitivity."""
        # FDA REQUIRED_NEGATIVE (contraindicated, not sensitivity)
        fda = MagicMock()
        fda.drug_name = "dabrafenib"
        fda.requirement = BiomarkerRequirement.REQUIRED_NEGATIVE
        fda.locus_match = "variant"
        mock_evidence.fda_biomarker_evidence = [fda]

        # VICC says resistant
        vicc = MagicMock()
        vicc.drugs = ["dabrafenib"]
        vicc.response_type = "RESISTANCE"
        vicc.source = None
        vicc.locus_match = "variant"
        mock_evidence.vicc_evidence = [vicc]

        high, significant = detect_discordant_evidence_internal(mock_evidence)

        # No conflict since FDA REQUIRED_NEGATIVE is excluded
        assert len(high) == 0
        assert len(significant) == 0


# =============================================================================
# TEST check_discordant_evidence
# =============================================================================


class TestCheckDiscordantEvidence:
    """Tests for check_discordant_evidence function."""

    def test_adds_high_severity_gaps(self, mock_evidence, base_context):
        """HIGH conflicts should add HIGH severity gaps."""
        # Set up a cross-source conflict
        cgi = MagicMock()
        cgi.drug = "vemurafenib"
        cgi.association = "Sensitivity"
        cgi.locus_match = "variant"
        mock_evidence.cgi_biomarkers = [cgi]

        civic = MagicMock()
        civic.drugs = ["vemurafenib"]
        civic.is_resistance = True
        civic.is_sensitivity = False
        civic.locus_match = "variant"
        mock_evidence.civic_evidence = [civic]

        check_discordant_evidence(mock_evidence, base_context)

        discordant_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.DISCORDANT
        ]
        assert len(discordant_gaps) >= 1
        assert discordant_gaps[0].severity == GapSeverity.HIGH

    def test_adds_significant_severity_gaps(self, mock_evidence, base_context):
        """FDA vs VICC conflicts should add SIGNIFICANT severity gaps."""
        # Set up FDA vs VICC conflict (using a non-sensitizing variant)
        mock_evidence.identifiers.gene = "EGFR"
        mock_evidence.identifiers.variant = "G12C"  # Not a known sensitizing variant

        fda = MagicMock()
        fda.drug_name = "sotorasib"
        fda.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        fda.locus_match = "variant"
        mock_evidence.fda_biomarker_evidence = [fda]

        vicc = MagicMock()
        vicc.drugs = ["sotorasib"]
        vicc.response_type = "RESISTANCE"
        vicc.source = None
        vicc.locus_match = "variant"
        mock_evidence.vicc_evidence = [vicc]

        base_context.gene = "EGFR"
        base_context.variant = "G12C"

        check_discordant_evidence(mock_evidence, base_context)

        discordant_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.DISCORDANT
        ]
        significant_gaps = [
            g for g in discordant_gaps if g.severity == GapSeverity.SIGNIFICANT
        ]
        assert len(significant_gaps) >= 1

    def test_adds_poorly_characterized(self, mock_evidence, base_context):
        """Conflicts should add 'conflicting evidence' to poorly_characterized."""
        # Set up a conflict
        cgi = MagicMock()
        cgi.drug = "vemurafenib"
        cgi.association = "Sensitivity"
        cgi.locus_match = "variant"
        mock_evidence.cgi_biomarkers = [cgi]

        civic = MagicMock()
        civic.drugs = ["vemurafenib"]
        civic.is_resistance = True
        civic.is_sensitivity = False
        civic.locus_match = "variant"
        mock_evidence.civic_evidence = [civic]

        check_discordant_evidence(mock_evidence, base_context)

        assert "conflicting evidence" in base_context.poorly_characterized

    def test_no_conflict_no_gaps(self, mock_evidence, base_context):
        """No conflicts should not add any gaps."""
        check_discordant_evidence(mock_evidence, base_context)

        discordant_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.DISCORDANT
        ]
        assert len(discordant_gaps) == 0

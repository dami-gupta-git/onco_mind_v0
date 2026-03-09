"""Unit tests for gap_detection/context.py.

Tests GapDetectionContext dataclass and its methods.
"""

import pytest

from oncomind.insight_builder.gap_detection.context import GapDetectionContext
from oncomind.models.extracted.evidence_gaps import (
    GapCategory,
    GapSeverity,
    EvidenceGap,
    CharacterizedAspect,
)

# =============================================================================
# FIXTURES
# =============================================================================


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


@pytest.fixture
def context_no_tumor():
    """Create a GapDetectionContext without tumor type."""
    return GapDetectionContext(
        gene="TP53",
        variant="R175H",
        tumor_type=None,
        is_cancer_gene=True,
        has_pathogenic_signal=True,
    )


# =============================================================================
# TEST GapDetectionContext initialization
# =============================================================================


class TestGapDetectionContextInit:
    """Tests for GapDetectionContext initialization."""

    def test_basic_initialization(self, base_context):
        """Context should initialize with provided values."""
        assert base_context.gene == "BRAF"
        assert base_context.variant == "V600E"
        assert base_context.tumor_type == "Melanoma"
        assert base_context.is_cancer_gene is True
        assert base_context.has_pathogenic_signal is True

    def test_default_empty_lists(self, base_context):
        """Context should initialize with empty lists by default."""
        assert base_context.gaps == []
        assert base_context.well_characterized == []
        assert base_context.well_characterized_detailed == []
        assert base_context.poorly_characterized == []

    def test_default_flags(self, base_context):
        """Context should initialize with False flags by default."""
        assert base_context.has_clinical is False
        assert base_context.has_drug_data is False

    def test_tumor_type_optional(self, context_no_tumor):
        """Context should accept None for tumor_type."""
        assert context_no_tumor.tumor_type is None


# =============================================================================
# TEST add_well_characterized
# =============================================================================


class TestAddWellCharacterized:
    """Tests for add_well_characterized method."""

    def test_adds_to_well_characterized_list(self, base_context):
        """Should add aspect to well_characterized list."""
        base_context.add_well_characterized(
            aspect="FDA-approved therapy",
            basis="3 drugs approved",
            category=GapCategory.DRUG_RESPONSE,
        )

        assert len(base_context.well_characterized) == 1
        assert "Fda-Approved Therapy" in base_context.well_characterized

    def test_adds_to_well_characterized_detailed(self, base_context):
        """Should add CharacterizedAspect to well_characterized_detailed list."""
        base_context.add_well_characterized(
            aspect="FDA-approved therapy",
            basis="3 drugs approved",
            category=GapCategory.DRUG_RESPONSE,
        )

        assert len(base_context.well_characterized_detailed) == 1
        detail = base_context.well_characterized_detailed[0]
        assert isinstance(detail, CharacterizedAspect)
        assert detail.basis == "3 drugs approved"
        assert detail.category == GapCategory.DRUG_RESPONSE

    def test_aspect_is_title_cased(self, base_context):
        """Aspect should be converted to title case."""
        base_context.add_well_characterized(
            aspect="fda-approved therapy",
            basis="1 drug",
            category=GapCategory.DRUG_RESPONSE,
        )

        assert base_context.well_characterized[0] == "Fda-Approved Therapy"
        assert (
            base_context.well_characterized_detailed[0].aspect == "Fda-Approved Therapy"
        )

    def test_includes_matches_on(self, base_context):
        """Should include matches_on field when provided."""
        base_context.add_well_characterized(
            aspect="resistance mechanisms",
            basis="5 entries",
            category=GapCategory.RESISTANCE,
            matches_on="3 variant, 2 gene",
        )

        detail = base_context.well_characterized_detailed[0]
        assert detail.matches_on == "3 variant, 2 gene"

    def test_includes_tumor_match(self, base_context):
        """Should include tumor_match field when provided."""
        base_context.add_well_characterized(
            aspect="clinical evidence",
            basis="10 entries",
            category=GapCategory.CLINICAL,
            tumor_match="5 tumor, 3 other",
        )

        detail = base_context.well_characterized_detailed[0]
        assert detail.tumor_match == "5 tumor, 3 other"

    def test_category_optional(self, base_context):
        """Category should be optional (None allowed)."""
        base_context.add_well_characterized(
            aspect="some aspect",
            basis="basis text",
        )

        detail = base_context.well_characterized_detailed[0]
        assert detail.category is None

    def test_multiple_entries(self, base_context):
        """Should handle multiple entries correctly."""
        base_context.add_well_characterized("aspect one", "basis one")
        base_context.add_well_characterized("aspect two", "basis two")
        base_context.add_well_characterized("aspect three", "basis three")

        assert len(base_context.well_characterized) == 3
        assert len(base_context.well_characterized_detailed) == 3


# =============================================================================
# TEST add_gap
# =============================================================================


class TestAddGap:
    """Tests for add_gap method."""

    def test_adds_gap_to_list(self, base_context):
        """Should add EvidenceGap to gaps list."""
        base_context.add_gap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.HIGH,
            description="No clinical trials available",
        )

        assert len(base_context.gaps) == 1
        gap = base_context.gaps[0]
        assert isinstance(gap, EvidenceGap)
        assert gap.category == GapCategory.CLINICAL
        assert gap.severity == GapSeverity.HIGH
        assert gap.description == "No clinical trials available"

    def test_includes_suggested_studies(self, base_context):
        """Should include suggested_studies when provided."""
        base_context.add_gap(
            category=GapCategory.PRECLINICAL,
            severity=GapSeverity.MODERATE,
            description="Limited preclinical data",
            suggested_studies=["Cell line studies", "PDX models"],
        )

        gap = base_context.gaps[0]
        assert gap.suggested_studies == ["Cell line studies", "PDX models"]

    def test_includes_addressable_with(self, base_context):
        """Should include addressable_with when provided."""
        base_context.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description="No functional data",
            addressable_with=["DepMap", "COSMIC"],
        )

        gap = base_context.gaps[0]
        assert gap.addressable_with == ["DepMap", "COSMIC"]

    def test_default_empty_lists(self, base_context):
        """Should default to empty lists for optional fields."""
        base_context.add_gap(
            category=GapCategory.VALIDATION,
            severity=GapSeverity.CRITICAL,
            description="No validation",
        )

        gap = base_context.gaps[0]
        assert gap.suggested_studies == []
        assert gap.addressable_with == []

    def test_multiple_gaps(self, base_context):
        """Should handle multiple gaps correctly."""
        base_context.add_gap(GapCategory.CLINICAL, GapSeverity.HIGH, "Gap 1")
        base_context.add_gap(GapCategory.FUNCTIONAL, GapSeverity.MODERATE, "Gap 2")
        base_context.add_gap(GapCategory.PRECLINICAL, GapSeverity.MINOR, "Gap 3")

        assert len(base_context.gaps) == 3
        assert base_context.gaps[0].description == "Gap 1"
        assert base_context.gaps[1].description == "Gap 2"
        assert base_context.gaps[2].description == "Gap 3"

    def test_all_severities(self, base_context):
        """Should handle all severity levels."""
        severities = [
            GapSeverity.CRITICAL,
            GapSeverity.HIGH,
            GapSeverity.SIGNIFICANT,
            GapSeverity.MODERATE,
            GapSeverity.MINOR,
            GapSeverity.INFORMATIONAL,
        ]

        for sev in severities:
            base_context.add_gap(GapCategory.CLINICAL, sev, f"Gap with {sev.value}")

        assert len(base_context.gaps) == 6


# =============================================================================
# TEST add_poorly_characterized
# =============================================================================


class TestAddPoorlyCharacterized:
    """Tests for add_poorly_characterized method."""

    def test_adds_to_poorly_characterized_list(self, base_context):
        """Should add aspect to poorly_characterized list."""
        base_context.add_poorly_characterized("resistance mechanisms")

        assert len(base_context.poorly_characterized) == 1
        assert "resistance mechanisms" in base_context.poorly_characterized

    def test_preserves_original_case(self, base_context):
        """Should preserve the original case of the aspect."""
        base_context.add_poorly_characterized("NSCLC-specific data")

        assert base_context.poorly_characterized[0] == "NSCLC-specific data"

    def test_multiple_entries(self, base_context):
        """Should handle multiple entries correctly."""
        base_context.add_poorly_characterized("aspect one")
        base_context.add_poorly_characterized("aspect two")
        base_context.add_poorly_characterized("aspect three")

        assert len(base_context.poorly_characterized) == 3
        assert "aspect one" in base_context.poorly_characterized
        assert "aspect two" in base_context.poorly_characterized
        assert "aspect three" in base_context.poorly_characterized


# =============================================================================
# TEST context flags
# =============================================================================


class TestContextFlags:
    """Tests for context flags (has_clinical, has_drug_data)."""

    def test_flags_can_be_set(self, base_context):
        """Flags should be settable."""
        base_context.has_clinical = True
        base_context.has_drug_data = True

        assert base_context.has_clinical is True
        assert base_context.has_drug_data is True

    def test_flags_independent(self, base_context):
        """Flags should be independent of each other."""
        base_context.has_clinical = True
        base_context.has_drug_data = False

        assert base_context.has_clinical is True
        assert base_context.has_drug_data is False


# =============================================================================
# TEST combined usage
# =============================================================================


class TestCombinedUsage:
    """Tests for using multiple context methods together."""

    def test_full_workflow(self, base_context):
        """Should handle full gap detection workflow."""
        # Add well-characterized aspects
        base_context.add_well_characterized(
            "FDA-approved therapy",
            "2 drugs",
            GapCategory.DRUG_RESPONSE,
            matches_on="2 variant",
            tumor_match="2 tumor",
        )
        base_context.add_well_characterized(
            "clinical trials",
            "5 trials",
            GapCategory.CLINICAL,
        )

        # Add gaps
        base_context.add_gap(
            GapCategory.PRECLINICAL,
            GapSeverity.MODERATE,
            "Limited preclinical models",
            suggested_studies=["PDX models"],
        )

        # Add poorly characterized
        base_context.add_poorly_characterized("resistance mechanisms")

        # Set flags
        base_context.has_clinical = True
        base_context.has_drug_data = True

        # Verify all state
        assert len(base_context.well_characterized) == 2
        assert len(base_context.well_characterized_detailed) == 2
        assert len(base_context.gaps) == 1
        assert len(base_context.poorly_characterized) == 1
        assert base_context.has_clinical is True
        assert base_context.has_drug_data is True

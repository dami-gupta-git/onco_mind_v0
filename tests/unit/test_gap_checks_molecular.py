"""Unit tests for molecular gap check functions.

Tests: check_hotspot_context, check_functional_predictions,
check_gene_mechanism, check_validation_gap.
"""

from unittest.mock import MagicMock

from oncomind.insight_builder.gap_detector import GapDetectionContext
from oncomind.insight_builder.gap_detection.checks import (
    check_hotspot_context as _check_hotspot_context,
    check_functional_predictions as _check_functional_predictions,
    check_gene_mechanism as _check_gene_mechanism,
    check_validation_gap as _check_validation_gap,
)
from oncomind.models.extracted.evidence_gaps import GapCategory, GapSeverity

from tests.unit.conftest import create_mock_depmap_evidence

# =============================================================================
# TEST _check_hotspot_context
# =============================================================================


class TestCheckHotspotContext:
    """Tests for _check_hotspot_context function."""

    def test_known_hotspot_well_characterized(self, mock_evidence):
        """Known hotspot should be marked as well-characterized."""
        mock_evidence.identifiers.gene = "BRAF"
        mock_evidence.identifiers.variant = "V600E"

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
        )

        _check_hotspot_context(mock_evidence, ctx)

        assert any("hotspot" in w.lower() for w in ctx.well_characterized)
        # Should NOT add a gap for known hotspots
        hotspot_gaps = [g for g in ctx.gaps if "hotspot" in g.description.lower()]
        assert len(hotspot_gaps) == 0

    def test_hotspot_adjacent_flagged(self, mock_evidence):
        """Hotspot-adjacent variant should flag a functional gap."""
        mock_evidence.identifiers.gene = "BRAF"
        mock_evidence.identifiers.variant = "V598E"  # Near V600

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V598E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
        )

        _check_hotspot_context(mock_evidence, ctx)

        # Should be marked as near hotspot
        assert any("near hotspot" in w.lower() for w in ctx.well_characterized)
        # Should have a functional gap
        functional_gaps = [g for g in ctx.gaps if g.category == GapCategory.FUNCTIONAL]
        assert len(functional_gaps) >= 1

    def test_non_hotspot_no_special_marking(self, mock_evidence):
        """Non-hotspot variant should not get hotspot markings."""
        mock_evidence.identifiers.gene = "BRAF"
        mock_evidence.identifiers.variant = "V100E"  # Far from hotspots

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V100E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=False,
        )

        _check_hotspot_context(mock_evidence, ctx)

        assert not any("hotspot" in w.lower() for w in ctx.well_characterized)


# =============================================================================
# TEST _check_functional_predictions
# =============================================================================


class TestCheckFunctionalPredictions:
    """Tests for _check_functional_predictions function."""

    def test_with_alphamissense_score(self, mock_evidence, base_context):
        """AlphaMissense score should mark as well-characterized."""
        mock_evidence.functional.alphamissense_score = 0.95

        _check_functional_predictions(mock_evidence, base_context)

        assert any(
            "pathogenicity" in w.lower() for w in base_context.well_characterized
        )
        functional_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.FUNCTIONAL
        ]
        pathogenicity_gaps = [
            g for g in functional_gaps if "pathogenicity" in g.description.lower()
        ]
        assert len(pathogenicity_gaps) == 0

    def test_with_cadd_score(self, mock_evidence, base_context):
        """CADD score should mark as well-characterized."""
        mock_evidence.functional.cadd_score = 25.0

        _check_functional_predictions(mock_evidence, base_context)

        assert any(
            "pathogenicity" in w.lower() for w in base_context.well_characterized
        )

    def test_with_polyphen2(self, mock_evidence, base_context):
        """PolyPhen2 prediction should mark as well-characterized."""
        mock_evidence.functional.polyphen2_prediction = "probably_damaging"

        _check_functional_predictions(mock_evidence, base_context)

        assert any(
            "pathogenicity" in w.lower() for w in base_context.well_characterized
        )

    def test_no_predictions_adds_gap(self, mock_evidence, base_context):
        """Missing all predictions should add a SIGNIFICANT gap."""
        mock_evidence.functional.alphamissense_score = None
        mock_evidence.functional.alphamissense_prediction = None
        mock_evidence.functional.cadd_score = None
        mock_evidence.functional.polyphen2_prediction = None
        mock_evidence.functional.polyphen2_score = None
        mock_evidence.functional.sift_prediction = None
        mock_evidence.functional.sift_score = None

        _check_functional_predictions(mock_evidence, base_context)

        assert "pathogenicity predictions" in base_context.poorly_characterized
        functional_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.FUNCTIONAL
        ]
        assert len(functional_gaps) >= 1
        assert any(g.severity == GapSeverity.SIGNIFICANT for g in functional_gaps)


# =============================================================================
# TEST _check_gene_mechanism
# =============================================================================


class TestCheckGeneMechanism:
    """Tests for _check_gene_mechanism function."""

    def test_with_gene_role(self, mock_evidence, base_context):
        """Gene with known role should be well-characterized."""
        mock_evidence.context.gene_role = "oncogene"
        mock_evidence.context.pathway = "MAPK"

        _check_gene_mechanism(mock_evidence, base_context)

        assert any(
            "gene function" in w.lower() for w in base_context.well_characterized
        )

    def test_with_depmap_essentiality(self, mock_evidence, base_context):
        """DepMap essentiality data should be well-characterized."""
        mock_evidence.depmap_evidence = create_mock_depmap_evidence(
            gene="TESTGENE",
            variant="V100E",
            is_essential=True,
            mean_dependency_score=-0.8,
            dependency_pct=75.0,
            n_dependent_lines=150,
            n_total_lines=200,
        )

        _check_gene_mechanism(mock_evidence, base_context)

        assert any("essentiality" in w.lower() for w in base_context.well_characterized)

    def test_no_mechanism_data_adds_gap(self, mock_evidence, base_context):
        """Missing mechanism data should add a gap."""
        mock_evidence.context.gene_role = None
        mock_evidence.context.pathway = None
        mock_evidence.depmap_evidence = None

        _check_gene_mechanism(mock_evidence, base_context)

        assert "functional mechanism" in base_context.poorly_characterized
        functional_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.FUNCTIONAL
        ]
        assert len(functional_gaps) >= 1


# =============================================================================
# TEST _check_validation_gap
# =============================================================================


class TestCheckValidationGap:
    """Tests for _check_validation_gap function."""

    def test_strong_signal_no_validation_critical(self, mock_evidence):
        """Strong oncogenic signal without validation should be CRITICAL for cancer gene."""
        mock_evidence.depmap_evidence = create_mock_depmap_evidence(
            gene="BRAF",
            variant="V600E",
            is_essential=True,
            mean_dependency_score=-0.9,
            dependency_pct=85.0,
        )

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
        )

        _check_validation_gap(mock_evidence, ctx)

        validation_gaps = [g for g in ctx.gaps if g.category == GapCategory.VALIDATION]
        assert len(validation_gaps) >= 1
        assert validation_gaps[0].severity == GapSeverity.CRITICAL

    def test_strong_signal_with_validation_no_gap(self, mock_evidence):
        """Strong signal with therapeutic validation should NOT add gap."""
        mock_evidence.depmap_evidence = create_mock_depmap_evidence(
            gene="BRAF",
            variant="V600E",
            is_essential=True,
            mean_dependency_score=-0.9,
            dependency_pct=85.0,
        )
        mock_evidence.civic_assertions = [MagicMock()]

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
        )

        _check_validation_gap(mock_evidence, ctx)

        validation_gaps = [g for g in ctx.gaps if g.category == GapCategory.VALIDATION]
        assert len(validation_gaps) == 0

    def test_no_strong_signal_no_gap(self, mock_evidence, base_context):
        """Without strong oncogenic signal, no validation gap."""
        _check_validation_gap(mock_evidence, base_context)

        validation_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.VALIDATION
        ]
        assert len(validation_gaps) == 0

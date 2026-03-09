"""Unit tests for observational gap check functions.

Tests: check_prevalence, check_clinical_trials,
check_preclinical_models, check_literature_depth.
"""

from unittest.mock import MagicMock

from oncomind.insight_builder.gap_detector import GapDetectionContext
from oncomind.insight_builder.gap_detection.checks import (
    check_prevalence as _check_prevalence,
    check_clinical_trials as _check_clinical_trials,
    check_preclinical_models as _check_preclinical_models,
    check_literature_depth as _check_literature_depth,
)
from oncomind.models.extracted.evidence_gaps import GapCategory, GapSeverity

from tests.unit.conftest import (
    create_mock_depmap_evidence,
    create_mock_cbioportal_evidence,
    create_mock_cell_line,
)

# =============================================================================
# TEST _check_prevalence
# =============================================================================


class TestCheckPrevalence:
    """Tests for _check_prevalence function."""

    def test_with_cbioportal_data(self, mock_evidence, base_context):
        """cBioPortal data should mark as well-characterized."""
        mock_evidence.cbioportal_evidence = create_mock_cbioportal_evidence(
            gene="TESTGENE",
            variant="V100E",
            tumor_type="NSCLC",
            study_name="TCGA PanCancer Atlas",
            total_samples=1000,
            samples_with_gene_mutation=50,
            samples_with_exact_variant=10,
            has_data=True,
        )

        _check_prevalence(mock_evidence, base_context)

        assert any("observed" in w.lower() for w in base_context.well_characterized)

    def test_no_prevalence_informational_gap(self, mock_evidence, base_context):
        """Missing prevalence for non-cancer gene (in cBioPortal) is INFORMATIONAL."""
        mock_evidence.cbioportal_evidence = create_mock_cbioportal_evidence(
            gene="TESTGENE",
            variant="V100E",
            tumor_type="NSCLC",
            total_samples=1000,
            samples_with_gene_mutation=50,
            samples_with_exact_variant=0,  # Variant not seen
            has_data=True,
        )

        _check_prevalence(mock_evidence, base_context)

        prevalence_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.PREVALENCE
        ]
        assert len(prevalence_gaps) >= 1
        assert prevalence_gaps[0].severity == GapSeverity.INFORMATIONAL

    def test_no_prevalence_minor_for_cancer_gene_no_clinical(self, mock_evidence):
        """Missing prevalence for cancer gene without clinical evidence is MINOR."""
        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600X",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=False,
        )

        mock_evidence.cbioportal_evidence = create_mock_cbioportal_evidence(
            gene="BRAF",
            variant="V600X",
            tumor_type="Melanoma",
            total_samples=1000,
            samples_with_gene_mutation=100,
            samples_with_exact_variant=0,  # Variant not seen
            has_data=True,
        )

        _check_prevalence(mock_evidence, ctx)

        prevalence_gaps = [g for g in ctx.gaps if g.category == GapCategory.PREVALENCE]
        assert len(prevalence_gaps) >= 1
        assert prevalence_gaps[0].severity == GapSeverity.MINOR

    def test_no_prevalence_significant_for_clinical_cancer_gene(self, mock_evidence):
        """Missing prevalence for clinical cancer gene is SIGNIFICANT."""
        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
        )

        mock_evidence.cbioportal_evidence = create_mock_cbioportal_evidence(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            total_samples=1000,
            samples_with_gene_mutation=100,
            samples_with_exact_variant=0,  # Variant not seen
            has_data=True,
        )

        _check_prevalence(mock_evidence, ctx)

        prevalence_gaps = [g for g in ctx.gaps if g.category == GapCategory.PREVALENCE]
        assert len(prevalence_gaps) >= 1
        assert prevalence_gaps[0].severity == GapSeverity.SIGNIFICANT


# =============================================================================
# TEST _check_clinical_trials
# =============================================================================


class TestCheckClinicalTrials:
    """Tests for _check_clinical_trials function."""

    def test_with_trials(self, mock_evidence, base_context):
        """Active trials should mark as well-characterized."""
        mock_evidence.clinical_trials = [MagicMock(), MagicMock()]
        base_context.has_clinical = True

        _check_clinical_trials(mock_evidence, base_context)

        assert any("trial" in w.lower() for w in base_context.well_characterized)

    def test_no_trials_with_clinical_adds_informational_gap(self, mock_evidence):
        """Clinical variant without trials should add INFORMATIONAL gap."""
        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_clinical_trials(mock_evidence, ctx)

        clinical_gaps = [g for g in ctx.gaps if g.category == GapCategory.CLINICAL]
        trial_gaps = [g for g in clinical_gaps if "trial" in g.description.lower()]
        assert len(trial_gaps) >= 1
        assert trial_gaps[0].severity == GapSeverity.INFORMATIONAL

    def test_no_trials_without_clinical_no_gap(self, mock_evidence, base_context):
        """Non-clinical variant without trials should NOT add gap."""
        _check_clinical_trials(mock_evidence, base_context)

        clinical_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.CLINICAL
        ]
        trial_gaps = [g for g in clinical_gaps if "trial" in g.description.lower()]
        assert len(trial_gaps) == 0


# =============================================================================
# TEST _check_preclinical_models
# =============================================================================


class TestCheckPreclinicalModels:
    """Tests for _check_preclinical_models function."""

    def test_with_cell_line_models(self, mock_evidence, base_context):
        """Cell line models should mark as well-characterized."""
        lung_cell_line = create_mock_cell_line(
            name="A549",
            depmap_id="ACH-000681",
            primary_disease="Lung Cancer",
            subtype="Non-Small Cell Lung Cancer",
            has_mutation=True,
            mutation_details="V100E",
        )
        mock_evidence.depmap_evidence = create_mock_depmap_evidence(
            gene="TESTGENE",
            variant="V100E",
            cell_line_models=[lung_cell_line],
        )

        _check_preclinical_models(mock_evidence, base_context)

        assert any(
            "model" in w.lower() or "cell line" in w.lower()
            for w in base_context.well_characterized
        )

    def test_models_exist_but_wrong_tumor_type(self, mock_evidence):
        """Models with mutation but wrong tumor type should flag gap."""
        colon_cell_line = create_mock_cell_line(
            name="HCT116",
            depmap_id="ACH-000971",
            primary_disease="Colon Cancer",
            subtype="Colorectal Adenocarcinoma",
            has_mutation=True,
            mutation_details="G12D",
        )
        mock_evidence.depmap_evidence = create_mock_depmap_evidence(
            gene="KRAS",
            variant="G12D",
            cell_line_models=[colon_cell_line],
        )

        ctx = GapDetectionContext(
            gene="KRAS",
            variant="G12D",
            tumor_type="Pancreatic",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_preclinical_models(mock_evidence, ctx)

        preclinical_gaps = [
            g for g in ctx.gaps if g.category == GapCategory.PRECLINICAL
        ]
        assert len(preclinical_gaps) >= 1

    def test_no_models_with_clinical_adds_minor_gap(self, mock_evidence):
        """Clinical variant without models should add MINOR gap."""
        mock_evidence.context.gene_role = "oncogene"

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_preclinical_models(mock_evidence, ctx)

        preclinical_gaps = [
            g for g in ctx.gaps if g.category == GapCategory.PRECLINICAL
        ]
        assert len(preclinical_gaps) >= 1
        assert preclinical_gaps[0].severity == GapSeverity.MINOR

    def test_tumor_type_matching_with_alias(self, mock_evidence):
        """Test that tumor type matching uses aliases (e.g., Melanoma matches SKIN)."""
        skin_cell_line = create_mock_cell_line(
            name="A375",
            depmap_id="ACH-000219",
            primary_disease="SKIN",
            subtype="Cutaneous Melanoma",
            has_mutation=True,
            mutation_details="V600E",
        )
        mock_evidence.depmap_evidence = create_mock_depmap_evidence(
            gene="BRAF",
            variant="V600E",
            cell_line_models=[skin_cell_line],
        )

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_preclinical_models(mock_evidence, ctx)

        assert any("melanoma" in w.lower() for w in ctx.well_characterized)
        preclinical_gaps = [
            g for g in ctx.gaps if g.category == GapCategory.PRECLINICAL
        ]
        assert len(preclinical_gaps) == 0

    def test_tumor_type_matching_nsclc_lung(self, mock_evidence):
        """Test that NSCLC matches LUNG tissue."""
        lung_cell_line = create_mock_cell_line(
            name="H1975",
            depmap_id="ACH-000414",
            primary_disease="LUNG",
            subtype="Non-Small Cell Lung Cancer",
            has_mutation=True,
            mutation_details="L858R",
        )
        mock_evidence.depmap_evidence = create_mock_depmap_evidence(
            gene="EGFR",
            variant="L858R",
            cell_line_models=[lung_cell_line],
        )

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_preclinical_models(mock_evidence, ctx)

        assert any("nsclc" in w.lower() for w in ctx.well_characterized)
        preclinical_gaps = [
            g for g in ctx.gaps if g.category == GapCategory.PRECLINICAL
        ]
        assert len(preclinical_gaps) == 0


# =============================================================================
# TEST _check_literature_depth
# =============================================================================


class TestCheckLiteratureDepth:
    """Tests for _check_literature_depth function."""

    def test_literature_not_searched_no_gap(self, mock_evidence, base_context):
        """If literature wasn't searched, don't report gap."""
        mock_evidence.literature_searched = False

        _check_literature_depth(mock_evidence, base_context)

        functional_gaps = [
            g for g in base_context.gaps if "literature" in g.description.lower()
        ]
        assert len(functional_gaps) == 0

    def test_no_literature_cancer_gene_critical(self, mock_evidence):
        """Cancer gene with no literature should be CRITICAL gap."""
        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = []

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
        )

        _check_literature_depth(mock_evidence, ctx)

        literature_gaps = [g for g in ctx.gaps if "literature" in g.description.lower()]
        assert len(literature_gaps) >= 1
        assert literature_gaps[0].severity == GapSeverity.CRITICAL

    def test_no_literature_non_cancer_gene_significant(
        self, mock_evidence, base_context
    ):
        """Non-cancer gene with no literature should be SIGNIFICANT gap."""
        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = []

        _check_literature_depth(mock_evidence, base_context)

        literature_gaps = [
            g for g in base_context.gaps if "literature" in g.description.lower()
        ]
        assert len(literature_gaps) >= 1
        assert literature_gaps[0].severity == GapSeverity.SIGNIFICANT

    def test_limited_literature_poorly_characterized(self, mock_evidence, base_context):
        """Few articles should mark as poorly-characterized."""
        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = [MagicMock(), MagicMock()]  # Only 2

        _check_literature_depth(mock_evidence, base_context)

        assert any("literature" in p.lower() for p in base_context.poorly_characterized)

    def test_sufficient_literature_well_characterized(
        self, mock_evidence, base_context
    ):
        """Many articles should mark as well-characterized."""
        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = [MagicMock() for _ in range(10)]

        _check_literature_depth(mock_evidence, base_context)

        assert any("literature" in w.lower() for w in base_context.well_characterized)

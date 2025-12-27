"""Unit tests for individual gap detection functions.

Tests each _check_* function in gap_detector.py independently using
mock Evidence objects to verify correct gap identification.
"""

import pytest
from unittest.mock import MagicMock, PropertyMock

from oncomind.insight_builder.gap_detector import (
    GapDetectionContext,
    detect_evidence_gaps,
    _check_hotspot_context,
    _check_functional_predictions,
    _check_gene_mechanism,
    _check_clinical_evidence,
    _check_tumor_type_evidence,
    _check_drug_response,
    _check_resistance_mechanisms,
    _check_discordant_evidence,
    _check_prevalence,
    _check_clinical_trials,
    _check_preclinical_models,
    _check_literature_depth,
    _check_validation_gap,
    _has_pathogenic_signal,
    _check_tumor_specific_evidence,
)
from oncomind.models.evidence.evidence_gaps import GapCategory, GapSeverity
from oncomind.models.evidence.literature_knowledge import LiteratureKnowledge, DrugResistance


# =============================================================================
# FIXTURES
# =============================================================================

@pytest.fixture
def mock_evidence():
    """Create a mock Evidence object with default values."""
    evidence = MagicMock()

    # Identifiers
    evidence.identifiers.gene = "TESTGENE"
    evidence.identifiers.variant = "V100E"

    # Context
    evidence.context.tumor_type = "NSCLC"
    evidence.context.gene_role = None
    evidence.context.pathway = None

    # Functional scores (all None by default)
    evidence.functional.alphamissense_score = None
    evidence.functional.alphamissense_prediction = None
    evidence.functional.cadd_score = None
    evidence.functional.polyphen2_prediction = None
    evidence.functional.gnomad_exome_af = None
    evidence.functional.snpeff_effect = None

    # Evidence lists (empty by default)
    evidence.civic_assertions = []
    evidence.civic_evidence = []
    evidence.fda_approvals = []
    evidence.vicc_evidence = []
    evidence.cgi_biomarkers = []
    evidence.preclinical_biomarkers = []
    evidence.pubmed_articles = []
    evidence.clinical_trials = []
    evidence.clinvar_entries = []
    evidence.clinvar_significance = None

    # DepMap evidence
    evidence.depmap_evidence = None

    # cBioPortal evidence
    evidence.cbioportal_evidence = None

    # Literature flag
    evidence.literature_searched = False

    # LLM-extracted literature knowledge
    evidence.literature_knowledge = None

    return evidence


@pytest.fixture
def base_context():
    """Create a base GapDetectionContext for testing."""
    return GapDetectionContext(
        gene="TESTGENE",
        variant="V100E",
        tumor_type="NSCLC",
        is_cancer_gene=False,
        has_pathogenic_signal=False,
    )


# =============================================================================
# TEST GapDetectionContext
# =============================================================================

class TestGapDetectionContext:
    """Tests for the GapDetectionContext dataclass."""

    def test_add_well_characterized(self, base_context):
        """Test adding well-characterized aspects."""
        base_context.add_well_characterized("test aspect", "test basis")

        assert "test aspect" in base_context.well_characterized
        assert len(base_context.well_characterized_detailed) == 1
        assert base_context.well_characterized_detailed[0].aspect == "test aspect"
        assert base_context.well_characterized_detailed[0].basis == "test basis"

    def test_add_gap(self, base_context):
        """Test adding a gap."""
        base_context.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description="Test gap",
            suggested_studies=["Study 1"],
            addressable_with=["Tool 1"],
        )

        assert len(base_context.gaps) == 1
        assert base_context.gaps[0].category == GapCategory.FUNCTIONAL
        assert base_context.gaps[0].severity == GapSeverity.SIGNIFICANT

    def test_add_poorly_characterized(self, base_context):
        """Test adding poorly-characterized aspects."""
        base_context.add_poorly_characterized("test aspect")

        assert "test aspect" in base_context.poorly_characterized


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

        assert any("pathogenicity" in w.lower() for w in base_context.well_characterized)
        functional_gaps = [g for g in base_context.gaps if g.category == GapCategory.FUNCTIONAL]
        pathogenicity_gaps = [g for g in functional_gaps if "pathogenicity" in g.description.lower()]
        assert len(pathogenicity_gaps) == 0

    def test_with_cadd_score(self, mock_evidence, base_context):
        """CADD score should mark as well-characterized."""
        mock_evidence.functional.cadd_score = 25.0

        _check_functional_predictions(mock_evidence, base_context)

        assert any("pathogenicity" in w.lower() for w in base_context.well_characterized)

    def test_with_polyphen2(self, mock_evidence, base_context):
        """PolyPhen2 prediction should mark as well-characterized."""
        mock_evidence.functional.polyphen2_prediction = "probably_damaging"

        _check_functional_predictions(mock_evidence, base_context)

        assert any("pathogenicity" in w.lower() for w in base_context.well_characterized)

    def test_no_predictions_adds_gap(self, mock_evidence, base_context):
        """Missing all predictions should add a SIGNIFICANT gap."""
        _check_functional_predictions(mock_evidence, base_context)

        assert "pathogenicity predictions" in base_context.poorly_characterized
        functional_gaps = [g for g in base_context.gaps if g.category == GapCategory.FUNCTIONAL]
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

        assert any("gene function" in w.lower() for w in base_context.well_characterized)

    def test_with_depmap_essentiality(self, mock_evidence, base_context):
        """DepMap essentiality data should be well-characterized."""
        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.gene_dependency = MagicMock()
        mock_evidence.depmap_evidence.gene_dependency.mean_dependency_score = -0.8
        mock_evidence.depmap_evidence.gene_dependency.dependency_pct = 75.0
        mock_evidence.depmap_evidence.is_essential.return_value = True

        _check_gene_mechanism(mock_evidence, base_context)

        assert any("essentiality" in w.lower() for w in base_context.well_characterized)

    def test_no_mechanism_data_adds_gap(self, mock_evidence, base_context):
        """Missing mechanism data should add a gap."""
        _check_gene_mechanism(mock_evidence, base_context)

        assert "functional mechanism" in base_context.poorly_characterized
        functional_gaps = [g for g in base_context.gaps if g.category == GapCategory.FUNCTIONAL]
        assert len(functional_gaps) >= 1


# =============================================================================
# TEST _check_clinical_evidence
# =============================================================================

class TestCheckClinicalEvidence:
    """Tests for _check_clinical_evidence function."""

    def test_with_fda_approvals(self, mock_evidence, base_context):
        """FDA approvals should mark as well-characterized."""
        mock_evidence.fda_approvals = [MagicMock()]

        _check_clinical_evidence(mock_evidence, base_context)

        assert any("clinical" in w.lower() for w in base_context.well_characterized)
        assert base_context.has_clinical is True

    def test_with_civic_assertions(self, mock_evidence, base_context):
        """CIViC assertions should mark as well-characterized."""
        mock_evidence.civic_assertions = [MagicMock(), MagicMock()]

        _check_clinical_evidence(mock_evidence, base_context)

        assert any("clinical" in w.lower() for w in base_context.well_characterized)
        assert base_context.has_clinical is True

    def test_no_clinical_evidence_adds_critical_gap(self, mock_evidence, base_context):
        """Missing clinical evidence should add a CRITICAL gap."""
        _check_clinical_evidence(mock_evidence, base_context)

        assert "clinical evidence" in base_context.poorly_characterized
        assert base_context.has_clinical is False
        clinical_gaps = [g for g in base_context.gaps if g.category == GapCategory.CLINICAL]
        assert len(clinical_gaps) >= 1
        assert any(g.severity == GapSeverity.CRITICAL for g in clinical_gaps)


# =============================================================================
# TEST _check_tumor_type_evidence
# =============================================================================

class TestCheckTumorTypeEvidence:
    """Tests for _check_tumor_type_evidence function."""

    def test_no_tumor_type_skips(self, mock_evidence, base_context):
        """No tumor type should skip check entirely."""
        base_context.tumor_type = None

        _check_tumor_type_evidence(mock_evidence, base_context)

        assert len(base_context.gaps) == 0

    def test_with_tumor_specific_civic(self, mock_evidence):
        """Tumor-specific CIViC data should be well-characterized."""
        assertion = MagicMock()
        assertion.disease = "Non-Small Cell Lung Cancer"
        mock_evidence.civic_assertions = [assertion]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="L858R",
            tumor_type="Non-Small Cell Lung Cancer",  # Match the disease exactly
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
        )

        _check_tumor_type_evidence(mock_evidence, ctx)

        assert any("evidence in" in w.lower() for w in ctx.well_characterized)

    def test_no_tumor_specific_evidence_cancer_gene_critical(self, mock_evidence):
        """Cancer gene without tumor-specific data should get CRITICAL gap."""
        ctx = GapDetectionContext(
            gene="KRAS",
            variant="G12D",
            tumor_type="Pancreatic",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=False,  # No clinical evidence
        )

        _check_tumor_type_evidence(mock_evidence, ctx)

        tumor_gaps = [g for g in ctx.gaps if g.category == GapCategory.TUMOR_TYPE]
        assert len(tumor_gaps) >= 1
        assert tumor_gaps[0].severity == GapSeverity.CRITICAL


# =============================================================================
# TEST _check_drug_response
# =============================================================================

class TestCheckDrugResponse:
    """Tests for _check_drug_response function."""

    def test_with_cgi_biomarkers(self, mock_evidence, base_context):
        """CGI biomarkers should mark as well-characterized."""
        mock_evidence.cgi_biomarkers = [MagicMock()]

        _check_drug_response(mock_evidence, base_context)

        assert any("drug response" in w.lower() for w in base_context.well_characterized)
        assert base_context.has_drug_data is True

    def test_with_vicc_evidence(self, mock_evidence, base_context):
        """VICC evidence should mark as well-characterized."""
        mock_evidence.vicc_evidence = [MagicMock()]

        _check_drug_response(mock_evidence, base_context)

        assert base_context.has_drug_data is True

    def test_with_depmap_drug_sensitivities(self, mock_evidence, base_context):
        """DepMap drug sensitivities should mark as well-characterized."""
        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.drug_sensitivities = [MagicMock()]

        _check_drug_response(mock_evidence, base_context)

        assert any("depmap" in w.lower() for w in base_context.well_characterized)

    def test_no_drug_data_adds_gap(self, mock_evidence, base_context):
        """Missing drug data should add a SIGNIFICANT gap."""
        _check_drug_response(mock_evidence, base_context)

        assert "drug response" in base_context.poorly_characterized
        drug_gaps = [g for g in base_context.gaps if g.category == GapCategory.DRUG_RESPONSE]
        assert len(drug_gaps) >= 1
        assert drug_gaps[0].severity == GapSeverity.SIGNIFICANT


# =============================================================================
# TEST _check_resistance_mechanisms
# =============================================================================

class TestCheckResistanceMechanisms:
    """Tests for _check_resistance_mechanisms function."""

    def test_with_resistance_articles(self, mock_evidence):
        """Resistance articles should mark as well-characterized."""
        article = MagicMock()
        article.is_resistance_evidence.return_value = True
        mock_evidence.pubmed_articles = [article]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        assert any("resistance" in w.lower() for w in ctx.well_characterized)

    def test_with_cgi_resistance(self, mock_evidence):
        """CGI resistance data should mark as well-characterized."""
        biomarker = MagicMock()
        biomarker.association = "Resistance"
        mock_evidence.cgi_biomarkers = [biomarker]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        assert any("resistance" in w.lower() for w in ctx.well_characterized)

    def test_with_civic_assertion_resistance(self, mock_evidence):
        """CIViC assertion with is_resistance=True should mark as well-characterized."""
        assertion = MagicMock()
        assertion.is_resistance = True
        mock_evidence.civic_assertions = [assertion]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        assert any("resistance" in w.lower() for w in ctx.well_characterized)
        # Check that CIViC is mentioned in the basis
        well_char_detail = [w for w in ctx.well_characterized_detailed if "resistance" in w.aspect.lower()]
        assert len(well_char_detail) > 0
        assert "CIViC" in well_char_detail[0].basis

    def test_with_vicc_resistance(self, mock_evidence):
        """VICC evidence with resistance response type should mark as well-characterized."""
        vicc = MagicMock()
        vicc.response_type = "RESISTANCE"
        mock_evidence.vicc_evidence = [vicc]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        assert any("resistance" in w.lower() for w in ctx.well_characterized)
        # Check that VICC is mentioned in the basis
        well_char_detail = [w for w in ctx.well_characterized_detailed if "resistance" in w.aspect.lower()]
        assert len(well_char_detail) > 0
        assert "VICC" in well_char_detail[0].basis

    def test_with_llm_literature_knowledge_resistance(self, mock_evidence):
        """LLM literature knowledge with resistant_to should mark as well-characterized."""
        mock_evidence.literature_knowledge = LiteratureKnowledge(
            resistant_to=[
                DrugResistance(drug="Gefitinib", evidence="clinical", is_predictive=True),
                DrugResistance(drug="Erlotinib", evidence="clinical", is_predictive=True),
            ],
            mutation_type="secondary",
        )

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        assert any("resistance" in w.lower() for w in ctx.well_characterized)
        # Check that LLM is mentioned in the basis with drug count
        well_char_detail = [w for w in ctx.well_characterized_detailed if "resistance" in w.aspect.lower()]
        assert len(well_char_detail) > 0
        assert "LLM literature" in well_char_detail[0].basis
        assert "2 drugs" in well_char_detail[0].basis

    def test_with_llm_literature_non_predictive_resistance(self, mock_evidence):
        """Non-predictive resistance in literature should NOT mark as well-characterized."""
        mock_evidence.literature_knowledge = LiteratureKnowledge(
            resistant_to=[
                DrugResistance(drug="Gefitinib", evidence="clinical", is_predictive=False),
            ],
            mutation_type="secondary",
        )

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        # Non-predictive resistance should NOT be counted (only prognostic)
        # So this should add a gap, not well-characterized
        assert any("resistance" in p.lower() for p in ctx.poorly_characterized)

    def test_multiple_resistance_sources(self, mock_evidence):
        """Multiple resistance sources should all be included in basis."""
        # Add multiple sources
        article = MagicMock()
        article.is_resistance_evidence.return_value = True
        mock_evidence.pubmed_articles = [article]

        assertion = MagicMock()
        assertion.is_resistance = True
        mock_evidence.civic_assertions = [assertion]

        vicc = MagicMock()
        vicc.response_type = "RESISTANCE"
        mock_evidence.vicc_evidence = [vicc]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        assert any("resistance" in w.lower() for w in ctx.well_characterized)
        # Check that multiple sources are mentioned
        well_char_detail = [w for w in ctx.well_characterized_detailed if "resistance" in w.aspect.lower()]
        assert len(well_char_detail) > 0
        basis = well_char_detail[0].basis
        assert "PubMed" in basis
        assert "CIViC" in basis
        assert "VICC" in basis

    def test_no_resistance_with_clinical_adds_gap(self, mock_evidence):
        """Clinical variant without resistance data should add gap."""
        ctx = GapDetectionContext(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        resistance_gaps = [g for g in ctx.gaps if g.category == GapCategory.RESISTANCE]
        assert len(resistance_gaps) >= 1

    def test_no_resistance_without_clinical_no_gap(self, mock_evidence, base_context):
        """Non-clinical variant without resistance data should NOT add gap."""
        base_context.has_clinical = False
        base_context.has_drug_data = False

        _check_resistance_mechanisms(mock_evidence, base_context)

        resistance_gaps = [g for g in base_context.gaps if g.category == GapCategory.RESISTANCE]
        assert len(resistance_gaps) == 0


# =============================================================================
# TEST _check_prevalence
# =============================================================================

class TestCheckPrevalence:
    """Tests for _check_prevalence function."""

    def test_with_cbioportal_data(self, mock_evidence, base_context):
        """cBioPortal data should mark as well-characterized."""
        mock_evidence.cbioportal_evidence = MagicMock()
        mock_evidence.cbioportal_evidence.has_data.return_value = True
        mock_evidence.cbioportal_evidence.study_name = "TCGA"
        mock_evidence.cbioportal_evidence.variant_prevalence_pct = 5.2

        _check_prevalence(mock_evidence, base_context)

        assert any("prevalence" in w.lower() for w in base_context.well_characterized)

    def test_no_prevalence_minor_gap(self, mock_evidence, base_context):
        """Missing prevalence for non-cancer gene should be MINOR gap."""
        _check_prevalence(mock_evidence, base_context)

        prevalence_gaps = [g for g in base_context.gaps if g.category == GapCategory.PREVALENCE]
        assert len(prevalence_gaps) >= 1
        assert prevalence_gaps[0].severity == GapSeverity.MINOR

    def test_no_prevalence_significant_for_clinical_cancer_gene(self, mock_evidence):
        """Missing prevalence for clinical cancer gene should be SIGNIFICANT."""
        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
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

    def test_no_trials_with_clinical_adds_minor_gap(self, mock_evidence):
        """Clinical variant without trials should add MINOR gap."""
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
        assert trial_gaps[0].severity == GapSeverity.MINOR

    def test_no_trials_without_clinical_no_gap(self, mock_evidence, base_context):
        """Non-clinical variant without trials should NOT add gap."""
        _check_clinical_trials(mock_evidence, base_context)

        clinical_gaps = [g for g in base_context.gaps if g.category == GapCategory.CLINICAL]
        trial_gaps = [g for g in clinical_gaps if "trial" in g.description.lower()]
        assert len(trial_gaps) == 0


# =============================================================================
# TEST _check_preclinical_models
# =============================================================================

class TestCheckPreclinicalModels:
    """Tests for _check_preclinical_models function."""

    def test_with_cell_line_models(self, mock_evidence, base_context):
        """Cell line models should mark as well-characterized."""
        model = MagicMock()
        model.has_mutation = True
        model.primary_disease = "Lung Cancer"
        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.cell_line_models = [model]

        _check_preclinical_models(mock_evidence, base_context)

        assert any("model" in w.lower() or "cell line" in w.lower() for w in base_context.well_characterized)

    def test_models_exist_but_wrong_tumor_type(self, mock_evidence):
        """Models with mutation but wrong tumor type should flag gap."""
        model = MagicMock()
        model.has_mutation = True
        model.primary_disease = "Colon Cancer"
        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.cell_line_models = [model]

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

        preclinical_gaps = [g for g in ctx.gaps if g.category == GapCategory.PRECLINICAL]
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

        preclinical_gaps = [g for g in ctx.gaps if g.category == GapCategory.PRECLINICAL]
        assert len(preclinical_gaps) >= 1
        assert preclinical_gaps[0].severity == GapSeverity.MINOR

    def test_tumor_type_matching_with_alias(self, mock_evidence):
        """Test that tumor type matching uses aliases (e.g., Melanoma matches SKIN)."""
        model = MagicMock()
        model.has_mutation = True
        model.primary_disease = "SKIN"  # cBioPortal returns tissue as "SKIN"
        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.cell_line_models = [model]

        ctx = GapDetectionContext(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",  # User specifies "Melanoma"
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_preclinical_models(mock_evidence, ctx)

        # Should match because Melanoma and SKIN are related via TUMOR_TYPE_MAPPINGS
        assert any("melanoma" in w.lower() for w in ctx.well_characterized)
        # Should NOT have a gap for wrong tumor type
        preclinical_gaps = [g for g in ctx.gaps if g.category == GapCategory.PRECLINICAL]
        assert len(preclinical_gaps) == 0

    def test_tumor_type_matching_nsclc_lung(self, mock_evidence):
        """Test that NSCLC matches LUNG tissue."""
        model = MagicMock()
        model.has_mutation = True
        model.primary_disease = "LUNG"
        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.cell_line_models = [model]

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

        # Should match because NSCLC and LUNG are related
        assert any("nsclc" in w.lower() for w in ctx.well_characterized)
        preclinical_gaps = [g for g in ctx.gaps if g.category == GapCategory.PRECLINICAL]
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

        functional_gaps = [g for g in base_context.gaps if "literature" in g.description.lower()]
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

    def test_no_literature_non_cancer_gene_significant(self, mock_evidence, base_context):
        """Non-cancer gene with no literature should be SIGNIFICANT gap."""
        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = []

        _check_literature_depth(mock_evidence, base_context)

        literature_gaps = [g for g in base_context.gaps if "literature" in g.description.lower()]
        assert len(literature_gaps) >= 1
        assert literature_gaps[0].severity == GapSeverity.SIGNIFICANT

    def test_limited_literature_poorly_characterized(self, mock_evidence, base_context):
        """Few articles should mark as poorly-characterized."""
        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = [MagicMock(), MagicMock()]  # Only 2

        _check_literature_depth(mock_evidence, base_context)

        assert any("literature" in p.lower() for p in base_context.poorly_characterized)

    def test_sufficient_literature_well_characterized(self, mock_evidence, base_context):
        """Many articles should mark as well-characterized."""
        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = [MagicMock() for _ in range(10)]

        _check_literature_depth(mock_evidence, base_context)

        assert any("literature" in w.lower() for w in base_context.well_characterized)


# =============================================================================
# TEST _check_validation_gap
# =============================================================================

class TestCheckValidationGap:
    """Tests for _check_validation_gap function."""

    def test_strong_signal_no_validation_critical(self, mock_evidence):
        """Strong oncogenic signal without validation should be CRITICAL for cancer gene."""
        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.is_essential.return_value = True

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
        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.is_essential.return_value = True
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

        validation_gaps = [g for g in base_context.gaps if g.category == GapCategory.VALIDATION]
        assert len(validation_gaps) == 0


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
        mock_evidence.functional.cadd_score = 15.0

        assert _has_pathogenic_signal(mock_evidence) is False

    def test_polyphen2_damaging(self, mock_evidence):
        """PolyPhen2 damaging should return True."""
        mock_evidence.functional.polyphen2_prediction = "probably_damaging"

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_civic_assertions(self, mock_evidence):
        """CIViC assertions should return True."""
        mock_evidence.civic_assertions = [MagicMock()]

        assert _has_pathogenic_signal(mock_evidence) is True

    def test_fda_approvals(self, mock_evidence):
        """FDA approvals should return True."""
        mock_evidence.fda_approvals = [MagicMock()]

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
        assert _has_pathogenic_signal(mock_evidence) is False


# =============================================================================
# TEST detect_evidence_gaps (integration of all checks)
# =============================================================================

class TestDetectEvidenceGapsIntegration:
    """Integration tests for the main detect_evidence_gaps function."""

    def test_well_studied_variant_few_gaps(self, mock_evidence):
        """Well-studied variant should have few gaps."""
        # Set up as well-studied
        mock_evidence.identifiers.gene = "BRAF"
        mock_evidence.identifiers.variant = "V600E"
        mock_evidence.context.tumor_type = "Melanoma"
        mock_evidence.context.gene_role = "oncogene"
        mock_evidence.functional.alphamissense_score = 0.98
        mock_evidence.functional.alphamissense_prediction = "pathogenic"
        mock_evidence.functional.cadd_score = 30.0
        mock_evidence.civic_assertions = [MagicMock()]
        mock_evidence.fda_approvals = [MagicMock()]
        mock_evidence.cgi_biomarkers = [MagicMock()]
        mock_evidence.literature_searched = True
        mock_evidence.pubmed_articles = [MagicMock() for _ in range(20)]

        gaps = detect_evidence_gaps(mock_evidence)

        # Should have good quality
        assert gaps.overall_evidence_quality in ("comprehensive", "moderate")
        assert len(gaps.well_characterized) > len(gaps.gaps)

    def test_unknown_variant_many_gaps(self, mock_evidence):
        """Unknown variant should have many gaps."""
        mock_evidence.identifiers.gene = "UNKNOWNGENE"
        mock_evidence.identifiers.variant = "X999Y"
        mock_evidence.literature_searched = True

        gaps = detect_evidence_gaps(mock_evidence)

        # Should have poor quality
        assert gaps.overall_evidence_quality in ("limited", "minimal")
        assert len(gaps.gaps) > 0

    def test_returns_evidence_gaps_object(self, mock_evidence):
        """Should return EvidenceGaps object with all fields."""
        gaps = detect_evidence_gaps(mock_evidence)

        assert hasattr(gaps, 'gaps')
        assert hasattr(gaps, 'overall_evidence_quality')
        assert hasattr(gaps, 'well_characterized')
        assert hasattr(gaps, 'well_characterized_detailed')
        assert hasattr(gaps, 'poorly_characterized')
        assert hasattr(gaps, 'research_priority')

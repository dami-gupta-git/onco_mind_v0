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
    _compute_overall_quality,
)
from oncomind.models.evidence.evidence_gaps import GapCategory, GapSeverity, EvidenceGap
from oncomind.models.evidence.literature_knowledge import LiteratureKnowledge, LitDrugResistance


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
    evidence.functional.polyphen2_score = None
    evidence.functional.sift_prediction = None
    evidence.functional.sift_score = None
    evidence.functional.gnomad_exome_af = None
    evidence.functional.gnomad_genome_af = None
    evidence.functional.snpeff_effect = None

    # Evidence lists (empty by default)
    evidence.civic_assertions = []
    evidence.civic_evidence = []
    evidence.fda_approvals = []
    evidence.vicc_evidence = []
    evidence.cgi_biomarkers = []
    evidence.preclinical_biomarkers = []
    evidence.early_phase_biomarkers = []
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
        """Test adding well-characterized aspects with title casing and category."""
        base_context.add_well_characterized(
            "test aspect",
            "test basis",
            category=GapCategory.FUNCTIONAL
        )

        # Aspect should be title-cased
        assert "Test Aspect" in base_context.well_characterized
        assert len(base_context.well_characterized_detailed) == 1
        assert base_context.well_characterized_detailed[0].aspect == "Test Aspect"
        assert base_context.well_characterized_detailed[0].basis == "test basis"
        assert base_context.well_characterized_detailed[0].category == GapCategory.FUNCTIONAL

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
        # Mock get_fda_approved_therapies to return empty list
        mock_evidence.get_fda_approved_therapies.return_value = []

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
# TEST _check_tumor_specific_evidence
# =============================================================================

class TestCheckTumorSpecificEvidence:
    """Tests for _check_tumor_specific_evidence helper function."""

    def test_tumor_evidence_match_sources_string(self, mock_evidence):
        """TumorEvidenceMatch should build correct sources string."""
        # Create CIViC assertion matching tumor - use tumor type that will match
        assertion = MagicMock()
        assertion.disease = "Lung Cancer"  # Use "Lung Cancer" so it matches "Lung"
        assertion.variant_level = None
        mock_evidence.civic_assertions = [assertion]

        # Create VICC evidence matching tumor
        vicc = MagicMock()
        vicc.disease = "Lung adenocarcinoma"  # Will match "Lung"
        vicc.variant_level = None
        mock_evidence.vicc_evidence = [vicc, vicc]  # 2 VICC items

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "Lung")

        assert tumor_match.has_tumor_evidence is True
        sources = tumor_match.sources_str
        assert "1 CIViC Assertions" in sources
        assert "2 VICC" in sources

    def test_tumor_evidence_match_locus_levels(self, mock_evidence):
        """TumorEvidenceMatch should track variant/codon/gene match levels."""
        from oncomind.models.evidence.base import EvidenceLevel

        # Create CIViC assertion with variant-level match
        assertion_variant = MagicMock()
        assertion_variant.disease = "Lung Cancer"
        assertion_variant.variant_level = EvidenceLevel(level="variant", scope="specific")

        # Create CIViC assertion with gene-level match
        assertion_gene = MagicMock()
        assertion_gene.disease = "Lung Cancer"
        assertion_gene.variant_level = EvidenceLevel(level="gene", scope="unspecified")

        mock_evidence.civic_assertions = [assertion_variant, assertion_gene]

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "Lung Cancer")

        assert tumor_match.total_matches == 2
        assert tumor_match.total_variant_level == 1
        assert tumor_match.total_gene_level == 1
        matches_on = tumor_match.matches_on_str
        assert "1 variant" in matches_on
        assert "1 gene" in matches_on

    def test_tumor_evidence_match_no_matches(self, mock_evidence):
        """TumorEvidenceMatch should handle no tumor matches correctly."""
        # Create evidence that doesn't match the tumor type
        assertion = MagicMock()
        assertion.disease = "Melanoma"
        assertion.variant_level = None
        mock_evidence.civic_assertions = [assertion]

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "Pancreatic")

        assert tumor_match.has_tumor_evidence is False
        assert tumor_match.total_matches == 0
        assert tumor_match.sources_str == "No tumor-specific evidence"
        assert tumor_match.matches_on_str is None

    def test_tumor_evidence_includes_all_cgi_tiers(self, mock_evidence):
        """TumorEvidenceMatch should include all CGI tiers (FDA, preclinical, early phase)."""
        # Create CGI FDA-approved biomarker
        cgi_fda = MagicMock()
        cgi_fda.tumor_type = "Lung Cancer"
        cgi_fda.variant_level = None
        mock_evidence.cgi_biomarkers = [cgi_fda]

        # Create preclinical biomarker
        preclin = MagicMock()
        preclin.tumor_type = "Lung Cancer"
        preclin.variant_level = None
        mock_evidence.preclinical_biomarkers = [preclin]

        # Create early phase biomarker
        early = MagicMock()
        early.tumor_type = "Lung Cancer"
        early.variant_level = None
        mock_evidence.early_phase_biomarkers = [early]

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "Lung Cancer")

        assert tumor_match.has_tumor_evidence is True
        assert tumor_match.total_matches == 3
        # All CGI tiers combined under "CGI"
        assert tumor_match.cgi_biomarkers is not None
        assert tumor_match.cgi_biomarkers.count == 3

    def test_tumor_evidence_fda_parsing(self, mock_evidence):
        """TumorEvidenceMatch should use FDA approval's parse_indication_for_tumor method."""
        # Create FDA approval with indication parsing
        fda = MagicMock()
        fda.indication = "Non-small cell lung cancer"
        fda.parse_indication_for_tumor = MagicMock(return_value={"tumor_match": True})
        fda.variant_level = None
        mock_evidence.fda_approvals = [fda]

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "NSCLC")

        assert tumor_match.has_tumor_evidence is True
        fda.parse_indication_for_tumor.assert_called_once_with("NSCLC")


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
        """DepMap drug sensitivities should mark as well-characterized only with tumor-matched cell lines."""
        # Create mock cell line matching the tumor type
        mock_cell_line = MagicMock()
        mock_cell_line.has_mutation = True
        mock_cell_line.primary_disease = "Lung"  # Matches NSCLC

        mock_evidence.depmap_evidence = MagicMock()
        mock_evidence.depmap_evidence.drug_sensitivities = [MagicMock()]
        mock_evidence.depmap_evidence.cell_line_models = [mock_cell_line]

        _check_drug_response(mock_evidence, base_context)

        assert any("depmap" in w.lower() for w in base_context.well_characterized)

    def test_no_drug_data_adds_gap(self, mock_evidence, base_context):
        """Missing drug data should add a SIGNIFICANT gap."""
        _check_drug_response(mock_evidence, base_context)

        assert "drug response" in base_context.poorly_characterized
        drug_gaps = [g for g in base_context.gaps if g.category == GapCategory.DRUG_RESPONSE]
        assert len(drug_gaps) >= 1
        assert drug_gaps[0].severity == GapSeverity.SIGNIFICANT

    def test_drug_response_counts_sources_correctly(self, mock_evidence, base_context):
        """Drug response should count CGI, VICC, and FDA sources correctly."""
        # Create 2 CGI biomarkers
        cgi1 = MagicMock()
        cgi1.variant_level = None
        cgi1.tumor_type = "Lung Cancer"
        cgi2 = MagicMock()
        cgi2.variant_level = None
        cgi2.tumor_type = "Melanoma"
        mock_evidence.cgi_biomarkers = [cgi1, cgi2]

        # Create 3 VICC evidence items
        mock_evidence.vicc_evidence = [MagicMock(), MagicMock(), MagicMock()]
        for v in mock_evidence.vicc_evidence:
            v.variant_level = None
            v.disease = "Other Cancer"

        # Create 1 FDA approval
        fda = MagicMock()
        fda.variant_level = None
        fda.match_level = "variant"
        fda.indication = "Melanoma"
        fda.parse_indication_for_tumor = MagicMock(return_value={"tumor_match": False})
        mock_evidence.fda_approvals = [fda]

        _check_drug_response(mock_evidence, base_context)

        # Check the basis string contains correct counts
        drug_resp = [w for w in base_context.well_characterized_detailed if "drug response" in w.aspect.lower()]
        assert len(drug_resp) == 1
        basis = drug_resp[0].basis
        assert "2 CGI" in basis
        assert "3 VICC" in basis
        assert "1 FDA" in basis

    def test_drug_response_tumor_match_counting(self, mock_evidence):
        """Drug response should track tumor matches vs others."""
        # Create CGI biomarkers - 1 matching lung, 1 not
        # Tumor matching uses: tumor_lower in tumor_type.lower()
        # So "lung" in "lung cancer" = True
        cgi_match = MagicMock()
        cgi_match.variant_level = None
        cgi_match.tumor_type = "Lung Cancer"  # Contains "lung"
        cgi_no_match = MagicMock()
        cgi_no_match.variant_level = None
        cgi_no_match.tumor_type = "Melanoma"  # Does NOT contain "lung"
        mock_evidence.cgi_biomarkers = [cgi_match, cgi_no_match]

        # Create VICC evidence - 2 matching
        vicc1 = MagicMock()
        vicc1.variant_level = None
        vicc1.disease = "Lung adenocarcinoma"  # Contains "lung"
        vicc2 = MagicMock()
        vicc2.variant_level = None
        vicc2.disease = "Lung squamous"  # Contains "lung"
        mock_evidence.vicc_evidence = [vicc1, vicc2]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="L858R",
            tumor_type="Lung",  # Use "Lung" so it matches via substring
            is_cancer_gene=True,
            has_pathogenic_signal=True,
        )

        _check_drug_response(mock_evidence, ctx)

        # Check tumor_match field
        drug_resp = [w for w in ctx.well_characterized_detailed if "drug response" in w.aspect.lower()]
        assert len(drug_resp) == 1
        tumor_match = drug_resp[0].tumor_match
        assert tumor_match is not None
        # Should have 3 tumor matches (1 CGI + 2 VICC) and 1 other
        assert "3 tumor" in tumor_match
        assert "1 other" in tumor_match

    def test_drug_response_locus_match_levels(self, mock_evidence, base_context):
        """Drug response should track match levels (variant, codon, gene)."""
        from oncomind.models.evidence.base import EvidenceLevel

        # Create CGI biomarker with variant-level match
        cgi_variant = MagicMock()
        cgi_variant.variant_level = EvidenceLevel(level="variant", scope="specific")
        cgi_variant.tumor_type = None
        mock_evidence.cgi_biomarkers = [cgi_variant]

        # Create VICC evidence with gene-level match
        vicc_gene = MagicMock()
        vicc_gene.variant_level = EvidenceLevel(level="gene", scope="unspecified")
        vicc_gene.disease = None
        mock_evidence.vicc_evidence = [vicc_gene]

        # Create FDA approval with codon-level match
        fda_codon = MagicMock()
        fda_codon.variant_level = EvidenceLevel(level="codon", scope="specific")
        fda_codon.match_level = None  # Will fall back to variant_level
        fda_codon.indication = None
        mock_evidence.fda_approvals = [fda_codon]

        _check_drug_response(mock_evidence, base_context)

        # Check matches_on field
        drug_resp = [w for w in base_context.well_characterized_detailed if "drug response" in w.aspect.lower()]
        assert len(drug_resp) == 1
        matches_on = drug_resp[0].matches_on
        assert matches_on is not None
        assert "1 variant" in matches_on
        assert "1 codon" in matches_on
        assert "1 gene" in matches_on

    def test_drug_response_excludes_preclinical_biomarkers(self, mock_evidence, base_context):
        """Drug response should NOT include preclinical or early phase biomarkers."""
        # Set up preclinical biomarkers (should be excluded from Drug Response)
        mock_evidence.preclinical_biomarkers = [MagicMock(), MagicMock()]
        mock_evidence.early_phase_biomarkers = [MagicMock()]

        _check_drug_response(mock_evidence, base_context)

        # Drug response should NOT be well-characterized (no FDA-tier evidence)
        assert "drug response" in base_context.poorly_characterized
        # But preclinical row should exist
        preclin = [w for w in base_context.well_characterized_detailed
                   if "preclinical" in w.aspect.lower()]
        assert len(preclin) == 1
        assert "2 preclinical" in preclin[0].basis
        assert "1 early phase" in preclin[0].basis


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
                LitDrugResistance(drug="Gefitinib", evidence="clinical", is_predictive=True),
                LitDrugResistance(drug="Erlotinib", evidence="clinical", is_predictive=True),
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
                LitDrugResistance(drug="Gefitinib", evidence="clinical", is_predictive=False),
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

    def test_resistance_mechanisms_locus_match_levels(self, mock_evidence):
        """Resistance mechanisms should track match levels (variant, codon, gene)."""
        from oncomind.models.evidence.base import EvidenceLevel

        # Create CGI biomarker with resistance at variant level
        cgi_resistance = MagicMock()
        cgi_resistance.association = "Resistance"
        cgi_resistance.variant_level = EvidenceLevel(level="variant", scope="specific")
        cgi_resistance.tumor_type = None
        mock_evidence.cgi_biomarkers = [cgi_resistance]

        # Create VICC evidence with resistance at gene level
        vicc_resistance = MagicMock()
        vicc_resistance.response_type = "RESISTANCE"
        vicc_resistance.variant_level = EvidenceLevel(level="gene", scope="unspecified")
        vicc_resistance.disease = None
        mock_evidence.vicc_evidence = [vicc_resistance]

        # Create CIViC assertion with resistance at codon level
        civic_resistance = MagicMock()
        civic_resistance.is_resistance = True
        civic_resistance.variant_level = EvidenceLevel(level="codon", scope="specific")
        civic_resistance.disease = None
        mock_evidence.civic_assertions = [civic_resistance]

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

        # Check matches_on field
        resist_detail = [w for w in ctx.well_characterized_detailed if "resistance" in w.aspect.lower()]
        assert len(resist_detail) == 1
        matches_on = resist_detail[0].matches_on
        assert matches_on is not None
        assert "1 variant" in matches_on
        assert "1 codon" in matches_on
        assert "1 gene" in matches_on

    def test_resistance_mechanisms_tumor_match_counting(self, mock_evidence):
        """Resistance mechanisms should track tumor matches vs others."""
        # Tumor matching uses: tumor_lower in tumor_type/disease.lower()
        # So "lung" in "lung cancer" = True

        # Create CGI biomarker with resistance - matching tumor
        cgi_resistance = MagicMock()
        cgi_resistance.association = "Resistance"
        cgi_resistance.variant_level = None
        cgi_resistance.tumor_type = "Lung Cancer"  # Contains "lung"
        mock_evidence.cgi_biomarkers = [cgi_resistance]

        # Create VICC evidence with resistance - not matching tumor
        vicc_resistance = MagicMock()
        vicc_resistance.response_type = "RESISTANCE"
        vicc_resistance.variant_level = None
        vicc_resistance.disease = "Melanoma"  # Does NOT contain "lung"
        mock_evidence.vicc_evidence = [vicc_resistance]

        # Create CIViC assertion with resistance - matching tumor
        civic_resistance = MagicMock()
        civic_resistance.is_resistance = True
        civic_resistance.variant_level = None
        civic_resistance.disease = "Lung adenocarcinoma"  # Contains "lung"
        mock_evidence.civic_assertions = [civic_resistance]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="Lung",  # Use "Lung" so it matches via substring
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        # Check tumor_match field
        resist_detail = [w for w in ctx.well_characterized_detailed if "resistance" in w.aspect.lower()]
        assert len(resist_detail) == 1
        tumor_match = resist_detail[0].tumor_match
        assert tumor_match is not None
        # Should have 2 tumor matches (1 CGI + 1 CIViC) and 1 other (VICC)
        assert "2 tumor" in tumor_match
        assert "1 other" in tumor_match

    def test_resistance_mechanisms_counts_sources_correctly(self, mock_evidence):
        """Resistance mechanisms should count all source types correctly."""
        # 1 PubMed article with resistance
        article = MagicMock()
        article.is_resistance_evidence.return_value = True
        mock_evidence.pubmed_articles = [article]

        # 2 CGI biomarkers with resistance
        cgi1 = MagicMock()
        cgi1.association = "Resistance"
        cgi1.variant_level = None
        cgi1.tumor_type = None
        cgi2 = MagicMock()
        cgi2.association = "RESISTANCE"
        cgi2.variant_level = None
        cgi2.tumor_type = None
        mock_evidence.cgi_biomarkers = [cgi1, cgi2]

        # 1 CIViC assertion with resistance
        civic = MagicMock()
        civic.is_resistance = True
        civic.variant_level = None
        civic.disease = None
        mock_evidence.civic_assertions = [civic]

        # 3 VICC evidence with resistance
        mock_evidence.vicc_evidence = []
        for _ in range(3):
            vicc = MagicMock()
            vicc.response_type = "RESISTANCE"
            vicc.variant_level = None
            vicc.disease = None
            mock_evidence.vicc_evidence.append(vicc)

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

        # Check the basis string contains correct counts
        resist_detail = [w for w in ctx.well_characterized_detailed if "resistance" in w.aspect.lower()]
        assert len(resist_detail) == 1
        basis = resist_detail[0].basis
        assert "1 PubMed" in basis
        assert "2 CGI" in basis
        assert "1 CIViC" in basis
        assert "3 VICC" in basis


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
        mock_evidence.cbioportal_evidence.samples_with_exact_variant = 10

        _check_prevalence(mock_evidence, base_context)

        assert any("observed" in w.lower() for w in base_context.well_characterized)

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
        # 1 minor gap with 0 well-characterized = moderate (small positive score)
        result_no_credit = _compute_overall_quality([minor_gap], 0)

        # Same gap with 5 well-characterized = comprehensive (negative net score)
        result_with_credit = _compute_overall_quality([minor_gap], 5)

        # With credit should be equal or better than without
        quality_order = ["comprehensive", "moderate", "limited", "minimal"]
        assert quality_order.index(result_with_credit) <= quality_order.index(result_no_credit)

    def test_many_well_characterized_yields_comprehensive(self):
        """Many well-characterized aspects should yield comprehensive even with gaps."""
        significant_gap = EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.SIGNIFICANT,
            description="Test gap",
        )
        # 1 significant gap but 10 well-characterized aspects
        result = _compute_overall_quality([significant_gap], 10)
        assert result == "comprehensive"

    def test_critical_gaps_hard_to_overcome(self):
        """Critical gaps in high-weight categories should be hard to overcome."""
        # VALIDATION has weight 3.5, critical multiplier 3.0 = 10.5 points
        critical_gap = EvidenceGap(
            category=GapCategory.VALIDATION,
            severity=GapSeverity.CRITICAL,
            description="Critical validation gap",
        )
        # 3 well-characterized = 4.5 credit, net = 6.0 -> limited
        result = _compute_overall_quality([critical_gap], 3)
        assert result in ("moderate", "limited", "minimal")

    def test_many_gaps_yields_poor_quality(self):
        """Many gaps should yield poor quality even with some well-characterized."""
        gaps = [
            EvidenceGap(category=GapCategory.CLINICAL, severity=GapSeverity.CRITICAL, description="Gap 1"),
            EvidenceGap(category=GapCategory.FUNCTIONAL, severity=GapSeverity.SIGNIFICANT, description="Gap 2"),
            EvidenceGap(category=GapCategory.PRECLINICAL, severity=GapSeverity.SIGNIFICANT, description="Gap 3"),
            EvidenceGap(category=GapCategory.VALIDATION, severity=GapSeverity.CRITICAL, description="Gap 4"),
        ]
        result = _compute_overall_quality(gaps, 3)
        assert result in ("limited", "minimal")

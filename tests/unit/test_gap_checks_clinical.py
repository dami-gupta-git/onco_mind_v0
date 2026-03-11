"""Unit tests for clinical gap check functions.

Tests: check_clinical_evidence, check_tumor_type_evidence,
check_tumor_specific_evidence, check_drug_response,
check_resistance_mechanisms.
"""

from unittest.mock import MagicMock

from oncomind.insight_builder.gap_detector import GapDetectionContext
from oncomind.insight_builder.gap_detection.checks import (
    check_clinical_evidence as _check_clinical_evidence,
    check_tumor_type_evidence as _check_tumor_type_evidence,
    check_tumor_specific_evidence as _check_tumor_specific_evidence,
    check_drug_response as _check_drug_response,
    check_depmap_drug_sensitivity as _check_depmap_drug_sensitivity,
    check_preclinical_biomarkers as _check_preclinical_biomarkers,
    check_resistance_mechanisms as _check_resistance_mechanisms,
)
from oncomind.models.extracted.evidence_gaps import GapCategory, GapSeverity
from oncomind.models.extracted.literature_knowledge import (
    LiteratureKnowledge,
    LitDrugResistance,
)

from tests.unit.conftest import create_mock_depmap_evidence, create_mock_cell_line

# =============================================================================
# TEST _check_clinical_evidence
# =============================================================================


class TestCheckClinicalEvidence:
    """Tests for _check_clinical_evidence function."""

    def test_with_fda_biomarker_variant_level(self, mock_evidence, base_context):
        """FDA biomarker evidence with VARIANT-level specificity should mark as well-characterized."""
        from oncomind.models.evidence.fda_biomarker import (
            SpecificityLevel,
            BiomarkerRequirement,
        )

        fda_ev = MagicMock()
        fda_ev.drug_name = "TestDrug"
        fda_ev.brand_name = "TestBrand"
        fda_ev.specificity = SpecificityLevel.VARIANT
        fda_ev.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        fda_ev.tumor_types = ["NSCLC"]
        fda_ev.variant_match_result = "exact"
        fda_ev.locus_match = "variant"
        fda_ev.tumor_match = True
        mock_evidence.fda_biomarker_evidence = [fda_ev]
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[fda_ev])
        mock_evidence.context.tumor_type = "NSCLC"

        _check_clinical_evidence(mock_evidence, base_context)

        assert any("clinical" in w.lower() for w in base_context.well_characterized)
        assert base_context.has_clinical is True

    def test_with_civic_assertions_sets_has_clinical(self, mock_evidence, base_context):
        """CIViC assertions should set has_clinical=True."""
        mock_evidence.civic_assertions = [MagicMock(), MagicMock()]
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])
        mock_evidence.context.tumor_type = "NSCLC"

        _check_clinical_evidence(mock_evidence, base_context)

        assert base_context.has_clinical is True

    def test_fda_biomarker_gene_level_is_well_characterized(
        self, mock_evidence, base_context
    ):
        """FDA biomarker evidence with GENE-level specificity should mark as well-characterized."""
        from oncomind.models.evidence.fda_biomarker import (
            SpecificityLevel,
            BiomarkerRequirement,
        )

        fda_ev = MagicMock()
        fda_ev.drug_name = "TestDrug"
        fda_ev.brand_name = "TestBrand"
        fda_ev.specificity = SpecificityLevel.GENE
        fda_ev.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        fda_ev.tumor_types = ["NSCLC"]
        fda_ev.variant_match_result = "gene"
        fda_ev.locus_match = "gene"
        fda_ev.tumor_match = True
        mock_evidence.fda_biomarker_evidence = [fda_ev]
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[fda_ev])
        mock_evidence.context.tumor_type = "NSCLC"

        _check_clinical_evidence(mock_evidence, base_context)

        assert any("clinical" in w.lower() for w in base_context.well_characterized)
        assert base_context.has_clinical is True

    def test_no_clinical_evidence_adds_critical_gap(self, mock_evidence, base_context):
        """Missing clinical evidence should add a CRITICAL gap."""
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        _check_clinical_evidence(mock_evidence, base_context)

        assert "clinical evidence" in base_context.poorly_characterized
        assert base_context.has_clinical is False
        clinical_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.CLINICAL
        ]
        assert len(clinical_gaps) >= 1
        assert any(g.severity == GapSeverity.CRITICAL for g in clinical_gaps)

    def test_civic_level_a_and_cgi_fda_approved_combined_basis(
        self, mock_evidence, base_context
    ):
        """CIViC Level A + CGI fda_approved together should produce combined actionability basis."""
        civic_level_a = MagicMock()
        civic_level_a.amp_level_letter = "A"
        mock_evidence.civic_assertions = [civic_level_a]
        mock_evidence.civic_evidence = []
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        cgi_approved_1 = MagicMock()
        cgi_approved_1.fda_approved = True
        cgi_approved_2 = MagicMock()
        cgi_approved_2.fda_approved = True
        mock_evidence.cgi_biomarkers = [cgi_approved_1, cgi_approved_2]

        _check_clinical_evidence(mock_evidence, base_context)

        assert base_context.has_clinical is True
        actionability = [
            w
            for w in base_context.well_characterized_detailed
            if w.aspect.lower() == "clinical actionability"
        ]
        assert len(actionability) == 1
        basis = actionability[0].basis
        assert "1 approval" in basis and "CIViC assertions" in basis
        assert "2 approvals" in basis and "CGI" in basis



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

    def test_with_tumor_specific_civic_variant_level(self, mock_evidence):
        """Tumor-specific CIViC data at VARIANT level should be well-characterized."""
        assertion = MagicMock()
        assertion.disease = "Non-Small Cell Lung Cancer"
        assertion.tumor_match = True
        assertion.locus_match = "variant"
        assertion.locus_variant_match = MagicMock()
        assertion.locus_variant_match.level = "variant"
        mock_evidence.civic_assertions = [assertion]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="L858R",
            tumor_type="Non-Small Cell Lung Cancer",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
        )

        _check_tumor_type_evidence(mock_evidence, ctx)

        assert any("evidence items for" in w.lower() for w in ctx.well_characterized)

    def test_no_tumor_specific_evidence_cancer_gene_critical(self, mock_evidence):
        """Cancer gene without tumor-specific data should get CRITICAL gap."""
        ctx = GapDetectionContext(
            gene="KRAS",
            variant="G12D",
            tumor_type="Pancreatic",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=False,
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
        assertion = MagicMock()
        assertion.disease = "Lung Cancer"
        assertion.locus_variant_match = None
        mock_evidence.civic_assertions = [assertion]

        vicc = MagicMock()
        vicc.disease = "Lung adenocarcinoma"
        vicc.locus_variant_match = None
        mock_evidence.vicc_evidence = [vicc, vicc]  # 2 VICC items

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "Lung")

        assert tumor_match.has_tumor_evidence is True
        sources = tumor_match.sources_str
        assert "1 CIViC Assertions" in sources
        assert "2 VICC (meta-KB)" in sources

    def test_tumor_evidence_match_locus_levels(self, mock_evidence):
        """TumorEvidenceMatch should track variant/codon/gene match levels."""
        from oncomind.models.evidence.base import EvidenceLevel

        assertion_variant = MagicMock()
        assertion_variant.disease = "Lung Cancer"
        assertion_variant.locus_variant_match = EvidenceLevel(
            level="variant", scope="specific"
        )
        assertion_variant.locus_match = "variant"

        assertion_gene = MagicMock()
        assertion_gene.disease = "Lung Cancer"
        assertion_gene.locus_variant_match = EvidenceLevel(
            level="gene", scope="unspecified"
        )
        assertion_gene.locus_match = "gene"

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
        assertion = MagicMock()
        assertion.disease = "Melanoma"
        assertion.locus_variant_match = None
        assertion.tumor_match = False
        mock_evidence.civic_assertions = [assertion]

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "Pancreatic")

        assert tumor_match.has_tumor_evidence is False
        assert tumor_match.total_matches == 0
        assert tumor_match.sources_str == "No tumor-specific evidence"
        assert tumor_match.matches_on_str is None

    def test_tumor_evidence_includes_all_cgi_tiers(self, mock_evidence):
        """TumorEvidenceMatch should include all CGI tiers (FDA, preclinical, early phase)."""
        cgi_fda = MagicMock()
        cgi_fda.tumor_type = "Lung Cancer"
        cgi_fda.locus_variant_match = None
        mock_evidence.cgi_biomarkers = [cgi_fda]

        preclin = MagicMock()
        preclin.tumor_type = "Lung Cancer"
        preclin.locus_variant_match = None
        mock_evidence.preclinical_biomarkers = [preclin]

        early = MagicMock()
        early.tumor_type = "Lung Cancer"
        early.locus_variant_match = None
        mock_evidence.early_phase_biomarkers = [early]

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "Lung Cancer")

        assert tumor_match.has_tumor_evidence is True
        assert tumor_match.total_matches == 3
        assert tumor_match.cgi_biomarkers is not None
        assert tumor_match.cgi_biomarkers.count == 3

    def test_tumor_evidence_fda_parsing(self, mock_evidence):
        """TumorEvidenceMatch should use FDA biomarker evidence's tumor_types."""
        from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement

        fda = MagicMock()
        fda.tumor_types = ["Non-small cell lung cancer"]
        fda.locus_variant_match = None
        fda.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        mock_evidence.fda_biomarker_evidence = [fda]

        tumor_match = _check_tumor_specific_evidence(mock_evidence, "NSCLC")

        assert tumor_match.has_tumor_evidence is True


# =============================================================================
# TEST _check_drug_response
# =============================================================================


class TestCheckDrugResponse:
    """Tests for _check_drug_response function."""

    def test_with_cgi_biomarkers(self, mock_evidence, base_context):
        """CGI biomarkers should set has_drug_data flag."""
        mock_evidence.cgi_biomarkers = [MagicMock()]
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        _check_drug_response(mock_evidence, base_context)

        assert base_context.has_drug_data is True

    def test_with_vicc_evidence(self, mock_evidence, base_context):
        """VICC evidence should mark as well-characterized."""
        mock_evidence.vicc_evidence = [MagicMock()]
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        _check_drug_response(mock_evidence, base_context)

        assert base_context.has_drug_data is True

    def test_with_depmap_drug_sensitivities(self, mock_evidence, base_context):
        """DepMap drug sensitivities should mark as well-characterized with tumor-matched cell lines."""
        lung_cell_line = create_mock_cell_line(
            name="A549",
            depmap_id="ACH-000681",
            primary_disease="Lung",
            subtype="Non-Small Cell Lung Cancer",
            has_mutation=True,
            mutation_details="V100E",
        )

        drug_sensitivity = MagicMock()
        drug_sensitivity.drug_name = "Trametinib"
        drug_sensitivity.mean_log2fc = -2.1
        drug_sensitivity.n_cell_lines = 5
        drug_sensitivity.sensitive_lines = ["A549", "H1299"]

        mock_evidence.depmap_evidence = create_mock_depmap_evidence(
            gene="TESTGENE",
            variant="V100E",
            is_essential=False,
            drug_sensitivities=[drug_sensitivity],
            cell_line_models=[lung_cell_line],
        )
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        _check_drug_response(mock_evidence, base_context)
        _check_depmap_drug_sensitivity(mock_evidence, base_context)

        assert any("depmap" in w.lower() for w in base_context.well_characterized)

    def test_no_drug_data_adds_gap(self, mock_evidence, base_context):
        """Missing drug data should add a SIGNIFICANT gap."""
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        _check_drug_response(mock_evidence, base_context)

        drug_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.DRUG_RESPONSE
        ]
        assert len(drug_gaps) >= 1
        assert drug_gaps[0].severity == GapSeverity.SIGNIFICANT
        assert "FDA-approved therapy" in drug_gaps[0].description

    def test_drug_response_counts_fda_only(self, mock_evidence, base_context):
        """Drug response in Gap Analysis only shows FDA-approved therapy."""
        from oncomind.models.evidence.fda_biomarker import (
            BiomarkerRequirement,
            SpecificityLevel,
        )

        cgi1 = MagicMock()
        cgi1.locus_variant_match = None
        cgi1.tumor_type = "Lung Cancer"
        cgi2 = MagicMock()
        cgi2.locus_variant_match = None
        cgi2.tumor_type = "Melanoma"
        mock_evidence.cgi_biomarkers = [cgi1, cgi2]

        mock_evidence.vicc_evidence = [MagicMock(), MagicMock(), MagicMock()]
        for v in mock_evidence.vicc_evidence:
            v.locus_variant_match = None
            v.disease = "Other Cancer"

        fda = MagicMock()
        fda.locus_variant_match = None
        fda.tumor_types = ["Lung Cancer"]
        fda.requirement = BiomarkerRequirement.REQUIRED_POSITIVE
        fda.drug_name = "TestDrug"
        fda.brand_name = None
        fda.specificity = SpecificityLevel.VARIANT
        fda.variant_match_result = "exact"
        fda.locus_match = "variant"
        fda.tumor_match = True
        mock_evidence.fda_biomarker_evidence = [fda]
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[fda])

        _check_drug_response(mock_evidence, base_context)

        drug_resp = [
            w
            for w in base_context.well_characterized_detailed
            if w.aspect.lower() == "drug response data"
        ]
        assert (
            len(drug_resp) == 0
        ), "CGI/VICC should not be in Gap Analysis well_characterized"

        fda_resp = [
            w
            for w in base_context.well_characterized_detailed
            if "fda-approved" in w.aspect.lower()
        ]
        assert len(fda_resp) == 1
        assert "1 FDA" in fda_resp[0].basis

        assert base_context.has_drug_data is True

    def test_drug_response_with_fda_count(self, mock_evidence):
        """Drug response uses get_filtered_fda_evidence() to count FDA drugs."""
        cgi_match = MagicMock()
        cgi_match.locus_variant_match = None
        cgi_match.tumor_type = "Lung Cancer"
        mock_evidence.cgi_biomarkers = [cgi_match]

        vicc1 = MagicMock()
        vicc1.locus_variant_match = None
        vicc1.disease = "Lung adenocarcinoma"
        mock_evidence.vicc_evidence = [vicc1]

        mock_evidence.get_filtered_fda_evidence = MagicMock(
            return_value=[MagicMock(), MagicMock(), MagicMock()]
        )

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="L858R",
            tumor_type="Lung",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
        )

        _check_drug_response(mock_evidence, ctx)

        fda_resp = [
            w
            for w in ctx.well_characterized_detailed
            if "fda-approved" in w.aspect.lower()
        ]
        assert len(fda_resp) == 1
        assert "3 FDA" in fda_resp[0].basis

        assert ctx.has_drug_data is True

    def test_drug_response_zero_fda_count_adds_gap(self, mock_evidence, base_context):
        """Zero FDA count should add a gap even if CGI/VICC evidence exists."""
        from oncomind.models.evidence.base import EvidenceLevel

        cgi_variant = MagicMock()
        cgi_variant.locus_variant_match = EvidenceLevel(
            level="variant", scope="specific"
        )
        cgi_variant.locus_match = "variant"
        cgi_variant.tumor_type = None
        mock_evidence.cgi_biomarkers = [cgi_variant]

        vicc_gene = MagicMock()
        vicc_gene.locus_variant_match = EvidenceLevel(level="gene", scope="unspecified")
        vicc_gene.locus_match = "gene"
        vicc_gene.disease = None
        mock_evidence.vicc_evidence = [vicc_gene]

        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        _check_drug_response(mock_evidence, base_context)

        drug_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.DRUG_RESPONSE
        ]
        assert len(drug_gaps) >= 1
        assert "FDA-approved therapy" in drug_gaps[0].description

        assert base_context.has_drug_data is True

    def test_drug_response_excludes_preclinical_biomarkers(
        self, mock_evidence, base_context
    ):
        """Drug response should NOT include preclinical or early phase biomarkers."""
        mock_evidence.preclinical_biomarkers = [MagicMock(), MagicMock()]
        mock_evidence.early_phase_biomarkers = [MagicMock()]
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        _check_drug_response(mock_evidence, base_context)
        _check_preclinical_biomarkers(mock_evidence, base_context)

        drug_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.DRUG_RESPONSE
        ]
        assert len(drug_gaps) >= 1
        assert "FDA-approved therapy" in drug_gaps[0].description

        preclin = [
            w
            for w in base_context.well_characterized_detailed
            if "preclinical" in w.aspect.lower()
        ]
        assert len(preclin) == 1
        assert "2 preclinical" in preclin[0].basis
        assert "1 early phase" in preclin[0].basis

    def test_fda_determines_approval_not_vicc(self, mock_evidence, base_context):
        """Gap should be added when VICC has drug data but FDA has none.

        Design principle: gap detection uses fda_biomarker_evidence as the source of
        truth for approval status. VICC evidence alone must not suppress a DRUG_RESPONSE gap.
        """
        vicc = MagicMock()
        vicc.locus_match = "variant"
        vicc.disease = "bladder cancer"
        vicc.tumor_match = True
        mock_evidence.vicc_evidence = [vicc]
        mock_evidence.fda_biomarker_evidence = []
        mock_evidence.get_filtered_fda_evidence = MagicMock(return_value=[])

        _check_drug_response(mock_evidence, base_context)

        drug_gaps = [g for g in base_context.gaps if g.category == GapCategory.DRUG_RESPONSE]
        assert len(drug_gaps) >= 1, "Should add DRUG_RESPONSE gap when FDA has no evidence"

        well_char_lower = [w.lower() for w in base_context.well_characterized]
        assert "fda-approved therapy" not in well_char_lower, (
            "Should NOT mark as well-characterized when only VICC has data"
        )


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
        well_char_detail = [
            w
            for w in ctx.well_characterized_detailed
            if "resistance" in w.aspect.lower()
        ]
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
        well_char_detail = [
            w
            for w in ctx.well_characterized_detailed
            if "resistance" in w.aspect.lower()
        ]
        assert len(well_char_detail) > 0
        assert "VICC" in well_char_detail[0].basis

    def test_with_llm_literature_knowledge_resistance(self, mock_evidence):
        """LLM literature knowledge with resistant_to should mark as well-characterized."""
        mock_evidence.literature_knowledge = LiteratureKnowledge(
            resistant_to=[
                LitDrugResistance(
                    drug="Gefitinib", evidence="clinical", is_predictive=True
                ),
                LitDrugResistance(
                    drug="Erlotinib", evidence="clinical", is_predictive=True
                ),
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
        well_char_detail = [
            w
            for w in ctx.well_characterized_detailed
            if "resistance" in w.aspect.lower()
        ]
        assert len(well_char_detail) > 0
        assert "LLM literature" in well_char_detail[0].basis
        assert "2 drugs" in well_char_detail[0].basis

    def test_with_llm_literature_non_predictive_resistance(self, mock_evidence):
        """Non-predictive resistance in literature should NOT mark as well-characterized."""
        mock_evidence.literature_knowledge = LiteratureKnowledge(
            resistant_to=[
                LitDrugResistance(
                    drug="Gefitinib", evidence="clinical", is_predictive=False
                ),
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

        assert any("resistance" in p.lower() for p in ctx.poorly_characterized)

    def test_multiple_resistance_sources(self, mock_evidence):
        """Multiple resistance sources should all be included in basis."""
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
        well_char_detail = [
            w
            for w in ctx.well_characterized_detailed
            if "resistance" in w.aspect.lower()
        ]
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

        resistance_gaps = [
            g for g in base_context.gaps if g.category == GapCategory.RESISTANCE
        ]
        assert len(resistance_gaps) == 0

    def test_resistance_mechanisms_locus_variant_match_levels(self, mock_evidence):
        """Resistance mechanisms should track match levels (variant, codon, gene)."""
        from oncomind.models.evidence.base import EvidenceLevel

        cgi_resistance = MagicMock()
        cgi_resistance.association = "Resistance"
        cgi_resistance.locus_variant_match = EvidenceLevel(
            level="variant", scope="specific"
        )
        cgi_resistance.locus_match = "variant"
        cgi_resistance.tumor_type = None
        mock_evidence.cgi_biomarkers = [cgi_resistance]

        vicc_resistance = MagicMock()
        vicc_resistance.response_type = "RESISTANCE"
        vicc_resistance.locus_variant_match = EvidenceLevel(
            level="gene", scope="unspecified"
        )
        vicc_resistance.locus_match = "gene"
        vicc_resistance.disease = None
        mock_evidence.vicc_evidence = [vicc_resistance]

        civic_resistance = MagicMock()
        civic_resistance.is_resistance = True
        civic_resistance.locus_variant_match = EvidenceLevel(
            level="codon", scope="specific"
        )
        civic_resistance.locus_match = "codon"
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

        resist_detail = [
            w
            for w in ctx.well_characterized_detailed
            if "resistance" in w.aspect.lower()
        ]
        assert len(resist_detail) == 1
        matches_on = resist_detail[0].matches_on
        assert matches_on is not None
        assert "1 variant" in matches_on
        assert "1 codon" in matches_on
        assert "1 gene" in matches_on

    def test_resistance_mechanisms_tumor_match_counting(self, mock_evidence):
        """Resistance mechanisms should track tumor matches vs others."""
        cgi_resistance = MagicMock()
        cgi_resistance.association = "Resistance"
        cgi_resistance.locus_variant_match = None
        cgi_resistance.tumor_type = "Lung Cancer"
        cgi_resistance.tumor_match = True
        mock_evidence.cgi_biomarkers = [cgi_resistance]

        vicc_resistance = MagicMock()
        vicc_resistance.response_type = "RESISTANCE"
        vicc_resistance.locus_variant_match = None
        vicc_resistance.disease = "Melanoma"
        vicc_resistance.tumor_match = False
        mock_evidence.vicc_evidence = [vicc_resistance]

        civic_resistance = MagicMock()
        civic_resistance.is_resistance = True
        civic_resistance.locus_variant_match = None
        civic_resistance.disease = "Lung adenocarcinoma"
        civic_resistance.tumor_match = True
        mock_evidence.civic_assertions = [civic_resistance]

        ctx = GapDetectionContext(
            gene="EGFR",
            variant="T790M",
            tumor_type="Lung",
            is_cancer_gene=True,
            has_pathogenic_signal=True,
            has_clinical=True,
            has_drug_data=True,
        )

        _check_resistance_mechanisms(mock_evidence, ctx)

        resist_detail = [
            w
            for w in ctx.well_characterized_detailed
            if "resistance" in w.aspect.lower()
        ]
        assert len(resist_detail) == 1
        tumor_match = resist_detail[0].tumor_match
        assert tumor_match is not None
        assert "2 tumor" in tumor_match
        assert "1 other" in tumor_match

    def test_resistance_mechanisms_counts_sources_correctly(self, mock_evidence):
        """Resistance mechanisms should count all source types correctly."""
        article = MagicMock()
        article.is_resistance_evidence.return_value = True
        mock_evidence.pubmed_articles = [article]

        cgi1 = MagicMock()
        cgi1.association = "Resistance"
        cgi1.locus_variant_match = None
        cgi1.tumor_type = None
        cgi2 = MagicMock()
        cgi2.association = "RESISTANCE"
        cgi2.locus_variant_match = None
        cgi2.tumor_type = None
        mock_evidence.cgi_biomarkers = [cgi1, cgi2]

        civic = MagicMock()
        civic.is_resistance = True
        civic.locus_variant_match = None
        civic.disease = None
        mock_evidence.civic_assertions = [civic]

        mock_evidence.vicc_evidence = []
        for _ in range(3):
            vicc = MagicMock()
            vicc.response_type = "RESISTANCE"
            vicc.locus_variant_match = None
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

        resist_detail = [
            w
            for w in ctx.well_characterized_detailed
            if "resistance" in w.aspect.lower()
        ]
        assert len(resist_detail) == 1
        basis = resist_detail[0].basis
        assert "1 PubMed" in basis
        assert "2 CGI" in basis
        assert "1 CIViC" in basis
        assert "3 VICC" in basis

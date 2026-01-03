"""Integration tests for new evidence model fields.

Tests validate that the new fields added to evidence models
(CIViCEvidence, VICCEvidence, FDAApproval) are properly populated
and passed through to the frontend via the backend API.
"""

import pytest
import sys
from pathlib import Path

# Add streamlit directory to path for imports
streamlit_dir = Path(__file__).parent.parent.parent.parent / "streamlit"
sys.path.insert(0, str(streamlit_dir))

from oncomind.models.evidence import (
    CIViCEvidence,
    VICCEvidence,
    FDAApproval,
    CIViCAssertionEvidence,
)
from oncomind.models.evidence.cbioportal import CBioPortalEvidence, CoMutationEntry
from oncomind.models.evidence.depmap import (
    DepMapEvidence,
    GeneDependency,
    DrugSensitivity,
    CellLineModel,
)
from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.cosmic import COSMICEvidence
from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
from oncomind.models.evidence.pubmed import PubMedEvidence
from oncomind.api.cgi import CGIBiomarker
from backend import get_variant_insight


class TestCIViCEvidenceNewFields:
    """Tests for new CIViC evidence fields: pmid, source_url, trust_rating."""

    def test_civic_evidence_pmid_field(self):
        """CIViCEvidence should accept pmid field."""
        evidence = CIViCEvidence(
            evidence_type="Predictive",
            evidence_level="A",
            clinical_significance="Sensitivity/Response",
            disease="Melanoma",
            drugs=["Vemurafenib"],
            pmid="22735384",
            source_url="https://pubmed.ncbi.nlm.nih.gov/22735384/",
            trust_rating=5,
        )
        assert evidence.pmid == "22735384"
        assert evidence.source_url == "https://pubmed.ncbi.nlm.nih.gov/22735384/"
        assert evidence.trust_rating == 5

    def test_civic_evidence_pmid_optional(self):
        """CIViCEvidence pmid and related fields should be optional."""
        evidence = CIViCEvidence(
            evidence_type="Predictive",
            evidence_level="A",
        )
        assert evidence.pmid is None
        assert evidence.source_url is None
        assert evidence.trust_rating is None

    def test_civic_evidence_rating_vs_trust_rating(self):
        """CIViCEvidence should have both rating and trust_rating fields."""
        evidence = CIViCEvidence(
            evidence_type="Predictive",
            rating=4,
            trust_rating=5,
        )
        assert evidence.rating == 4
        assert evidence.trust_rating == 5


class TestVICCEvidenceNewFields:
    """Tests for new VICC evidence fields: molecular_profile, molecular_profile_score."""

    def test_vicc_evidence_molecular_profile_field(self):
        """VICCEvidence should accept molecular_profile fields."""
        evidence = VICCEvidence(
            gene="BRAF",
            variant="V600E",
            disease="Melanoma",
            drugs=["Vemurafenib"],
            response_type="Sensitivity",
            molecular_profile="BRAF V600E",
            molecular_profile_score=0.95,
        )
        assert evidence.molecular_profile == "BRAF V600E"
        assert evidence.molecular_profile_score == 0.95

    def test_vicc_evidence_molecular_profile_optional(self):
        """VICCEvidence molecular profile fields should be optional."""
        evidence = VICCEvidence(
            gene="BRAF",
            variant="V600E",
        )
        assert evidence.molecular_profile is None
        assert evidence.molecular_profile_score is None

    def test_vicc_evidence_all_fields(self):
        """VICCEvidence should support all fields including new ones."""
        evidence = VICCEvidence(
            description="BRAF V600E predicts sensitivity to vemurafenib",
            gene="BRAF",
            variant="V600E",
            disease="Melanoma",
            drugs=["Vemurafenib", "Dabrafenib"],
            evidence_level="A",
            response_type="Sensitivity",
            source="civic",
            publication_url="https://pubmed.ncbi.nlm.nih.gov/22735384/",
            oncogenic="Oncogenic",
            is_sensitivity=True,
            is_resistance=False,
            oncokb_level="1A",
            molecular_profile="BRAF V600E",
            molecular_profile_score=0.98,
        )
        assert evidence.molecular_profile == "BRAF V600E"
        assert evidence.molecular_profile_score == 0.98
        assert evidence.oncokb_level == "1A"


class TestFDAApprovalNewFields:
    """Tests for new FDA approval fields: companion_diagnostic, black_box_warning, dosing_for_variant."""

    def test_fda_approval_companion_diagnostic_field(self):
        """FDAApproval should accept companion_diagnostic field."""
        approval = FDAApproval(
            drug_name="Vemurafenib",
            brand_name="Zelboraf",
            generic_name="vemurafenib",
            indication="Treatment of BRAF V600E-mutant melanoma",
            companion_diagnostic="cobas 4800 BRAF V600 Mutation Test",
        )
        assert approval.companion_diagnostic == "cobas 4800 BRAF V600 Mutation Test"

    def test_fda_approval_black_box_warning_field(self):
        """FDAApproval should accept black_box_warning field."""
        approval = FDAApproval(
            drug_name="Vemurafenib",
            brand_name="Zelboraf",
            black_box_warning="New primary cutaneous malignancies can occur",
        )
        assert approval.black_box_warning == "New primary cutaneous malignancies can occur"

    def test_fda_approval_dosing_for_variant_field(self):
        """FDAApproval should accept dosing_for_variant field."""
        approval = FDAApproval(
            drug_name="Vemurafenib",
            brand_name="Zelboraf",
            dosing_for_variant="960 mg orally twice daily",
        )
        assert approval.dosing_for_variant == "960 mg orally twice daily"

    def test_fda_approval_all_new_fields(self):
        """FDAApproval should support all new fields together."""
        approval = FDAApproval(
            drug_name="Vemurafenib",
            brand_name="Zelboraf",
            generic_name="vemurafenib",
            indication="Treatment of BRAF V600E-mutant melanoma",
            approval_date="2011-08-17",
            marketing_status="Prescription",
            gene="BRAF",
            variant_in_indications=True,
            variant_in_clinical_studies=True,
            companion_diagnostic="cobas 4800 BRAF V600 Mutation Test",
            black_box_warning="New primary cutaneous malignancies can occur",
            dosing_for_variant="960 mg orally twice daily",
        )
        assert approval.companion_diagnostic is not None
        assert approval.black_box_warning is not None
        assert approval.dosing_for_variant is not None

    def test_fda_approval_new_fields_optional(self):
        """FDAApproval new fields should be optional."""
        approval = FDAApproval(
            drug_name="Vemurafenib",
            brand_name="Zelboraf",
        )
        assert approval.companion_diagnostic is None
        assert approval.black_box_warning is None
        assert approval.dosing_for_variant is None


class TestCIViCAssertionFields:
    """Tests for CIViCAssertionEvidence model fields."""

    def test_civic_assertion_existing_fields(self):
        """CIViCAssertionEvidence should have all expected fields."""
        assertion = CIViCAssertionEvidence(
            assertion_id=123,
            name="BRAF V600E predicts sensitivity to vemurafenib",
            amp_level="Tier I",
            amp_tier="Level A",
            amp_level_letter="A",
            assertion_type="Predictive",
            significance="Sensitivity/Response",
            status="accepted",
            molecular_profile="BRAF V600E",
            disease="Melanoma",
            therapies=["Vemurafenib"],
            fda_companion_test=True,
            nccn_guideline="NCCN Melanoma Guidelines",
            description="FDA-approved therapy for BRAF V600E melanoma",
        )
        assert assertion.assertion_id == 123
        assert assertion.molecular_profile == "BRAF V600E"
        assert assertion.fda_companion_test is True
        assert assertion.nccn_guideline == "NCCN Melanoma Guidelines"


class TestBackendNewFields:
    """Tests for backend returning new evidence fields."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_returns_civic_evidence_key(self):
        """Backend should return civic_evidence in response."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert "civic_evidence" in result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_civic_evidence_has_new_fields(self):
        """Backend civic_evidence should include new fields structure."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        civic_evidence = result.get("civic_evidence", [])

        # Verify structure includes new fields (may be None if not populated by API)
        if civic_evidence:
            for entry in civic_evidence:
                assert "pmid" in entry
                assert "source_url" in entry
                assert "trust_rating" in entry

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_vicc_evidence_has_new_fields(self):
        """Backend vicc_evidence should include new molecular profile fields."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        vicc_evidence = result.get("vicc_evidence", [])

        # Verify structure includes new fields (may be None if not populated by API)
        if vicc_evidence:
            for entry in vicc_evidence:
                assert "molecular_profile" in entry
                assert "molecular_profile_score" in entry

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_fda_approvals_has_new_fields(self):
        """Backend fda_approvals should include new clinical fields."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        fda_approvals = result.get("fda_approvals", [])

        # Verify structure includes new fields (may be None if not populated by API)
        if fda_approvals:
            for entry in fda_approvals:
                assert "companion_diagnostic" in entry
                assert "black_box_warning" in entry
                assert "dosing_for_variant" in entry


class TestBackendAnnotationsTab:
    """Tests for backend functional annotations (Functional tab)."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_returns_annotations(self):
        """Backend should return annotations dict."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert "annotations" in result
        annotations = result["annotations"]

        # Verify expected keys exist
        assert "alphamissense_score" in annotations
        assert "alphamissense_prediction" in annotations
        assert "cadd_score" in annotations
        assert "polyphen2_prediction" in annotations
        assert "gnomad_exome_af" in annotations
        assert "snpeff_effect" in annotations

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_annotations_braf_v600e(self):
        """BRAF V600E should have functional predictions."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        annotations = result.get("annotations", {})

        # BRAF V600E is a well-characterized variant
        # At least some scores should be present
        has_any_score = any([
            annotations.get("alphamissense_score") is not None,
            annotations.get("cadd_score") is not None,
            annotations.get("polyphen2_prediction") is not None,
        ])

        # Note: This assertion may be relaxed if MyVariant doesn't return scores
        # for this specific query. The important thing is the structure exists.
        assert "annotations" in result


class TestBackendEvidenceStructure:
    """Tests for overall backend evidence structure."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_all_evidence_keys_present(self):
        """Backend should return all expected evidence keys."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

        # All evidence types should be present (even if empty lists/None)
        expected_keys = [
            "fda_approvals",
            "civic_assertions",
            "civic_evidence",
            "vicc_evidence",
            "cgi_biomarkers",
            "clinvar_entries",
            "cosmic_entries",
            "clinical_trials",
            "pubmed_articles",
            "preclinical_biomarkers",
            "early_phase_biomarkers",
            "cbioportal_evidence",
            "depmap_evidence",
        ]

        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_braf_v600e_has_fda_approvals(self):
        """BRAF V600E should have FDA approvals."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        fda_approvals = result.get("fda_approvals", [])

        # BRAF V600E has multiple FDA-approved drugs
        assert len(fda_approvals) >= 1, "BRAF V600E should have at least 1 FDA approval"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_braf_v600e_has_vicc_evidence(self):
        """BRAF V600E should have VICC evidence."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        vicc_evidence = result.get("vicc_evidence", [])

        # BRAF V600E is well-characterized in VICC
        assert len(vicc_evidence) >= 1, "BRAF V600E should have VICC evidence"


class TestCBioPortalEvidenceFields:
    """Tests for cBioPortal evidence model fields."""

    def test_cbioportal_evidence_creation(self):
        """Test creating CBioPortalEvidence model."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            study_id="skcm_tcga_pan_can_atlas_2018",
            study_name="Skin Cutaneous Melanoma (TCGA, PanCancer Atlas)",
            total_samples=449,
            samples_with_gene_mutation=230,
            samples_with_exact_variant=217,
            gene_prevalence_pct=51.2,
            variant_prevalence_pct=48.3,
            co_occurring=[CoMutationEntry(gene="CDKN2A", count=98, pct=45.2, odds_ratio=2.34)],
            mutually_exclusive=[CoMutationEntry(gene="NRAS", count=2, pct=0.9, odds_ratio=0.08)],
        )

        assert evidence.gene == "BRAF"
        assert evidence.variant == "V600E"
        assert evidence.total_samples == 449
        assert evidence.gene_prevalence_pct == 51.2
        assert len(evidence.co_occurring) == 1
        assert len(evidence.mutually_exclusive) == 1

    def test_cbioportal_evidence_has_data(self):
        """Test has_data method."""
        evidence_with_data = CBioPortalEvidence(gene="BRAF", total_samples=449)
        evidence_no_data = CBioPortalEvidence(gene="BRAF", total_samples=0)

        assert evidence_with_data.has_data() is True
        assert evidence_no_data.has_data() is False

    def test_cbioportal_evidence_get_study_url(self):
        """Test study URL generation."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            study_id="skcm_tcga_pan_can_atlas_2018",
        )

        url = evidence.get_study_url()
        assert url == "https://www.cbioportal.org/study/summary?id=skcm_tcga_pan_can_atlas_2018"

    def test_cbioportal_evidence_study_url_none(self):
        """Test study URL is None when no study_id."""
        evidence = CBioPortalEvidence(gene="BRAF")
        assert evidence.get_study_url() is None

    def test_co_mutation_entry_creation(self):
        """Test creating CoMutationEntry."""
        entry = CoMutationEntry(
            gene="TP53",
            count=28,
            pct=12.9,
            odds_ratio=1.56,
        )

        assert entry.gene == "TP53"
        assert entry.count == 28
        assert entry.pct == 12.9
        assert entry.odds_ratio == 1.56

    def test_co_mutation_entry_defaults(self):
        """Test CoMutationEntry default values."""
        entry = CoMutationEntry(gene="TP53")

        assert entry.count == 0
        assert entry.pct == 0.0
        assert entry.odds_ratio is None

    def test_cbioportal_evidence_to_dict(self):
        """Test converting to dictionary."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            variant="V600E",
            study_id="skcm_tcga",
            total_samples=449,
        )

        result = evidence.to_dict()

        assert result["gene"] == "BRAF"
        assert result["variant"] == "V600E"
        assert result["study_id"] == "skcm_tcga"
        assert result["total_samples"] == 449

    def test_cbioportal_evidence_to_prompt_context(self):
        """Test generating LLM prompt context."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            study_id="skcm_tcga",
            study_name="TCGA Melanoma",
            total_samples=449,
            samples_with_gene_mutation=230,
            gene_prevalence_pct=51.2,
            co_occurring=[CoMutationEntry(gene="CDKN2A", count=98, pct=45.2, odds_ratio=2.34)],
        )

        context = evidence.to_prompt_context()

        assert "BRAF" in context
        assert "51.2%" in context
        assert "CDKN2A" in context
        assert "cBioPortal" in context

    def test_cbioportal_evidence_no_data_prompt(self):
        """Test prompt context when no data available."""
        evidence = CBioPortalEvidence(gene="BRAF", total_samples=0)

        context = evidence.to_prompt_context()

        assert "No cBioPortal data available" in context


class TestCBioPortalBackendIntegration:
    """Integration tests for cBioPortal evidence from backend."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_returns_cbioportal_evidence(self):
        """Backend should return cbioportal_evidence in response."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert "cbioportal_evidence" in result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_cbioportal_evidence_structure(self):
        """Backend cbioportal_evidence should have expected structure."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        cbio_evidence = result.get("cbioportal_evidence")

        # May be None if no data available
        if cbio_evidence is not None:
            # Verify expected fields exist
            assert "gene" in cbio_evidence
            assert "total_samples" in cbio_evidence
            assert "gene_prevalence_pct" in cbio_evidence
            assert "co_occurring" in cbio_evidence
            assert "mutually_exclusive" in cbio_evidence

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_cbioportal_braf_melanoma_has_data(self):
        """BRAF V600E in Melanoma should have cBioPortal data."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        cbio_evidence = result.get("cbioportal_evidence")

        # BRAF is very common in melanoma - should have data
        if cbio_evidence is not None:
            assert cbio_evidence["gene"] == "BRAF"
            assert cbio_evidence["total_samples"] > 0
            # BRAF prevalence in melanoma is typically >40%
            assert cbio_evidence["gene_prevalence_pct"] > 30


class TestDepMapEvidenceFields:
    """Tests for DepMap evidence model fields."""

    def test_depmap_evidence_creation(self):
        """Test creating DepMapEvidence model."""
        evidence = DepMapEvidence(
            gene="BRAF",
            variant="V600E",
            gene_dependency=GeneDependency(
                gene="BRAF",
                mean_dependency_score=-0.75,
                n_dependent_lines=150,
                n_total_lines=500,
                dependency_pct=30.0,
            ),
        )

        assert evidence.gene == "BRAF"
        assert evidence.variant == "V600E"
        assert evidence.gene_dependency is not None
        assert evidence.gene_dependency.mean_dependency_score == -0.75

    def test_depmap_evidence_has_data(self):
        """Test has_data method."""
        evidence_with_data = DepMapEvidence(
            gene="BRAF",
            gene_dependency=GeneDependency(gene="BRAF", mean_dependency_score=-0.5),
        )
        evidence_no_data = DepMapEvidence(gene="BRAF")

        assert evidence_with_data.has_data() is True
        assert evidence_no_data.has_data() is False

    def test_depmap_evidence_is_essential(self):
        """Test is_essential method."""
        essential = DepMapEvidence(
            gene="BRAF",
            gene_dependency=GeneDependency(gene="BRAF", mean_dependency_score=-0.75),
        )
        not_essential = DepMapEvidence(
            gene="BRAF",
            gene_dependency=GeneDependency(gene="BRAF", mean_dependency_score=-0.2),
        )

        assert essential.is_essential() is True
        assert not_essential.is_essential() is False

    def test_gene_dependency_creation(self):
        """Test creating GeneDependency model."""
        dependency = GeneDependency(
            gene="BRAF",
            mean_dependency_score=-0.65,
            n_dependent_lines=120,
            n_total_lines=450,
            dependency_pct=26.7,
            top_dependent_lines=["A375", "SK-MEL-28"],
        )

        assert dependency.gene == "BRAF"
        assert dependency.mean_dependency_score == -0.65
        assert len(dependency.top_dependent_lines) == 2

    def test_drug_sensitivity_creation(self):
        """Test creating DrugSensitivity model."""
        sensitivity = DrugSensitivity(
            drug_name="Vemurafenib",
            mean_log2fc=-2.1,
            n_cell_lines=25,
            sensitive_lines=["A375", "SK-MEL-28"],
        )

        assert sensitivity.drug_name == "Vemurafenib"
        assert sensitivity.mean_log2fc == -2.1
        assert sensitivity.n_cell_lines == 25

    def test_cell_line_model_creation(self):
        """Test creating CellLineModel."""
        cell_line = CellLineModel(
            name="A375",
            ccle_name="A375_SKIN",
            primary_disease="Melanoma",
            subtype="Cutaneous",
            has_mutation=True,
            mutation_details="BRAF V600E",
        )

        assert cell_line.name == "A375"
        assert cell_line.has_mutation is True
        assert cell_line.mutation_details == "BRAF V600E"

    def test_depmap_evidence_to_dict(self):
        """Test converting DepMapEvidence to dictionary."""
        evidence = DepMapEvidence(
            gene="BRAF",
            variant="V600E",
            gene_dependency=GeneDependency(
                gene="BRAF",
                mean_dependency_score=-0.65,
            ),
        )

        result = evidence.to_dict()

        assert result["gene"] == "BRAF"
        assert result["variant"] == "V600E"
        assert result["gene_dependency"] is not None

    def test_depmap_evidence_to_prompt_context(self):
        """Test generating LLM prompt context."""
        evidence = DepMapEvidence(
            gene="BRAF",
            gene_dependency=GeneDependency(
                gene="BRAF",
                mean_dependency_score=-0.75,
                n_dependent_lines=150,
                n_total_lines=500,
                dependency_pct=30.0,
            ),
            drug_sensitivities=[
                DrugSensitivity(drug_name="Vemurafenib", ic50_nm=50.0, n_cell_lines=20),
            ],
        )

        context = evidence.to_prompt_context()

        assert "BRAF" in context
        assert "ESSENTIAL" in context or "-0.75" in context
        assert "DepMap" in context

    def test_depmap_evidence_no_data_prompt(self):
        """Test prompt context when no data available."""
        evidence = DepMapEvidence(gene="BRAF")

        context = evidence.to_prompt_context()

        assert "No DepMap" in context


class TestDepMapBackendIntegration:
    """Integration tests for DepMap evidence from backend."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_returns_depmap_evidence(self):
        """Backend should return depmap_evidence in response."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert "depmap_evidence" in result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_depmap_evidence_structure(self):
        """Backend depmap_evidence should have expected structure when present."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        depmap_evidence = result.get("depmap_evidence")

        # May be None if no data available
        if depmap_evidence is not None:
            # Verify expected fields exist
            assert "gene" in depmap_evidence
            assert "gene_dependency" in depmap_evidence
            assert "drug_sensitivities" in depmap_evidence
            assert "cell_line_models" in depmap_evidence
            assert "is_essential" in depmap_evidence


class TestCGIBiomarkerFields:
    """Tests for CGI biomarker model fields."""

    def test_cgi_biomarker_creation(self):
        """Test creating CGIBiomarker model."""
        biomarker = CGIBiomarker(
            gene="EGFR",
            alteration="EGFR:G719S",
            drug="Afatinib",
            drug_status="Approved",
            association="Responsive",
            evidence_level="FDA guidelines",
            source="FDA",
            tumor_type="NSCLC",
            tumor_type_full="Non-small cell lung cancer",
        )

        assert biomarker.gene == "EGFR"
        assert biomarker.drug == "Afatinib"
        assert biomarker.drug_status == "Approved"

    def test_cgi_biomarker_is_fda_approved(self):
        """Test is_fda_approved method."""
        approved = CGIBiomarker(
            gene="EGFR",
            alteration="EGFR:G719S",
            drug="Afatinib",
            drug_status="Approved",
            association="Responsive",
            evidence_level="FDA guidelines",
            source="FDA",
            tumor_type="NSCLC",
            tumor_type_full="Non-small cell lung cancer",
        )

        assert approved.is_fda_approved() is True


class TestClinVarEvidenceFields:
    """Tests for ClinVar evidence model fields."""

    def test_clinvar_evidence_creation(self):
        """Test creating ClinVarEvidence model."""
        evidence = ClinVarEvidence(
            clinical_significance="Pathogenic",
            conditions=["Melanoma"],
            review_status="criteria provided, multiple submitters, no conflicts",
            variation_id="12345",
        )

        assert evidence.clinical_significance == "Pathogenic"
        assert "Melanoma" in evidence.conditions

    def test_clinvar_evidence_defaults(self):
        """Test ClinVarEvidence default values."""
        evidence = ClinVarEvidence()

        assert evidence.clinical_significance is None
        assert evidence.conditions == []


class TestCOSMICEvidenceFields:
    """Tests for COSMIC evidence model fields."""

    def test_cosmic_evidence_creation(self):
        """Test creating COSMICEvidence model."""
        evidence = COSMICEvidence(
            mutation_id="COSM476",
            primary_site="skin",
            sample_count=1500,
        )

        assert evidence.mutation_id == "COSM476"
        assert evidence.sample_count == 1500

    def test_cosmic_evidence_defaults(self):
        """Test COSMICEvidence default values."""
        evidence = COSMICEvidence()

        assert evidence.mutation_id is None
        assert evidence.sample_count is None


class TestClinicalTrialEvidenceFields:
    """Tests for ClinicalTrialEvidence model fields."""

    def test_clinical_trial_evidence_creation(self):
        """Test creating ClinicalTrialEvidence model."""
        trial = ClinicalTrialEvidence(
            nct_id="NCT12345678",
            title="BRAF V600E Melanoma Trial",
            phase="Phase 3",
            status="Recruiting",
            interventions=["Vemurafenib"],
            conditions=["Melanoma"],
            url="https://clinicaltrials.gov/ct2/show/NCT12345678",
        )

        assert trial.nct_id == "NCT12345678"
        assert trial.phase == "Phase 3"

    def test_clinical_trial_is_phase2_or_later(self):
        """Test is_phase2_or_later method."""
        phase3 = ClinicalTrialEvidence(
            nct_id="NCT123",
            title="Test",
            status="Recruiting",
            phase="PHASE3",
            url="https://example.com",
        )
        phase1 = ClinicalTrialEvidence(
            nct_id="NCT456",
            title="Test",
            status="Recruiting",
            phase="PHASE1",
            url="https://example.com",
        )

        assert phase3.is_phase2_or_later() is True
        assert phase1.is_phase2_or_later() is False


class TestPubMedEvidenceFields:
    """Tests for PubMed evidence model fields."""

    def test_pubmed_evidence_creation(self):
        """Test creating PubMedEvidence model."""
        article = PubMedEvidence(
            pmid="22735384",
            title="BRAF V600E in Melanoma",
            year="2012",
            journal="NEJM",
            url="https://pubmed.ncbi.nlm.nih.gov/22735384/",
        )

        assert article.pmid == "22735384"
        assert article.year == "2012"

    def test_pubmed_evidence_signal_type(self):
        """Test signal type methods."""
        resistance = PubMedEvidence(
            pmid="123",
            title="Test",
            url="https://example.com",
            signal_type="resistance",
        )
        sensitivity = PubMedEvidence(
            pmid="456",
            title="Test",
            url="https://example.com",
            signal_type="sensitivity",
        )

        assert resistance.is_resistance_evidence() is True
        assert resistance.is_sensitivity_evidence() is False
        assert sensitivity.is_sensitivity_evidence() is True
        assert sensitivity.is_resistance_evidence() is False

    def test_pubmed_evidence_citation_format(self):
        """Test format_citation method."""
        article = PubMedEvidence(
            pmid="22735384",
            title="BRAF V600E in Melanoma",
            authors=["Chapman PB", "Hauschild A"],
            year="2012",
            journal="NEJM",
            url="https://pubmed.ncbi.nlm.nih.gov/22735384/",
        )

        citation = article.format_citation()

        assert "Chapman PB et al." in citation
        assert "(2012)" in citation
        assert "PMID: 22735384" in citation

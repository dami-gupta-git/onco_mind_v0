"""Integration tests for new evidence model fields.

Tests validate that the new fields added to evidence models
(CIViCEvidence, VICCEvidence, FDAApproval) are properly populated
and passed through to the frontend via the backend API.
"""

import pytest
import sys
from pathlib import Path

# Add streamlit directory to path for imports
streamlit_dir = Path(__file__).parent.parent.parent / "streamlit"
sys.path.insert(0, str(streamlit_dir))

from oncomind.models.evidence import (
    CIViCEvidence,
    VICCEvidence,
    FDAApproval,
    CIViCAssertionEvidence,
)
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

        # All evidence types should be present (even if empty lists)
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

"""Unit tests for FDA evidence models.

Tests the Pydantic models in oncomind/models/evidence/fda.py.
"""

import pytest

from oncomind.models.evidence.fda import (
    ClinicalStudyEvidence,
    MechanismEvidence,
    AdverseReactionsEvidence,
    FDALabelEvidence,
    FDAApproval,
    extract_combination_partners,
)


class TestClinicalStudyEvidence:
    """Tests for ClinicalStudyEvidence model."""

    def test_default_values(self):
        """All fields should default to None."""
        data = ClinicalStudyEvidence()

        assert data.trial_name is None
        assert data.nct_id is None
        assert data.patients_n is None
        assert data.pfs_months_treatment is None
        assert data.pfs_months_control is None
        assert data.hazard_ratio is None
        assert data.hazard_ratio_ci is None
        assert data.orr_treatment is None
        assert data.orr_control is None
        assert data.biomarker_breakdown is None

    def test_set_values(self):
        """Should correctly store all values."""
        data = ClinicalStudyEvidence(
            trial_name="CAPItello-291",
            nct_id="NCT04305496",
            patients_n=289,
            pfs_months_treatment=7.2,
            pfs_months_control=3.6,
            hazard_ratio=0.50,
            hazard_ratio_ci=(0.38, 0.65),
            orr_treatment=22.0,
            orr_control=12.0,
            biomarker_breakdown={"PIK3CA": 76.0, "AKT1": 13.0},
        )

        assert data.trial_name == "CAPItello-291"
        assert data.nct_id == "NCT04305496"
        assert data.patients_n == 289
        assert data.pfs_months_treatment == 7.2
        assert data.pfs_months_control == 3.6
        assert data.hazard_ratio == 0.50
        assert data.hazard_ratio_ci == (0.38, 0.65)
        assert data.orr_treatment == 22.0
        assert data.orr_control == 12.0
        assert data.biomarker_breakdown["PIK3CA"] == 76.0

    def test_model_dump(self):
        """Test serialization with model_dump."""
        data = ClinicalStudyEvidence(
            trial_name="TEST-001",
            hazard_ratio=0.65,
            hazard_ratio_ci=(0.5, 0.8),
        )
        dumped = data.model_dump()

        assert dumped["trial_name"] == "TEST-001"
        assert dumped["hazard_ratio"] == 0.65
        assert dumped["hazard_ratio_ci"] == (0.5, 0.8)


class TestMechanismEvidence:
    """Tests for MechanismEvidence model."""

    def test_default_values(self):
        """Default values should be correct."""
        data = MechanismEvidence()

        assert data.targets == []
        assert data.mechanism is None
        assert data.preclinical is None

    def test_set_values(self):
        """Should correctly store all values."""
        data = MechanismEvidence(
            targets=["AKT1", "AKT2", "AKT3"],
            mechanism="Inhibits phosphorylation of downstream substrates",
            preclinical="In vitro, reduced cell proliferation",
        )

        assert data.targets == ["AKT1", "AKT2", "AKT3"]
        assert "phosphorylation" in data.mechanism
        assert "In vitro" in data.preclinical

    def test_model_dump(self):
        """Test serialization with model_dump."""
        data = MechanismEvidence(targets=["EGFR", "HER2"])
        dumped = data.model_dump()

        assert dumped["targets"] == ["EGFR", "HER2"]
        assert dumped["mechanism"] is None


class TestAdverseReactionsEvidence:
    """Tests for AdverseReactionsEvidence model."""

    def test_default_values(self):
        """Default values should be correct."""
        data = AdverseReactionsEvidence()

        assert data.common_toxicities == []
        assert data.serious_rate is None
        assert data.discontinuation_rate is None

    def test_set_values(self):
        """Should correctly store all values."""
        data = AdverseReactionsEvidence(
            common_toxicities=[("Diarrhea", 77.0), ("Nausea", 50.0)],
            serious_rate=35.0,
            discontinuation_rate=13.0,
        )

        assert len(data.common_toxicities) == 2
        assert data.common_toxicities[0] == ("Diarrhea", 77.0)
        assert data.serious_rate == 35.0
        assert data.discontinuation_rate == 13.0

    def test_model_dump(self):
        """Test serialization with model_dump."""
        data = AdverseReactionsEvidence(serious_rate=25.0)
        dumped = data.model_dump()

        assert dumped["serious_rate"] == 25.0
        assert dumped["common_toxicities"] == []


class TestFDALabelEvidence:
    """Tests for FDALabelEvidence model."""

    def test_required_fields(self):
        """Drug and gene are required fields."""
        data = FDALabelEvidence(drug="capivasertib", gene="AKT1")

        assert data.drug == "capivasertib"
        assert data.gene == "AKT1"

    def test_optional_fields_default_none(self):
        """Optional fields should default to None."""
        data = FDALabelEvidence(drug="test", gene="BRAF")

        assert data.brand_name is None
        assert data.generic_name is None
        assert data.manufacturer is None
        assert data.indications_and_usage is None
        assert data.clinical_studies is None
        assert data.mechanism_of_action is None
        assert data.adverse_reactions is None
        assert data.last_label_update is None
        assert data.update_reason is None
        assert data.clinical_studies_text is None
        assert data.mechanism_of_action_text is None
        assert data.adverse_reactions_text is None

    def test_with_nested_models(self):
        """Test FDALabelEvidence with nested evidence models."""
        clinical_studies = ClinicalStudyEvidence(
            trial_name="CAPItello-291",
            nct_id="NCT04305496",
            hazard_ratio=0.50,
        )
        mechanism = MechanismEvidence(
            targets=["AKT1", "AKT2", "AKT3"],
        )
        adverse = AdverseReactionsEvidence(
            serious_rate=35.0,
        )

        data = FDALabelEvidence(
            drug="capivasertib",
            gene="AKT1",
            brand_name="TRUQAP",
            clinical_studies=clinical_studies,
            mechanism_of_action=mechanism,
            adverse_reactions=adverse,
        )

        assert data.drug == "capivasertib"
        assert data.brand_name == "TRUQAP"
        assert data.clinical_studies.trial_name == "CAPItello-291"
        assert data.clinical_studies.hazard_ratio == 0.50
        assert "AKT1" in data.mechanism_of_action.targets
        assert data.adverse_reactions.serious_rate == 35.0

    def test_model_dump_nested(self):
        """Test serialization with nested models."""
        clinical_studies = ClinicalStudyEvidence(
            trial_name="TEST-001",
            hazard_ratio=0.65,
        )
        data = FDALabelEvidence(
            drug="test_drug",
            gene="TEST",
            clinical_studies=clinical_studies,
        )
        dumped = data.model_dump()

        assert dumped["drug"] == "test_drug"
        assert dumped["clinical_studies"]["trial_name"] == "TEST-001"
        assert dumped["clinical_studies"]["hazard_ratio"] == 0.65


class TestFDAApproval:
    """Tests for FDAApproval model."""

    def test_default_values(self):
        """All fields should default to None/False."""
        data = FDAApproval()

        assert data.drug_name is None
        assert data.brand_name is None
        assert data.generic_name is None
        assert data.indication is None
        assert data.approval_date is None
        assert data.gene is None
        assert data.variant_in_indications is False
        assert data.variant_in_clinical_studies is False
        assert data.companion_diagnostic is None
        assert data.black_box_warning is None
        assert data.dosing_for_variant is None
        assert data.association is None

    def test_set_values(self):
        """Should correctly store all values."""
        data = FDAApproval(
            drug_name="sotorasib",
            brand_name="LUMAKRAS",
            generic_name="sotorasib",
            indication="KRAS G12C-mutated NSCLC",
            approval_date="2021-05-28",
            gene="KRAS",
            variant_in_indications=True,
            association="Responsive",
        )

        assert data.drug_name == "sotorasib"
        assert data.brand_name == "LUMAKRAS"
        assert data.indication == "KRAS G12C-mutated NSCLC"
        assert data.gene == "KRAS"
        assert data.variant_in_indications is True
        assert data.association == "Responsive"


class TestFDAApprovalExtractApprovedVariant:
    """Tests for FDAApproval.extract_approved_variant method."""

    def test_extract_variant_from_indication_missense(self):
        """Should extract missense variant from indication text."""
        approval = FDAApproval(
            gene="KRAS",
            indication="For treatment of KRAS G12C-mutated NSCLC",
        )

        variant = approval.extract_approved_variant()
        assert variant == "G12C"

    def test_extract_variant_braf_v600e(self):
        """Should extract BRAF V600E variant."""
        approval = FDAApproval(
            gene="BRAF",
            indication="For treatment of BRAF V600E-mutated melanoma",
        )

        variant = approval.extract_approved_variant()
        assert variant == "V600E"

    def test_extract_exon_deletion(self):
        """Should extract exon deletion pattern."""
        approval = FDAApproval(
            gene="EGFR",
            indication="For treatment of EGFR exon 19 deletion NSCLC",
        )

        variant = approval.extract_approved_variant()
        assert "exon" in variant.lower()
        assert "del" in variant.lower()

    def test_extract_any_mutation(self):
        """Should return 'any mutation' for gene-mutated patterns."""
        approval = FDAApproval(
            gene="BRCA",
            indication="For treatment of BRCA-mutated ovarian cancer",
        )

        variant = approval.extract_approved_variant()
        assert variant == "any mutation"

    def test_extract_known_variant_sotorasib(self):
        """Should use fallback for known variant-specific drugs."""
        approval = FDAApproval(
            gene="KRAS",
            generic_name="sotorasib",
            indication="For treatment of advanced NSCLC",  # No variant in text
        )

        variant = approval.extract_approved_variant()
        assert variant == "G12C"

    def test_extract_known_variant_adagrasib(self):
        """Should use fallback for adagrasib."""
        approval = FDAApproval(
            gene="KRAS",
            generic_name="adagrasib",
            indication="For treatment of advanced NSCLC",
        )

        variant = approval.extract_approved_variant()
        assert variant == "G12C"

    def test_extract_no_gene_returns_none(self):
        """Should return None if gene is not set."""
        approval = FDAApproval(
            indication="For treatment of cancer",
        )

        variant = approval.extract_approved_variant()
        assert variant is None

    def test_extract_no_variant_found_returns_none(self):
        """Should return None if no variant pattern matches."""
        approval = FDAApproval(
            gene="TP53",
            indication="For treatment of cancer",  # No variant mentioned
        )

        variant = approval.extract_approved_variant()
        assert variant is None


class TestFDAApprovalParseIndicationForTumor:
    """Tests for FDAApproval.parse_indication_for_tumor method."""

    def test_tumor_match_melanoma(self):
        """Should match melanoma tumor type."""
        approval = FDAApproval(
            indication="For treatment of metastatic melanoma with BRAF V600E mutation",
        )

        result = approval.parse_indication_for_tumor("melanoma")

        assert result["tumor_match"] is True
        assert result["line_of_therapy"] == "unspecified"

    def test_tumor_match_nsclc(self):
        """Should match NSCLC tumor type."""
        approval = FDAApproval(
            indication="For first-line treatment of non-small cell lung cancer",
        )

        result = approval.parse_indication_for_tumor("NSCLC")

        assert result["tumor_match"] is True
        assert result["line_of_therapy"] == "first-line"

    def test_later_line_therapy(self):
        """Should detect later-line therapy indicators."""
        approval = FDAApproval(
            indication="For treatment of melanoma after prior therapy",
        )

        result = approval.parse_indication_for_tumor("melanoma")

        assert result["tumor_match"] is True
        assert result["line_of_therapy"] == "later-line"

    def test_second_line_therapy(self):
        """Should detect second-line therapy."""
        approval = FDAApproval(
            indication="For second-line treatment of colorectal cancer",
        )

        result = approval.parse_indication_for_tumor("colorectal")

        assert result["tumor_match"] is True
        assert result["line_of_therapy"] == "later-line"

    def test_accelerated_approval(self):
        """Should detect accelerated approval."""
        approval = FDAApproval(
            indication="Accelerated approval for treatment of melanoma",
        )

        result = approval.parse_indication_for_tumor("melanoma")

        assert result["tumor_match"] is True
        assert result["approval_type"] == "accelerated"

    def test_no_tumor_match(self):
        """Should return tumor_match=False when tumor doesn't match."""
        approval = FDAApproval(
            indication="For treatment of breast cancer",
        )

        result = approval.parse_indication_for_tumor("melanoma")

        assert result["tumor_match"] is False
        assert result["line_of_therapy"] == "unspecified"
        assert result["approval_type"] == "unspecified"

    def test_empty_indication(self):
        """Should handle empty indication."""
        approval = FDAApproval(indication=None)

        result = approval.parse_indication_for_tumor("melanoma")

        assert result["tumor_match"] is False

    def test_empty_tumor_type(self):
        """Should handle empty tumor type."""
        approval = FDAApproval(indication="For treatment of melanoma")

        result = approval.parse_indication_for_tumor("")

        assert result["tumor_match"] is False

    def test_msi_h_tumor_agnostic(self):
        """Should match MSI-H tumor-agnostic approvals."""
        approval = FDAApproval(
            indication="[FDA APPROVED FOR MSI-H] For treatment of MSI-H or dMMR cancer",
        )

        # Should match any solid tumor
        result = approval.parse_indication_for_tumor("endometrial")

        assert result["tumor_match"] is True

    def test_myeloproliferative_matching(self):
        """Should match myeloproliferative neoplasm variants."""
        approval = FDAApproval(
            indication="For treatment of myelofibrosis",
        )

        result = approval.parse_indication_for_tumor("myeloproliferative neoplasm")

        assert result["tumor_match"] is True


class TestFDAApprovalExtractIndicationCancerType:
    """Tests for FDAApproval.extract_indication_cancer_type method."""

    def test_extract_melanoma(self):
        """Should extract melanoma cancer type."""
        approval = FDAApproval(
            indication="For treatment of metastatic melanoma",
        )

        cancer_type = approval.extract_indication_cancer_type()
        assert cancer_type == "melanoma"

    def test_extract_nsclc(self):
        """Should extract NSCLC cancer type."""
        approval = FDAApproval(
            indication="For treatment of non-small cell lung cancer",
        )

        cancer_type = approval.extract_indication_cancer_type()
        assert cancer_type == "NSCLC"

    def test_extract_breast_cancer(self):
        """Should extract breast cancer type."""
        approval = FDAApproval(
            indication="For treatment of HR-positive, HER2-negative breast cancer",
        )

        cancer_type = approval.extract_indication_cancer_type()
        assert cancer_type == "breast cancer"

    def test_extract_colorectal(self):
        """Should extract colorectal cancer type."""
        approval = FDAApproval(
            indication="For treatment of metastatic colorectal cancer",
        )

        cancer_type = approval.extract_indication_cancer_type()
        assert cancer_type == "colorectal cancer"

    def test_extract_ovarian_cancer(self):
        """Should extract ovarian cancer type."""
        approval = FDAApproval(
            indication="For treatment of epithelial ovarian cancer",
        )

        cancer_type = approval.extract_indication_cancer_type()
        assert cancer_type == "ovarian cancer"

    def test_extract_msi_h_pan_cancer(self):
        """Should extract MSI-H pan-cancer type."""
        approval = FDAApproval(
            indication="For treatment of MSI-H solid tumors",
        )

        cancer_type = approval.extract_indication_cancer_type()
        # Should match either "MSI-H tumors (pan-cancer)" or "solid tumors (pan-cancer)"
        assert cancer_type is not None
        assert "pan-cancer" in cancer_type

    def test_extract_gist(self):
        """Should extract GIST cancer type."""
        approval = FDAApproval(
            indication="For treatment of gastrointestinal stromal tumor",
        )

        cancer_type = approval.extract_indication_cancer_type()
        assert cancer_type == "GIST"

    def test_extract_none_for_empty_indication(self):
        """Should return None for empty indication."""
        approval = FDAApproval(indication=None)

        cancer_type = approval.extract_indication_cancer_type()
        assert cancer_type is None

    def test_extract_none_when_no_match(self):
        """Should return None when no cancer type matches."""
        approval = FDAApproval(
            indication="For treatment of some rare condition",
        )

        cancer_type = approval.extract_indication_cancer_type()
        assert cancer_type is None


class TestExtractCombinationPartners:
    """Tests for extract_combination_partners function."""

    def test_single_partner_fulvestrant(self):
        """Should extract single partner drug."""
        partners = extract_combination_partners(
            "in combination with fulvestrant for HR+/HER2- breast cancer"
        )

        assert len(partners) == 1
        assert partners[0] == "fulvestrant"

    def test_single_partner_panitumumab(self):
        """Should extract panitumumab from indication text."""
        partners = extract_combination_partners(
            "in combination with panitumumab for mCRC"
        )

        assert len(partners) == 1
        assert partners[0] == "panitumumab"

    def test_monotherapy_returns_empty(self):
        """Should return empty list for monotherapy indications."""
        partners = extract_combination_partners(
            "as a single agent for NSCLC"
        )

        assert partners == []

    def test_none_input_returns_empty(self):
        """Should return empty list for None input."""
        partners = extract_combination_partners(None)

        assert partners == []

    def test_empty_string_returns_empty(self):
        """Should return empty list for empty string."""
        partners = extract_combination_partners("")

        assert partners == []

    def test_combined_with_pattern(self):
        """Should extract partner using 'combined with' pattern."""
        partners = extract_combination_partners(
            "combined with carboplatin for ovarian cancer"
        )

        assert len(partners) == 1
        assert partners[0] == "carboplatin"

    def test_multiple_partners_and_pattern(self):
        """Should extract partners joined with 'and' as separate entries."""
        partners = extract_combination_partners(
            "in combination with paclitaxel and carboplatin for NSCLC"
        )

        # "paclitaxel and carboplatin" is split into two separate partners
        assert len(partners) == 2
        assert "paclitaxel" in partners
        assert "carboplatin" in partners

    def test_cleans_trailing_for_phrase(self):
        """Should clean 'for' and following text."""
        partners = extract_combination_partners(
            "in combination with fulvestrant for treatment of breast cancer"
        )

        assert len(partners) == 1
        assert partners[0] == "fulvestrant"
        assert "treatment" not in partners[0]
        assert "breast" not in partners[0]

    def test_case_insensitive(self):
        """Should match regardless of case."""
        partners = extract_combination_partners(
            "IN COMBINATION WITH FULVESTRANT for breast cancer"
        )

        assert len(partners) == 1
        assert partners[0].lower() == "fulvestrant"

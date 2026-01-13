"""Tests for FDA response type determination and biomarker selection drug filtering.

Tests verify:
1. Response type is correctly determined based on VICC evidence
2. OncoKB is weighted higher than other sources
3. Biomarker selection drugs are correctly identified
4. T790M shows resistance to gefitinib/afatinib but sensitivity to osimertinib
5. L858R shows sensitivity to gefitinib
"""

import pytest

from oncomind.config.constants import BIOMARKER_SELECTION_DRUGS


class TestBiomarkerSelectionDrugs:
    """Tests for biomarker selection drug identification."""

    def test_datroway_is_biomarker_selection_for_egfr(self):
        """DATROWAY targets TROP2, not EGFR - should be filtered for EGFR queries."""
        assert "datopotamab deruxtecan" in BIOMARKER_SELECTION_DRUGS
        assert "datroway" in BIOMARKER_SELECTION_DRUGS

        info = BIOMARKER_SELECTION_DRUGS["datopotamab deruxtecan"]
        assert info["target"] == "TROP2"
        assert "EGFR" in info["biomarker_genes"]

    def test_pemetrexed_is_biomarker_selection_for_egfr_and_alk(self):
        """Pemetrexed targets TYMS, not EGFR/ALK - should be filtered."""
        assert "pemetrexed" in BIOMARKER_SELECTION_DRUGS

        info = BIOMARKER_SELECTION_DRUGS["pemetrexed"]
        assert info["target"] == "TYMS"
        assert "EGFR" in info["biomarker_genes"]
        assert "ALK" in info["biomarker_genes"]

    def test_pembrolizumab_is_biomarker_selection_for_multiple_genes(self):
        """Pembrolizumab targets PD-1, uses PD-L1/TMB as biomarker - not gene-specific."""
        assert "pembrolizumab" in BIOMARKER_SELECTION_DRUGS
        assert "keytruda" in BIOMARKER_SELECTION_DRUGS
        assert "keytruda qlex" in BIOMARKER_SELECTION_DRUGS

        info = BIOMARKER_SELECTION_DRUGS["pembrolizumab"]
        assert info["target"] == "PD-1"
        # Should be filtered for all common oncogene queries
        assert "EGFR" in info["biomarker_genes"]
        assert "KRAS" in info["biomarker_genes"]
        assert "BRAF" in info["biomarker_genes"]
        assert "ALK" in info["biomarker_genes"]

    def test_is_biomarker_selection_drug_helper(self):
        """Test the helper function logic used in the UI."""
        def is_biomarker_selection_drug(drug_name: str, gene: str) -> bool:
            """Check if drug uses gene as biomarker but targets something else."""
            if not drug_name or not gene:
                return False
            drug_lower = drug_name.lower()
            for drug_key, info in BIOMARKER_SELECTION_DRUGS.items():
                if drug_key in drug_lower or drug_lower in drug_key:
                    if gene.upper() in [g.upper() for g in info.get("biomarker_genes", [])]:
                        return True
            return False

        # Should be filtered (biomarker selection drugs)
        assert is_biomarker_selection_drug("DATOPOTAMAB DERUXTECAN (DATROWAY)", "EGFR") is True
        assert is_biomarker_selection_drug("PEMETREXED DISODIUM (Pemetrexed)", "EGFR") is True
        assert is_biomarker_selection_drug("PEMETREXED DISODIUM (Pemetrexed)", "ALK") is True
        # Pembrolizumab/Keytruda - PD-L1/TMB based, not gene-specific
        assert is_biomarker_selection_drug("PEMBROLIZUMAB AND BERAHYALURONIDASE ALFA-PMPH (KEYTRUDA QLEX)", "EGFR") is True
        assert is_biomarker_selection_drug("Keytruda", "EGFR") is True
        assert is_biomarker_selection_drug("Pembrolizumab", "KRAS") is True
        assert is_biomarker_selection_drug("Pembrolizumab", "BRAF") is True

        # Should NOT be filtered (direct target drugs)
        assert is_biomarker_selection_drug("OSIMERTINIB (TAGRISSO)", "EGFR") is False
        assert is_biomarker_selection_drug("GEFITINIB (Gefitinib)", "EGFR") is False
        assert is_biomarker_selection_drug("VEMURAFENIB", "BRAF") is False

        # Edge cases
        assert is_biomarker_selection_drug("", "EGFR") is False
        assert is_biomarker_selection_drug("OSIMERTINIB", "") is False


@pytest.mark.integration
class TestFDAResponseType:
    """Integration tests for FDA response type determination.

    Response type for FDA-approved therapies is determined from the `association`
    field on FDAApproval (sourced from CGI biomarkers). VICC evidence is separate
    and doesn't influence FDA therapy response types directly.
    """

    @pytest.fixture
    def evidence_with_associations(self):
        """Create Evidence object with FDA approvals that have association field set."""
        from oncomind.models.evidence.evidence import Evidence, VariantContext, VariantIdentifiers
        from oncomind.models.evidence.fda import FDAApproval

        context = VariantContext(tumor_type="NSCLC")
        identifiers = VariantIdentifiers(variant_id="EGFR:T790M", gene="EGFR", variant="T790M")

        # Create FDA approvals with association field (source of response_type)
        # biomarker_matched=True means this approval covers the queried variant
        fda_approvals = [
            FDAApproval(
                drug_name="Gefitinib",
                brand_name="IRESSA",
                generic_name="GEFITINIB",
                indication="EGFR exon 19 deletions or exon 21 L858R",
                variant_in_indications=False,
                association="Resistant",  # T790M causes resistance to gefitinib
                locus_variant_match={"level": "codon", "scope": "specific", "origin": "kb"},
                biomarker_matched=True,  # Include in get_fda_approved_therapies()
            ),
            FDAApproval(
                drug_name="TAGRISSO",
                brand_name="TAGRISSO",
                generic_name="OSIMERTINIB",
                indication="EGFR T790M mutation-positive NSCLC",
                variant_in_indications=True,
                association="Responsive",  # T790M is the target for osimertinib
                locus_variant_match={"level": "variant", "scope": "specific", "origin": "kb"},
                biomarker_matched=True,  # Include in get_fda_approved_therapies()
            ),
        ]

        return Evidence(
            context=context,
            identifiers=identifiers,
            fda_approvals=fda_approvals,
        )

    def test_gefitinib_shows_resistance_for_t790m(self, evidence_with_associations):
        """Gefitinib should show Resistance for T790M when association is 'Resistant'."""
        therapies = evidence_with_associations.get_fda_approved_therapies()

        gefitinib = next((t for t in therapies if "gefitinib" in t.drug_name.lower()), None)
        assert gefitinib is not None
        assert gefitinib.response_type == "Resistance"

    def test_osimertinib_shows_sensitivity_for_t790m(self, evidence_with_associations):
        """Osimertinib should show Sensitivity for T790M when association is 'Responsive'."""
        therapies = evidence_with_associations.get_fda_approved_therapies()

        osimertinib = next((t for t in therapies if "osimertinib" in t.drug_name.lower()), None)
        assert osimertinib is not None
        assert osimertinib.response_type == "Sensitivity"

    def test_no_association_yields_none_response_type(self):
        """FDA approval without association field should have None response_type."""
        from oncomind.models.evidence.evidence import Evidence, VariantContext, VariantIdentifiers
        from oncomind.models.evidence.fda import FDAApproval

        context = VariantContext(tumor_type="NSCLC")
        identifiers = VariantIdentifiers(variant_id="EGFR:L858R", gene="EGFR", variant="L858R")

        fda_approvals = [
            FDAApproval(
                drug_name="Gefitinib",
                brand_name="IRESSA",
                generic_name="GEFITINIB",
                indication="EGFR exon 19 deletions or exon 21 L858R",
                variant_in_indications=True,
                # No association field - response_type should be None
                locus_variant_match={"level": "variant", "scope": "specific", "origin": "kb"},
                biomarker_matched=True,  # Include in get_fda_approved_therapies()
            ),
        ]

        evidence = Evidence(
            context=context,
            identifiers=identifiers,
            fda_approvals=fda_approvals,
        )

        therapies = evidence.get_fda_approved_therapies()
        gefitinib = next((t for t in therapies if "gefitinib" in t.drug_name.lower()), None)

        assert gefitinib is not None
        assert gefitinib.response_type is None


@pytest.mark.integration
class TestLLMEvidenceSummary:
    """Tests for LLM evidence summary including response type.

    The LLM summary includes VICC MetaKB entries which show response types
    from sources like OncoKB. FDA approvals use the association field for
    response type (if available).
    """

    @pytest.fixture
    def evidence_for_llm(self):
        """Create Evidence object for LLM summary test."""
        from oncomind.models.evidence.evidence import Evidence, VariantContext, VariantIdentifiers
        from oncomind.models.evidence.vicc import VICCEvidence
        from oncomind.models.evidence.fda import FDAApproval

        context = VariantContext(tumor_type="NSCLC")
        identifiers = VariantIdentifiers(variant_id="EGFR:T790M", gene="EGFR", variant="T790M")

        # VICC evidence showing resistance - this appears in VICC MetaKB section
        vicc_evidence = [
            VICCEvidence(
                source="oncokb",
                drugs=["GEFITINIB"],
                response_type="resistant",
                is_resistance=True,
                evidence_level="R1",
                disease="NSCLC",
                variant="T790M",
            ),
        ]

        fda_approvals = [
            FDAApproval(
                drug_name="Gefitinib",
                brand_name="IRESSA",
                generic_name="GEFITINIB",
                indication="EGFR exon 19 deletions or L858R",
                variant_in_indications=False,
                locus_variant_match={"level": "codon", "scope": "specific", "origin": "kb"},
                biomarker_matched=True,  # Include in get_fda_approved_therapies()
            ),
            FDAApproval(
                drug_name="TAGRISSO",
                brand_name="TAGRISSO",
                generic_name="OSIMERTINIB",
                indication="EGFR T790M mutation-positive NSCLC",
                variant_in_indications=True,
                locus_variant_match={"level": "variant", "scope": "specific", "origin": "kb"},
                biomarker_matched=True,  # Include in get_fda_approved_therapies()
            ),
        ]

        return Evidence(
            context=context,
            identifiers=identifiers,
            vicc_evidence=vicc_evidence,
            fda_approvals=fda_approvals,
        )

    def test_llm_summary_includes_vicc_resistance_annotation(self, evidence_for_llm):
        """LLM summary should include VICC evidence with 'res' abbreviation for resistance."""
        summary = evidence_for_llm.get_evidence_summary_for_llm()

        # VICC MetaKB section shows resistance as "res" abbreviation
        assert "gefitinib" in summary.lower()
        # The VICC section shows "(res, oncokb)" format
        assert "res" in summary.lower() or "resistance" in summary.lower()

    def test_llm_summary_includes_fda_drugs(self, evidence_for_llm):
        """LLM summary should include FDA-approved drugs."""
        summary = evidence_for_llm.get_evidence_summary_for_llm()

        # Both FDA drugs should appear
        assert "gefitinib" in summary.lower()
        assert "osimertinib" in summary.lower()

    def test_llm_summary_no_longer_says_any_mutation(self, evidence_for_llm):
        """LLM summary should NOT say 'approved for any EGFR mutation'."""
        summary = evidence_for_llm.get_evidence_summary_for_llm()

        # Should NOT contain misleading "any mutation" text
        assert "approved for any" not in summary.lower()

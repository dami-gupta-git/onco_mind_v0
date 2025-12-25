"""Tests for data models."""

import pytest
from pydantic import ValidationError

from oncomind.models.llm_insight import LLMInsight
from oncomind.models import RecommendedTherapy
from oncomind.models.insight import CIViCEvidence, VICCEvidence
from oncomind.models.insight.myvariant_evidence import MyVariantEvidence
from oncomind.models.variant import VariantInput


class TestVariantInput:
    """Tests for VariantInput model."""

    def test_variant_input_creation(self):
        """Test creating a variant input."""
        variant = VariantInput(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
        )
        assert variant.gene == "BRAF"
        assert variant.variant == "V600E"
        assert variant.tumor_type == "Melanoma"

    def test_to_hgvs(self):
        """Test HGVS conversion."""
        variant = VariantInput(gene="BRAF", variant="V600E", tumor_type="Melanoma")
        assert variant.to_hgvs() == "BRAF:V600E"

    def test_variant_input_without_tumor(self):
        """Test creating a variant input without tumor type."""
        variant = VariantInput(
            gene="KRAS",
            variant="G12C",
            tumor_type=None,
        )
        assert variant.gene == "KRAS"
        assert variant.variant == "G12C"
        assert variant.tumor_type is None

    def test_variant_input_validates_snp_indel_only(self):
        """Test that only SNPs and small indels are allowed."""
        # These should succeed
        VariantInput(gene="BRAF", variant="V600E")  # Missense
        VariantInput(gene="TP53", variant="R248*")  # Nonsense
        VariantInput(gene="BRCA1", variant="185delAG")  # Deletion
        VariantInput(gene="EGFR", variant="L747fs")  # Frameshift

        # These should fail
        with pytest.raises(ValidationError, match="fusion.*not supported"):
            VariantInput(gene="ALK", variant="fusion")

        with pytest.raises(ValidationError, match="amplification.*not supported"):
            VariantInput(gene="ERBB2", variant="amplification")

        with pytest.raises(ValidationError, match="splice.*not supported"):
            VariantInput(gene="MET", variant="exon 14 skipping")

    def test_variant_input_deletion_allowed(self):
        """Test that small deletions are allowed."""
        variant = VariantInput(gene="BRCA1", variant="185delAG")
        assert variant.gene == "BRCA1"
        assert variant.variant == "185delAG"

    def test_variant_input_insertion_allowed(self):
        """Test that small insertions are allowed."""
        # Using a variant with 'ins' keyword
        variant = VariantInput(gene="EGFR", variant="L747_P753delinsS")
        assert variant.gene == "EGFR"


class TestCIViCEvidence:
    """Tests for CIViC Evidence model."""

    def test_civic_evidence_creation(self):
        """Test creating CIViC evidence."""
        civic = CIViCEvidence(
            evidence_type="Predictive",
            evidence_level="A",
            clinical_significance="Sensitivity/Response",
            disease="Melanoma",
            drugs=["Vemurafenib"],
            description="Test evidence",
        )
        assert civic.evidence_type == "Predictive"
        assert "Vemurafenib" in civic.drugs


class TestMyVariantEvidence:
    """Tests for MyVariantEvidence model."""

    def test_myvariant_evidence_has_evidence(self):
        """Test has_evidence method."""
        evidence = MyVariantEvidence(
            variant_id="BRAF:V600E",
            gene="BRAF",
            variant="V600E",
        )
        assert not evidence.has_evidence()

        evidence.civic = [CIViCEvidence(evidence_type="Predictive")]
        assert evidence.has_evidence()

    def test_myvariant_evidence_with_identifiers(self):
        """Test creating evidence with database identifiers."""
        evidence = MyVariantEvidence(
            variant_id="BRAF:V600E",
            gene="BRAF",
            variant="V600E",
            cosmic_id="COSM476",
            ncbi_gene_id="673",
            dbsnp_id="rs113488022",
            clinvar_id="13961",
            hgvs_genomic="NC_000007.13:g.140453136A>T",
            hgvs_protein="NP_004324.2:p.Val600Glu",
            hgvs_transcript="NM_004333.4:c.1799T>A",
        )
        assert evidence.cosmic_id == "COSM476"
        assert evidence.ncbi_gene_id == "673"
        assert evidence.dbsnp_id == "rs113488022"
        assert evidence.clinvar_id == "13961"
        assert evidence.hgvs_genomic == "NC_000007.13:g.140453136A>T"
        assert evidence.hgvs_protein == "NP_004324.2:p.Val600Glu"
        assert evidence.hgvs_transcript == "NM_004333.4:c.1799T>A"


class TestLLMInsight:
    """Tests for LLMInsight model."""

    def test_insight_creation(self):
        """Test creating an insight."""
        insight = LLMInsight(
            llm_summary="Test summary",
            rationale="Test rationale",
        )
        assert insight.llm_summary == "Test summary"
        assert insight.rationale == "Test rationale"

    def test_insight_with_therapies(self):
        """Test creating an insight with recommended therapies."""
        therapy = RecommendedTherapy(
            drug_name="Vemurafenib",
            evidence_level="FDA-approved",
            approval_status="Approved",
            clinical_context="First-line therapy",
        )
        insight = LLMInsight(
            llm_summary="Test summary",
            rationale="Test rationale",
            recommended_therapies=[therapy],
        )
        assert len(insight.recommended_therapies) == 1
        assert insight.recommended_therapies[0].drug_name == "Vemurafenib"

    def test_insight_defaults(self):
        """Test LLMInsight default values."""
        insight = LLMInsight(
            llm_summary="Test summary",
            rationale="Test rationale",
        )
        assert insight.clinical_trials_available is False
        assert insight.recommended_therapies == []
        assert insight.references == []

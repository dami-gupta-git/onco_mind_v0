"""Tests for cBioPortal API client."""

import pytest
from unittest.mock import AsyncMock, patch, MagicMock

from oncomind.api.cbioportal import CBioPortalClient, CoMutationData, CBioPortalError
from oncomind.models.evidence.cbioportal import CBioPortalEvidence, CoMutationEntry


class TestCoMutationData:
    """Tests for CoMutationData dataclass."""

    def test_creation(self):
        """Test creating CoMutationData."""
        data = CoMutationData(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            study_id="skcm_tcga_pan_can_atlas_2018",
            total_samples=449,
            samples_with_gene_mutation=230,
            samples_with_exact_variant=217,
            gene_prevalence_pct=51.2,
            variant_prevalence_pct=48.3,
            co_occurring=[{"gene": "CDKN2A", "count": 98, "pct": 45.2, "odds_ratio": 2.34}],
            mutually_exclusive=[{"gene": "NRAS", "count": 2, "pct": 0.9, "odds_ratio": 0.08}],
            study_name="Skin Cutaneous Melanoma (TCGA, PanCancer Atlas)",
        )

        assert data.gene == "BRAF"
        assert data.variant == "V600E"
        assert data.total_samples == 449
        assert data.gene_prevalence_pct == 51.2
        assert len(data.co_occurring) == 1
        assert len(data.mutually_exclusive) == 1

    def test_to_dict(self):
        """Test converting CoMutationData to dictionary."""
        data = CoMutationData(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            study_id="skcm_tcga",
            total_samples=449,
            samples_with_gene_mutation=230,
            samples_with_exact_variant=217,
            gene_prevalence_pct=51.2,
            variant_prevalence_pct=48.3,
            co_occurring=[],
            mutually_exclusive=[],
            study_name="TCGA Melanoma",
        )

        result = data.to_dict()

        assert result["gene"] == "BRAF"
        assert result["study_name"] == "TCGA Melanoma"
        assert result["total_samples"] == 449


class TestCBioPortalEvidence:
    """Tests for CBioPortalEvidence model."""

    def test_has_data_true(self):
        """Test has_data returns True when data exists."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            total_samples=449,
        )
        assert evidence.has_data() is True

    def test_has_data_false(self):
        """Test has_data returns False when no data."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            total_samples=0,
        )
        assert evidence.has_data() is False

    def test_get_study_url(self):
        """Test generating study URL."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            study_id="skcm_tcga_pan_can_atlas_2018",
        )

        url = evidence.get_study_url()

        assert url == "https://www.cbioportal.org/study/summary?id=skcm_tcga_pan_can_atlas_2018"

    def test_get_study_url_none(self):
        """Test study URL is None when no study_id."""
        evidence = CBioPortalEvidence(gene="BRAF")

        assert evidence.get_study_url() is None

    def test_to_dict(self):
        """Test converting evidence to dictionary."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            study_id="skcm_tcga",
            study_name="TCGA Melanoma",
            total_samples=449,
            samples_with_gene_mutation=230,
            samples_with_exact_variant=217,
            gene_prevalence_pct=51.2,
            variant_prevalence_pct=48.3,
            co_occurring=[CoMutationEntry(gene="TP53", count=28, pct=12.9, odds_ratio=1.56)],
            mutually_exclusive=[],
        )

        result = evidence.to_dict()

        assert result["gene"] == "BRAF"
        assert result["study_name"] == "TCGA Melanoma"
        assert len(result["co_occurring"]) == 1
        assert result["co_occurring"][0]["gene"] == "TP53"

    def test_to_prompt_context_with_data(self):
        """Test generating LLM prompt context."""
        evidence = CBioPortalEvidence(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            study_id="skcm_tcga",
            study_name="TCGA Melanoma",
            total_samples=449,
            samples_with_gene_mutation=230,
            samples_with_exact_variant=217,
            gene_prevalence_pct=51.2,
            variant_prevalence_pct=48.3,
            co_occurring=[CoMutationEntry(gene="CDKN2A", count=98, pct=45.2, odds_ratio=2.34)],
            mutually_exclusive=[CoMutationEntry(gene="NRAS", count=2, pct=0.9, odds_ratio=0.08)],
        )

        context = evidence.to_prompt_context()

        assert "BRAF" in context
        assert "51.2%" in context
        assert "CDKN2A" in context
        assert "NRAS" in context
        assert "cBioPortal" in context

    def test_to_prompt_context_no_data(self):
        """Test prompt context when no data available."""
        evidence = CBioPortalEvidence(gene="BRAF", total_samples=0)

        context = evidence.to_prompt_context()

        assert "No cBioPortal data available" in context


class TestCoMutationEntry:
    """Tests for CoMutationEntry model."""

    def test_creation(self):
        """Test creating a CoMutationEntry."""
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

    def test_defaults(self):
        """Test default values."""
        entry = CoMutationEntry(gene="TP53")

        assert entry.count == 0
        assert entry.pct == 0.0
        assert entry.odds_ratio is None


class TestCBioPortalClient:
    """Tests for CBioPortalClient class."""

    def test_get_study_ids_melanoma(self):
        """Test study ID lookup for melanoma."""
        client = CBioPortalClient()

        study_ids = client._get_study_ids("melanoma")

        assert "skcm_tcga_pan_can_atlas_2018" in study_ids

    def test_get_study_ids_nsclc(self):
        """Test study ID lookup for NSCLC."""
        client = CBioPortalClient()

        study_ids = client._get_study_ids("nsclc")

        # Should contain NSCLC-specific studies
        assert any("nsclc" in s or "luad" in s for s in study_ids)

    def test_get_study_ids_colorectal(self):
        """Test study ID lookup for colorectal."""
        client = CBioPortalClient()

        study_ids = client._get_study_ids("colorectal")

        assert "coadread_tcga_pan_can_atlas_2018" in study_ids

    def test_get_study_ids_default(self):
        """Test fallback to default study for unknown tumor type."""
        client = CBioPortalClient()

        study_ids = client._get_study_ids("unknown_tumor_type")

        assert "msk_impact_2017" in study_ids

    def test_get_study_ids_none(self):
        """Test fallback to default study when no tumor type provided."""
        client = CBioPortalClient()

        study_ids = client._get_study_ids(None)

        assert "msk_impact_2017" in study_ids

    def test_get_study_ids_case_insensitive(self):
        """Test tumor type lookup is case insensitive."""
        client = CBioPortalClient()

        study_ids_lower = client._get_study_ids("melanoma")
        study_ids_upper = client._get_study_ids("MELANOMA")
        study_ids_mixed = client._get_study_ids("Melanoma")

        assert study_ids_lower == study_ids_upper == study_ids_mixed

    def test_cancer_genes_co_occurrence_list(self):
        """Test that CANCER_GENES_CO_OCCURRENCE contains expected genes."""
        from oncomind.config.constants import CANCER_GENES_CO_OCCURRENCE

        # Check key tumor suppressors
        assert "TP53" in CANCER_GENES_CO_OCCURRENCE
        assert "PTEN" in CANCER_GENES_CO_OCCURRENCE
        assert "CDKN2A" in CANCER_GENES_CO_OCCURRENCE

        # Check key oncogenes
        assert "KRAS" in CANCER_GENES_CO_OCCURRENCE
        assert "BRAF" in CANCER_GENES_CO_OCCURRENCE
        assert "PIK3CA" in CANCER_GENES_CO_OCCURRENCE

        # Check DDR genes
        assert "ATM" in CANCER_GENES_CO_OCCURRENCE
        assert "BRCA1" in CANCER_GENES_CO_OCCURRENCE
        assert "BRCA2" in CANCER_GENES_CO_OCCURRENCE

    def test_cancer_genes_co_occurrence_count(self):
        """Test CANCER_GENES_CO_OCCURRENCE has reasonable count (~35 genes)."""
        from oncomind.config.constants import CANCER_GENES_CO_OCCURRENCE

        # Should have between 30-40 genes
        assert 30 <= len(CANCER_GENES_CO_OCCURRENCE) <= 40

    @pytest.mark.asyncio
    async def test_get_entrez_id_cached(self):
        """Test that Entrez ID is cached after first lookup."""
        client = CBioPortalClient()
        client._gene_cache["BRAF"] = 673

        result = await client._get_entrez_id("BRAF")

        assert result == 673

    @pytest.mark.asyncio
    async def test_fetch_co_mutation_data_no_study(self):
        """Test fetch returns None when no study found."""
        client = CBioPortalClient()

        with patch.object(client, "_get_study_ids", return_value=[]):
            result = await client.fetch_co_mutation_data("BRAF", "V600E", "Melanoma")

        assert result is None

    @pytest.mark.asyncio
    async def test_fetch_co_mutation_data_no_entrez_id(self):
        """Test fetch returns None when gene not found."""
        client = CBioPortalClient()

        with patch.object(client, "_get_entrez_id", return_value=None):
            result = await client.fetch_co_mutation_data("UNKNOWNGENE", "X123Y", "Melanoma")

        assert result is None

    @pytest.mark.asyncio
    async def test_context_manager(self):
        """Test async context manager protocol."""
        async with CBioPortalClient() as client:
            assert client._client is not None

        # Client should be closed after exiting context
        assert client._client is None

    @pytest.mark.asyncio
    async def test_close(self):
        """Test closing the client."""
        client = CBioPortalClient()
        # Force client creation
        client._get_client()
        assert client._client is not None

        await client.close()

        assert client._client is None



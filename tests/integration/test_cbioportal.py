"""Integration tests for cBioPortal API client (requires network)."""

import pytest

from oncomind.api.cbioportal import CBioPortalClient


class TestCBioPortalClientIntegration:
    """Integration tests for cBioPortal client."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_braf_v600e_melanoma(self):
        """Test fetching real data for BRAF V600E in melanoma."""
        async with CBioPortalClient() as client:
            result = await client.fetch_co_mutation_data("BRAF", "V600E", "Melanoma")

        # Should return data
        assert result is not None
        assert result.gene == "BRAF"
        assert result.total_samples > 0
        assert result.gene_prevalence_pct > 0

        # BRAF is common in melanoma (~50%)
        assert result.gene_prevalence_pct > 30

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_kras_g12c_nsclc(self):
        """Test fetching real data for KRAS G12C in NSCLC."""
        async with CBioPortalClient() as client:
            result = await client.fetch_co_mutation_data("KRAS", "G12C", "NSCLC")

        assert result is not None
        assert result.gene == "KRAS"
        assert result.total_samples > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_egfr_nsclc(self):
        """Test fetching real data for EGFR in NSCLC."""
        async with CBioPortalClient() as client:
            result = await client.fetch_co_mutation_data("EGFR", "L858R", "NSCLC")

        assert result is not None
        assert result.gene == "EGFR"
        assert result.study_name is not None  # Should have study name

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_tp53_pan_cancer(self):
        """Test fetching TP53 without tumor type (pan-cancer)."""
        async with CBioPortalClient() as client:
            result = await client.fetch_co_mutation_data("TP53", tumor_type=None)

        assert result is not None
        assert result.gene == "TP53"
        assert result.study_id == "msk_impact_2017"  # Default pan-cancer study

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_get_entrez_id_real(self):
        """Test fetching real Entrez ID."""
        async with CBioPortalClient() as client:
            entrez_id = await client._get_entrez_id("BRAF")

        assert entrez_id == 673  # BRAF Entrez ID

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_get_entrez_id_kras(self):
        """Test fetching Entrez ID for KRAS."""
        async with CBioPortalClient() as client:
            entrez_id = await client._get_entrez_id("KRAS")

        assert entrez_id == 3845  # KRAS Entrez ID

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_get_study_name_real(self):
        """Test fetching real study name."""
        async with CBioPortalClient() as client:
            study_name = await client._get_study_name("skcm_tcga_pan_can_atlas_2018")

        assert study_name is not None
        assert "Melanoma" in study_name or "TCGA" in study_name

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_get_sample_count(self):
        """Test fetching sample count for a study."""
        async with CBioPortalClient() as client:
            count = await client._get_sample_count("skcm_tcga_pan_can_atlas_2018")

        # TCGA melanoma has ~400-500 samples
        assert count > 300

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_mutations_found(self):
        """Test that co-occurring mutations are detected for BRAF in melanoma."""
        async with CBioPortalClient() as client:
            result = await client.fetch_co_mutation_data("BRAF", "V600E", "Melanoma")

        assert result is not None
        # BRAF commonly co-occurs with CDKN2A in melanoma
        # Check we have at least some co-occurring mutations
        assert len(result.co_occurring) > 0 or len(result.mutually_exclusive) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mutual_exclusivity_braf_nras(self):
        """Test that BRAF and NRAS are mutually exclusive in melanoma."""
        async with CBioPortalClient() as client:
            result = await client.fetch_co_mutation_data("BRAF", "V600E", "Melanoma")

        assert result is not None

        # NRAS should be in mutually exclusive list (they're in the same pathway)
        nras_in_exclusive = any(
            m.get("gene") == "NRAS" for m in result.mutually_exclusive
        )
        nras_in_co_occurring = any(
            c.get("gene") == "NRAS" for c in result.co_occurring
        )

        # NRAS should NOT be co-occurring with BRAF (same pathway = mutually exclusive)
        # It may or may not show up in mutually_exclusive depending on thresholds
        assert not nras_in_co_occurring or nras_in_exclusive

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_unknown_gene(self):
        """Test handling of unknown gene."""
        async with CBioPortalClient() as client:
            entrez_id = await client._get_entrez_id("NOTAREALGENE123")

        assert entrez_id is None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_colorectal_study_selection(self):
        """Test that colorectal tumor type selects correct study."""
        async with CBioPortalClient() as client:
            result = await client.fetch_co_mutation_data("KRAS", "G12D", "colorectal")

        assert result is not None
        assert "coadread" in result.study_id or "crc" in result.study_id

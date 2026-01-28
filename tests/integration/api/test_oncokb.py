"""Integration tests for OncoKB API client."""

import pytest

from oncomind.api.oncokb import fetch_cancer_gene_list, _cancer_gene_cache
import oncomind.api.oncokb as oncokb_module


class TestFetchCancerGeneList:
    """Tests for fetch_cancer_gene_list()."""

    @pytest.fixture(autouse=True)
    def clear_cache(self):
        """Clear the module-level cache before each test."""
        oncokb_module._cancer_gene_cache = None
        yield
        oncokb_module._cancer_gene_cache = None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_non_empty_set(self):
        """Should return a non-empty set of gene symbols."""
        genes = await fetch_cancer_gene_list()
        assert isinstance(genes, set)
        assert len(genes) > 0, "OncoKB cancer gene list should not be empty"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_contains_well_known_cancer_genes(self):
        """Should contain well-known cancer genes."""
        genes = await fetch_cancer_gene_list()
        expected_genes = {"BRAF", "EGFR", "KRAS", "TP53", "PIK3CA", "ALK", "BRCA1", "BRCA2"}
        for gene in expected_genes:
            assert gene in genes, f"Expected {gene} in OncoKB cancer gene list"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_contains_pdgfra_and_kit(self):
        """PDGFRA and KIT should be in the OncoKB cancer gene list."""
        genes = await fetch_cancer_gene_list()
        assert "PDGFRA" in genes, "PDGFRA should be in OncoKB cancer gene list"
        assert "KIT" in genes, "KIT should be in OncoKB cancer gene list"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_genes_are_uppercase(self):
        """All gene symbols should be uppercase."""
        genes = await fetch_cancer_gene_list()
        for gene in genes:
            assert gene == gene.upper(), f"Gene {gene} should be uppercase"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_caches_result(self):
        """Second call should return cached result."""
        genes1 = await fetch_cancer_gene_list()
        genes2 = await fetch_cancer_gene_list()
        assert genes1 is genes2, "Second call should return same cached object"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_reasonable_gene_count(self):
        """Should return a reasonable number of cancer genes (OncoKB has ~700+)."""
        genes = await fetch_cancer_gene_list()
        assert len(genes) > 500, f"Expected >500 cancer genes, got {len(genes)}"

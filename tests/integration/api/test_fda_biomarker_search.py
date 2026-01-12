"""Integration tests for FDA biomarker search API.

These tests make REAL API calls to the OpenFDA endpoint.
Run with: pytest tests/integration/api/test_fda_biomarker_search.py -v -s

Tests verify that search_fda_by_biomarker correctly finds FDA-approved drugs
for various gene/variant combinations, including recent approvals that may
not be in curated databases.
"""

import pytest
from oncomind.api.fda_drugs import (
    OpenFDAClient,
    search_fda_by_biomarker,
    FDALabelInfo,
)


class TestFDABiomarkerSearchLive:
    """Live OpenFDA biomarker search integration tests.

    These tests hit the real OpenFDA API to verify search and parsing logic.
    """

    @pytest.fixture
    def client(self):
        """Create an OpenFDA client for testing."""
        return OpenFDAClient()

    # =========================================================================
    # IDH1/IDH2 Tests - Including Vorasidenib (2024 approval)
    # =========================================================================

    @pytest.mark.asyncio
    async def test_idh1_finds_vorasidenib(self, client):
        """Test that IDH1 search finds Vorasidenib (approved Aug 2024 for glioma).

        This tests detection of recent FDA approvals that may not be in
        curated databases like CGI or VICC yet.
        """
        results = await client.search_by_biomarker("IDH1")

        brand_names = [r.brand_name.upper() for r in results if r.brand_name]
        generic_names = [r.generic_name.upper() for r in results if r.generic_name]

        print(f"\nFound {len(results)} drugs for IDH1:")
        for r in results:
            print(f"  - {r.brand_name} ({r.generic_name})")

        # Vorasidenib (brand: VORANIGO) should be found
        assert "VORANIGO" in brand_names or "VORASIDENIB" in generic_names, (
            f"Expected Vorasidenib/VORANIGO for IDH1, got brands: {brand_names}"
        )

    @pytest.mark.asyncio
    async def test_idh1_finds_ivosidenib(self, client):
        """Test that IDH1 search finds Ivosidenib (TIBSOVO)."""
        results = await client.search_by_biomarker("IDH1")

        brand_names = [r.brand_name.upper() for r in results if r.brand_name]

        assert "TIBSOVO" in brand_names, f"Expected TIBSOVO for IDH1, got: {brand_names}"

    @pytest.mark.asyncio
    async def test_idh1_finds_olutasidenib(self, client):
        """Test that IDH1 search finds Olutasidenib (REZLIDHIA)."""
        results = await client.search_by_biomarker("IDH1")

        brand_names = [r.brand_name.upper() for r in results if r.brand_name]

        assert "REZLIDHIA" in brand_names, f"Expected REZLIDHIA for IDH1, got: {brand_names}"

    @pytest.mark.asyncio
    async def test_idh1_biomarker_specificity_is_gene_level(self, client):
        """Test that IDH1 drugs have gene-level biomarker specificity.

        IDH inhibitors are approved for IDH1-mutated cancers, not specific variants.
        """
        results = await client.search_by_biomarker("IDH1")

        for r in results:
            if r.brand_name and r.brand_name.upper() in ["VORANIGO", "TIBSOVO", "REZLIDHIA"]:
                assert r.biomarker_specificity == "gene", (
                    f"{r.brand_name} should have gene-level specificity, got: {r.biomarker_specificity}"
                )
                assert "IDH1" in r.target_genes, (
                    f"{r.brand_name} should have IDH1 in target_genes, got: {r.target_genes}"
                )

    # =========================================================================
    # BRAF Tests - Variant-specific approvals (V600E, V600K)
    # =========================================================================

    @pytest.mark.asyncio
    async def test_braf_finds_vemurafenib(self, client):
        """Test that BRAF search finds Vemurafenib (ZELBORAF)."""
        results = await client.search_by_biomarker("BRAF")

        brand_names = [r.brand_name.upper() for r in results if r.brand_name]

        print(f"\nFound {len(results)} drugs for BRAF:")
        for r in results[:10]:
            print(f"  - {r.brand_name} ({r.generic_name})")

        assert "ZELBORAF" in brand_names, f"Expected ZELBORAF for BRAF, got: {brand_names}"

    @pytest.mark.asyncio
    async def test_braf_v600e_filters_correctly(self, client):
        """Test that BRAF V600E variant filter works."""
        results = await client.search_by_biomarker("BRAF", variant="V600E")

        print(f"\nFound {len(results)} drugs for BRAF V600E:")
        for r in results:
            print(f"  - {r.brand_name}: specificity={r.biomarker_specificity}, variants={r.specific_variants}")

        # Should find BRAF V600E-specific drugs
        assert len(results) > 0, "Should find at least one BRAF V600E drug"

        # All results should either have V600E in specific_variants or be gene-level
        for r in results:
            has_v600e = any("V600E" in v.upper() for v in r.specific_variants) if r.specific_variants else False
            is_gene_level = r.biomarker_specificity == "gene"
            has_v600_mention = "V600" in (r.indications_and_usage or "").upper()

            assert has_v600e or is_gene_level or has_v600_mention, (
                f"{r.brand_name} should be V600E-specific or gene-level"
            )

    @pytest.mark.asyncio
    async def test_braf_variant_specific_approvals(self, client):
        """Test that BRAF drugs have variant-level specificity detected."""
        results = await client.search_by_biomarker("BRAF")

        # Find vemurafenib/dabrafenib/encorafenib - these are V600-specific
        variant_specific_drugs = []
        for r in results:
            if r.biomarker_specificity == "variant" and r.specific_variants:
                variant_specific_drugs.append((r.brand_name, r.specific_variants))

        print(f"\nVariant-specific BRAF drugs:")
        for name, variants in variant_specific_drugs:
            print(f"  - {name}: {variants}")

        assert len(variant_specific_drugs) > 0, "Should find variant-specific BRAF approvals"

    # =========================================================================
    # EGFR Tests
    # =========================================================================

    @pytest.mark.asyncio
    async def test_egfr_finds_osimertinib(self, client):
        """Test that EGFR search finds Osimertinib (TAGRISSO)."""
        results = await client.search_by_biomarker("EGFR")

        brand_names = [r.brand_name.upper() for r in results if r.brand_name]

        assert "TAGRISSO" in brand_names, f"Expected TAGRISSO for EGFR, got: {brand_names}"

    @pytest.mark.asyncio
    async def test_egfr_exon19_del_filter(self, client):
        """Test that EGFR exon 19 deletion filter works."""
        # Use "exon 19" as variant to test exon-level filtering
        results = await client.search_by_biomarker("EGFR", variant="exon 19")

        print(f"\nFound {len(results)} drugs for EGFR exon 19:")
        for r in results:
            print(f"  - {r.brand_name}")

        # Should find EGFR TKIs
        assert len(results) > 0, "Should find EGFR exon 19 drugs"

    # =========================================================================
    # ALK Tests (Fusion genes)
    # =========================================================================

    @pytest.mark.asyncio
    async def test_alk_finds_crizotinib(self, client):
        """Test that ALK search finds Crizotinib (XALKORI)."""
        results = await client.search_by_biomarker("ALK")

        brand_names = [r.brand_name.upper() for r in results if r.brand_name]

        print(f"\nFound {len(results)} drugs for ALK:")
        for r in results:
            print(f"  - {r.brand_name}")

        assert "XALKORI" in brand_names, f"Expected XALKORI for ALK, got: {brand_names}"

    # =========================================================================
    # Convenience Function Tests
    # =========================================================================

    @pytest.mark.asyncio
    async def test_convenience_function_search_fda_by_biomarker(self):
        """Test the module-level convenience function."""
        results = await search_fda_by_biomarker("IDH1")

        assert len(results) > 0, "Convenience function should return results"
        assert all(isinstance(r, FDALabelInfo) for r in results), "Should return FDALabelInfo objects"

    @pytest.mark.asyncio
    async def test_convenience_function_with_variant(self):
        """Test convenience function with variant parameter."""
        results = await search_fda_by_biomarker("BRAF", variant="V600E")

        assert len(results) > 0, "Should find BRAF V600E drugs"

    # =========================================================================
    # Edge Cases and Error Handling
    # =========================================================================

    @pytest.mark.asyncio
    async def test_unknown_gene_returns_empty(self, client):
        """Test that searching for unknown gene returns empty list."""
        results = await client.search_by_biomarker("NOTAREALGENE12345")

        assert results == [], f"Unknown gene should return empty list, got: {results}"

    @pytest.mark.asyncio
    async def test_results_are_deduplicated(self, client):
        """Test that results don't contain duplicate drugs."""
        results = await client.search_by_biomarker("BRAF")

        seen_drugs = set()
        for r in results:
            key = (r.generic_name or r.brand_name or r.drug_name).lower()
            assert key not in seen_drugs, f"Duplicate drug found: {key}"
            seen_drugs.add(key)

    @pytest.mark.asyncio
    async def test_all_results_have_required_fields(self, client):
        """Test that all results have essential fields populated."""
        results = await client.search_by_biomarker("IDH1")

        for r in results:
            assert r.drug_name, "drug_name should be set"
            assert r.api_success, "api_success should be True"
            # Note: Some FDA API results may have missing brand/generic names
            # in the openfda section, so we only require drug_name to be set
            # (which defaults to "Unknown" if both are missing)

    # =========================================================================
    # Indications Text Tests
    # =========================================================================

    @pytest.mark.asyncio
    async def test_indications_text_captured(self, client):
        """Test that indications text is captured for results."""
        results = await client.search_by_biomarker("IDH1")

        for r in results:
            if r.brand_name and r.brand_name.upper() == "VORANIGO":
                assert r.indications_and_usage, "VORANIGO should have indications text"
                assert "IDH1" in r.indications_and_usage.upper(), "Indications should mention IDH1"
                # Check for glioma indication
                assert "GLIOMA" in r.indications_and_usage.upper() or "ASTROCYTOMA" in r.indications_and_usage.upper(), (
                    "VORANIGO indications should mention glioma or astrocytoma"
                )

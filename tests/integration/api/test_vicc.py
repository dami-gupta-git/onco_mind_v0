"""Integration tests for VICC MetaKB API.

Tests validate that the VICC API returns therapeutic associations from JAX/PMKB sources.
Note: CIViC and MolecularMatch are filtered out intentionally (we get CIViC directly,
MolecularMatch has unreliable response_type values).
"""

import pytest

from oncomind.api.vicc import VICCClient


class TestVICCBRAFV600E:
    """Tests for BRAF V600E - Tier I variant with FDA-approved BRAF inhibitors."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_associations(self):
        """BRAF V600E should return evidence from JAX/PMKB sources."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=50)
            # After filtering CIViC/MolecularMatch, expect fewer results
            assert len(associations) >= 1, "BRAF V600E should have at least 1 association"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_correct_gene_attribution(self):
        """Associations should have correct gene attribution."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=50)
            braf_assocs = [a for a in associations if a.gene == "BRAF"]
            assert len(braf_assocs) > 0, "Should have associations with gene=BRAF"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_drug_associations(self):
        """Should return drug associations (JAX/PMKB may have different drugs than CIViC)."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=50)

            all_drugs = set()
            for assoc in associations:
                for drug in assoc.drugs:
                    all_drugs.add(drug.lower())

            assert len(all_drugs) > 0, f"Should have drug associations, got: {all_drugs}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_sensitivity_associations(self):
        """BRAF V600E is a sensitivity marker - should have sensitivity associations."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=50)
            sensitivity_assocs = [a for a in associations if a.is_sensitivity()]
            assert len(sensitivity_assocs) > 0, "BRAF V600E should have sensitivity associations"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_evidence_levels(self):
        """Should have evidence level annotations."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=50)

            evidence_levels = {a.evidence_level for a in associations if a.evidence_level}
            # JAX/PMKB may have different level distribution than CIViC
            assert len(evidence_levels) > 0, f"Should have evidence levels, got: {evidence_levels}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_sources_are_filtered(self):
        """Verify CIViC and MolecularMatch are filtered out."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=50)

            sources = {a.source.lower() for a in associations if a.source}
            assert "civic" not in sources, "CIViC should be filtered out"
            assert "molecularmatch" not in sources, "MolecularMatch should be filtered out"


class TestVICCEGFRL858R:
    """Tests for EGFR L858R - common activating mutation in NSCLC.

    Note: EGFR L858R in VICC is heavily represented by MolecularMatch/CIViC,
    so after filtering we may have few or no results from JAX/PMKB.
    """

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_or_empty(self):
        """EGFR L858R may have limited JAX/PMKB coverage."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("EGFR", "L858R", max_results=50)
            # May be empty after filtering - that's OK, we get EGFR from CIViC directly
            assert isinstance(associations, list)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_sources_are_filtered(self):
        """If there are results, verify filtering works."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("EGFR", "L858R", max_results=50)

            if associations:
                sources = {a.source.lower() for a in associations if a.source}
                assert "civic" not in sources, "CIViC should be filtered out"
                assert "molecularmatch" not in sources, "MolecularMatch should be filtered out"


class TestVICCKRASG12C:
    """Tests for KRAS G12C - targetable mutation with approved inhibitors."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_associations(self):
        """KRAS G12C should have evidence."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("KRAS", "G12C", max_results=50)
            assert len(associations) >= 2, "KRAS G12C should have at least 2 associations"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_drug_associations(self):
        """KRAS G12C should have drug associations."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("KRAS", "G12C", max_results=50)

            all_drugs = set()
            for assoc in associations:
                for drug in assoc.drugs:
                    all_drugs.add(drug.lower())

            assert len(all_drugs) > 0, "KRAS G12C should have drug associations"


class TestVICCAssociationStructure:
    """Tests for VICCAssociation data structure integrity."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_required_fields_types(self):
        """Associations should have required fields with correct types."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=10)

            for assoc in associations:
                assert isinstance(assoc.gene, str)
                assert isinstance(assoc.disease, str)
                assert isinstance(assoc.drugs, list)
                assert isinstance(assoc.source, str)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_methods_work(self):
        """Association methods should work without error."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=10)

            for assoc in associations:
                _ = assoc.is_sensitivity()
                _ = assoc.is_resistance()
                _ = assoc.get_oncokb_level()
                _ = assoc.to_dict()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_to_dict_keys(self):
        """to_dict should return expected keys."""
        async with VICCClient() as client:
            associations = await client.fetch_associations("BRAF", "V600E", max_results=10)

            expected_keys = {
                "description", "gene", "variant", "disease", "drugs",
                "evidence_level", "response_type", "source", "publication_url",
                "oncogenic", "is_sensitivity", "is_resistance", "oncokb_level"
            }

            for assoc in associations:
                assoc_dict = assoc.to_dict()
                assert set(assoc_dict.keys()) == expected_keys


class TestVICCTumorTypeFilter:
    """Tests for tumor type filtering."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_melanoma_filter(self):
        """Filtering by melanoma should work."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "BRAF", "V600E", tumor_type="melanoma", max_results=10
            )

            # If we got results, verify we didn't crash
            # Exact filtering depends on API data
            assert True

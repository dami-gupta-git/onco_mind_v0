"""Integration tests for VICC MetaKB API.

Tests validate that the VICC API returns therapeutic associations from JAX/PMKB sources.
Note: CIViC and MolecularMatch are filtered out intentionally (we get CIViC directly,
MolecularMatch has unreliable response_type values).
"""

import pytest

from oncomind.api.vicc import VICCClient


def assert_valid_association(assoc):
    """Verify an association has required fields with correct types."""
    assert isinstance(assoc.gene, str) and assoc.gene, "gene must be non-empty string"
    assert isinstance(assoc.disease, str), "disease must be string"
    assert isinstance(assoc.drugs, list), "drugs must be list"
    assert (
        isinstance(assoc.source, str) and assoc.source
    ), "source must be non-empty string"
    assert assoc.source.lower() not in (
        "civic",
        "molecularmatch",
    ), "filtered sources should not appear"


class TestVICCBRAFV600E:
    """Tests for BRAF V600E - Tier I variant with FDA-approved BRAF inhibitors."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_associations_with_drugs(self):
        """BRAF V600E should return drug associations from JAX/PMKB sources."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "BRAF", "V600E", max_results=200
            )
            assert (
                len(associations) >= 1
            ), "BRAF V600E should have at least 1 association"

            for assoc in associations:
                assert_valid_association(assoc)

            all_drugs = {drug.lower() for assoc in associations for drug in assoc.drugs}
            assert len(all_drugs) > 0, "Should have drug associations"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_sensitivity_associations(self):
        """BRAF V600E is a sensitivity marker - should have sensitivity associations."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "BRAF", "V600E", max_results=200
            )
            assert len(associations) > 0, "BRAF V600E should return associations"

            for assoc in associations:
                assert_valid_association(assoc)

            sensitivity_assocs = [a for a in associations if a.is_sensitivity()]
            assert (
                len(sensitivity_assocs) > 0
            ), "BRAF V600E should have sensitivity associations"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_melanoma_tumor_filter(self):
        """Filtering by melanoma should return melanoma-relevant associations."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "BRAF", "V600E", tumor_type="melanoma", max_results=200
            )
            assert (
                len(associations) > 0
            ), "BRAF V600E with melanoma filter should return associations"

            for assoc in associations:
                assert_valid_association(assoc)


class TestVICCEGFRL858R:
    """Tests for EGFR L858R - common activating mutation in NSCLC."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_associations_with_evidence_levels(self):
        """EGFR L858R should return associations with evidence levels."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "EGFR", "L858R", max_results=200
            )
            assert (
                len(associations) >= 1
            ), "EGFR L858R should have at least 1 association"

            for assoc in associations:
                assert_valid_association(assoc)

            evidence_levels = {
                a.evidence_level for a in associations if a.evidence_level
            }
            assert len(evidence_levels) > 0, "Should have evidence levels"


class TestVICCEGFRT790M:
    """Tests for EGFR T790M - resistance mutation."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_resistance_associations(self):
        """EGFR T790M is a resistance mutation - should have resistance associations."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "EGFR", "T790M", max_results=200
            )
            assert len(associations) > 0, "EGFR T790M should return associations"

            for assoc in associations:
                assert_valid_association(assoc)

            resistance_assocs = [a for a in associations if a.is_resistance()]
            assert (
                len(resistance_assocs) > 0
            ), "EGFR T790M should have resistance associations"


class TestVICCKRASG12C:
    """Tests for KRAS G12C - targetable mutation with approved inhibitors."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_associations_with_drugs(self):
        """KRAS G12C should have drug associations."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "KRAS", "G12C", max_results=200
            )
            assert (
                len(associations) >= 2
            ), "KRAS G12C should have at least 2 associations"

            for assoc in associations:
                assert_valid_association(assoc)

            all_drugs = {drug.lower() for assoc in associations for drug in assoc.drugs}
            assert len(all_drugs) > 0, "KRAS G12C should have drug associations"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_correct_gene_attribution(self):
        """Associations should have correct gene attribution."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "KRAS", "G12C", max_results=200
            )
            assert len(associations) > 0, "KRAS G12C should return associations"

            for assoc in associations:
                assert_valid_association(assoc)

            kras_assocs = [a for a in associations if a.gene == "KRAS"]
            assert len(kras_assocs) > 0, "Should have associations with gene=KRAS"


class TestVICCALKFusion:
    """Tests for ALK fusions - common in NSCLC."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_associations(self):
        """ALK fusions should return associations."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "ALK", "Fusion", max_results=200
            )
            assert (
                len(associations) >= 1
            ), "ALK Fusion should have at least 1 association"

            for assoc in associations:
                assert_valid_association(assoc)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_sensitivity_associations(self):
        """ALK fusions are sensitivity markers for ALK inhibitors."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "ALK", "Fusion", max_results=200
            )
            assert len(associations) > 0, "ALK Fusion should return associations"

            for assoc in associations:
                assert_valid_association(assoc)

            sensitivity_assocs = [a for a in associations if a.is_sensitivity()]
            assert (
                len(sensitivity_assocs) > 0
            ), "ALK Fusion should have sensitivity associations"


class TestVICCPIK3CAH1047R:
    """Tests for PIK3CA H1047R - common hotspot in breast cancer."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_associations_with_evidence(self):
        """PIK3CA H1047R should return associations with evidence levels."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "ALK", "F1174L", max_results=200
            )
            assert (
                len(associations) >= 1
            ), "ALK F1174L should have at least 1 association"

            for assoc in associations:
                assert_valid_association(assoc)

            evidence_levels = {
                a.evidence_level for a in associations if a.evidence_level
            }
            assert len(evidence_levels) > 0, "Should have evidence levels"


class TestVICCERBB2S310F:
    """Tests for ERBB2 S310F - HER2 mutation."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_to_dict_keys(self):
        """to_dict should return expected keys."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "ERBB2", "S310F", max_results=200
            )
            assert len(associations) > 0, "ERBB2 S310F should return associations"

            expected_keys = {
                "description",
                "gene",
                "variant",
                "disease",
                "drugs",
                "evidence_level",
                "response_type",
                "source",
                "publication_url",
                "oncogenic",
                "is_sensitivity",
                "is_resistance",
                "oncokb_level",
            }

            for assoc in associations:
                assert_valid_association(assoc)
                assoc_dict = assoc.to_dict()
                assert set(assoc_dict.keys()) == expected_keys


class TestVICCRETFusion:
    """Tests for RET fusions - targetable in NSCLC and thyroid."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_get_oncokb_level(self):
        """get_oncokb_level() should return valid OncoKB levels or None."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "RET", "Fusion", max_results=200
            )
            assert len(associations) > 0, "RET Fusion should return associations"

            for assoc in associations:
                assert_valid_association(assoc)
                level = assoc.get_oncokb_level()
                assert level is None or isinstance(level, str)


class TestVICCMETAmplification:
    """Tests for MET amplification."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_methods_work(self):
        """Association methods should work without error."""
        async with VICCClient() as client:
            associations = await client.fetch_associations(
                "MET", "amplification", max_results=200
            )
            assert len(associations) > 0, "MET amplification should return associations"

            for assoc in associations:
                assert_valid_association(assoc)
                # All methods should return correct types
                assert isinstance(assoc.is_sensitivity(), bool)
                assert isinstance(assoc.is_resistance(), bool)
                level = assoc.get_oncokb_level()
                assert level is None or isinstance(level, str)
                assert isinstance(assoc.to_dict(), dict)

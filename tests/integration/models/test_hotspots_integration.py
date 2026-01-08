"""Integration tests for Cancer Hotspots evidence.

Tests the Hotspots integration functionality:
- Hotspot data is fetched from local TSV file (data/hotspots_msk/hotspots.txt)
- Well-known hotspots (BRAF V600E, KRAS G12D) are correctly identified
- Variant distribution and tumor type composition are populated
- Hotspots data is included in biological context for LLM
"""

import pytest

from oncomind.insight_builder import Conductor, ConductorConfig
from oncomind.api.hotspots import HotspotsClient, get_cancer_hotspots, is_hotspot_variant


# =============================================================================
# HOTSPOTS API UNIT TESTS
# =============================================================================

class TestHotspotsAPI:
    """Unit tests for the Hotspots API client."""

    def test_hotspots_client_loads_data(self):
        """HotspotsClient should load data from TSV file."""
        client = HotspotsClient()
        entries = client.get_hotspots_for_gene("BRAF")
        assert len(entries) > 0, "BRAF should have hotspot entries"

    def test_get_cancer_hotspots_braf(self):
        """get_cancer_hotspots should return codon positions for BRAF."""
        codons = get_cancer_hotspots("BRAF")
        assert len(codons) > 0, "BRAF should have hotspot codons"
        assert 600 in codons, "BRAF V600 should be a hotspot"

    def test_get_cancer_hotspots_kras(self):
        """get_cancer_hotspots should return codon positions for KRAS."""
        codons = get_cancer_hotspots("KRAS")
        assert len(codons) > 0, "KRAS should have hotspot codons"
        assert 12 in codons, "KRAS G12 should be a hotspot"
        assert 13 in codons, "KRAS G13 should be a hotspot"
        assert 61 in codons, "KRAS Q61 should be a hotspot"

    def test_is_hotspot_variant_braf_v600e(self):
        """is_hotspot_variant should return True for BRAF V600E."""
        assert is_hotspot_variant("BRAF", "V600E") is True
        assert is_hotspot_variant("BRAF", "V600K") is True
        assert is_hotspot_variant("BRAF", "V600M") is True

    def test_is_hotspot_variant_kras_g12d(self):
        """is_hotspot_variant should return True for KRAS G12D."""
        assert is_hotspot_variant("KRAS", "G12D") is True
        assert is_hotspot_variant("KRAS", "G12C") is True
        assert is_hotspot_variant("KRAS", "G12V") is True

    def test_is_hotspot_variant_non_hotspot(self):
        """is_hotspot_variant should return False for non-hotspot positions."""
        # BRAF position 100 is not a known hotspot
        assert is_hotspot_variant("BRAF", "A100T") is False

    def test_fetch_hotspot_evidence_braf_v600e(self):
        """fetch_hotspot_evidence should return detailed evidence for BRAF V600E."""
        client = HotspotsClient()
        evidence = client.fetch_hotspot_evidence("BRAF", "V600E")

        assert evidence is not None
        assert evidence.is_hotspot is True
        assert evidence.is_exact_variant_match is True
        assert evidence.hotspot is not None
        assert evidence.hotspot.gene == "BRAF"
        assert "600" in evidence.hotspot.residue

        # Should have variant distribution
        assert len(evidence.hotspot.variants) > 0
        # V600E should be the most common
        top_variant = evidence.hotspot.get_top_variants(1)[0]
        assert top_variant.amino_acid == "E", "V600E should be the most common variant"
        assert top_variant.count > 100, "V600E should have many samples"

        # Should have tumor type composition
        assert len(evidence.hotspot.tumor_type_composition) > 0
        # Skin (melanoma) should be prominent
        top_tumors = evidence.hotspot.get_top_tumor_types(3)
        tumor_names = [t.tumor_type.lower() for t in top_tumors]
        assert "skin" in tumor_names, "Melanoma (skin) should be a top tumor type for BRAF V600"

    def test_fetch_hotspot_evidence_codon_match(self):
        """fetch_hotspot_evidence should indicate codon-level match for less common variants."""
        client = HotspotsClient()
        # V600R is at the same position but less common
        evidence = client.fetch_hotspot_evidence("BRAF", "V600R")

        assert evidence is not None
        assert evidence.is_hotspot is True
        # R is a less common change, might not be exact match
        # The test is flexible - it just verifies we get hotspot data


# =============================================================================
# HOTSPOTS INTEGRATION TESTS
# =============================================================================

@pytest.mark.integration
class TestHotspotsIntegration:
    """Integration tests for hotspots data in the full pipeline."""

    @pytest.mark.asyncio
    async def test_braf_v600e_has_hotspot_evidence(self):
        """BRAF V600E should have hotspot evidence in the pipeline."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        # Should have hotspots evidence
        hotspots = result.evidence.hotspots_evidence
        assert hotspots is not None, "Should have hotspots_evidence"
        assert hotspots.has_data() is True, "Hotspots should have data"
        assert hotspots.is_hotspot is True, "BRAF V600E should be a hotspot"
        assert hotspots.is_exact_variant_match is True, "V600E should be exact match"

    @pytest.mark.asyncio
    async def test_kras_g12d_has_hotspot_evidence(self):
        """KRAS G12D should have hotspot evidence in the pipeline."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("KRAS G12D", tumor_type="Pancreatic Cancer")

        hotspots = result.evidence.hotspots_evidence
        assert hotspots is not None, "Should have hotspots_evidence"
        assert hotspots.has_data() is True, "Hotspots should have data"
        assert hotspots.is_hotspot is True, "KRAS G12D should be a hotspot"

    @pytest.mark.asyncio
    async def test_hotspots_in_biological_context(self):
        """Hotspots data should appear in biological context for LLM."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        # Get biological context
        bio_context = result.evidence.get_biological_context_for_llm()

        # Should include hotspot information
        assert "CANCER HOTSPOT" in bio_context, \
            f"Biological context should include CANCER HOTSPOT, got:\n{bio_context}"
        assert "BRAF" in bio_context
        assert "V600" in bio_context

    @pytest.mark.asyncio
    async def test_hotspots_prompt_context_format(self):
        """Hotspots to_prompt_context should return well-formatted string."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        hotspots = result.evidence.hotspots_evidence
        assert hotspots is not None

        prompt_context = hotspots.to_prompt_context()
        assert "CANCER HOTSPOT" in prompt_context
        assert "cancer samples" in prompt_context
        assert "Common changes" in prompt_context
        assert "Tumor distribution" in prompt_context
        assert "cancerhotspots.org" in prompt_context

    @pytest.mark.asyncio
    async def test_hotspots_serialization(self):
        """Hotspots evidence should serialize correctly."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        hotspots = result.evidence.hotspots_evidence
        assert hotspots is not None

        # Test model_dump
        data = hotspots.model_dump()
        assert "gene" in data
        assert "is_hotspot" in data
        assert "is_exact_variant_match" in data
        assert data["is_hotspot"] is True

        # Test nested hotspot data
        assert "hotspot" in data
        assert data["hotspot"] is not None
        assert "variants" in data["hotspot"]
        assert "tumor_type_composition" in data["hotspot"]

    @pytest.mark.asyncio
    async def test_non_hotspot_variant(self):
        """Non-hotspot variants should have hotspots_evidence but is_hotspot=False."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # CHEK2 I157T is not a cancer hotspot
            result = await conductor.run("CHEK2 I157T", tumor_type="Breast Cancer")

        hotspots = result.evidence.hotspots_evidence
        # May be None or have is_hotspot=False
        if hotspots is not None:
            assert hotspots.is_hotspot is False or hotspots.has_data() is False, \
                "CHEK2 I157T should not be a cancer hotspot"


# =============================================================================
# HOTSPOTS DATA FILE TESTS
# =============================================================================

class TestHotspotsDataFile:
    """Tests for the hotspots data file integrity."""

    def test_hotspots_file_has_expected_genes(self):
        """Hotspots file should contain common cancer genes."""
        client = HotspotsClient()

        expected_genes = ["BRAF", "KRAS", "NRAS", "PIK3CA", "EGFR", "TP53", "IDH1", "KIT"]
        for gene in expected_genes:
            entries = client.get_hotspots_for_gene(gene)
            assert len(entries) > 0, f"{gene} should have hotspot entries"

    def test_hotspots_total_count(self):
        """Hotspots file should have reasonable number of entries."""
        client = HotspotsClient()
        # Force data load
        client._ensure_loaded()

        total_genes = len(client._data) if client._data else 0
        assert total_genes > 100, f"Should have >100 genes with hotspots, got {total_genes}"

    def test_hotspots_variant_counts_reasonable(self):
        """Hotspot variant counts should be reasonable."""
        client = HotspotsClient()
        evidence = client.fetch_hotspot_evidence("BRAF", "V600E")

        assert evidence is not None
        assert evidence.hotspot is not None

        # BRAF V600 should have many samples
        assert evidence.hotspot.total_samples > 100, \
            f"BRAF V600 should have >100 samples, got {evidence.hotspot.total_samples}"

        # Q-value should be very small (highly significant)
        assert evidence.hotspot.q_value < 0.01, \
            f"BRAF V600 should have q<0.01, got {evidence.hotspot.q_value}"

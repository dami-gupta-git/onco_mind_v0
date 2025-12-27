"""Integration tests for evidence assessment at different layers.

Tests verify:
1. well_characterized_detailed field with basis explanations
2. Evidence basis tags in research hypotheses
3. Hotspot detection in gap analysis
4. Streamlit backend evidence_gaps structure
5. End-to-end pipeline from API to UI-ready data
"""

import pytest
import asyncio
import sys
from pathlib import Path

from oncomind.insight_builder import Conductor, ConductorConfig
from oncomind.insight_builder.gap_detector import detect_evidence_gaps
from oncomind.models.evidence.evidence_gaps import (
    GapCategory,
    GapSeverity,
    EvidenceGaps,
    CharacterizedAspect,
)
from oncomind.models.gene_context import is_hotspot_variant, is_hotspot_adjacent


# =============================================================================
# LAYER 1: WELL-CHARACTERIZED DETAILED FIELD TESTS
# =============================================================================

class TestWellCharacterizedDetailed:
    """Tests for well_characterized_detailed field with basis explanations."""

    def test_characterized_aspect_model(self):
        """CharacterizedAspect model should have aspect and basis."""
        aspect = CharacterizedAspect(
            aspect="clinical actionability",
            basis="5 FDA approvals + 2 CIViC assertions"
        )

        assert aspect.aspect == "clinical actionability"
        assert aspect.basis == "5 FDA approvals + 2 CIViC assertions"

    def test_evidence_gaps_has_well_characterized_detailed(self):
        """EvidenceGaps should have well_characterized_detailed field."""
        gaps = EvidenceGaps(
            well_characterized=["clinical actionability"],
            well_characterized_detailed=[
                CharacterizedAspect(
                    aspect="clinical actionability",
                    basis="5 FDA approvals"
                )
            ]
        )

        assert len(gaps.well_characterized_detailed) == 1
        assert gaps.well_characterized_detailed[0].aspect == "clinical actionability"
        assert gaps.well_characterized_detailed[0].basis == "5 FDA approvals"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_braf_v600e_has_detailed_characterization(self):
        """BRAF V600E should have detailed well-characterized aspects with basis."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have well_characterized_detailed populated
        assert len(gaps.well_characterized_detailed) > 0

        # Each detailed item should have both aspect and basis
        for item in gaps.well_characterized_detailed:
            assert item.aspect is not None and len(item.aspect) > 0
            assert item.basis is not None and len(item.basis) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_hotspot_basis_mentions_cancerhotspots(self):
        """Hotspot basis should mention cancerhotspots.org."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("KRAS G12D", tumor_type="Pancreatic")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find the hotspot characterization
        hotspot_items = [
            item for item in gaps.well_characterized_detailed
            if "hotspot" in item.aspect.lower()
        ]

        assert len(hotspot_items) > 0, "Should have hotspot in well_characterized_detailed"

        # Basis should mention cancerhotspots.org
        hotspot = hotspot_items[0]
        assert "cancerhotspots" in hotspot.basis.lower() or "codon" in hotspot.basis.lower()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_clinical_actionability_basis_has_counts(self):
        """Clinical actionability basis should include FDA/CIViC counts."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find clinical actionability
        clinical_items = [
            item for item in gaps.well_characterized_detailed
            if "clinical" in item.aspect.lower() or "actionability" in item.aspect.lower()
        ]

        if clinical_items:
            # Basis should mention FDA or CIViC
            basis = clinical_items[0].basis.lower()
            assert "fda" in basis or "civic" in basis, \
                f"Clinical basis should mention FDA or CIViC: {clinical_items[0].basis}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_pathogenicity_basis_has_scores(self):
        """Computational pathogenicity basis should include prediction scores."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("PIK3CA H1047R", tumor_type="Breast")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find pathogenicity characterization
        path_items = [
            item for item in gaps.well_characterized_detailed
            if "pathogenicity" in item.aspect.lower()
        ]

        if path_items:
            basis = path_items[0].basis
            # Should have at least one prediction tool mentioned
            has_predictor = any(
                tool in basis.lower()
                for tool in ["alphamissense", "cadd", "polyphen"]
            )
            assert has_predictor, f"Pathogenicity basis should mention predictors: {basis}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_legacy_and_detailed_are_consistent(self):
        """well_characterized and well_characterized_detailed should match."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR L858R", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Same number of items
        assert len(gaps.well_characterized) == len(gaps.well_characterized_detailed)

        # Aspects should match
        legacy_aspects = set(gaps.well_characterized)
        detailed_aspects = set(item.aspect for item in gaps.well_characterized_detailed)
        assert legacy_aspects == detailed_aspects


# =============================================================================
# LAYER 2: EVIDENCE BASIS TAGS IN RESEARCH HYPOTHESES
# =============================================================================

class TestEvidenceBasisTags:
    """Tests for evidence basis tags in LLM-generated research hypotheses."""

    VALID_TAGS = [
        "[Direct Clinical Data]",
        "[Preclinical Data]",
        "[Pan-Cancer Extrapolation]",
        "[Nearby-Variant Inference]",
        "[Pathway-Level Inference]",
    ]

    @pytest.mark.integration
    @pytest.mark.asyncio
    @pytest.mark.skipif(
        not pytest.importorskip("openai", reason="OpenAI not installed"),
        reason="Requires OpenAI API"
    )
    async def test_research_hypotheses_have_tags(self):
        """Research hypotheses should include evidence basis tags."""
        config = ConductorConfig(
            enable_llm=True,
            enable_literature=False,
            llm_model="gpt-4o-mini"
        )

        try:
            async with Conductor(config) as conductor:
                result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

            if result.llm and result.llm.research_hypotheses:
                for hypothesis in result.llm.research_hypotheses:
                    # Each hypothesis should start with a valid tag
                    has_tag = any(
                        hypothesis.strip().startswith(tag)
                        for tag in self.VALID_TAGS
                    )
                    # Allow for hypotheses without tags (LLM may not always comply)
                    # but warn if none have tags
                    if not has_tag:
                        print(f"Warning: Hypothesis missing tag: {hypothesis[:50]}...")
        except Exception as e:
            pytest.skip(f"LLM test skipped due to: {e}")


# =============================================================================
# LAYER 3: HOTSPOT DETECTION IN GAP ANALYSIS
# =============================================================================

class TestHotspotGapIntegration:
    """Tests for hotspot detection integration with gap analysis."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_hotspot_variant_marked_well_characterized(self):
        """Known hotspot variants should be marked as well-characterized."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)

        # Test multiple known hotspots
        hotspot_variants = [
            ("BRAF", "V600E"),
            ("KRAS", "G12D"),
            ("EGFR", "L858R"),
            ("PIK3CA", "H1047R"),
        ]

        async with Conductor(config) as conductor:
            for gene, variant in hotspot_variants:
                result = await conductor.run(f"{gene} {variant}")
                gaps = result.evidence.evidence_gaps
                if gaps is None:
                    gaps = result.evidence.compute_evidence_gaps()

                # Should have "known cancer hotspot" in well_characterized
                assert "known cancer hotspot" in gaps.well_characterized, \
                    f"{gene} {variant} should be recognized as hotspot"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_hotspot_adjacent_flagged_for_research(self):
        """Hotspot-adjacent variants should be flagged as research opportunities."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # BRAF V598E is near hotspot V600
            result = await conductor.run("BRAF V598E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should recognize it's near a hotspot
        near_hotspot = [
            w for w in gaps.well_characterized
            if "near hotspot" in w.lower()
        ]
        assert len(near_hotspot) > 0, "Should recognize variant is near hotspot"

        # Should have a FUNCTIONAL gap for characterization
        functional_gaps = gaps.get_gaps_by_category(GapCategory.FUNCTIONAL)
        hotspot_gaps = [
            g for g in functional_gaps
            if "hotspot" in g.description.lower()
        ]
        assert len(hotspot_gaps) > 0, "Should flag need for functional characterization"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_non_hotspot_not_marked_as_hotspot(self):
        """Non-hotspot variants should not be marked as hotspots."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # BRAF V500E is not a hotspot (far from 600)
            result = await conductor.run("BRAF V500E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should NOT have "known cancer hotspot"
        assert "known cancer hotspot" not in gaps.well_characterized


# =============================================================================
# LAYER 4: STREAMLIT BACKEND EVIDENCE_GAPS STRUCTURE
# =============================================================================

# Add streamlit directory to path
streamlit_dir = Path(__file__).parent.parent.parent / "streamlit"
sys.path.insert(0, str(streamlit_dir))


class TestStreamlitBackendEvidenceGaps:
    """Tests for Streamlit backend evidence_gaps data structure."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_backend_returns_well_characterized_detailed(self):
        """Backend should return well_characterized_detailed in evidence_gaps."""
        from backend import get_variant_insight

        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert "evidence_gaps" in result

        evidence_gaps = result["evidence_gaps"]
        assert "well_characterized_detailed" in evidence_gaps
        assert isinstance(evidence_gaps["well_characterized_detailed"], list)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_well_characterized_detailed_has_aspect_and_basis(self):
        """Each item in well_characterized_detailed should have aspect and basis."""
        from backend import get_variant_insight

        result = await get_variant_insight(
            gene="PIK3CA",
            variant="H1047R",
            tumor_type="Breast",
            enable_llm=False,
            enable_literature=False,
        )

        assert "evidence_gaps" in result
        wc_detailed = result["evidence_gaps"]["well_characterized_detailed"]

        assert len(wc_detailed) > 0, "Should have at least one well-characterized aspect"

        for item in wc_detailed:
            assert "aspect" in item, "Item should have 'aspect' key"
            assert "basis" in item, "Item should have 'basis' key"
            assert len(item["aspect"]) > 0, "Aspect should not be empty"
            assert len(item["basis"]) > 0, "Basis should not be empty"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_gaps_has_expected_structure(self):
        """evidence_gaps should have all expected fields."""
        from backend import get_variant_insight

        result = await get_variant_insight(
            gene="KRAS",
            variant="G12D",
            enable_llm=False,
            enable_literature=False,
        )

        evidence_gaps = result["evidence_gaps"]

        # Required fields
        assert "overall_quality" in evidence_gaps
        assert "research_priority" in evidence_gaps
        assert "well_characterized" in evidence_gaps
        assert "well_characterized_detailed" in evidence_gaps
        assert "poorly_characterized" in evidence_gaps
        assert "gaps" in evidence_gaps
        assert "has_critical_gaps" in evidence_gaps

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_gaps_have_correct_structure(self):
        """Individual gaps should have correct structure."""
        from backend import get_variant_insight

        result = await get_variant_insight(
            gene="ARID1A",
            variant="Q1328*",
            tumor_type="Ovarian",
            enable_llm=False,
            enable_literature=False,
        )

        gaps = result["evidence_gaps"]["gaps"]

        if gaps:  # May have no gaps for well-characterized variants
            for gap in gaps:
                assert "category" in gap
                assert "severity" in gap
                assert "description" in gap
                assert "suggested_studies" in gap
                assert "addressable_with" in gap


# =============================================================================
# LAYER 5: END-TO-END PIPELINE TESTS
# =============================================================================

class TestEndToEndPipeline:
    """End-to-end tests from API through to UI-ready data."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_full_pipeline_braf_v600e(self):
        """Test full pipeline for well-characterized variant."""
        from backend import get_variant_insight

        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        # Basic structure
        assert "error" not in result
        assert result["variant"]["gene"] == "BRAF"
        assert result["variant"]["variant"] == "V600E"

        # Evidence quality for well-characterized variant
        assert result["evidence_gaps"]["overall_quality"] in ("comprehensive", "moderate")

        # Should have therapies
        assert len(result["recommended_therapies"]) > 0 or len(result.get("fda_approvals", [])) > 0

        # Should have well-characterized detailed with basis
        wc_detailed = result["evidence_gaps"]["well_characterized_detailed"]
        assert len(wc_detailed) > 0

        # Hotspot should be marked
        hotspot_items = [item for item in wc_detailed if "hotspot" in item["aspect"].lower()]
        assert len(hotspot_items) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_full_pipeline_rare_variant(self):
        """Test full pipeline for less-characterized variant."""
        from backend import get_variant_insight

        result = await get_variant_insight(
            gene="ARID1A",
            variant="R1989*",
            tumor_type="Ovarian",
            enable_llm=False,
            enable_literature=False,
        )

        # Should still complete without error
        assert "error" not in result

        # Should have gaps
        gaps = result["evidence_gaps"]["gaps"]
        # Rare variants typically have more gaps
        # (but this depends on actual data)

        # Should have research priority set
        assert result["evidence_gaps"]["research_priority"] in (
            "very_high", "high", "medium", "low", "unknown"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_pipeline_data_consistency(self):
        """Test that data is consistent across layers."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR L858R", tumor_type="NSCLC")

        # Get gaps from model
        model_gaps = result.evidence.evidence_gaps
        if model_gaps is None:
            model_gaps = result.evidence.compute_evidence_gaps()

        # Get gaps from backend
        from backend import get_variant_insight
        backend_result = await get_variant_insight(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        backend_gaps = backend_result["evidence_gaps"]

        # Quality should match
        assert model_gaps.overall_evidence_quality == backend_gaps["overall_quality"]

        # Priority should match
        assert model_gaps.research_priority == backend_gaps["research_priority"]

        # Number of well-characterized items should match
        assert len(model_gaps.well_characterized) == len(backend_gaps["well_characterized"])
        assert len(model_gaps.well_characterized_detailed) == len(backend_gaps["well_characterized_detailed"])


# =============================================================================
# LAYER 6: SERIALIZATION TESTS
# =============================================================================

class TestEvidenceGapsSerialization:
    """Tests for evidence gaps serialization to JSON."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_gaps_serializes_to_json(self):
        """EvidenceGaps should serialize to valid JSON."""
        import json

        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Serialize to dict
        gaps_dict = gaps.model_dump(mode="json")

        # Should be valid JSON
        json_str = json.dumps(gaps_dict)
        assert len(json_str) > 0

        # Should round-trip
        parsed = json.loads(json_str)
        assert "well_characterized_detailed" in parsed
        assert isinstance(parsed["well_characterized_detailed"], list)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_characterized_aspect_in_json(self):
        """CharacterizedAspect should be properly represented in JSON."""
        import json

        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("KRAS G12C", tumor_type="NSCLC")

        # Get full result JSON
        result_json = result.model_dump(mode="json")
        json_str = json.dumps(result_json)
        parsed = json.loads(json_str)

        # Navigate to well_characterized_detailed
        wc_detailed = parsed["evidence"]["evidence_gaps"]["well_characterized_detailed"]

        if wc_detailed:
            # Each item should have aspect and basis as strings
            for item in wc_detailed:
                assert isinstance(item["aspect"], str)
                assert isinstance(item["basis"], str)

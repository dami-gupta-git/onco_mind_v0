"""Integration tests for evidence gap detection.

Tests the gap detection functionality with real API calls to verify:
1. Functional data is fetched and gaps are correctly identified
2. Hotspot detection works correctly
3. Discordant evidence detection filters intra-source noise
4. Research priority is computed correctly
5. Well-characterized aspects are properly identified
"""

import pytest
import asyncio

from oncomind.insight_builder import Conductor, ConductorConfig
from oncomind.insight_builder.gap_detector import (
    detect_evidence_gaps,
    _detect_discordant_evidence,
    _has_pathogenic_signal,
    _normalize_source,
)
from oncomind.models.gene_context import (
    is_hotspot_variant,
    is_hotspot_adjacent,
    CANCER_HOTSPOTS,
)
from oncomind.models.evidence.evidence_gaps import (
    GapCategory,
    GapSeverity,
    EvidenceGaps,
)


# =============================================================================
# HOTSPOT DETECTION UNIT TESTS
# =============================================================================

class TestHotspotDetection:
    """Unit tests for hotspot detection functions."""

    def test_is_hotspot_variant_kras_g12d(self):
        """KRAS G12D is a known hotspot."""
        assert is_hotspot_variant("KRAS", "G12D") is True

    def test_is_hotspot_variant_braf_v600e(self):
        """BRAF V600E is a known hotspot."""
        assert is_hotspot_variant("BRAF", "V600E") is True

    def test_is_hotspot_variant_tp53_r248w(self):
        """TP53 R248W is a known hotspot."""
        assert is_hotspot_variant("TP53", "R248W") is True

    def test_is_hotspot_variant_egfr_l858r(self):
        """EGFR L858R is a known hotspot."""
        assert is_hotspot_variant("EGFR", "L858R") is True

    def test_is_hotspot_variant_non_hotspot(self):
        """BRAF V500E is NOT a hotspot."""
        assert is_hotspot_variant("BRAF", "V500E") is False

    def test_is_hotspot_variant_unknown_gene(self):
        """Unknown gene returns False."""
        assert is_hotspot_variant("FAKEGENE", "V600E") is False

    def test_is_hotspot_adjacent_near_braf_600(self):
        """BRAF V598E is adjacent to hotspot V600."""
        is_adj, hotspot = is_hotspot_adjacent("BRAF", "V598E", window=5)
        assert is_adj is True
        assert hotspot == 600

    def test_is_hotspot_adjacent_near_kras_12(self):
        """KRAS G14D is adjacent to hotspot G12."""
        is_adj, hotspot = is_hotspot_adjacent("KRAS", "G14D", window=5)
        assert is_adj is True
        assert hotspot == 12

    def test_is_hotspot_adjacent_at_hotspot(self):
        """Variant AT a hotspot is not 'adjacent'."""
        is_adj, hotspot = is_hotspot_adjacent("BRAF", "V600E", window=5)
        assert is_adj is False
        assert hotspot is None

    def test_is_hotspot_adjacent_far_from_hotspot(self):
        """BRAF V500E is far from hotspots."""
        is_adj, hotspot = is_hotspot_adjacent("BRAF", "V500E", window=5)
        assert is_adj is False
        assert hotspot is None

    def test_hotspot_case_insensitivity(self):
        """Hotspot detection should be case-insensitive for gene."""
        assert is_hotspot_variant("kras", "G12D") is True
        assert is_hotspot_variant("Kras", "g12d") is True

    def test_hotspot_with_p_prefix(self):
        """Hotspot detection works with p. prefix."""
        assert is_hotspot_variant("BRAF", "p.V600E") is True


# =============================================================================
# SOURCE NORMALIZATION TESTS
# =============================================================================

class TestSourceNormalization:
    """Tests for source normalization to prevent duplicate counting."""

    def test_normalize_civic(self):
        """CIViC normalizes to CIViC."""
        assert _normalize_source("CIViC") == "CIViC"
        assert _normalize_source("civic") == "CIViC"

    def test_normalize_vicc_civic(self):
        """VICC/civic normalizes to CIViC."""
        assert _normalize_source("VICC/civic") == "CIViC"
        assert _normalize_source("vicc/civic") == "CIViC"

    def test_normalize_cgi(self):
        """CGI normalizes correctly."""
        assert _normalize_source("CGI") == "CGI"
        assert _normalize_source("vicc/cgi") == "CGI"

    def test_normalize_oncokb(self):
        """OncoKB normalizes correctly."""
        assert _normalize_source("oncokb") == "OncoKB"
        assert _normalize_source("VICC/oncokb") == "OncoKB"

    def test_normalize_molecularmatch(self):
        """MolecularMatch normalizes correctly."""
        assert _normalize_source("molecularmatch") == "MolecularMatch"
        assert _normalize_source("VICC/molecularmatch") == "MolecularMatch"


# =============================================================================
# GAP DETECTION INTEGRATION TESTS (LIVE API)
# =============================================================================

@pytest.mark.integration
class TestGapDetectionIntegration:
    """Integration tests for gap detection with live API calls."""

    @pytest.fixture
    def event_loop(self):
        """Create event loop for async tests."""
        loop = asyncio.new_event_loop()
        yield loop
        loop.close()

    @pytest.mark.asyncio
    async def test_kras_g12d_has_functional_data(self):
        """KRAS G12D should have functional data after MyVariant query fix."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("KRAS G12D", tumor_type="NSCLC")

        # Check functional data is present
        func = result.functional
        assert func.alphamissense_score is not None, "AlphaMissense score should be present"
        assert func.alphamissense_prediction is not None, "AlphaMissense prediction should be present"
        assert func.snpeff_effect is not None, "SnpEff effect should be present"

    @pytest.mark.asyncio
    async def test_kras_g12d_is_hotspot(self):
        """KRAS G12D should be recognized as a hotspot."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("KRAS G12D", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        assert "known cancer hotspot" in gaps.well_characterized

    @pytest.mark.asyncio
    async def test_kras_g12d_no_discordant_noise(self):
        """KRAS G12D should NOT have noisy intra-source discordant gaps."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("KRAS G12D", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Get discordant gaps
        discordant = gaps.get_gaps_by_category(GapCategory.DISCORDANT)

        # Should have 0 or very few discordant gaps (not 10+ noisy ones)
        assert len(discordant) <= 2, f"Too many discordant gaps: {len(discordant)}"

        # Check no gap has repeated "CIViC, CIViC, CIViC..." pattern
        for gap in discordant:
            assert "CIViC, CIViC" not in gap.description, \
                f"Duplicate CIViC sources in gap: {gap.description}"

    @pytest.mark.asyncio
    async def test_braf_v600e_comprehensive_evidence(self):
        """BRAF V600E in Melanoma should have comprehensive evidence."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # BRAF V600E is well-studied
        assert gaps.overall_evidence_quality in ("comprehensive", "moderate")
        assert "known cancer hotspot" in gaps.well_characterized
        assert "computational pathogenicity" in gaps.well_characterized

    @pytest.mark.asyncio
    async def test_rare_variant_has_gaps(self):
        """A rare variant should have more gaps than well-studied ones."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # Use a less common variant
            result = await conductor.run("ARID1A Q1328*", tumor_type="Ovarian")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have at least some gaps
        assert len(gaps.gaps) > 0, "Rare variant should have evidence gaps"

    @pytest.mark.asyncio
    async def test_evidence_quality_levels(self):
        """Test that evidence quality is computed and is valid."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR L858R", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        valid_qualities = ("comprehensive", "moderate", "limited", "minimal", "unknown")
        assert gaps.overall_evidence_quality in valid_qualities

    @pytest.mark.asyncio
    async def test_research_priority_levels(self):
        """Test that research priority is computed and is valid."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("PIK3CA H1047R", tumor_type="Breast")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        valid_priorities = ("very_high", "high", "medium", "low", "unknown")
        assert gaps.research_priority in valid_priorities


# =============================================================================
# HOTSPOT-ADJACENT VARIANT TESTS
# =============================================================================

@pytest.mark.integration
class TestHotspotAdjacentVariants:
    """Tests for hotspot-adjacent variant detection and gap flagging."""

    @pytest.mark.asyncio
    async def test_hotspot_adjacent_flagged_as_research_opportunity(self):
        """Hotspot-adjacent variants should be flagged for functional study."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # KRAS G14D is near hotspot G12
            result = await conductor.run("KRAS G14D", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have "near hotspot" in well_characterized
        near_hotspot = [w for w in gaps.well_characterized if "near hotspot" in w.lower()]
        assert len(near_hotspot) > 0, "Should recognize variant is near hotspot"

        # Should have FUNCTIONAL gap for rare-near-hotspot
        functional_gaps = gaps.get_gaps_by_category(GapCategory.FUNCTIONAL)
        hotspot_related = [g for g in functional_gaps if "hotspot" in g.description.lower()]
        assert len(hotspot_related) > 0, "Should flag functional gap for hotspot-adjacent"


# =============================================================================
# EVIDENCE GAPS MODEL TESTS
# =============================================================================

class TestEvidenceGapsModel:
    """Unit tests for EvidenceGaps model methods."""

    def test_has_critical_gaps(self):
        """Test has_critical_gaps method."""
        from oncomind.models.evidence.evidence_gaps import EvidenceGaps, EvidenceGap

        gaps = EvidenceGaps(
            gaps=[
                EvidenceGap(
                    category=GapCategory.CLINICAL,
                    severity=GapSeverity.CRITICAL,
                    description="No clinical evidence"
                )
            ]
        )
        assert gaps.has_critical_gaps() is True

        gaps_minor = EvidenceGaps(
            gaps=[
                EvidenceGap(
                    category=GapCategory.PREVALENCE,
                    severity=GapSeverity.MINOR,
                    description="Prevalence unknown"
                )
            ]
        )
        assert gaps_minor.has_critical_gaps() is False

    def test_get_gaps_by_category(self):
        """Test filtering gaps by category."""
        from oncomind.models.evidence.evidence_gaps import EvidenceGaps, EvidenceGap

        gaps = EvidenceGaps(
            gaps=[
                EvidenceGap(
                    category=GapCategory.FUNCTIONAL,
                    severity=GapSeverity.SIGNIFICANT,
                    description="Functional gap"
                ),
                EvidenceGap(
                    category=GapCategory.CLINICAL,
                    severity=GapSeverity.CRITICAL,
                    description="Clinical gap"
                ),
                EvidenceGap(
                    category=GapCategory.FUNCTIONAL,
                    severity=GapSeverity.MINOR,
                    description="Another functional gap"
                ),
            ]
        )

        functional = gaps.get_gaps_by_category(GapCategory.FUNCTIONAL)
        assert len(functional) == 2

        clinical = gaps.get_gaps_by_category(GapCategory.CLINICAL)
        assert len(clinical) == 1

    def test_top_gaps_sorted_by_severity(self):
        """Test top_gaps returns gaps sorted by severity."""
        from oncomind.models.evidence.evidence_gaps import EvidenceGaps, EvidenceGap

        gaps = EvidenceGaps(
            gaps=[
                EvidenceGap(
                    category=GapCategory.PREVALENCE,
                    severity=GapSeverity.MINOR,
                    description="Minor gap"
                ),
                EvidenceGap(
                    category=GapCategory.CLINICAL,
                    severity=GapSeverity.CRITICAL,
                    description="Critical gap"
                ),
                EvidenceGap(
                    category=GapCategory.FUNCTIONAL,
                    severity=GapSeverity.SIGNIFICANT,
                    description="Significant gap"
                ),
            ]
        )

        top = gaps.top_gaps(n=3)
        assert top[0].severity == GapSeverity.CRITICAL
        assert top[1].severity == GapSeverity.SIGNIFICANT
        assert top[2].severity == GapSeverity.MINOR

    def test_to_dict_for_llm(self):
        """Test LLM dict conversion includes required fields."""
        from oncomind.models.evidence.evidence_gaps import EvidenceGaps, EvidenceGap

        gaps = EvidenceGaps(
            gaps=[
                EvidenceGap(
                    category=GapCategory.DISCORDANT,
                    severity=GapSeverity.SIGNIFICANT,
                    description="Conflicting drug response",
                    suggested_studies=["Meta-analysis"]
                )
            ],
            overall_evidence_quality="moderate",
            research_priority="medium",
            well_characterized=["hotspot"],
            poorly_characterized=["resistance"],
        )

        llm_dict = gaps.to_dict_for_llm()

        assert llm_dict["overall_quality"] == "moderate"
        assert llm_dict["research_priority"] == "medium"
        assert "hotspot" in llm_dict["well_characterized"]
        assert "resistance" in llm_dict["knowledge_gaps"]
        assert len(llm_dict["conflicting_evidence"]) == 1

    def test_to_summary(self):
        """Test human-readable summary generation."""
        from oncomind.models.evidence.evidence_gaps import EvidenceGaps, EvidenceGap

        gaps = EvidenceGaps(
            gaps=[
                EvidenceGap(
                    category=GapCategory.CLINICAL,
                    severity=GapSeverity.CRITICAL,
                    description="No clinical evidence"
                )
            ],
            overall_evidence_quality="limited",
            research_priority="high",
            well_characterized=["computational pathogenicity"],
            poorly_characterized=["clinical data"],
        )

        summary = gaps.to_summary()

        assert "LIMITED" in summary
        assert "HIGH" in summary
        assert "computational pathogenicity" in summary
        assert "No clinical evidence" in summary


# =============================================================================
# DISCORDANT EVIDENCE DETECTION TESTS
# =============================================================================

class TestDiscordantEvidenceDetection:
    """Tests for discordant evidence detection logic."""

    def test_no_conflict_when_same_source(self):
        """No conflict should be flagged when only CIViC data exists."""
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="KRAS:G12D",
                gene="KRAS",
                variant="G12D"
            ),
            civic_evidence=[
                CIViCEvidence(drugs=["Erlotinib"], clinical_significance="Resistance"),
                CIViCEvidence(drugs=["Erlotinib"], clinical_significance="Sensitivity"),
            ]
        )

        conflicts = _detect_discordant_evidence(evidence)

        # Should NOT flag intra-source conflict
        assert len(conflicts) == 0, f"Intra-source conflict should not be flagged: {conflicts}"

    def test_combination_therapy_not_flagged(self):
        """Combination therapies should not be mixed with monotherapy."""
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="KRAS:G12D",
                gene="KRAS",
                variant="G12D"
            ),
            civic_evidence=[
                # Monotherapy resistance
                CIViCEvidence(drugs=["Erlotinib"], clinical_significance="Resistance"),
                # Combination therapy sensitivity (different context)
                CIViCEvidence(drugs=["Erlotinib", "Sotorasib"], clinical_significance="Sensitivity"),
            ]
        )

        conflicts = _detect_discordant_evidence(evidence)

        # Should NOT flag conflict - different contexts (mono vs combo)
        assert len(conflicts) == 0, f"Combo therapy should not create conflict: {conflicts}"


# =============================================================================
# CANCER HOTSPOTS DATA TESTS
# =============================================================================

class TestCancerHotspotsData:
    """Tests to verify cancer hotspots data is correct."""

    def test_all_major_genes_have_hotspots(self):
        """Major cancer genes should have hotspot definitions."""
        major_genes = ["BRAF", "KRAS", "EGFR", "PIK3CA", "TP53", "IDH1", "KIT"]
        for gene in major_genes:
            assert gene in CANCER_HOTSPOTS, f"{gene} should have hotspots defined"

    def test_hotspot_values_are_positive_integers(self):
        """All hotspot codons should be positive integers."""
        for gene, codons in CANCER_HOTSPOTS.items():
            for codon in codons:
                assert isinstance(codon, int), f"{gene} codon {codon} should be int"
                assert codon > 0, f"{gene} codon {codon} should be positive"

    def test_no_duplicate_hotspots_per_gene(self):
        """Each gene should not have duplicate hotspot codons."""
        for gene, codons in CANCER_HOTSPOTS.items():
            assert len(codons) == len(set(codons)), \
                f"{gene} has duplicate hotspot codons"

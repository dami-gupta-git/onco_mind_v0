"""Integration tests for core evidence gap detection.

Tests the core gap detection functionality:
1. Hotspot detection
2. Source normalization
3. Gap detection integration with live APIs
4. Hotspot-adjacent variant detection
5. EvidenceGaps model methods
6. Cancer hotspots data validation

For specialized tests, see:
- test_discordant_evidence.py: Conflict detection tests
- test_tumor_matching.py: Tumor type matching tests
- test_civic_integration.py: CIViC identifier tests
"""

import pytest
import asyncio

from oncomind.insight_builder import Conductor, ConductorConfig
from oncomind.insight_builder.gap_detector import (
    detect_evidence_gaps,
    _normalize_source,
)
from oncomind.models.gene_context import (
    is_hotspot_variant,
    is_hotspot_adjacent,
)
from oncomind.api.hotspots import get_cancer_hotspots
from oncomind.models.extracted.evidence_gaps import (
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
        """BRAF V598E is adjacent to hotspot V597 (the nearest)."""
        is_adj, hotspot = is_hotspot_adjacent("BRAF", "V598E", window=5)
        assert is_adj is True
        # Returns the NEAREST hotspot (597 is 1 codon away, 600 is 2 away)
        assert hotspot == 597

    def test_is_hotspot_adjacent_near_kras_12(self):
        """KRAS G14D is adjacent to hotspot G13 (the nearest)."""
        is_adj, hotspot = is_hotspot_adjacent("KRAS", "G14D", window=5)
        assert is_adj is True
        # Returns the NEAREST hotspot (13 is 1 codon away, 12 is 2 away)
        assert hotspot == 13

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

        assert "Known Cancer Hotspot" in gaps.well_characterized

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
        """BRAF V600E in Melanoma should have comprehensive evidence.

        This variant is FDA-approved with multiple drugs and well-studied.
        FDA vs VICC resistance reports are NOT flagged as conflicts because
        VICC resistance typically represents acquired resistance (expected
        clinical behavior) rather than intrinsic resistance.
        """
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # BRAF V600E is well-studied - should be comprehensive or moderate
        assert gaps.overall_evidence_quality in ("comprehensive", "moderate")
        # Check case-insensitively
        well_char_lower = [w.lower() for w in gaps.well_characterized]
        assert "known cancer hotspot" in well_char_lower
        assert "computational pathogenicity" in well_char_lower

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
        from oncomind.models.extracted.evidence_gaps import EvidenceGaps, EvidenceGap

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
        from oncomind.models.extracted.evidence_gaps import EvidenceGaps, EvidenceGap

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
        from oncomind.models.extracted.evidence_gaps import EvidenceGaps, EvidenceGap

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
        from oncomind.models.extracted.evidence_gaps import EvidenceGaps, EvidenceGap

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
        from oncomind.models.extracted.evidence_gaps import EvidenceGaps, EvidenceGap

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
# CANCER HOTSPOTS DATA TESTS
# =============================================================================

class TestCancerHotspotsData:
    """Tests to verify cancer hotspots data is correct."""

    def test_all_major_genes_have_hotspots(self):
        """Major cancer genes should have hotspot definitions."""
        major_genes = ["BRAF", "KRAS", "EGFR", "PIK3CA", "TP53", "IDH1", "KIT"]
        for gene in major_genes:
            codons = get_cancer_hotspots(gene)
            assert len(codons) > 0, f"{gene} should have hotspots defined"

    def test_hotspot_values_are_positive_integers(self):
        """All hotspot codons should be positive integers."""
        # Test a subset of major genes
        major_genes = ["BRAF", "KRAS", "EGFR", "PIK3CA", "TP53"]
        for gene in major_genes:
            codons = get_cancer_hotspots(gene)
            for codon in codons:
                assert isinstance(codon, int), f"{gene} codon {codon} should be int"
                assert codon > 0, f"{gene} codon {codon} should be positive"

    def test_no_duplicate_hotspots_per_gene(self):
        """Each gene should not have duplicate hotspot codons."""
        # Test a subset of major genes
        major_genes = ["BRAF", "KRAS", "EGFR", "PIK3CA", "TP53"]
        for gene in major_genes:
            codons = get_cancer_hotspots(gene)
            assert len(codons) == len(set(codons)), \
                f"{gene} has duplicate hotspot codons"


# =============================================================================
# CIVIC SENSITIVITY/RESISTANCE INTEGRATION TESTS
# =============================================================================

@pytest.mark.integration
class TestCivicSensitivityResistance:
    """Integration tests for CIViC is_sensitivity and is_resistance in gap detection.

    These tests verify that:
    1. CIViC evidence with direction=SUPPORTS is correctly categorized
    2. CIViC evidence with direction=DOES_NOT_SUPPORT is excluded
    3. Gap detector correctly counts sensitive vs resistant drugs from CIViC
    """

    @pytest.mark.asyncio
    async def test_braf_v600e_has_sensitivity_evidence(self):
        """BRAF V600E should have CIViC sensitivity evidence for BRAF inhibitors."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        # Get CIViC evidence items
        civic_evidence = result.evidence.civic_evidence

        # Should have sensitivity evidence (vemurafenib, dabrafenib, etc.)
        sensitivity_evidence = [e for e in civic_evidence if e.is_sensitivity]
        assert len(sensitivity_evidence) > 0, "BRAF V600E should have CIViC sensitivity evidence"

        # Verify these are SUPPORTS direction
        for e in sensitivity_evidence:
            if e.evidence_direction:
                assert e.evidence_direction.upper() == "SUPPORTS", (
                    f"is_sensitivity=True but direction={e.evidence_direction}"
                )

    @pytest.mark.asyncio
    async def test_egfr_t790m_has_resistance_evidence(self):
        """EGFR T790M should have CIViC resistance evidence."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR T790M", tumor_type="NSCLC")

        # Get CIViC evidence items
        civic_evidence = result.evidence.civic_evidence

        # T790M is a resistance mutation - should have resistance evidence
        resistance_evidence = [e for e in civic_evidence if e.is_resistance]
        assert len(resistance_evidence) > 0, "EGFR T790M should have CIViC resistance evidence"

        # Verify these are SUPPORTS direction with RESIST in significance
        for e in resistance_evidence:
            assert "RESIST" in e.clinical_significance.upper(), (
                f"is_resistance=True but significance={e.clinical_significance}"
            )

    @pytest.mark.asyncio
    async def test_does_not_support_excluded_from_sensitivity(self):
        """Evidence with DOES_NOT_SUPPORT direction should not be counted as sensitivity."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        civic_evidence = result.evidence.civic_evidence

        # Find evidence with "Sensitivity/Response" significance
        sensitivity_topic_evidence = [
            e for e in civic_evidence
            if e.clinical_significance and (
                "SENSITIV" in e.clinical_significance.upper() or
                "RESPONSE" in e.clinical_significance.upper()
            )
        ]

        # For each sensitivity topic, verify is_sensitivity respects direction
        for e in sensitivity_topic_evidence:
            if e.evidence_direction:
                direction = e.evidence_direction.upper()
                if direction == "DOES_NOT_SUPPORT":
                    assert e.is_sensitivity is False, (
                        f"DOES_NOT_SUPPORT should give is_sensitivity=False, "
                        f"got is_sensitivity={e.is_sensitivity}"
                    )
                elif direction == "SUPPORTS":
                    assert e.is_sensitivity is True, (
                        f"SUPPORTS should give is_sensitivity=True"
                    )

    @pytest.mark.asyncio
    async def test_civic_evidence_direction_in_gap_detector(self):
        """Gap detector should only count SUPPORTS evidence for drug categorization."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        civic_evidence = result.evidence.civic_evidence

        # Count sensitivity evidence (only SUPPORTS)
        sensitivity_count = sum(1 for e in civic_evidence if e.is_sensitivity)

        # Count all sensitivity-topic evidence (regardless of direction)
        all_sensitivity_topic = sum(
            1 for e in civic_evidence
            if e.clinical_significance and (
                "SENSITIV" in e.clinical_significance.upper() or
                "RESPONSE" in e.clinical_significance.upper()
            )
        )

        # If there's DOES_NOT_SUPPORT evidence, counts should differ
        # Or at minimum, sensitivity_count <= all_sensitivity_topic
        assert sensitivity_count <= all_sensitivity_topic, (
            f"is_sensitivity count ({sensitivity_count}) should not exceed "
            f"all sensitivity topic count ({all_sensitivity_topic})"
        )

    @pytest.mark.asyncio
    async def test_sensitivity_resistance_mutually_exclusive(self):
        """is_sensitivity and is_resistance should be mutually exclusive."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        civic_evidence = result.evidence.civic_evidence

        # No evidence should have both is_sensitivity=True AND is_resistance=True
        for e in civic_evidence:
            assert not (e.is_sensitivity and e.is_resistance), (
                f"Evidence {e.evidence_id} has both is_sensitivity=True and is_resistance=True"
            )


# =============================================================================
# FDA APPROVAL GAP TESTS
# =============================================================================

@pytest.mark.integration
class TestFDAApprovalGap:
    """Integration tests for FDA-approved therapy gap detection.

    These tests verify that variants without FDA-approved therapies
    correctly show the "No FDA-approved therapy" gap.
    """

    @pytest.mark.asyncio
    async def test_variant_without_fda_approval_shows_gap(self):
        """ERBB2 S310F in Bladder Cancer has no FDA-approved therapy - should show gap."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("ERBB2 S310F", tumor_type="Bladder Cancer")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have "No FDA-approved therapy" gap
        drug_response_gaps = gaps.get_gaps_by_category(GapCategory.DRUG_RESPONSE)
        fda_gap = [g for g in drug_response_gaps if "no fda-approved therapy" in g.description.lower()]

        assert len(fda_gap) > 0, (
            f"Expected 'No FDA-approved therapy' gap for ERBB2 S310F. "
            f"Got gaps: {[g.description for g in drug_response_gaps]}"
        )

        # The gap should be SIGNIFICANT severity
        assert fda_gap[0].severity == GapSeverity.SIGNIFICANT

    @pytest.mark.asyncio
    async def test_variant_without_fda_approval_still_has_drug_response(self):
        """ERBB2 S310F has VICC/CGI drug response data but no FDA approval."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("ERBB2 S310F", tumor_type="Bladder Cancer")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have drug response data in well_characterized (from VICC/CGI)
        well_char_lower = [w.lower() for w in gaps.well_characterized]
        assert "drug response data" in well_char_lower, (
            "Should have 'drug response data' well-characterized from VICC/CGI data"
        )

        # But should NOT have "FDA-approved therapy" in well_characterized
        assert "fda-approved therapy" not in well_char_lower, (
            "Should NOT have 'FDA-approved therapy' since there's none"
        )

    @pytest.mark.asyncio
    async def test_variant_with_fda_approval_no_gap(self):
        """BRAF V600E in Melanoma has FDA-approved therapy - should NOT have FDA gap."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should NOT have "No FDA-approved therapy" gap
        drug_response_gaps = gaps.get_gaps_by_category(GapCategory.DRUG_RESPONSE)
        fda_gap = [g for g in drug_response_gaps if "no fda-approved therapy" in g.description.lower()]

        assert len(fda_gap) == 0, (
            f"BRAF V600E Melanoma should NOT have 'No FDA-approved therapy' gap. "
            f"Got: {[g.description for g in fda_gap]}"
        )

        # Should have "FDA-approved therapy" in well_characterized
        well_char_lower = [w.lower() for w in gaps.well_characterized]
        assert "fda-approved therapy" in well_char_lower, (
            "BRAF V600E Melanoma should have 'FDA-approved therapy' well-characterized"
        )

    @pytest.mark.asyncio
    async def test_fda_biomarker_evidence_determines_approval(self):
        """FDA approval is determined by fda_biomarker_evidence, not CGI/VICC."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("ERBB2 S310F", tumor_type="Bladder Cancer")

        # Verify ERBB2 S310F has VICC evidence but no FDA biomarker evidence
        assert len(result.evidence.vicc_evidence) > 0, "Should have VICC evidence"
        assert len(result.evidence.fda_biomarker_evidence) == 0, (
            "Should NOT have FDA biomarker evidence"
        )

        # This confirms the gap is correctly based on FDA data, not VICC

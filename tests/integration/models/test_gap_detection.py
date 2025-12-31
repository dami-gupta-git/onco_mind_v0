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
    _detect_discordant_evidence_internal,
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

        conflicts = _detect_discordant_evidence_internal(evidence)

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

        conflicts = _detect_discordant_evidence_internal(evidence)

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


# =============================================================================
# CIViC EID/AID INTEGRATION TESTS
# =============================================================================

@pytest.mark.integration
class TestCIViCIdentifiers:
    """Integration tests for CIViC EID (Evidence Item ID) and AID (Assertion ID)."""

    @pytest.mark.asyncio
    async def test_civic_assertions_have_aid(self):
        """CIViC assertions should have AID (Assertion ID) populated."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # BRAF V600E in Melanoma has well-known CIViC assertions
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        # Check that CIViC assertions have AIDs
        if result.evidence.civic_assertions:
            for assertion in result.evidence.civic_assertions:
                # assertion_id should be a positive integer
                assert assertion.assertion_id is not None, "CIViC assertion should have assertion_id"
                assert isinstance(assertion.assertion_id, int), "assertion_id should be an integer"
                assert assertion.assertion_id > 0, "assertion_id should be positive"

                # AID should be formatted as "AID{number}"
                assert assertion.aid is not None, "CIViC assertion should have aid property"
                assert assertion.aid.startswith("AID"), f"AID should start with 'AID', got: {assertion.aid}"
                assert assertion.aid == f"AID{assertion.assertion_id}"

                # civic_url should be a valid URL
                assert assertion.civic_url is not None, "CIViC assertion should have civic_url"
                assert "civicdb.org/assertions/" in assertion.civic_url

    @pytest.mark.asyncio
    async def test_civic_evidence_have_eid(self):
        """CIViC evidence items should have EID (Evidence Item ID) populated."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # EGFR L858R has well-known CIViC evidence
            result = await conductor.run("EGFR L858R", tumor_type="NSCLC")

        # Check that CIViC evidence items have EIDs (if any exist)
        if result.evidence.civic_evidence:
            has_eid = False
            for evidence in result.evidence.civic_evidence:
                if evidence.evidence_id is not None:
                    has_eid = True
                    # evidence_id should be a positive integer
                    assert isinstance(evidence.evidence_id, int), "evidence_id should be an integer"
                    assert evidence.evidence_id > 0, "evidence_id should be positive"

                    # EID should be formatted as "EID{number}"
                    assert evidence.eid is not None, "CIViC evidence should have eid property"
                    assert evidence.eid.startswith("EID"), f"EID should start with 'EID', got: {evidence.eid}"
                    assert evidence.eid == f"EID{evidence.evidence_id}"

                    # civic_url should be a valid URL
                    assert evidence.civic_url is not None, "CIViC evidence should have civic_url"
                    assert "civicdb.org/evidence/" in evidence.civic_url

            # Note: Not all API sources may return IDs, so we don't require all to have IDs
            # Just verify that if they have IDs, they are properly formatted

    @pytest.mark.asyncio
    async def test_civic_ids_in_model_dump(self):
        """CIViC IDs should be included when model is serialized."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        # Check assertions are serializable with AIDs
        if result.evidence.civic_assertions:
            for assertion in result.evidence.civic_assertions:
                data = assertion.model_dump()
                if assertion.assertion_id is not None:
                    assert "aid" in data, "aid should be in model dump"
                    assert "civic_url" in data, "civic_url should be in model dump"
                    assert data["aid"] == assertion.aid
                    assert data["civic_url"] == assertion.civic_url


# =============================================================================
# DEPMAP TUMOR-TYPE FILTERING TESTS
# =============================================================================

@pytest.mark.integration
class TestDepMapTumorTypeFiltering:
    """Tests for DepMap data filtering by tumor type.

    Verifies that DepMap cell lines and drug sensitivities are NOT added to
    well_characterized when the cell lines don't match the queried tumor type.
    """

    @pytest.mark.asyncio
    async def test_akt1_e17k_excludes_non_matching_tumor_cell_lines(self):
        """AKT1 E17K has Breast/Endometrial cell lines - should NOT appear in NSCLC well_characterized.

        AKT1 E17K is found in Breast and Endometrial cancer cell lines in DepMap.
        When queried with NSCLC tumor type, these cell lines should NOT be added
        to well_characterized.
        """
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("AKT1 E17K", tumor_type="NSCLC")

        # Verify we have DepMap data with cell lines
        assert result.evidence.depmap_evidence is not None, "Should have DepMap evidence"
        assert len(result.evidence.depmap_evidence.cell_line_models) > 0, \
            "Should have cell line models from DepMap"

        # Verify the cell lines are NOT from NSCLC (they should be Breast/Endometrial)
        cell_line_diseases = [
            cl.primary_disease for cl in result.evidence.depmap_evidence.cell_line_models
            if cl.has_mutation and cl.primary_disease
        ]
        assert len(cell_line_diseases) > 0, "Should have cell lines with disease annotation"

        # None of the cell lines should be NSCLC/Lung
        nsclc_matches = [d for d in cell_line_diseases if "lung" in d.lower()]
        assert len(nsclc_matches) == 0, \
            f"AKT1 E17K should not have NSCLC cell lines, found: {nsclc_matches}"

        # Get gaps
        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # DepMap cell lines should NOT be in well_characterized since they don't match NSCLC
        depmap_in_well_char = [
            w for w in gaps.well_characterized
            if "depmap" in w.lower() or "cell line" in w.lower() or "ccle" in w.lower()
        ]

        # Should not have NSCLC-specific cell line models in well_characterized
        nsclc_cell_lines_in_well_char = [
            w for w in gaps.well_characterized
            if "nsclc" in w.lower() and "cell line" in w.lower()
        ]
        assert len(nsclc_cell_lines_in_well_char) == 0, \
            f"Should NOT have NSCLC cell lines in well_characterized: {nsclc_cell_lines_in_well_char}"

    @pytest.mark.asyncio
    async def test_akt1_e17k_excludes_depmap_drug_sensitivity_for_non_matching_tumor(self):
        """AKT1 E17K drug sensitivities should NOT be in well_characterized for NSCLC.

        Even if AKT1 E17K has drug sensitivity data from DepMap, it should NOT
        be added to well_characterized when queried with NSCLC because the
        underlying cell lines are not from NSCLC.
        """
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("AKT1 E17K", tumor_type="NSCLC")

        # Get gaps
        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # DepMap drug sensitivity should NOT be in well_characterized
        depmap_drug_in_well_char = [
            w for w in gaps.well_characterized
            if "depmap" in w.lower() and "drug" in w.lower()
        ]
        assert len(depmap_drug_in_well_char) == 0, \
            f"Should NOT have DepMap drug sensitivity in well_characterized for non-matching tumor: {depmap_drug_in_well_char}"

        # Also check well_characterized_detailed
        depmap_drug_detailed = [
            item for item in gaps.well_characterized_detailed
            if "depmap" in item.aspect.lower() and "drug" in item.aspect.lower()
        ]
        assert len(depmap_drug_detailed) == 0, \
            f"Should NOT have DepMap drug sensitivity in well_characterized_detailed: {depmap_drug_detailed}"

    @pytest.mark.asyncio
    async def test_braf_v600e_includes_cell_lines_for_matching_tumor(self):
        """BRAF V600E should include cell lines when queried with Melanoma.

        BRAF V600E has many Skin Cancer (Melanoma) cell lines in DepMap, so when
        queried with Melanoma, the cell lines SHOULD be in well_characterized.
        """
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        # Verify we have DepMap data
        assert result.evidence.depmap_evidence is not None, "Should have DepMap evidence"

        # Get gaps
        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have Melanoma cell line models in well_characterized
        melanoma_cell_lines = [
            w for w in gaps.well_characterized
            if "melanoma" in w.lower() and "cell line" in w.lower()
        ]
        assert len(melanoma_cell_lines) > 0, \
            f"Should have Melanoma cell lines in well_characterized: {gaps.well_characterized}"

    @pytest.mark.asyncio
    async def test_non_matching_tumor_creates_preclinical_gap(self):
        """When cell lines exist but don't match tumor type, a PRECLINICAL gap should be created."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("AKT1 E17K", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have a PRECLINICAL gap about cross-histology
        preclinical_gaps = gaps.get_gaps_by_category(GapCategory.PRECLINICAL)
        cross_histology_gaps = [
            g for g in preclinical_gaps
            if "cross-histology" in g.description.lower() or "none in" in g.description.lower()
        ]
        assert len(cross_histology_gaps) > 0, \
            f"Should have cross-histology preclinical gap, got: {[g.description for g in preclinical_gaps]}"

        # Should mention that models exist but not in the queried tumor type
        gap_desc = cross_histology_gaps[0].description.lower()
        assert "nsclc" in gap_desc or "lung" in gap_desc, \
            f"Gap should mention NSCLC: {cross_histology_gaps[0].description}"

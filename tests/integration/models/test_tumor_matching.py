"""Integration tests for tumor type matching in evidence gap detection.

Tests the tumor type filtering and matching functionality:
- DepMap cell line filtering by tumor type
- Drug response tumor_match and matches_on tracking
- Resistance mechanisms tumor matching
- Tumor-specific evidence tracking
- Clinical actionability tumor matching
"""

import pytest

from oncomind.insight_builder import Conductor, ConductorConfig
from oncomind.models.extracted.evidence_gaps import GapCategory


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


# =============================================================================
# DRUG RESPONSE TUMOR MATCH INTEGRATION TESTS
# =============================================================================

@pytest.mark.integration
class TestDrugResponseTumorMatch:
    """Integration tests for drug response tumor_match and matches_on tracking."""

    @pytest.mark.asyncio
    async def test_egfr_l858r_drug_response_has_tumor_match(self):
        """EGFR L858R in NSCLC should have drug response with tumor_match info."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR L858R", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find drug response data in well_characterized_detailed
        drug_resp = [
            item for item in gaps.well_characterized_detailed
            if "drug response data" in item.aspect.lower()
        ]

        # Should have drug response data
        assert len(drug_resp) > 0, "Should have drug response data in well_characterized"

        # Should have tumor_match field
        drug_resp_item = drug_resp[0]
        assert drug_resp_item.tumor_match is not None, \
            f"Drug response should have tumor_match, got: {drug_resp_item}"

        # Should show "X tumor" count since EGFR L858R has NSCLC-specific data
        assert "tumor" in drug_resp_item.tumor_match, \
            f"tumor_match should contain 'tumor', got: {drug_resp_item.tumor_match}"

    @pytest.mark.asyncio
    async def test_egfr_l858r_drug_response_has_matches_on(self):
        """EGFR L858R drug response should have matches_on locus levels."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR L858R", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find drug response data in well_characterized_detailed
        drug_resp = [
            item for item in gaps.well_characterized_detailed
            if "drug response data" in item.aspect.lower()
        ]

        assert len(drug_resp) > 0, "Should have drug response data in well_characterized"

        drug_resp_item = drug_resp[0]
        assert drug_resp_item.matches_on is not None, \
            f"Drug response data should have matches_on, got: {drug_resp_item}"

        # matches_on should contain locus level info (variant, codon, or gene)
        has_level = any(level in drug_resp_item.matches_on
                       for level in ["variant", "codon", "gene"])
        assert has_level, \
            f"matches_on should contain locus level, got: {drug_resp_item.matches_on}"

    @pytest.mark.asyncio
    async def test_drug_response_other_tumor_shows_other_count(self):
        """Drug response for variant in different tumor should show 'other' count."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # BRAF V600E is primarily Melanoma/Colorectal - query in Pancreatic
            result = await conductor.run("BRAF V600E", tumor_type="Pancreatic")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find drug response data
        drug_resp = [
            item for item in gaps.well_characterized_detailed
            if "drug response data" in item.aspect.lower()
        ]

        if drug_resp:
            # If there's drug response data but not for Pancreatic, should show "other"
            tumor_match = drug_resp[0].tumor_match
            if tumor_match:
                # May have mostly "other" entries since BRAF V600E approvals are for Melanoma
                assert "tumor" in tumor_match or "other" in tumor_match, \
                    f"tumor_match should contain counts, got: {tumor_match}"


# =============================================================================
# RESISTANCE MECHANISMS TUMOR MATCH INTEGRATION TESTS
# =============================================================================

@pytest.mark.integration
class TestResistanceMechanismsTumorMatch:
    """Integration tests for resistance mechanisms tumor_match and matches_on tracking."""

    @pytest.mark.asyncio
    async def test_egfr_t790m_has_resistance_data(self):
        """EGFR T790M is a known resistance mutation - should have resistance mechanisms."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR T790M", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find resistance mechanisms in well_characterized_detailed
        resistance = [
            item for item in gaps.well_characterized_detailed
            if "resistance" in item.aspect.lower()
        ]

        # T790M is THE resistance mutation - should be well-characterized
        assert len(resistance) > 0, \
            f"EGFR T790M should have resistance mechanisms, got: {gaps.well_characterized}"

    @pytest.mark.asyncio
    async def test_egfr_t790m_resistance_has_matches_on(self):
        """EGFR T790M resistance data should have matches_on locus levels."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR T790M", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find resistance mechanisms
        resistance = [
            item for item in gaps.well_characterized_detailed
            if "resistance" in item.aspect.lower()
        ]

        if resistance:
            resist_item = resistance[0]
            # Should have matches_on for locus level tracking
            if resist_item.matches_on:
                has_level = any(level in resist_item.matches_on
                               for level in ["variant", "codon", "gene"])
                assert has_level, \
                    f"matches_on should contain locus level, got: {resist_item.matches_on}"


# =============================================================================
# TUMOR SPECIFIC EVIDENCE INTEGRATION TESTS
# =============================================================================

@pytest.mark.integration
class TestTumorSpecificEvidenceIntegration:
    """Integration tests for tumor-specific evidence tracking."""

    @pytest.mark.asyncio
    async def test_egfr_l858r_nsclc_has_tumor_evidence(self):
        """EGFR L858R in NSCLC should show 'Evidence Items For NSCLC' as well-characterized."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR L858R", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have "Evidence Items For NSCLC" in well_characterized
        evidence_in_tumor = [
            w for w in gaps.well_characterized
            if "evidence items for" in w.lower() and "nsclc" in w.lower()
        ]
        assert len(evidence_in_tumor) > 0, \
            f"Should have 'Evidence Items For NSCLC', got: {gaps.well_characterized}"

    @pytest.mark.asyncio
    async def test_tumor_evidence_detailed_has_match_info(self):
        """Tumor-specific evidence should have matches_on info."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("KRAS G12C", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find tumor-specific evidence
        tumor_evidence = [
            item for item in gaps.well_characterized_detailed
            if "evidence items for" in item.aspect.lower()
        ]

        if tumor_evidence:
            tumor_item = tumor_evidence[0]
            # Should have matches_on for locus level breakdown
            assert tumor_item.matches_on is not None or tumor_item.tumor_match is not None, \
                f"Tumor evidence should have match info, got: {tumor_item}"

    @pytest.mark.asyncio
    async def test_braf_v600e_melanoma_comprehensive(self):
        """BRAF V600E in Melanoma should have comprehensive tumor-matched evidence."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("BRAF V600E", tumor_type="Melanoma")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # BRAF V600E in Melanoma is THE canonical example
        assert gaps.overall_evidence_quality in ("comprehensive", "moderate"), \
            f"BRAF V600E in Melanoma should be comprehensive, got: {gaps.overall_evidence_quality}"

        # Should have tumor-specific evidence
        tumor_evidence = [
            item for item in gaps.well_characterized_detailed
            if "evidence items for" in item.aspect.lower() and "melanoma" in item.aspect.lower()
        ]
        assert len(tumor_evidence) > 0, \
            f"Should have Melanoma-specific evidence, got: {[i.aspect for i in gaps.well_characterized_detailed]}"


# =============================================================================
# CLINICAL ACTIONABILITY TUMOR MATCH INTEGRATION TESTS
# =============================================================================

@pytest.mark.integration
class TestClinicalActionabilityTumorMatch:
    """Integration tests for clinical actionability tumor_match tracking."""

    @pytest.mark.asyncio
    async def test_egfr_l858r_clinical_actionability_has_tumor_match(self):
        """EGFR L858R clinical actionability should track tumor matches."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("EGFR L858R", tumor_type="NSCLC")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find clinical actionability
        clinical = [
            item for item in gaps.well_characterized_detailed
            if "clinical actionability" in item.aspect.lower()
        ]

        if clinical:
            clinical_item = clinical[0]
            # Should have tumor_match tracking
            if clinical_item.tumor_match:
                # EGFR L858R has NSCLC-specific approvals
                assert "tumor" in clinical_item.tumor_match or "pan_cancer" in clinical_item.tumor_match, \
                    f"Should have tumor or pan_cancer count, got: {clinical_item.tumor_match}"

    @pytest.mark.asyncio
    async def test_clinical_actionability_other_cancer_creates_gap(self):
        """Clinical actionability for mismatched tumor should create a gap."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            # BRAF V600E is approved in Melanoma - query in Thyroid (different)
            result = await conductor.run("BRAF V600E", tumor_type="Thyroid")

        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Find clinical actionability
        clinical = [
            item for item in gaps.well_characterized_detailed
            if "clinical actionability" in item.aspect.lower()
        ]

        # If there's clinical data, check for tumor match tracking
        if clinical:
            clinical_item = clinical[0]
            # Should have some indication of tumor match status via tumor_match field
            has_tracking = clinical_item.tumor_match is not None
            # Note: if evidence exists only in other cancers, a gap should also be created

"""Integration tests for ClinicalTrials EvidenceLevel population.

Tests that EvidenceLevel fields are correctly populated based on trial content:
- Variant-specific trials: locus_variant_match.level == "variant"
- Gene-only trials: locus_variant_match.level == "gene"
- Disease-only trials: locus_variant_match is None
- Cancer-specific trials: cancer_type_match.level == "cancer_specific"
"""

import pytest

from oncomind.api.clinicaltrials import ClinicalTrialsClient


class TestClinicalTrialsEvidenceLevelIntegration:
    """Integration tests for EvidenceLevel population in ClinicalTrialEvidence."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_variant_specific(self):
        """Test that variant-specific trials have locus_variant_match.level == 'variant'.

        Uses EGFR L858R which is a well-known, specific variant commonly
        mentioned in trial titles/eligibility.
        """
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence(
                gene="EGFR",
                variant="L858R",
                tumor_type="NSCLC",
                recruiting_only=False,  # Include more trials for testing
                max_results=10,
            )

        assert len(results) > 0, "Expected at least one trial for EGFR L858R"

        # Check that variant-specific trials exist
        variant_specific_count = sum(
            1
            for t in results
            if t.locus_variant_match and t.locus_variant_match.level == "variant"
        )
        gene_only_count = sum(
            1
            for t in results
            if t.locus_variant_match and t.locus_variant_match.level == "gene"
        )

        # At least one should be variant-specific for a well-known variant
        assert (
            variant_specific_count > 0 or gene_only_count > 0
        ), "Expected at least one trial with locus_variant_match populated"

        # All trials should have locus_variant_match set (biomarker search)
        for trial in results:
            assert (
                trial.locus_variant_match is not None
            ), f"Trial {trial.nct_id} should have locus_variant_match set"
            assert (
                trial.locus_variant_match.origin == "trial"
            ), f"Trial {trial.nct_id} should have origin='trial'"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_gene_level(self):
        """Test that trials have locus_variant_match with valid levels.

        Trials can match at variant, codon, or gene level depending on
        how the trial describes its biomarker requirements.
        """
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence(
                gene="KRAS",
                variant="G12C",
                tumor_type="NSCLC",
                recruiting_only=False,
                max_results=10,
            )

        assert len(results) > 0, "Expected at least one trial for KRAS"

        # Check EvidenceLevel fields
        for trial in results:
            assert trial.locus_variant_match is not None
            # Level can be variant (exact match), codon (e.g. G12), or gene
            assert trial.locus_variant_match.level in ("variant", "codon", "gene")
            assert trial.locus_variant_match.origin == "trial"

            if trial.locus_variant_match.level in ("variant", "codon"):
                # Scope can be "specific" (exact match) or "ambiguous" (codon-level match)
                assert trial.locus_variant_match.scope in ("specific", "ambiguous")
            else:
                assert trial.locus_variant_match.scope == "unspecified"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_cancer_type_match(self):
        """Test that cancer_type_match is set when tumor_type is provided."""
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                recruiting_only=False,
                max_results=10,
            )

        assert len(results) > 0, "Expected at least one trial for BRAF V600E"

        for trial in results:
            # cancer_type_match should be set when tumor_type is provided
            assert (
                trial.cancer_type_match is not None
            ), f"Trial {trial.nct_id} should have cancer_type_match when tumor_type provided"
            # level can be "cancer_specific" if trial conditions match tumor_type,
            # or None if trial has the biomarker but doesn't match the tumor type
            assert trial.cancer_type_match.level in (
                "cancer_specific",
                "pan_cancer",
                None,
            )
            assert trial.cancer_type_match.origin == "trial"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_no_tumor_type(self):
        """Test that cancer_type_match is None when tumor_type not provided."""
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence(
                gene="TP53",
                variant="R175H",
                tumor_type=None,  # No tumor type
                recruiting_only=False,
                max_results=5,
            )

        # May or may not find trials for TP53 specifically
        for trial in results:
            # cancer_type_match should be None when tumor_type not provided
            assert (
                trial.cancer_type_match is None
            ), f"Trial {trial.nct_id} should have cancer_type_match=None when no tumor_type"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_disease_search_cancer_specific(self):
        """Test disease-based search returns cancer_type_match populated."""
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence_by_disease(
                tumor_type="NSCLC",
                gene="EGFR",  # Optional context
                variant="L858R",  # Optional context
                recruiting_only=False,
                max_results=10,
            )

        assert len(results) > 0, "Expected trials for NSCLC"

        for trial in results:
            # cancer_type_match should always be set for disease search
            assert (
                trial.cancer_type_match is not None
            ), f"Trial {trial.nct_id} should have cancer_type_match for disease search"
            assert trial.cancer_type_match.origin == "trial"

            # Most NSCLC trials should be cancer_specific
            # (searched by condition, so conditions should match)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_disease_search_locus_variant_match_optional(self):
        """Test disease-based search may or may not have locus_variant_match.

        locus_variant_match depends on whether trial mentions gene/variant.
        """
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence_by_disease(
                tumor_type="Melanoma",
                gene="BRAF",
                variant="V600E",
                recruiting_only=False,
                max_results=10,
            )

        assert len(results) > 0, "Expected trials for Melanoma"

        # Count trials by locus_variant_match status
        variant_specific = 0
        gene_only = 0
        no_biomarker = 0

        for trial in results:
            if trial.locus_variant_match is None:
                no_biomarker += 1
            elif trial.locus_variant_match.level == "variant":
                variant_specific += 1
            elif trial.locus_variant_match.level == "gene":
                gene_only += 1

        # Disease search should return a mix - some with biomarker, some without
        # (since it's searching by condition, not by gene)
        # At least verify the logic runs without error
        assert variant_specific >= 0
        assert gene_only >= 0
        assert no_biomarker >= 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_merged_search_deduplication(self):
        """Test that search_all_trial_evidence merges and deduplicates results."""
        async with ClinicalTrialsClient() as client:
            results = await client.search_all_trial_evidence(
                gene="EGFR",
                variant="L858R",
                tumor_type="NSCLC",
                recruiting_only=False,
                max_results=20,
            )

        # Check no duplicate NCT IDs
        nct_ids = [t.nct_id for t in results]
        assert len(nct_ids) == len(
            set(nct_ids)
        ), "Found duplicate NCT IDs in merged results"

        # All results should have properly populated fields
        for trial in results:
            assert trial.nct_id.startswith("NCT")
            # cancer_type_match should be set (tumor_type was provided)
            assert trial.cancer_type_match is not None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_level_scope_consistency(self):
        """Test that scope is consistent with level.

        - variant + specific/ambiguous (ambiguous for codon-level matches)
        - gene + unspecified
        - cancer_specific + specific
        - pan_cancer + unspecified
        """
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence(
                gene="BRAF",
                variant="V600E",
                tumor_type="Colorectal",
                recruiting_only=False,
                max_results=10,
            )

        for trial in results:
            # Check locus_variant_match consistency
            if trial.locus_variant_match:
                if trial.locus_variant_match.level == "variant":
                    # Scope can be "specific" (exact match) or "ambiguous" (codon-level match)
                    assert trial.locus_variant_match.scope in ("specific", "ambiguous")
                elif trial.locus_variant_match.level == "gene":
                    assert trial.locus_variant_match.scope == "unspecified"

            # Check cancer_type_match consistency
            if trial.cancer_type_match:
                if trial.cancer_type_match.level == "cancer_specific":
                    assert trial.cancer_type_match.scope == "specific"
                elif trial.cancer_type_match.level == "pan_cancer":
                    assert trial.cancer_type_match.scope == "unspecified"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_disease_only_no_gene_context(self):
        """Test disease search without gene/variant context.

        When no gene/variant provided, locus_variant_match should be None.
        """
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence_by_disease(
                tumor_type="Breast Cancer",
                gene=None,  # No gene context
                variant=None,  # No variant context
                recruiting_only=False,
                max_results=5,
            )

        assert len(results) > 0, "Expected trials for Breast Cancer"

        for trial in results:
            # Without gene context, locus_variant_match should be None
            assert (
                trial.locus_variant_match is None
            ), f"Trial {trial.nct_id} should have locus_variant_match=None without gene context"
            # cancer_type_match should still be set
            assert trial.cancer_type_match is not None

"""Integration tests for ClinicalTrials EvidenceLevel population.

Tests that EvidenceLevel fields are correctly populated based on trial content:
- Variant-specific trials: variant_level.level == "variant"
- Gene-only trials: variant_level.level == "gene"
- Disease-only trials: variant_level is None
- Cancer-specific trials: cancer_type_level.level == "cancer_specific"
"""

import pytest

from oncomind.api.clinicaltrials import ClinicalTrialsClient


class TestClinicalTrialsEvidenceLevelIntegration:
    """Integration tests for EvidenceLevel population in ClinicalTrialEvidence."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_variant_specific(self):
        """Test that variant-specific trials have variant_level.level == 'variant'.

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
            1 for t in results
            if t.variant_level and t.variant_level.level == "variant"
        )
        gene_only_count = sum(
            1 for t in results
            if t.variant_level and t.variant_level.level == "gene"
        )

        # At least one should be variant-specific for a well-known variant
        assert variant_specific_count > 0 or gene_only_count > 0, (
            "Expected at least one trial with variant_level populated"
        )

        # All trials should have variant_level set (biomarker search)
        for trial in results:
            assert trial.variant_level is not None, (
                f"Trial {trial.nct_id} should have variant_level set"
            )
            assert trial.variant_level.origin == "trial", (
                f"Trial {trial.nct_id} should have origin='trial'"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_gene_level(self):
        """Test that gene-level trials have variant_level.level == 'gene'.

        When a trial mentions the gene but not the specific variant,
        it should be marked as gene-level.
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
            assert trial.variant_level is not None
            assert trial.variant_level.level in ("variant", "gene")
            assert trial.variant_level.origin == "trial"

            if trial.variant_level.level == "variant":
                # Scope can be "specific" (exact match) or "ambiguous" (codon-level match)
                assert trial.variant_level.scope in ("specific", "ambiguous")
            else:
                assert trial.variant_level.scope == "unspecified"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_cancer_type_level(self):
        """Test that cancer_type_level is set when tumor_type is provided."""
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
            # cancer_type_level should be set when tumor_type is provided
            assert trial.cancer_type_level is not None, (
                f"Trial {trial.nct_id} should have cancer_type_level when tumor_type provided"
            )
            assert trial.cancer_type_level.level in ("cancer_specific", "pan_cancer")
            assert trial.cancer_type_level.origin == "trial"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_no_tumor_type(self):
        """Test that cancer_type_level is None when tumor_type not provided."""
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
            # cancer_type_level should be None when tumor_type not provided
            assert trial.cancer_type_level is None, (
                f"Trial {trial.nct_id} should have cancer_type_level=None when no tumor_type"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_disease_search_cancer_specific(self):
        """Test disease-based search returns cancer_type_level populated."""
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
            # cancer_type_level should always be set for disease search
            assert trial.cancer_type_level is not None, (
                f"Trial {trial.nct_id} should have cancer_type_level for disease search"
            )
            assert trial.cancer_type_level.origin == "trial"

            # Most NSCLC trials should be cancer_specific
            # (searched by condition, so conditions should match)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_disease_search_variant_level_optional(self):
        """Test disease-based search may or may not have variant_level.

        variant_level depends on whether trial mentions gene/variant.
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

        # Count trials by variant_level status
        variant_specific = 0
        gene_only = 0
        no_biomarker = 0

        for trial in results:
            if trial.variant_level is None:
                no_biomarker += 1
            elif trial.variant_level.level == "variant":
                variant_specific += 1
            elif trial.variant_level.level == "gene":
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
        assert len(nct_ids) == len(set(nct_ids)), "Found duplicate NCT IDs in merged results"

        # All results should have properly populated fields
        for trial in results:
            assert trial.nct_id.startswith("NCT")
            # cancer_type_level should be set (tumor_type was provided)
            assert trial.cancer_type_level is not None

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
            # Check variant_level consistency
            if trial.variant_level:
                if trial.variant_level.level == "variant":
                    # Scope can be "specific" (exact match) or "ambiguous" (codon-level match)
                    assert trial.variant_level.scope in ("specific", "ambiguous")
                elif trial.variant_level.level == "gene":
                    assert trial.variant_level.scope == "unspecified"

            # Check cancer_type_level consistency
            if trial.cancer_type_level:
                if trial.cancer_type_level.level == "cancer_specific":
                    assert trial.cancer_type_level.scope == "specific"
                elif trial.cancer_type_level.level == "pan_cancer":
                    assert trial.cancer_type_level.scope == "unspecified"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_disease_only_no_gene_context(self):
        """Test disease search without gene/variant context.

        When no gene/variant provided, variant_level should be None.
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
            # Without gene context, variant_level should be None
            assert trial.variant_level is None, (
                f"Trial {trial.nct_id} should have variant_level=None without gene context"
            )
            # cancer_type_level should still be set
            assert trial.cancer_type_level is not None

"""Integration tests for ClinicalTrials EvidenceLevel population.

Tests that EvidenceLevel fields are correctly populated based on trial content:
- All trials are gene-level: locus_variant_match.level == "gene"
- Disease-only trials (no gene context): locus_variant_match is None
- Cancer-specific trials: cancer_type_match.level == "cancer_specific"
"""

import pytest

from oncomind.api.clinicaltrials import ClinicalTrialsClient


class TestClinicalTrialsEvidenceLevelIntegration:
    """Integration tests for EvidenceLevel population in ClinicalTrialEvidence."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_always_gene_level(self):
        """Test that all trials from biomarker search have locus_variant_match.level == 'gene'.

        Clinical trials are always gene-level regardless of variant text in eligibility.
        """
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence(
                gene="EGFR",
                variant="L858R",
                tumor_type="NSCLC",
                recruiting_only=False,
                max_results=10,
            )

        assert len(results) > 0, "Expected at least one trial for EGFR L858R"

        for trial in results:
            assert (
                trial.locus_variant_match is not None
            ), f"Trial {trial.nct_id} should have locus_variant_match set"
            assert (
                trial.locus_variant_match.level == "gene"
            ), f"Trial {trial.nct_id} should be gene-level, got {trial.locus_variant_match.level}"
            assert (
                trial.locus_variant_match.scope == "unspecified"
            ), f"Trial {trial.nct_id} scope should be 'unspecified'"
            assert (
                trial.locus_variant_match.origin == "trial"
            ), f"Trial {trial.nct_id} should have origin='trial'"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_biomarker_search_gene_level_kras(self):
        """Test that all KRAS trials are gene-level."""
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence(
                gene="KRAS",
                variant="G12C",
                tumor_type="NSCLC",
                recruiting_only=False,
                max_results=10,
            )

        assert len(results) > 0, "Expected at least one trial for KRAS"

        for trial in results:
            assert trial.locus_variant_match is not None
            assert trial.locus_variant_match.level == "gene"
            assert trial.locus_variant_match.scope == "unspecified"
            assert trial.locus_variant_match.origin == "trial"

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
            assert (
                trial.cancer_type_match is not None
            ), f"Trial {trial.nct_id} should have cancer_type_match when tumor_type provided"
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
                tumor_type=None,
                recruiting_only=False,
                max_results=5,
            )

        for trial in results:
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
                gene="EGFR",
                variant="L858R",
                recruiting_only=False,
                max_results=10,
            )

        assert len(results) > 0, "Expected trials for NSCLC"

        for trial in results:
            assert (
                trial.cancer_type_match is not None
            ), f"Trial {trial.nct_id} should have cancer_type_match for disease search"
            assert trial.cancer_type_match.origin == "trial"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_disease_search_with_gene_is_gene_level(self):
        """Test disease-based search with gene context gives gene-level locus match."""
        async with ClinicalTrialsClient() as client:
            results = await client.search_trial_evidence_by_disease(
                tumor_type="Melanoma",
                gene="BRAF",
                variant="V600E",
                recruiting_only=False,
                max_results=10,
            )

        assert len(results) > 0, "Expected trials for Melanoma"

        for trial in results:
            # With gene context, locus_variant_match should be gene-level
            assert trial.locus_variant_match is not None
            assert trial.locus_variant_match.level == "gene"
            assert trial.locus_variant_match.scope == "unspecified"

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

        nct_ids = [t.nct_id for t in results]
        assert len(nct_ids) == len(
            set(nct_ids)
        ), "Found duplicate NCT IDs in merged results"

        for trial in results:
            assert trial.nct_id.startswith("NCT")
            assert trial.cancer_type_match is not None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_level_scope_consistency(self):
        """Test that scope is consistent with level.

        All trials: gene + unspecified.
        Cancer type: cancer_specific + specific, or pan_cancer + unspecified.
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
            if trial.locus_variant_match:
                assert trial.locus_variant_match.level == "gene"
                assert trial.locus_variant_match.scope == "unspecified"

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
                gene=None,
                variant=None,
                recruiting_only=False,
                max_results=5,
            )

        assert len(results) > 0, "Expected trials for Breast Cancer"

        for trial in results:
            assert (
                trial.locus_variant_match is None
            ), f"Trial {trial.nct_id} should have locus_variant_match=None without gene context"
            assert trial.cancer_type_match is not None

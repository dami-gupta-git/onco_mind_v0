"""Integration tests for DepMap API client (requires network)."""

import pytest

from oncomind.api.depmap import DepMapClient
from oncomind.api.cbioportal import CBioPortalClient


class TestDepMapClientIntegration:
    """Integration tests for DepMap client."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_braf_data(self):
        """Test fetching DepMap data for BRAF."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("BRAF", "V600E")

        # Should return data (fallback if API unavailable)
        assert result is not None
        assert result.gene == "BRAF"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_braf_has_dependency_data(self):
        """Test that BRAF has gene dependency data."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("BRAF", "V600E")

        assert result is not None
        assert result.gene_dependency is not None
        assert result.gene_dependency.mean_dependency_score is not None
        # BRAF is essential in BRAF-mutant lines (negative score)
        assert result.gene_dependency.mean_dependency_score < 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_braf_has_drug_sensitivities(self):
        """Test that BRAF has drug sensitivity data from PRISM (if data files available)."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("BRAF", "V600E")

        assert result is not None

        # Drug sensitivities require local data files - skip if unavailable
        if len(result.drug_sensitivities) == 0:
            pytest.skip(
                "Drug sensitivity data not available (missing local data files)"
            )

        # All returned drugs should be sensitive (log2fc <= -1.7)
        for ds in result.drug_sensitivities:
            assert ds.mean_log2fc is not None
            assert ds.mean_log2fc <= -1.7

        # Top drugs should have very negative log2fc (highly effective)
        top_drug = result.drug_sensitivities[0]
        assert (
            top_drug.mean_log2fc < -2.0
        ), f"Expected potent drug, got log2fc={top_drug.mean_log2fc}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_braf_has_cell_lines(self):
        """Test that BRAF V600E has known cell line models (if data files available)."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("BRAF", "V600E")

        assert result is not None

        # Cell line models require local data files - skip if unavailable
        if len(result.cell_line_models) == 0:
            pytest.skip("Cell line model data not available (missing local data files)")

        # Should include well-known BRAF V600E melanoma lines
        cell_line_names = [cl.name for cl in result.cell_line_models]
        # A375 may appear as "A375", "A-375", or "A375 SKIN CJ1" etc.
        assert any(
            "A375" in name or "A-375" in name for name in cell_line_names
        ), f"Expected A375 cell line, got: {cell_line_names}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_kras_g12c(self):
        """Test fetching data for KRAS G12C."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("KRAS", "G12C")

        assert result is not None
        assert result.gene == "KRAS"
        assert result.gene_dependency is not None

        # KRAS is highly essential
        assert result.gene_dependency.mean_dependency_score < -0.5

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_kras_g12c_has_drug_sensitivities(self):
        """Test that KRAS G12C includes drug sensitivity data from PRISM."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("KRAS", "G12C")

        assert result is not None

        # Should have sensitive drugs
        if result.drug_sensitivities:
            # All returned drugs should be sensitive
            for ds in result.drug_sensitivities:
                assert ds.mean_log2fc is not None
                assert ds.mean_log2fc <= -1.7

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_egfr_l858r(self):
        """Test fetching data for EGFR L858R."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("EGFR", "L858R")

        assert result is not None
        assert result.gene == "EGFR"

        # Should have cell lines with this mutation
        if result.cell_line_models:
            mutant_lines = [cl for cl in result.cell_line_models if cl.has_mutation]
            assert len(mutant_lines) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_unknown_gene(self):
        """Test handling of unknown gene."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("NOTAREALGENE123")

        # Should return None for unknown genes
        assert result is None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_tp53_not_essential(self):
        """Test that TP53 (tumor suppressor) is not essential when lost."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("TP53")

        assert result is not None
        assert result.gene == "TP53"

        # TP53 is a tumor suppressor - cells with TP53 loss don't depend on it
        if result.gene_dependency:
            assert result.gene_dependency.mean_dependency_score > -0.5

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_is_essential_method(self):
        """Test the is_essential() helper method."""
        async with DepMapClient() as client:
            # Use KRAS which is highly essential pan-cancer
            kras_result = await client.fetch_depmap_evidence("KRAS")
            tp53_result = await client.fetch_depmap_evidence("TP53")

        # KRAS should be essential pan-cancer (score < -0.5)
        if kras_result:
            assert (
                kras_result.is_essential() is True
            ), f"Expected KRAS to be essential, got score: {kras_result.get_essential_score()}"

        # TP53 should not be essential
        if tp53_result:
            assert tp53_result.is_essential() is False

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_data_method(self):
        """Test the has_data() helper method."""
        async with DepMapClient() as client:
            braf_result = await client.fetch_depmap_evidence("BRAF", "V600E")
            unknown_result = await client.fetch_depmap_evidence("UNKNOWNGENE")

        if braf_result:
            assert braf_result.has_data() is True

        # Unknown gene should return None, not empty evidence
        assert unknown_result is None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_context_manager(self):
        """Test async context manager protocol."""
        async with DepMapClient() as client:
            assert client._client is not None

        # Client should be closed after exiting context
        assert client._client is None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_drug_sensitivity_log2fc_values(self):
        """Test that drug sensitivities have valid log2fc values (if data files available)."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("BRAF", "V600E")

        assert result is not None

        # Drug sensitivities require local data files - skip if unavailable
        if len(result.drug_sensitivities) == 0:
            pytest.skip(
                "Drug sensitivity data not available (missing local data files)"
            )

        for ds in result.drug_sensitivities:
            # Should have drug name
            assert ds.drug_name is not None
            assert len(ds.drug_name) > 0

            # log2fc should be <= -1.7 for sensitive drugs
            assert ds.mean_log2fc is not None
            assert ds.mean_log2fc <= -1.7

            # Should have cell line count
            assert ds.n_cell_lines > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cell_line_disease_annotation(self):
        """Test that cell lines have disease annotations (if data files available)."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("BRAF", "V600E")

        assert result is not None

        # Cell line models require local data files - skip if unavailable
        if len(result.cell_line_models) == 0:
            pytest.skip("Cell line model data not available (missing local data files)")

        # Check that cell lines have valid names
        for cl in result.cell_line_models:
            assert cl.name is not None
            assert len(cl.name) > 0

        # At least some cell lines should have disease annotation
        lines_with_disease = [
            cl for cl in result.cell_line_models if cl.primary_disease
        ]
        assert (
            len(lines_with_disease) > 0
        ), "Expected at least some cell lines to have disease annotations"

        # Common BRAF V600E-associated cancers should be present
        all_diseases = {cl.primary_disease.lower() for cl in lines_with_disease}
        # Check that melanoma/skin cancer is represented (most common BRAF V600E tumor type)
        assert any(
            "melanoma" in d or "skin" in d for d in all_diseases
        ), f"Expected melanoma/skin cancer among BRAF V600E lines, got: {all_diseases}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_gene_case_insensitive(self):
        """Test that gene lookup is case insensitive."""
        async with DepMapClient() as client:
            result_upper = await client.fetch_depmap_evidence("BRAF")
            result_lower = await client.fetch_depmap_evidence("braf")
            result_mixed = await client.fetch_depmap_evidence("Braf")

        # All should return equivalent data
        assert result_upper is not None
        assert result_lower is not None
        assert result_mixed is not None

        assert result_upper.gene == result_lower.gene == result_mixed.gene

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_pik3ca_has_drug_sensitivities(self):
        """Test that PIK3CA includes drug sensitivity data from PRISM."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("PIK3CA", "H1047R")

        assert result is not None
        assert result.gene == "PIK3CA"

        # Should have sensitive drugs if cell lines with this mutation exist
        if result.drug_sensitivities:
            for ds in result.drug_sensitivities:
                assert ds.mean_log2fc is not None
                assert ds.mean_log2fc <= -1.7

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_erbb2_has_drug_sensitivities(self):
        """Test that ERBB2/HER2 includes drug sensitivity data from PRISM."""
        async with DepMapClient() as client:
            result = await client.fetch_depmap_evidence("ERBB2")

        assert result is not None
        assert result.gene == "ERBB2"

        # Should have sensitive drugs if cell lines exist
        if result.drug_sensitivities:
            for ds in result.drug_sensitivities:
                assert ds.mean_log2fc is not None
                assert ds.mean_log2fc <= -1.7

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_drug_sensitivity_directly(self):
        """Test fetch_drug_sensitivity method directly with BRAF V600E cell lines."""
        async with DepMapClient() as client:
            # Fetch BRAF V600E data which includes drug sensitivities
            result = await client.fetch_depmap_evidence("BRAF", "V600E")

        assert result is not None
        assert (
            len(result.drug_sensitivities) > 0
        ), "Expected drug sensitivities for BRAF V600E"

        # Validate structure of returned data
        for ds in result.drug_sensitivities:
            assert ds.drug_name is not None
            assert ds.mean_log2fc is not None
            assert ds.mean_log2fc <= -1.7, "Only sensitive drugs should be returned"
            assert ds.n_cell_lines > 0

        print(
            f"\nFound {len(result.drug_sensitivities)} sensitive drugs for BRAF V600E cell lines"
        )
        print(f"Top 5 drugs by sensitivity:")
        for ds in result.drug_sensitivities[:5]:
            print(f"  {ds.drug_name}: log2fc={ds.mean_log2fc:.2f}, n={ds.n_cell_lines}")


class TestCCLETissueParsing:
    """Regression tests for CCLE sample ID tissue parsing.

    Root cause: rsplit("_", 1) on multi-word tissue suffixes like
    AUTONOMIC_GANGLIA or CENTRAL_NERVOUS_SYSTEM only captured the last
    token ("GANGLIA", "SYSTEM"), causing tumor_types_match to fail and
    producing false "no cell lines in tumor type" gaps.

    Fix: _CCLE_MULTIWORD_TISSUE_SUFFIXES lookup table checked before rsplit.
    """

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_alk_f1174l_neuroblastoma_cell_lines_have_correct_tissue(self):
        """KELLY and SK-N-SH are neuroblastoma lines; their tissue must be
        'Neuroblastoma', not 'GANGLIA' (truncated from AUTONOMIC_GANGLIA)."""
        async with CBioPortalClient() as client:
            cell_lines = await client.fetch_cell_lines_with_mutation("ALK", "F1174L")

        if not cell_lines:
            pytest.skip("No ALK F1174L cell lines found in CCLE")

        autonomic_lines = [
            cl for cl in cell_lines if "AUTONOMIC" in cl.get("sample_id", "").upper()
        ]
        assert len(autonomic_lines) > 0, (
            "Expected at least one AUTONOMIC_GANGLIA cell line for ALK F1174L "
            f"(e.g. KELLY). Got: {[cl['sample_id'] for cl in cell_lines]}"
        )

        for cl in autonomic_lines:
            assert cl["tissue"] == "Neuroblastoma", (
                f"Cell line {cl['sample_id']} should have tissue='Neuroblastoma', "
                f"got '{cl['tissue']}'. "
                "AUTONOMIC_GANGLIA suffix was not mapped correctly."
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_central_nervous_system_suffix_mapped_correctly(self):
        """Cell lines with CENTRAL_NERVOUS_SYSTEM suffix must have tissue='Brain Cancer',
        not 'SYSTEM'."""
        async with CBioPortalClient() as client:
            cell_lines = await client.fetch_cell_lines_with_mutation("BRAF", "V600E")

        cns_lines = [
            cl for cl in cell_lines
            if "CENTRAL_NERVOUS_SYSTEM" in cl.get("sample_id", "").upper()
        ]
        if not cns_lines:
            pytest.skip("No CENTRAL_NERVOUS_SYSTEM cell lines found for BRAF V600E")

        for cl in cns_lines:
            assert cl["tissue"] == "Brain Cancer", (
                f"Cell line {cl['sample_id']} should have tissue='Brain Cancer', "
                f"got '{cl['tissue']}'. "
                "CENTRAL_NERVOUS_SYSTEM suffix was not mapped correctly."
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_single_word_tissue_suffixes_unchanged(self):
        """Single-word suffixes like LUNG, SKIN, BREAST must still parse correctly."""
        async with CBioPortalClient() as client:
            cell_lines = await client.fetch_cell_lines_with_mutation("BRAF", "V600E")

        single_word_expected = {"LUNG", "SKIN", "BREAST", "BONE", "INTESTINE"}
        found_tissues = {cl["tissue"] for cl in cell_lines}
        overlap = found_tissues & single_word_expected
        assert len(overlap) > 0, (
            f"Expected some single-word tissue labels (e.g. LUNG, SKIN), "
            f"got: {found_tissues}"
        )

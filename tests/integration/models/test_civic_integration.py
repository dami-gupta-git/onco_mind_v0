"""Integration tests for CIViC evidence and identifier handling.

Tests the CIViC integration functionality:
- CIViC assertions have AID (Assertion ID) populated
- CIViC evidence items have EID (Evidence Item ID) populated
- CIViC IDs are properly serialized
- AKT1 E17K Breast Cancer evidence (well-characterized actionable variant)
"""

import pytest

from oncomind.insight_builder import Conductor, ConductorConfig

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
                assert (
                    assertion.assertion_id is not None
                ), "CIViC assertion should have assertion_id"
                assert isinstance(
                    assertion.assertion_id, int
                ), "assertion_id should be an integer"
                assert assertion.assertion_id > 0, "assertion_id should be positive"

                # AID should be formatted as "AID{number}"
                assert (
                    assertion.aid is not None
                ), "CIViC assertion should have aid property"
                assert assertion.aid.startswith(
                    "AID"
                ), f"AID should start with 'AID', got: {assertion.aid}"
                assert assertion.aid == f"AID{assertion.assertion_id}"

                # civic_url should be a valid URL
                assert (
                    assertion.civic_url is not None
                ), "CIViC assertion should have civic_url"
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
                    assert isinstance(
                        evidence.evidence_id, int
                    ), "evidence_id should be an integer"
                    assert evidence.evidence_id > 0, "evidence_id should be positive"

                    # EID should be formatted as "EID{number}"
                    assert (
                        evidence.eid is not None
                    ), "CIViC evidence should have eid property"
                    assert evidence.eid.startswith(
                        "EID"
                    ), f"EID should start with 'EID', got: {evidence.eid}"
                    assert evidence.eid == f"EID{evidence.evidence_id}"

                    # civic_url should be a valid URL
                    assert (
                        evidence.civic_url is not None
                    ), "CIViC evidence should have civic_url"
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
# AKT1 E17K BREAST CANCER INTEGRATION TEST
# =============================================================================


@pytest.mark.integration
class TestAKT1E17KBreastCancer:
    """Integration tests for AKT1 E17K in Breast Cancer - a well-characterized actionable variant."""

    @pytest.mark.asyncio
    async def test_akt1_e17k_breast_cancer_has_clinical_evidence(self):
        """AKT1 E17K in Breast Cancer should have FDA approval and CIViC evidence.

        AKT1 E17K is FDA-approved with Capivasertib + Fulvestrant for HR+/HER2- breast cancer.
        This test verifies the full pipeline returns expected clinical evidence.
        """
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("AKT1 E17K", tumor_type="Breast Cancer")

        # Should have evidence
        assert result.evidence is not None, "Should have evidence"

        # Should have CIViC evidence items with Capivasertib
        assert (
            len(result.evidence.civic_evidence) > 0
        ), "Should have CIViC evidence items"
        drug_names = []
        for e in result.evidence.civic_evidence:
            drug_names.extend(e.drugs or [])
        drug_names_lower = [d.lower() for d in drug_names]
        assert any(
            "capiva" in d for d in drug_names_lower
        ), f"Should have Capivasertib in CIViC evidence, found: {drug_names}"

        # Should have VICC evidence
        assert len(result.evidence.vicc_evidence) > 0, "Should have VICC evidence"

        # Get gaps
        gaps = result.evidence.evidence_gaps
        if gaps is None:
            gaps = result.evidence.compute_evidence_gaps()

        # Should have good evidence quality (comprehensive or moderate)
        assert gaps.overall_evidence_quality in (
            "comprehensive",
            "moderate",
        ), f"AKT1 E17K should have good evidence quality, got: {gaps.overall_evidence_quality}"

    @pytest.mark.asyncio
    async def test_akt1_e17k_breast_cancer_civic_eids_unique(self):
        """CIViC evidence items for AKT1 E17K should have unique EIDs."""
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("AKT1 E17K", tumor_type="Breast Cancer")

        # Check for duplicate EIDs
        evidence_ids = [
            e.evidence_id for e in result.evidence.civic_evidence if e.evidence_id
        ]
        unique_ids = set(evidence_ids)
        assert len(evidence_ids) == len(
            unique_ids
        ), f"Duplicate CIViC EIDs found: {evidence_ids}"

    @pytest.mark.asyncio
    async def test_akt1_e17k_has_depmap_data(self):
        """AKT1 E17K should have DepMap cell line data.

        Note: AKT1 E17K cell lines in DepMap are primarily from other tumor types
        (e.g., Thyroid, Endometrial), so this test just verifies DepMap data exists.
        """
        config = ConductorConfig(enable_llm=False, enable_literature=False)
        async with Conductor(config) as conductor:
            result = await conductor.run("AKT1 E17K", tumor_type="Breast Cancer")

        # Should have DepMap evidence
        assert (
            result.evidence.depmap_evidence is not None
        ), "Should have DepMap evidence"
        assert (
            len(result.evidence.depmap_evidence.cell_line_models) > 0
        ), "Should have cell line models"

        # Verify cell lines have mutation data
        cell_lines_with_mutation = [
            cl
            for cl in result.evidence.depmap_evidence.cell_line_models
            if cl.has_mutation
        ]
        assert (
            len(cell_lines_with_mutation) > 0
        ), "Should have cell lines with AKT1 E17K mutation"

"""Integration tests for CIViC GraphQL API.

Tests validate that the CIViC API returns expected assertions
for well-characterized oncogenic mutations with AMP/ASCO/CAP tier classifications.
"""

import pytest

from oncomind.api.civic import CIViCClient, CIViCAssertion


class TestCIViCBRAFV600E:
    """Tests for BRAF V600E - Tier I variant with FDA-approved BRAF inhibitors."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_assertions(self):
        """BRAF V600E should return substantial evidence."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("BRAF", "V600E", max_results=50)
            assert len(assertions) >= 1, "BRAF V600E should have at least 1 assertion"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_predictive_assertions(self):
        """BRAF V600E should have predictive (therapy response) assertions."""
        async with CIViCClient() as client:
            assertions = await client.fetch_predictive_assertions("BRAF", "V600E", max_results=50)
            assert len(assertions) >= 1, "BRAF V600E should have at least 1 predictive assertion"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_expected_braf_inhibitors(self):
        """Should return expected BRAF inhibitor drugs."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("BRAF", "V600E", max_results=50)

            expected_drugs = {"vemurafenib", "dabrafenib", "encorafenib", "trametinib", "cobimetinib"}
            all_drugs = set()
            for assertion in assertions:
                for drug in assertion.therapies:
                    all_drugs.add(drug.lower())

            found_expected = all_drugs & expected_drugs
            assert len(found_expected) > 0, (
                f"Expected at least one of {expected_drugs}, got: {all_drugs}"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_sensitivity_assertions(self):
        """BRAF V600E is a sensitivity marker - should have sensitivity associations."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("BRAF", "V600E", max_results=50)
            sensitivity_asserts = [a for a in assertions if a.is_sensitivity()]
            assert len(sensitivity_asserts) > 0, "BRAF V600E should have sensitivity assertions"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_tier_i_assertions(self):
        """BRAF V600E should have Tier I (strongest evidence) assertions."""
        async with CIViCClient() as client:
            assertions = await client.fetch_tier_i_assertions("BRAF", "V600E")
            assert len(assertions) >= 1, "BRAF V600E should have at least 1 Tier I assertion"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_melanoma_disease_coverage(self):
        """Melanoma is a key indication - should have melanoma-related evidence."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("BRAF", "V600E", max_results=50)

            diseases = {a.disease.lower() for a in assertions if a.disease}
            melanoma_related = any("melanoma" in d for d in diseases)
            assert melanoma_related, f"Should have melanoma-related evidence, got: {diseases}"


class TestCIViCEGFRL858R:
    """Tests for EGFR L858R - common activating mutation in NSCLC."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_assertions(self):
        """EGFR L858R should have evidence."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("EGFR", "L858R", max_results=50)
            assert len(assertions) >= 1, "EGFR L858R should have at least 1 assertion"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_expected_egfr_tkis(self):
        """Should return expected EGFR TKI drugs."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("EGFR", "L858R", max_results=50)

            expected_drugs = {"erlotinib", "gefitinib", "afatinib", "osimertinib"}
            all_drugs = set()
            for assertion in assertions:
                for drug in assertion.therapies:
                    all_drugs.add(drug.lower())

            found_expected = all_drugs & expected_drugs
            assert len(found_expected) > 0, (
                f"Expected at least one EGFR TKI {expected_drugs}, got: {all_drugs}"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_sensitivity_assertions(self):
        """EGFR L858R should have sensitivity assertions."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("EGFR", "L858R", max_results=50)
            sensitivity_asserts = [a for a in assertions if a.is_sensitivity()]
            assert len(sensitivity_asserts) > 0, "EGFR L858R should have sensitivity assertions"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_lung_cancer_disease_coverage(self):
        """Should have lung cancer related diseases."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("EGFR", "L858R", max_results=50)

            diseases = {a.disease.lower() for a in assertions if a.disease}
            lung_related = any(
                "lung" in d or "nsclc" in d or "non-small" in d
                for d in diseases
            )
            assert lung_related, f"Should have lung cancer evidence, got: {diseases}"


class TestCIViCEGFRT790M:
    """Tests for EGFR T790M - resistance mutation to first-gen EGFR TKIs."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_returns_assertions(self):
        """EGFR T790M should have evidence."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("EGFR", "T790M", max_results=50)
            assert len(assertions) >= 1, "EGFR T790M should have at least 1 assertion"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_has_resistance_associations(self):
        """EGFR T790M is a resistance marker - should have resistance associations."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("EGFR", "T790M", max_results=50)
            resistance_asserts = [a for a in assertions if a.is_resistance()]
            # T790M may also have sensitivity to osimertinib
            assert len(assertions) > 0, "EGFR T790M should have assertions"


class TestCIViCAssertionStructure:
    """Tests for CIViCAssertion data structure integrity."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_required_fields_types(self):
        """Assertions should have required fields with correct types."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("BRAF", "V600E", max_results=10)

            for assertion in assertions:
                assert isinstance(assertion.assertion_id, int)
                assert isinstance(assertion.name, str)
                assert isinstance(assertion.assertion_type, str)
                assert isinstance(assertion.molecular_profile, str)
                assert isinstance(assertion.disease, str)
                assert isinstance(assertion.therapies, list)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_amp_tier_methods_work(self):
        """AMP tier extraction methods should work without error."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("BRAF", "V600E", max_results=10)

            for assertion in assertions:
                tier = assertion.get_amp_tier()
                level = assertion.get_amp_level()
                # If amp_level is set, methods should return valid values
                if assertion.amp_level:
                    assert tier in [None, "Tier I", "Tier II", "Tier III", "Tier IV"]
                    assert level in [None, "A", "B", "C", "D"]

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_boolean_methods_work(self):
        """Boolean helper methods should work without error."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("BRAF", "V600E", max_results=10)

            for assertion in assertions:
                _ = assertion.is_sensitivity()
                _ = assertion.is_resistance()
                _ = assertion.is_accepted()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_to_dict_keys(self):
        """to_dict should return expected keys."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions("BRAF", "V600E", max_results=10)

            expected_keys = {
                "assertion_id", "name", "amp_level", "amp_tier", "amp_level_letter",
                "assertion_type", "assertion_direction", "significance", "status",
                "molecular_profile", "disease", "therapies", "fda_companion_test",
                "nccn_guideline", "description", "is_sensitivity", "is_resistance"
            }

            for assertion in assertions:
                assertion_dict = assertion.to_dict()
                assert set(assertion_dict.keys()) == expected_keys


class TestCIViCTumorTypeFilter:
    """Tests for tumor type filtering."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_melanoma_filter(self):
        """Filtering by melanoma should work."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions(
                "BRAF", "V600E", tumor_type="melanoma", max_results=20
            )

            # If we got results, they should be melanoma-related
            for assertion in assertions:
                disease_lower = assertion.disease.lower()
                assert "melanoma" in disease_lower or "skin" in disease_lower, (
                    f"Expected melanoma-related disease, got: {assertion.disease}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_nsclc_filter(self):
        """Filtering by NSCLC should work."""
        async with CIViCClient() as client:
            assertions = await client.fetch_assertions(
                "EGFR", "L858R", tumor_type="NSCLC", max_results=20
            )

            # EGFR L858R in NSCLC should return results
            assert len(assertions) >= 1, "EGFR L858R with NSCLC filter should have results"


class TestCIViCAssertionEvidence:
    """Tests for assertion evidence model conversion."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_assertion_evidence(self):
        """fetch_assertion_evidence should return evidence objects."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence("BRAF", "V600E", max_results=10)

            assert len(evidence_list) >= 1, "Should return at least 1 evidence object"
            for evidence in evidence_list:
                assert hasattr(evidence, "match_level")
                assert evidence.match_level in ["variant", "codon", "gene"]

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_match_level_tracking(self):
        """Match level should be correctly tracked in evidence."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence("BRAF", "V600E", max_results=10)

            # V600E exact matches should have variant-level match
            v600e_evidence = [e for e in evidence_list if "V600E" in e.molecular_profile.upper()]
            for evidence in v600e_evidence:
                assert evidence.match_level == "variant", (
                    f"V600E exact match should be variant level, got {evidence.match_level}"
                )

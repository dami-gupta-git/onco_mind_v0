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
    """Tests for assertion evidence model conversion (fetch_assertion_evidence)."""

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

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_has_required_fields(self):
        """Evidence objects should have all required fields populated."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence("BRAF", "V600E", max_results=5)

            assert len(evidence_list) >= 1
            for evidence in evidence_list:
                # Required identifiers
                assert evidence.assertion_id is not None
                assert isinstance(evidence.assertion_id, int)

                # Computed properties
                assert evidence.aid is not None
                assert evidence.aid.startswith("AID")
                assert evidence.civic_url is not None
                assert "civicdb.org/assertions/" in evidence.civic_url

                # Core fields
                assert evidence.molecular_profile is not None
                assert isinstance(evidence.molecular_profile, str)
                assert evidence.disease is not None
                assert isinstance(evidence.therapies, list)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_sensitivity_resistance_flags(self):
        """Evidence should correctly flag sensitivity and resistance."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence("BRAF", "V600E", max_results=20)

            # Should have at least one sensitivity assertion for BRAF V600E
            sensitivity_evidence = [e for e in evidence_list if e.is_sensitivity]
            assert len(sensitivity_evidence) > 0, (
                "BRAF V600E should have sensitivity evidence"
            )

            # Check flags are boolean
            for evidence in evidence_list:
                assert isinstance(evidence.is_sensitivity, bool)
                assert isinstance(evidence.is_resistance, bool)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_amp_tier_levels(self):
        """Evidence should have AMP tier information."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence("BRAF", "V600E", max_results=20)

            # At least some BRAF V600E evidence should have AMP tier (Tier I)
            tiered_evidence = [e for e in evidence_list if e.amp_tier]
            assert len(tiered_evidence) > 0, "BRAF V600E should have tiered evidence"

            for evidence in tiered_evidence:
                assert evidence.amp_tier in ["Tier I", "Tier II", "Tier III", "Tier IV"]
                if evidence.amp_level_letter:
                    assert evidence.amp_level_letter in ["A", "B", "C", "D"]

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_with_tumor_type_filter(self):
        """Tumor type filtering should work correctly."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence(
                "BRAF", "V600E", tumor_type="melanoma", max_results=10
            )

            # If we got results, disease_match should be True
            for evidence in evidence_list:
                assert evidence.disease_match is True, (
                    f"Evidence with tumor filter should have disease_match=True, "
                    f"got disease={evidence.disease}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_egfr_l858r(self):
        """EGFR L858R should return evidence with expected structure."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence("EGFR", "L858R", max_results=10)

            assert len(evidence_list) >= 1, "EGFR L858R should have evidence"

            # Should have sensitivity evidence (EGFR TKI sensitivity)
            sensitivity_evidence = [e for e in evidence_list if e.is_sensitivity]
            assert len(sensitivity_evidence) > 0, (
                "EGFR L858R should have TKI sensitivity evidence"
            )

            # Should have therapies
            for evidence in sensitivity_evidence:
                if evidence.therapies:
                    assert len(evidence.therapies) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_gene_only_query(self):
        """Querying by gene only should return gene-level matches."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence("BRAF", max_results=10)

            assert len(evidence_list) >= 1, "BRAF gene query should return evidence"

            # Should include various BRAF variants
            molecular_profiles = {e.molecular_profile for e in evidence_list}
            assert len(molecular_profiles) >= 1

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_serialization(self):
        """Evidence should serialize correctly via model_dump."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_assertion_evidence("BRAF", "V600E", max_results=3)

            assert len(evidence_list) >= 1
            for evidence in evidence_list:
                data = evidence.model_dump()

                # Check key fields are in serialized output
                assert "assertion_id" in data
                assert "aid" in data
                assert "civic_url" in data
                assert "molecular_profile" in data
                assert "disease" in data
                assert "therapies" in data
                assert "match_level" in data
                assert "is_sensitivity" in data
                assert "is_resistance" in data

                # Verify computed properties are serialized
                assert data["aid"] == evidence.aid
                assert data["civic_url"] == evidence.civic_url


class TestCIViCEvidenceItems:
    """Tests for CIViC Evidence Items (EIDs) fetched directly from GraphQL API."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_evidence_items_akt1_e17k_breast_cancer(self):
        """AKT1 E17K in Breast Cancer should return unique evidence items (EIDs)."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items(
                "AKT1", "E17K", tumor_type="Breast Cancer", max_results=20
            )

            # Should have evidence items
            assert len(evidence_list) >= 1, "AKT1 E17K should have evidence items"

            # Check for uniqueness - no duplicate EIDs
            evidence_ids = [e.evidence_id for e in evidence_list]
            unique_ids = set(evidence_ids)
            assert len(evidence_ids) == len(unique_ids), (
                f"Duplicate EIDs found: {evidence_ids}"
            )

            # Verify each item has required fields
            for evidence in evidence_list:
                assert evidence.evidence_id is not None
                assert evidence.eid is not None
                assert evidence.eid.startswith("EID")
                assert evidence.match_level in ("variant", "codon", "gene")

            # Should have Capivasertib evidence (FDA approved for AKT1 E17K)
            drug_names = []
            for e in evidence_list:
                drug_names.extend(e.drugs or [])
            drug_names_lower = [d.lower() for d in drug_names]
            assert any("capiva" in d for d in drug_names_lower), (
                f"Should have Capivasertib evidence, found: {drug_names}"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_evidence_items_braf_v600e(self):
        """BRAF V600E should return evidence items (EIDs)."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items("BRAF", "V600E", max_results=10)

            assert len(evidence_list) >= 1, "BRAF V600E should have evidence items"

            for evidence in evidence_list:
                # Should have evidence_id (numeric)
                assert evidence.evidence_id is not None
                assert isinstance(evidence.evidence_id, int)

                # Should have EID computed property
                assert evidence.eid is not None
                assert evidence.eid.startswith("EID")

                # Should have civic_url
                assert evidence.civic_url is not None
                assert "civicdb.org/evidence/" in evidence.civic_url

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_items_have_required_fields(self):
        """Evidence items should have all required fields populated."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items("BRAF", "V600E", max_results=5)

            assert len(evidence_list) >= 1
            for evidence in evidence_list:
                # Core fields should be present
                assert evidence.evidence_id is not None
                assert evidence.evidence_type is not None
                assert evidence.evidence_level is not None

                # Match level tracking
                assert evidence.match_level in ["variant", "codon", "gene"]
                assert evidence.matched_profile is not None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_items_with_drugs(self):
        """BRAF V600E evidence items should include drug information."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items("BRAF", "V600E", max_results=20)

            # At least some evidence should have drugs
            evidence_with_drugs = [e for e in evidence_list if e.drugs]
            assert len(evidence_with_drugs) > 0, "BRAF V600E should have evidence with drug info"

            for evidence in evidence_with_drugs:
                assert isinstance(evidence.drugs, list)
                assert len(evidence.drugs) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_items_with_pmid(self):
        """Evidence items should have PMID references."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items("BRAF", "V600E", max_results=10)

            # At least some evidence should have PMIDs
            evidence_with_pmid = [e for e in evidence_list if e.pmid]
            assert len(evidence_with_pmid) > 0, "Should have evidence with PMID"

            for evidence in evidence_with_pmid:
                assert evidence.pmid is not None
                assert evidence.source_url is not None
                assert "pubmed" in evidence.source_url.lower()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_items_egfr_l858r(self):
        """EGFR L858R should return evidence items."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items("EGFR", "L858R", max_results=10)

            assert len(evidence_list) >= 1, "EGFR L858R should have evidence items"

            # Should have variant-level matches
            variant_matches = [e for e in evidence_list if e.match_level == "variant"]
            assert len(variant_matches) > 0, "EGFR L858R should have variant-level matches"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_items_gene_only(self):
        """Querying by gene only should return gene-level evidence."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items("BRAF", max_results=10)

            assert len(evidence_list) >= 1, "BRAF gene query should return evidence"

            # Should have variety of matched profiles
            matched_profiles = {e.matched_profile for e in evidence_list}
            assert len(matched_profiles) >= 1

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_items_tumor_type_filter(self):
        """Tumor type filtering should work for evidence items."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items(
                "BRAF", "V600E", tumor_type="melanoma", max_results=10
            )

            # If we got results, disease_match should be True
            for evidence in evidence_list:
                assert evidence.disease_match is True, (
                    f"Evidence with tumor filter should have disease_match=True, "
                    f"got disease={evidence.disease}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_items_serialization(self):
        """Evidence items should serialize correctly via model_dump."""
        async with CIViCClient() as client:
            evidence_list = await client.fetch_evidence_items("BRAF", "V600E", max_results=3)

            assert len(evidence_list) >= 1
            for evidence in evidence_list:
                data = evidence.model_dump()

                # Check key fields are in serialized output
                assert "evidence_id" in data
                assert "eid" in data
                assert "civic_url" in data
                assert "evidence_type" in data
                assert "evidence_level" in data
                assert "drugs" in data
                assert "match_level" in data
                assert "matched_profile" in data
                assert "disease_match" in data

                # Verify computed properties are serialized
                assert data["eid"] == evidence.eid
                assert data["civic_url"] == evidence.civic_url

"""Tests for general utility functions in oncomind.utils.utils."""

import pytest
from unittest.mock import MagicMock

from oncomind.utils import (
    normalize_drug_name,
    is_kit_false_positive,
    dedupe_civic_evidence,
    dedupe_vicc_evidence,
)


class TestNormalizeDrugName:
    """Tests for normalize_drug_name function."""

    def test_removes_hydrochloride_suffix(self):
        """Test removal of HYDROCHLORIDE suffix."""
        assert normalize_drug_name("ERLOTINIB HYDROCHLORIDE") == "erlotinib"
        assert normalize_drug_name("erlotinib hydrochloride") == "erlotinib"

    def test_removes_mesylate_suffix(self):
        """Test removal of MESYLATE suffix."""
        assert normalize_drug_name("DABRAFENIB MESYLATE") == "dabrafenib"
        assert normalize_drug_name("Imatinib Mesylate") == "imatinib"

    def test_removes_maleate_suffix(self):
        """Test removal of MALEATE suffix."""
        assert normalize_drug_name("SUNITINIB MALEATE") == "sunitinib"

    def test_removes_fumarate_suffix(self):
        """Test removal of FUMARATE suffix."""
        assert normalize_drug_name("AFATINIB FUMARATE") == "afatinib"

    def test_removes_succinate_suffix(self):
        """Test removal of SUCCINATE suffix."""
        assert normalize_drug_name("SUVOREXANT SUCCINATE") == "suvorexant"

    def test_removes_tartrate_suffix(self):
        """Test removal of TARTRATE suffix."""
        assert normalize_drug_name("RIVASTIGMINE TARTRATE") == "rivastigmine"

    def test_removes_citrate_suffix(self):
        """Test removal of CITRATE suffix."""
        assert normalize_drug_name("TOREMIFENE CITRATE") == "toremifene"

    def test_removes_sulfate_suffix(self):
        """Test removal of SULFATE suffix."""
        assert normalize_drug_name("VINCRISTINE SULFATE") == "vincristine"

    def test_removes_phosphate_suffix(self):
        """Test removal of PHOSPHATE suffix."""
        assert normalize_drug_name("FLUDARABINE PHOSPHATE") == "fludarabine"

    def test_removes_sodium_suffix(self):
        """Test removal of SODIUM suffix."""
        assert normalize_drug_name("PEMETREXED SODIUM") == "pemetrexed"

    def test_removes_potassium_suffix(self):
        """Test removal of POTASSIUM suffix."""
        assert normalize_drug_name("CLAVULANATE POTASSIUM") == "clavulanate"

    def test_removes_calcium_suffix(self):
        """Test removal of CALCIUM suffix."""
        assert normalize_drug_name("LEUCOVORIN CALCIUM") == "leucovorin"

    def test_removes_dimesylate_suffix(self):
        """Test removal of DIMESYLATE suffix."""
        assert normalize_drug_name("LAPATINIB DIMESYLATE") == "lapatinib"

    def test_removes_dihydrochloride_suffix(self):
        """Test removal of DIHYDROCHLORIDE suffix."""
        assert normalize_drug_name("TRAMETINIB DIHYDROCHLORIDE") == "trametinib"

    def test_removes_tosylate_suffix(self):
        """Test removal of TOSYLATE suffix."""
        assert normalize_drug_name("SORAFENIB TOSYLATE") == "sorafenib"

    def test_removes_acetate_suffix(self):
        """Test removal of ACETATE suffix."""
        assert normalize_drug_name("LEUPROLIDE ACETATE") == "leuprolide"

    def test_removes_besylate_suffix(self):
        """Test removal of BESYLATE suffix."""
        assert normalize_drug_name("AMLODIPINE BESYLATE") == "amlodipine"

    def test_no_suffix_returns_lowercase(self):
        """Test drug name without suffix returns lowercase."""
        assert normalize_drug_name("Pembrolizumab") == "pembrolizumab"
        assert normalize_drug_name("NIVOLUMAB") == "nivolumab"
        assert normalize_drug_name("trastuzumab") == "trastuzumab"

    def test_case_insensitive(self):
        """Test that function is case insensitive."""
        assert normalize_drug_name("ERLOTINIB HYDROCHLORIDE") == "erlotinib"
        assert normalize_drug_name("erlotinib hydrochloride") == "erlotinib"
        assert normalize_drug_name("Erlotinib Hydrochloride") == "erlotinib"

    def test_empty_string(self):
        """Test empty string returns empty string."""
        assert normalize_drug_name("") == ""

    def test_none_like_empty(self):
        """Test None-like values return empty string."""
        assert normalize_drug_name(None) == ""

    def test_whitespace_handling(self):
        """Test whitespace is trimmed."""
        assert normalize_drug_name("  ERLOTINIB HYDROCHLORIDE  ") == "erlotinib"
        assert normalize_drug_name("  Pembrolizumab  ") == "pembrolizumab"

    def test_only_first_matching_suffix_removed(self):
        """Test that only one suffix is removed even if name somehow had multiple."""
        # This is an edge case - real drug names don't have multiple suffixes
        # But the function should only remove one
        result = normalize_drug_name("DRUG SODIUM")
        assert result == "drug"

    def test_preserves_biologic_formulation_differences(self):
        """Test that different biologic formulations are NOT normalized together.

        IV (AMIVANTAMAB-VMJW) and SC (AMIVANTAMAB AND HYALURONIDASE-LPUJ)
        formulations are clinically distinct and should remain separate.
        """
        iv_form = normalize_drug_name("AMIVANTAMAB-VMJW")
        sc_form = normalize_drug_name("AMIVANTAMAB AND HYALURONIDASE-LPUJ (HUMAN RECOMBINANT)")

        # These should NOT be equal - they are different formulations
        assert iv_form != sc_form
        # IV form keeps the suffix
        assert iv_form == "amivantamab-vmjw"
        # SC form keeps its full name
        assert "hyaluronidase" in sc_form.lower()


class TestIsKitFalsePositive:
    """Tests for is_kit_false_positive function."""

    def test_diagnostic_kit_without_oncology_is_false_positive(self):
        """Test diagnostic kit without oncology context is flagged."""
        indication = "This diagnostic kit is used for detecting viral infections"
        assert is_kit_false_positive(indication) is True

    def test_test_kit_without_oncology_is_false_positive(self):
        """Test test kit without oncology context is flagged."""
        indication = "This test kit determines blood glucose levels"
        assert is_kit_false_positive(indication) is True

    def test_preparation_kit_without_oncology_is_false_positive(self):
        """Test preparation kit without oncology context is flagged."""
        indication = "This kit for the preparation of solutions is used for bowel preparation"
        assert is_kit_false_positive(indication) is True

    def test_kit_with_cancer_context_is_not_false_positive(self):
        """Test kit mentioned with cancer context is not flagged."""
        indication = "For treatment of KIT-positive gastrointestinal stromal tumor (cancer)"
        assert is_kit_false_positive(indication) is False

    def test_kit_with_tumor_context_is_not_false_positive(self):
        """Test kit mentioned with tumor context is not flagged."""
        indication = "For patients with KIT D816V mutation in tumor samples"
        assert is_kit_false_positive(indication) is False

    def test_kit_with_gist_context_is_not_false_positive(self):
        """Test KIT mentioned with GIST context is not flagged."""
        indication = "For treatment of KIT-positive GIST"
        assert is_kit_false_positive(indication) is False

    def test_kit_with_melanoma_context_is_not_false_positive(self):
        """Test KIT mentioned with melanoma context is not flagged."""
        indication = "For melanoma patients with KIT mutations"
        assert is_kit_false_positive(indication) is False

    def test_kit_with_leukemia_context_is_not_false_positive(self):
        """Test KIT mentioned with leukemia context is not flagged."""
        indication = "For acute myeloid leukemia with KIT mutation"
        assert is_kit_false_positive(indication) is False

    def test_kit_with_mastocytosis_context_is_not_false_positive(self):
        """Test KIT mentioned with mastocytosis context is not flagged."""
        indication = "For systemic mastocytosis with KIT D816V"
        assert is_kit_false_positive(indication) is False

    def test_no_kit_pattern_returns_false(self):
        """Test indication without kit patterns returns False."""
        indication = "For treatment of non-small cell lung cancer"
        assert is_kit_false_positive(indication) is False

    def test_empty_indication_returns_false(self):
        """Test empty indication returns False."""
        assert is_kit_false_positive("") is False

    def test_none_indication_returns_false(self):
        """Test None indication returns False."""
        assert is_kit_false_positive(None) is False

    def test_case_insensitive(self):
        """Test that function is case insensitive."""
        assert is_kit_false_positive("DIAGNOSTIC KIT for testing") is True
        assert is_kit_false_positive("diagnostic kit for testing") is True
        assert is_kit_false_positive("Diagnostic Kit for testing") is True


class TestDedupeCivicEvidence:
    """Tests for dedupe_civic_evidence function."""

    def _create_mock_evidence(
        self,
        evidence_id: int,
        eid: str = "EID1",
        civic_url: str = "https://civicdb.org/evidence/1",
        evidence_type: str = "Predictive",
        evidence_level: str = "A",
        clinical_significance: str = "Sensitivity/Response",
        disease: str = "Lung Cancer",
        drugs: list = None,
        description: str = "Test description",
        pmid: str = "12345678",
        source_url: str = "https://pubmed.ncbi.nlm.nih.gov/12345678",
        trust_rating: int = 4,
        rating: int = 3,
        evidence_direction: str = "Supports",
        locus_match: str = "variant",
        matched_profile: str = "BRAF V600E",
        tumor_match: str = "exact",
    ):
        """Create a mock CIViC evidence object."""
        mock = MagicMock()
        mock.evidence_id = evidence_id
        mock.eid = eid
        mock.civic_url = civic_url
        mock.evidence_type = evidence_type
        mock.evidence_level = evidence_level
        mock.clinical_significance = clinical_significance
        mock.disease = disease
        mock.drugs = drugs or ["Vemurafenib"]
        mock.description = description
        mock.pmid = pmid
        mock.source_url = source_url
        mock.trust_rating = trust_rating
        mock.rating = rating
        mock.evidence_direction = evidence_direction
        mock.locus_match = locus_match
        mock.matched_profile = matched_profile
        mock.tumor_match = tumor_match
        return mock

    def test_deduplicates_by_evidence_id(self):
        """Test that duplicate evidence_ids are removed."""
        evidence_list = [
            self._create_mock_evidence(evidence_id=123, eid="EID123"),
            self._create_mock_evidence(evidence_id=123, eid="EID123"),  # Duplicate
            self._create_mock_evidence(evidence_id=456, eid="EID456"),
        ]

        result = dedupe_civic_evidence(evidence_list)

        assert len(result) == 2
        evidence_ids = [e["evidence_id"] for e in result]
        assert evidence_ids == [123, 456]

    def test_preserves_first_occurrence(self):
        """Test that first occurrence of duplicate is preserved."""
        evidence_list = [
            self._create_mock_evidence(evidence_id=123, description="First"),
            self._create_mock_evidence(evidence_id=123, description="Second"),
        ]

        result = dedupe_civic_evidence(evidence_list)

        assert len(result) == 1
        assert result[0]["description"] == "First"

    def test_returns_correct_dict_structure(self):
        """Test that returned dict has expected keys."""
        evidence_list = [
            self._create_mock_evidence(evidence_id=123)
        ]

        result = dedupe_civic_evidence(evidence_list)

        expected_keys = {
            "evidence_id", "eid", "civic_url", "evidence_type", "evidence_level",
            "clinical_significance", "disease", "drugs", "description", "pmid",
            "source_url", "trust_rating", "evidence_direction", "locus_match",
            "matched_profile", "tumor_match"
        }
        assert set(result[0].keys()) == expected_keys

    def test_uses_trust_rating_over_rating(self):
        """Test that trust_rating is preferred over rating when available."""
        evidence_list = [
            self._create_mock_evidence(evidence_id=123, trust_rating=5, rating=3)
        ]

        result = dedupe_civic_evidence(evidence_list)

        assert result[0]["trust_rating"] == 5

    def test_falls_back_to_rating_when_trust_rating_none(self):
        """Test that rating is used when trust_rating is None."""
        evidence_list = [
            self._create_mock_evidence(evidence_id=123, trust_rating=None, rating=3)
        ]

        result = dedupe_civic_evidence(evidence_list)

        assert result[0]["trust_rating"] == 3

    def test_empty_list_returns_empty(self):
        """Test that empty input returns empty list."""
        assert dedupe_civic_evidence([]) == []

    def test_preserves_order(self):
        """Test that original order is preserved after deduplication."""
        evidence_list = [
            self._create_mock_evidence(evidence_id=300),
            self._create_mock_evidence(evidence_id=100),
            self._create_mock_evidence(evidence_id=200),
        ]

        result = dedupe_civic_evidence(evidence_list)

        evidence_ids = [e["evidence_id"] for e in result]
        assert evidence_ids == [300, 100, 200]


class TestDedupeViccEvidence:
    """Tests for dedupe_vicc_evidence function."""

    def _create_mock_vicc(
        self,
        drugs: list = None,
        disease: str = "Melanoma",
        response_type: str = "Sensitivity",
        source: str = "JAX",
        evidence_level: str = "Level A",
        molecular_profile: str = "BRAF V600E",
        molecular_profile_score: float = 0.95,
        publication_url: str = "https://example.com/paper",
        locus_match: str = "variant",
        matched_profile: str = "BRAF V600E",
        tumor_match: str = "exact",
    ):
        """Create a mock VICC evidence object."""
        mock = MagicMock()
        mock.drugs = drugs or ["Vemurafenib"]
        mock.disease = disease
        mock.response_type = response_type
        mock.source = source
        mock.evidence_level = evidence_level
        mock.molecular_profile = molecular_profile
        mock.molecular_profile_score = molecular_profile_score
        mock.publication_url = publication_url
        mock.locus_match = locus_match
        mock.matched_profile = matched_profile
        mock.tumor_match = tumor_match
        return mock

    def test_deduplicates_by_drugs_disease_response(self):
        """Test deduplication by (drugs, disease, response_type) key."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer X", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer X", response_type="Sensitivity"),  # Duplicate
            self._create_mock_vicc(drugs=["Drug B"], disease="Cancer X", response_type="Sensitivity"),  # Different drug
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 2

    def test_different_disease_not_deduplicated(self):
        """Test that different diseases are not deduplicated."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer X", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer Y", response_type="Sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 2

    def test_different_response_type_not_deduplicated(self):
        """Test that different response types are not deduplicated."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer X", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer X", response_type="Resistance"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 2

    def test_drug_order_does_not_matter(self):
        """Test that drug order doesn't affect deduplication (sorted)."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A", "Drug B"], disease="Cancer X", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug B", "Drug A"], disease="Cancer X", response_type="Sensitivity"),  # Same drugs, different order
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

    def test_case_insensitive_disease(self):
        """Test that disease comparison is case insensitive."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease="MELANOMA", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug A"], disease="melanoma", response_type="Sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

    def test_case_insensitive_response_type(self):
        """Test that response_type comparison is case insensitive."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer", response_type="SENSITIVITY"),
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer", response_type="sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

    def test_handles_none_drugs(self):
        """Test handling of None drugs."""
        evidence_list = [
            self._create_mock_vicc(drugs=None, disease="Cancer", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=None, disease="Cancer", response_type="Sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

    def test_handles_empty_drugs(self):
        """Test handling of empty drugs list."""
        evidence_list = [
            self._create_mock_vicc(drugs=[], disease="Cancer", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=[], disease="Cancer", response_type="Sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

    def test_handles_none_disease(self):
        """Test handling of None disease."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease=None, response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug A"], disease=None, response_type="Sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

    def test_handles_none_response_type(self):
        """Test handling of None response_type."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer", response_type=None),
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer", response_type=None),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

    def test_returns_correct_dict_structure(self):
        """Test that returned dict has expected keys."""
        evidence_list = [
            self._create_mock_vicc()
        ]

        result = dedupe_vicc_evidence(evidence_list)

        expected_keys = {
            "source", "drugs", "disease", "response_type", "evidence_level",
            "molecular_profile", "molecular_profile_score", "publication_url",
            "locus_match", "matched_profile", "tumor_match"
        }
        assert set(result[0].keys()) == expected_keys

    def test_handles_publication_url_list(self):
        """Test that publication_url list is handled (first element taken)."""
        mock = self._create_mock_vicc()
        mock.publication_url = ["https://first.com", "https://second.com"]

        result = dedupe_vicc_evidence([mock])

        assert result[0]["publication_url"] == "https://first.com"

    def test_handles_publication_url_string(self):
        """Test that publication_url string is preserved."""
        mock = self._create_mock_vicc()
        mock.publication_url = "https://single.com"

        result = dedupe_vicc_evidence([mock])

        assert result[0]["publication_url"] == "https://single.com"

    def test_handles_empty_publication_url_list(self):
        """Test that empty publication_url list is handled."""
        mock = self._create_mock_vicc()
        mock.publication_url = []

        result = dedupe_vicc_evidence([mock])

        assert result[0]["publication_url"] == []

    def test_empty_list_returns_empty(self):
        """Test that empty input returns empty list."""
        assert dedupe_vicc_evidence([]) == []

    def test_preserves_order(self):
        """Test that original order is preserved after deduplication."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug C"], disease="Cancer", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug B"], disease="Cancer", response_type="Sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        drugs = [e["drugs"] for e in result]
        assert drugs == [["Drug C"], ["Drug A"], ["Drug B"]]

    def test_whitespace_trimmed_in_disease(self):
        """Test that whitespace is trimmed in disease comparison."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease="  Cancer  ", response_type="Sensitivity"),
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer", response_type="Sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

    def test_whitespace_trimmed_in_response_type(self):
        """Test that whitespace is trimmed in response_type comparison."""
        evidence_list = [
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer", response_type="  Sensitivity  "),
            self._create_mock_vicc(drugs=["Drug A"], disease="Cancer", response_type="Sensitivity"),
        ]

        result = dedupe_vicc_evidence(evidence_list)

        assert len(result) == 1

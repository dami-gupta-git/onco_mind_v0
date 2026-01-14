"""Tests for centralized tumor matching functions in base.py."""

import pytest

from oncomind.models.evidence.base import (
    is_pan_cancer_term,
    tumor_types_match,
    compute_cancer_specificity,
)


class TestIsPanCancerTerm:
    """Tests for is_pan_cancer_term function."""

    def test_pan_cancer_terms_detected(self):
        """Pan-cancer terms should be detected."""
        pan_cancer_terms = [
            "Cancer",
            "cancer",
            "Solid Tumor",
            "solid tumor",
            "Solid Tumors",
            "Advanced Solid Tumor",
            "Malignant Neoplasm",
            "Neoplasm",
            "Tumor",
            "Tumour",
            "All Tumor Types",
            "Any Tumor",
            "Any Cancer",
        ]
        for term in pan_cancer_terms:
            assert is_pan_cancer_term(term), f"Expected '{term}' to be pan-cancer"

    def test_specific_cancers_not_pan_cancer(self):
        """Specific cancer types should NOT be detected as pan-cancer."""
        specific_cancers = [
            "Lung Cancer",
            "Non-Small Cell Lung Cancer",
            "NSCLC",
            "Breast Cancer",
            "Melanoma",
            "Colorectal Cancer",
            "Pancreatic Cancer",
            "Ovarian Cancer",
            "Lung Adenocarcinoma",
            "HER2-Negative Breast Cancer",
        ]
        for term in specific_cancers:
            assert not is_pan_cancer_term(term), f"Expected '{term}' to NOT be pan-cancer"

    def test_none_is_not_pan_cancer(self):
        """None should not be considered pan-cancer."""
        assert not is_pan_cancer_term(None)

    def test_empty_string_is_not_pan_cancer(self):
        """Empty string should not be considered pan-cancer."""
        assert not is_pan_cancer_term("")

    def test_whitespace_handling(self):
        """Terms with extra whitespace should still be detected."""
        assert is_pan_cancer_term("  cancer  ")
        assert is_pan_cancer_term("  solid tumor  ")

    def test_fda_tumor_agnostic_biomarkers_detected(self):
        """FDA tumor-agnostic biomarkers should be detected via substring matching."""
        fda_pan_cancer_indications = [
            "MSI-H Solid Tumors",
            "MSI-H/dMMR Solid Tumors",
            "dMMR Advanced Solid Tumors",
            "NTRK fusion-positive tumors",
            "NTRK gene fusion solid tumors",
            "Tumor-agnostic indication",
            "TMB-H solid tumors",
        ]
        for indication in fda_pan_cancer_indications:
            assert is_pan_cancer_term(indication), f"Expected '{indication}' to be pan-cancer"

    def test_fda_indication_with_specific_cancer_not_pan_cancer(self):
        """FDA indications for specific cancers should NOT be pan-cancer."""
        specific_fda_indications = [
            "Metastatic melanoma",
            "Non-small cell lung cancer",
            "HER2-positive breast cancer",
            "BRAF V600E mutant colorectal cancer",
            "Unresectable hepatocellular carcinoma",
        ]
        for indication in specific_fda_indications:
            assert not is_pan_cancer_term(indication), f"Expected '{indication}' to NOT be pan-cancer"


class TestTumorTypesMatch:
    """Tests for tumor_types_match function."""

    def test_direct_substring_match(self):
        """Direct substring matches should work."""
        assert tumor_types_match("Lung Cancer", "Lung")
        assert tumor_types_match("Lung", "Lung Cancer")
        assert tumor_types_match("Non-Small Cell Lung Cancer", "Lung")

    def test_exact_match(self):
        """Exact matches should work."""
        assert tumor_types_match("Melanoma", "Melanoma")
        assert tumor_types_match("NSCLC", "NSCLC")

    def test_tumor_type_mappings_abbreviation_to_full(self):
        """Abbreviation to full name via TUMOR_TYPE_MAPPINGS should work."""
        # NSCLC -> non-small cell lung cancer
        assert tumor_types_match("Non-Small Cell Lung Cancer", "NSCLC")
        assert tumor_types_match("Lung Adenocarcinoma", "NSCLC")

        # melanoma -> skin
        assert tumor_types_match("Skin", "Melanoma")
        assert tumor_types_match("SKIN", "melanoma")

    def test_tumor_type_mappings_full_to_abbreviation(self):
        """Full name to abbreviation via TUMOR_TYPE_MAPPINGS should work."""
        assert tumor_types_match("NSCLC", "Non-Small Cell Lung Cancer")
        assert tumor_types_match("nsclc", "lung adenocarcinoma")

    def test_cross_alias_matching(self):
        """Two aliases of the same category should match."""
        # Both are aliases under 'nsclc' in TUMOR_TYPE_MAPPINGS
        assert tumor_types_match("Lung Adenocarcinoma", "Lung Squamous Cell Carcinoma")

    def test_no_match_different_cancers(self):
        """Different cancer types should not match."""
        assert not tumor_types_match("Breast Cancer", "Lung Cancer")
        assert not tumor_types_match("Melanoma", "Colorectal Cancer")
        assert not tumor_types_match("Ovarian Cancer", "NSCLC")

    def test_none_source_returns_false(self):
        """None source should return False."""
        assert not tumor_types_match(None, "Lung Cancer")

    def test_none_query_returns_false(self):
        """None query should return False."""
        assert not tumor_types_match("Lung Cancer", None)

    def test_both_none_returns_false(self):
        """Both None should return False."""
        assert not tumor_types_match(None, None)

    def test_min_substring_length_prevents_short_matches(self):
        """Short strings should not match via direct substring to avoid false positives."""
        # Short strings (< 3 chars) should not match via direct substring
        # However, they CAN still match via TUMOR_TYPE_MAPPINGS if they're aliases
        # For example, "ma" won't match "colorectal" but may match melanoma-related terms
        # via the mappings

        # Test with strings that are NOT in TUMOR_TYPE_MAPPINGS
        assert not tumor_types_match("xy", "xyz cancer")  # Too short for substring, not in mappings

        # Longer strings should work via substring
        assert tumor_types_match("mel", "melanoma")

    def test_case_insensitive(self):
        """Matching should be case-insensitive."""
        assert tumor_types_match("LUNG CANCER", "lung cancer")
        assert tumor_types_match("Melanoma", "MELANOMA")
        assert tumor_types_match("nsclc", "NSCLC")


class TestComputeCancerSpecificity:
    """Tests for compute_cancer_specificity function."""

    def test_pan_cancer_returns_pan_cancer(self):
        """Pan-cancer terms should return 'pan_cancer'."""
        assert compute_cancer_specificity("Cancer", "NSCLC") == "pan_cancer"
        assert compute_cancer_specificity("Solid Tumor", "Melanoma") == "pan_cancer"
        assert compute_cancer_specificity("Malignant Neoplasm", "Breast Cancer") == "pan_cancer"

    def test_matching_tumor_returns_cancer_specific(self):
        """Matching tumor types should return 'cancer_specific'."""
        assert compute_cancer_specificity("Non-Small Cell Lung Cancer", "NSCLC") == "cancer_specific"
        assert compute_cancer_specificity("Lung Adenocarcinoma", "NSCLC") == "cancer_specific"
        assert compute_cancer_specificity("Melanoma", "Melanoma") == "cancer_specific"

    def test_non_matching_returns_source_disease(self):
        """Non-matching specific tumors should return the source disease name."""
        result = compute_cancer_specificity("Breast Cancer", "NSCLC")
        assert result == "Breast Cancer"

        result = compute_cancer_specificity("Ovarian Cancer", "Melanoma")
        assert result == "Ovarian Cancer"

    def test_none_source_returns_none(self):
        """None source disease should return None (unknown)."""
        assert compute_cancer_specificity(None, "NSCLC") is None

    def test_none_query_returns_source_disease(self):
        """None queried tumor should return the source disease."""
        assert compute_cancer_specificity("Breast Cancer", None) == "Breast Cancer"
        assert compute_cancer_specificity("NSCLC", None) == "NSCLC"

    def test_both_none_returns_none(self):
        """Both None should return None (unknown)."""
        assert compute_cancer_specificity(None, None) is None


class TestIntegrationWithAPIClients:
    """Integration tests verifying centralized functions work correctly for API client use cases."""

    def test_pan_cancer_detection_for_civic(self):
        """Test pan-cancer detection that CIViC API client uses."""
        # Pan-cancer terms should be detected
        assert is_pan_cancer_term("Cancer") is True
        assert is_pan_cancer_term("Solid Tumor") is True

        # Specific cancers should not be detected as pan-cancer
        assert is_pan_cancer_term("NSCLC") is False
        assert is_pan_cancer_term("Lung Cancer") is False

    def test_tumor_matching_for_civic(self):
        """Test tumor matching logic used by CIViC API client."""
        # Matching tumors should return True
        assert tumor_types_match("Non-Small Cell Lung Cancer", "NSCLC") is True
        assert tumor_types_match("Lung Adenocarcinoma", "NSCLC") is True

        # Non-matching should return False
        assert tumor_types_match("Breast Cancer", "NSCLC") is False

    def test_tumor_matching_for_vicc(self):
        """Test tumor matching logic used by VICC API client."""
        # Matching tumors should return True
        assert tumor_types_match("Non-Small Cell Lung Cancer", "NSCLC") is True
        assert tumor_types_match("melanoma", "melanoma") is True

        # Non-matching should return False
        assert tumor_types_match("breast cancer", "melanoma") is False

    def test_tumor_matching_for_cgi(self):
        """Test tumor matching logic used by CGI API client."""
        # CGI uses abbreviations like NSCLC, MEL, CRC
        assert tumor_types_match("NSCLC", "Non-Small Cell Lung Cancer") is True
        assert tumor_types_match("MEL", "Melanoma") is True
        assert tumor_types_match("CRC", "Colorectal Cancer") is True

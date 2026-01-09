"""Integration tests for FDA Label Service.

These tests verify the FDA label service correctly:
1. Collects drugs from multiple sources (CGI, CIViC, VICC, FDA biomarkers)
2. Looks up labels in fda_labels.json cache
3. Fetches missing drugs from OpenFDA and updates cache

Run with: pytest tests/integration/api/test_fda_label_service.py -v -s

Use -s to see print output for debugging.
"""

import pytest
from pathlib import Path

from oncomind.api.fda_label_service import (
    collect_all_drugs,
    get_drugs_from_fda_biomarkers,
    get_fda_labels_for_drugs,
    FDALabelData,
)


class TestFDALabelServiceDrugCollection:
    """Test drug collection from various sources."""

    def test_collect_drugs_from_fda_biomarkers_akt1(self):
        """Test collecting drugs for AKT1 from FDA biomarkers file."""
        drugs = get_drugs_from_fda_biomarkers("AKT1")

        print(f"\nDrugs from FDA biomarkers for AKT1: {drugs}")

        # AKT1 should have capivasertib (Truqap)
        # The exact drugs depend on the fda_oncology_biomarkers.xlsx content
        assert isinstance(drugs, set)
        # May be empty if fda_oncology_biomarkers.xlsx doesn't exist or doesn't have AKT1

    def test_collect_drugs_from_fda_biomarkers_braf(self):
        """Test collecting drugs for BRAF from FDA biomarkers file."""
        drugs = get_drugs_from_fda_biomarkers("BRAF")

        print(f"\nDrugs from FDA biomarkers for BRAF: {drugs}")

        # BRAF should have vemurafenib, dabrafenib, etc.
        assert isinstance(drugs, set)

    def test_collect_all_drugs_akt1_e17k(self):
        """Test collecting all drugs for AKT1 E17K from all sources."""
        # Simulate drug lists from CGI, CIViC, VICC
        cgi_drugs = ["Capivasertib", "MK-2206"]
        civic_drugs = ["AZD5363", "MK-2206"]  # AZD5363 is capivasertib
        vicc_drugs = ["Capivasertib"]

        all_drugs = collect_all_drugs(
            gene="AKT1",
            cgi_drugs=cgi_drugs,
            civic_drugs=civic_drugs,
            vicc_drugs=vicc_drugs,
        )

        print(f"\nAll drugs collected for AKT1: {all_drugs}")

        # Should have at least the drugs we passed in
        assert "Capivasertib" in all_drugs
        assert "MK-2206" in all_drugs
        assert "AZD5363" in all_drugs

        # Should also have any drugs from fda_oncology_biomarkers.xlsx
        assert len(all_drugs) >= 3


class TestFDALabelServiceIntegration:
    """Integration tests for the full FDA label service flow."""

    def test_get_fda_labels_for_akt1_drugs(self):
        """Test getting FDA labels for AKT1-associated drugs."""
        drugs = {"capivasertib", "MK-2206", "AZD5363"}

        labels = get_fda_labels_for_drugs(drugs, gene="AKT1", fetch_missing=False)

        print(f"\nFDA labels for AKT1 drugs (cache only, no fetch):")
        for label in labels:
            print(f"  - {label.drug}: {label.brand_name or 'no brand'}")

        # Should return FDALabelData objects
        assert isinstance(labels, list)
        for label in labels:
            assert isinstance(label, FDALabelData)

    @pytest.mark.skip(reason="Enable to test OpenFDA fetching - may be slow")
    def test_get_fda_labels_with_fetch(self):
        """Test getting FDA labels with OpenFDA fetch for missing drugs."""
        # Use a drug that's likely not in the cache
        drugs = {"osimertinib"}

        labels = get_fda_labels_for_drugs(drugs, gene="EGFR", fetch_missing=True)

        print(f"\nFDA labels with fetch enabled:")
        for label in labels:
            print(f"  - {label.drug}: {label.brand_name or 'no brand'}")
            print(f"    indications: {label.indications_and_usage[:150] if label.indications_and_usage else 'None'}...")

        # Should have fetched the drug
        assert len(labels) > 0


class TestFDALabelServiceAKT1E17KBreastCancer:
    """Specific test for AKT1 E17K in Breast Cancer (as requested by user)."""

    def test_akt1_e17k_breast_cancer_full_flow(self):
        """Test the full FDA label service flow for AKT1 E17K in Breast Cancer.

        This simulates what happens when the evidence_aggregator processes
        AKT1 E17K in Breast Cancer and calls _ensure_fda_labels_for_evidence.
        """
        # Simulate drugs that would be collected from CGI/CIViC/VICC for AKT1 E17K
        # These are actual drugs associated with AKT1 E17K in breast cancer
        cgi_drugs = ["Capivasertib"]  # FDA-approved for PIK3CA/AKT1/PTEN-altered HR+ breast cancer
        civic_drugs = ["MK-2206", "AZD5363", "Ipatasertib"]  # AKT inhibitors in trials
        vicc_drugs = ["Capivasertib"]

        # Collect all drugs (also reads from fda_oncology_biomarkers.xlsx)
        all_drugs = collect_all_drugs(
            gene="AKT1",
            cgi_drugs=cgi_drugs,
            civic_drugs=civic_drugs,
            vicc_drugs=vicc_drugs,
        )

        print(f"\n=== AKT1 E17K Breast Cancer FDA Label Test ===")
        print(f"Collected {len(all_drugs)} drugs for AKT1:")
        for drug in sorted(all_drugs):
            print(f"  - {drug}")

        # Get FDA labels (cache lookup only, no API fetch to keep test fast)
        labels = get_fda_labels_for_drugs(all_drugs, gene="AKT1", fetch_missing=False)

        print(f"\nFound {len(labels)} FDA labels in cache:")
        for label in labels:
            brand = label.brand_name or "no brand"
            print(f"  - {label.drug} ({brand})")

        # Assertions
        assert len(all_drugs) >= 3, f"Expected at least 3 drugs, got {len(all_drugs)}"

        # Capivasertib should be in the drug list (FDA-approved for AKT1-altered breast cancer)
        drug_names_lower = {d.lower() for d in all_drugs}
        assert "capivasertib" in drug_names_lower, \
            f"Capivasertib should be in drugs for AKT1 E17K breast cancer"

        # If we have cached labels, verify they have expected fields
        for label in labels:
            assert isinstance(label, FDALabelData)
            assert label.drug is not None
            assert label.gene is not None

    def test_akt1_e17k_capivasertib_indication(self):
        """Test that capivasertib indication mentions AKT1 and breast cancer."""
        # Get labels for capivasertib specifically
        labels = get_fda_labels_for_drugs({"capivasertib"}, gene="AKT1", fetch_missing=False)

        if not labels:
            pytest.skip("Capivasertib not in cache - run prefetch or enable fetch_missing")

        label = labels[0]
        print(f"\n=== Capivasertib Label Analysis ===")
        print(f"Drug: {label.drug}")
        print(f"Brand: {label.brand_name}")
        print(f"Gene: {label.gene}")

        if label.indications_and_usage:
            print(f"\nIndications (first 500 chars):")
            print(f"  {label.indications_and_usage[:500]}...")

            # The indication should mention AKT1 and breast cancer
            indication_lower = label.indications_and_usage.lower()

            # Check for key terms (Truqap/capivasertib indication mentions AKT1)
            has_akt1 = "akt1" in indication_lower or "akt" in indication_lower
            has_breast = "breast" in indication_lower
            has_pik3ca = "pik3ca" in indication_lower
            has_pten = "pten" in indication_lower

            print(f"\nIndication analysis:")
            print(f"  Contains AKT1: {has_akt1}")
            print(f"  Contains breast: {has_breast}")
            print(f"  Contains PIK3CA: {has_pik3ca}")
            print(f"  Contains PTEN: {has_pten}")

            # Truqap is approved for PIK3CA/AKT1/PTEN-altered breast cancer
            assert has_breast, "Indication should mention breast cancer"
            assert has_akt1 or has_pik3ca or has_pten, \
                "Indication should mention AKT1, PIK3CA, or PTEN"
        else:
            print("\nNo indications text available in cache")


class TestFDALabelServiceEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_drug_set(self):
        """Test handling of empty drug set."""
        labels = get_fda_labels_for_drugs(set(), gene="BRAF", fetch_missing=False)

        assert labels == []

    def test_unknown_drug(self):
        """Test handling of unknown drugs."""
        labels = get_fda_labels_for_drugs(
            {"FAKE_DRUG_XYZ123"},
            gene="BRAF",
            fetch_missing=False,  # Don't hit API for fake drug
        )

        # Should return empty since drug not in cache
        assert labels == []

    def test_collect_drugs_unknown_gene(self):
        """Test collecting drugs for unknown gene."""
        drugs = get_drugs_from_fda_biomarkers("FAKE_GENE_XYZ")

        # Should return empty set, not error
        assert drugs == set()

    def test_collect_all_drugs_with_empty_lists(self):
        """Test collect_all_drugs with empty/None lists."""
        all_drugs = collect_all_drugs(
            gene="BRAF",
            cgi_drugs=[],
            civic_drugs=None,
            vicc_drugs=[],
        )

        # Should still return drugs from fda_oncology_biomarkers.xlsx if it has BRAF
        assert isinstance(all_drugs, set)

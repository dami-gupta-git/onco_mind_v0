"""Unit tests for biomarker specificity parsing."""

import pytest
from oncomind.insight_builder.fda_processor import (
    parse_biomarker_specificity,
    is_variant_covered,
    match_fda_approval,
)

class TestParseBiomarkerSpecificity:
    """Tests for parse_biomarker_specificity function."""

    def test_variant_level_kras_g12c(self):
        """Test parsing variant-level specificity like KRAS G12C."""
        text = "indicated for KRAS G12C mutated non-small cell lung cancer"
        result = parse_biomarker_specificity(text, "KRAS")

        assert result is not None
        assert result["level"] == "variant"
        assert result["codon"] == "G12"
        assert result["specified_variant"] == "G12C"

    def test_truqap_true(self):
        indication = """TRUQAP, in combination with fulvestrant, is indicated for 
        the treatment of adult patients with hormone receptor (HR)-positive, 
        HER2-negative breast cancer with one or more PIK3CA/AKT1/PTEN -alteration"""

        result = match_fda_approval("AKT1", "E17K", "breast cancer", indication)
        assert result["matched"]
        assert result["match_level"]=='gene'
        assert result["combination_partners"] ==['fulvestrant']

    def test_truqap_fail(self):
        indication = """TRUQAP, in combination with fulvestrant, is indicated for 
                the treatment of adult patients with hormone receptor (HR)-positive, 
                HER2-negative breast cancer with one or more PIK3CA/AKT1/PTEN -alteration"""
        result = match_fda_approval("AFGET1", "A424", "breast cancer", indication)
        assert not result["matched"]

    def test_codon_match_braf(self):
        # indication = """treatment of adult patients with KRAS G12C-mutated
        # locally advanced or metastatic NSCLC, as determined by an FDA-approved test, who have received at least one pri"""

        indication="""
        ZELBORAF is indicated for the treatment of patients with 
        unresectable or metastatic melanoma with BRAF V600 mutation as detected 
        by an FDA-approved test.
        """
        result = match_fda_approval("BRAF", "V600K", "melanoma", indication)
        assert result["match_level"]=='codon'

    def test_codon_match(self):
        indication = "treatment of adult patients with KRAS G12C-mutated NSCLC"

        # Exact match
        result = match_fda_approval("KRAS", "G12C", "NSCLC", indication)
        print(result)
        # {"matched": True, "match_level": "variant", ...}

        # Same codon, different variant
        result = match_fda_approval("KRAS", "G12D", "NSCLC", indication)
        print(result)
        # {"matched": False, "match_level": "codon", ...}

        # Different codon entirely
        result = match_fda_approval("KRAS", "Q61H", "NSCLC", indication)
        print(result)
        # {"matched": False, "match_level": None, ...}


    def test_cap(self):
        """Test Capivasertib - gene-level approval for PIK3CA/AKT1/PTEN."""
        text = "metastatic breast cancer with one or more PIK3CA/AKT1/PTEN -alteration as detected by an FDA-approved test following"
        specificity = parse_biomarker_specificity(text, "AKT1")

        assert specificity is not None
        assert specificity["level"] == "gene"

        covered, level = is_variant_covered("E17K", specificity)
        assert covered is True
        assert level == "gene"

    def test_variant_level_braf_v600e(self):
        """Test parsing variant-level specificity like BRAF V600E."""
        text = "for the treatment of patients with BRAF V600E mutation-positive melanoma"
        result = parse_biomarker_specificity(text, "BRAF")

        assert result is not None
        assert result["level"] == "variant"
        assert result["codon"] == "V600"
        assert result["specified_variant"] == "V600E"

    def test_variant_level_egfr_t790m(self):
        """Test parsing variant-level specificity like EGFR T790M."""
        text = "indicated for EGFR T790M mutation-positive NSCLC"
        result = parse_biomarker_specificity(text, "EGFR")

        assert result is not None
        assert result["level"] == "variant"
        assert result["codon"] == "T790"
        assert result["specified_variant"] == "T790M"

    def test_codon_level_braf_v600(self):
        """Test parsing codon-level specificity like BRAF V600."""
        text = "for patients with BRAF V600 mutation in melanoma"
        result = parse_biomarker_specificity(text, "BRAF")

        assert result is not None
        assert result["level"] == "codon"
        assert result["codon"] == "V600"

    def test_exon_level_egfr_exon_19(self):
        """Test parsing exon-level specificity like EGFR exon 19."""
        text = "indicated for EGFR exon 19 deletion or exon 21 L858R mutations"
        result = parse_biomarker_specificity(text, "EGFR")

        assert result is not None
        assert result["level"] == "exon"
        assert result["exon"] == "19"

    def test_exon_level_case_insensitive(self):
        """Test that exon matching is case insensitive."""
        text = "indicated for EGFR Exon 20 insertion mutations"
        result = parse_biomarker_specificity(text, "EGFR")

        assert result is not None
        assert result["level"] == "exon"
        assert result["exon"] == "20"

    def test_gene_level_alteration(self):
        """Test parsing gene-level specificity with 'alteration'."""
        text = "for patients with AKT1 alteration as detected by an FDA-approved test"
        result = parse_biomarker_specificity(text, "AKT1")

        assert result is not None
        assert result["level"] == "gene"

    def test_gene_level_mutation(self):
        """Test parsing gene-level specificity with 'mutation'."""
        text = "indicated for PIK3CA mutation positive breast cancer"
        result = parse_biomarker_specificity(text, "PIK3CA")

        assert result is not None
        assert result["level"] == "gene"

    def test_gene_level_with_slash_list(self):
        """Test parsing gene-level from slash-separated list."""
        text = "with PIK3CA/AKT1/PTEN -alteration as detected by an FDA-approved test"
        result = parse_biomarker_specificity(text, "AKT1")

        assert result is not None
        assert result["level"] == "gene"

    def test_no_match_returns_none(self):
        """Test that no match returns None."""
        text = "indicated for breast cancer"
        result = parse_biomarker_specificity(text, "BRAF")

        assert result is None

    def test_no_match_different_gene(self):
        """Test that wrong gene returns None."""
        text = "indicated for KRAS G12C mutated cancer"
        result = parse_biomarker_specificity(text, "BRAF")

        assert result is None

    def test_case_insensitive_gene(self):
        """Test that gene matching works with different cases in text."""
        text = "for patients with kras G12C mutation"
        result = parse_biomarker_specificity(text, "KRAS")

        assert result is not None
        assert result["level"] == "variant"
        assert result["codon"] == "G12"
        assert result["specified_variant"] == "G12C"

    def test_real_capivasertib_indication(self):
        """Test with real Capivasertib (TRUQAP) indication text."""
        text = """TRUQAP, in combination with fulvestrant, is indicated for the treatment
        of adult patients with hormone receptor (HR)-positive, human epidermal growth
        factor receptor 2 (HER2)-negative, locally advanced or metastatic breast cancer
        with one or more PIK3CA/AKT1/PTEN -alteration as detected by an FDA-approved test"""

        result = parse_biomarker_specificity(text, "AKT1")

        assert result is not None
        assert result["level"] == "gene"

    def test_real_osimertinib_indication(self):
        """Test with real Osimertinib (TAGRISSO) T790M indication text."""
        text = """TAGRISSO is indicated for the treatment of adult patients with
        metastatic EGFR T790M mutation-positive non-small cell lung cancer (NSCLC),
        as detected by an FDA-approved test, whose disease has progressed on or
        after EGFR TKI therapy"""

        result = parse_biomarker_specificity(text, "EGFR")

        assert result is not None
        assert result["level"] == "variant"
        assert result["codon"] == "T790"
        assert result["specified_variant"] == "T790M"

    def test_real_sotorasib_indication(self):
        """Test with real Sotorasib (LUMAKRAS) indication text."""
        text = """LUMAKRAS is indicated for the treatment of adult patients with
        KRAS G12C-mutated locally advanced or metastatic non-small cell lung cancer"""

        result = parse_biomarker_specificity(text, "KRAS")

        assert result is not None
        assert result["level"] == "variant"
        assert result["codon"] == "G12"
        assert result["specified_variant"] == "G12C"


class TestIsVariantCovered:
    """Tests for is_variant_covered function.

    is_variant_covered returns (covered: bool, match_level: str | None)
    """

    def test_gene_level_covers_any_variant(self):
        """Test that gene-level approval covers any variant."""
        specificity = {"level": "gene"}

        covered, level = is_variant_covered("E17K", specificity)
        assert covered is True
        assert level == "gene"

        covered, level = is_variant_covered("V600E", specificity)
        assert covered is True

        covered, level = is_variant_covered("G12C", specificity)
        assert covered is True

    def test_variant_level_exact_match(self):
        """Test variant-level approval requires exact match."""
        specificity = {"level": "variant", "codon": "G12", "specified_variant": "G12C"}

        covered, level = is_variant_covered("G12C", specificity)
        assert covered is True
        assert level == "variant"

        # Same codon, different variant - not covered but match_level is "codon"
        covered, level = is_variant_covered("G12D", specificity)
        assert covered is False
        assert level == "codon"

        covered, level = is_variant_covered("G12V", specificity)
        assert covered is False
        assert level == "codon"

    def test_variant_level_case_insensitive(self):
        """Test variant matching is case insensitive."""
        specificity = {"level": "variant", "codon": "G12", "specified_variant": "G12C"}

        covered, level = is_variant_covered("g12c", specificity)
        assert covered is True
        assert level == "variant"

        covered, level = is_variant_covered("G12c", specificity)
        assert covered is True

    def test_codon_level_covers_same_codon(self):
        """Test codon-level approval covers all variants at same codon."""
        specificity = {"level": "codon", "codon": "V600"}

        covered, level = is_variant_covered("V600E", specificity)
        assert covered is True
        assert level == "codon"

        covered, _ = is_variant_covered("V600K", specificity)
        assert covered is True

        covered, _ = is_variant_covered("V600R", specificity)
        assert covered is True

        covered, level = is_variant_covered("V601E", specificity)
        assert covered is False
        assert level is None

    def test_codon_level_different_codon(self):
        """Test codon-level approval doesn't cover different codons."""
        specificity = {"level": "codon", "codon": "G12"}

        covered, _ = is_variant_covered("G12C", specificity)
        assert covered is True

        covered, _ = is_variant_covered("G12D", specificity)
        assert covered is True

        covered, level = is_variant_covered("G13D", specificity)
        assert covered is False
        assert level is None

    def test_exon_level_returns_false(self):
        """Test exon-level returns False (needs exon mapping)."""
        specificity = {"level": "exon", "exon": "19"}

        # We don't have exon mapping yet, so this should return False
        covered, level = is_variant_covered("L858R", specificity)
        assert covered is False
        assert level is None

    def test_real_scenario_akt1_e17k_match(self):
        text = "metastatic breast cancer with one or more PIK3CA/AKT1/PTEN -alteration as detected by an FDA-approved test following"
        specificity = parse_biomarker_specificity(text, "AKT1")

        covered, level = is_variant_covered("E17K", specificity)
        assert covered is True
        assert level == "gene"

    def test_real_scenario_kras_g12c_variant_level(self):
        """Test KRAS G12C with variant-level approval (Sotorasib)."""
        # Sotorasib is approved specifically for KRAS G12C
        specificity = {"level": "variant", "codon": "G12", "specified_variant": "G12C"}

        covered, level = is_variant_covered("G12C", specificity)
        assert covered is True
        assert level == "variant"

        # G12D is same codon but not covered
        covered, level = is_variant_covered("G12D", specificity)
        assert covered is False
        assert level == "codon"

    def test_real_scenario_egfr_t790m_variant_level(self):
        """Test EGFR T790M with variant-level approval (Osimertinib)."""
        # Osimertinib is approved specifically for EGFR T790M
        specificity = {"level": "variant", "codon": "T790", "specified_variant": "T790M"}

        covered, level = is_variant_covered("T790M", specificity)
        assert covered is True
        assert level == "variant"

        covered, level = is_variant_covered("L858R", specificity)
        assert covered is False
        assert level is None

    def test_real_scenario_braf_v600_codon_level(self):
        """Test BRAF V600 with codon-level approval."""
        # Some drugs are approved for BRAF V600 mutations (V600E, V600K, etc.)
        specificity = {"level": "codon", "codon": "V600"}

        covered, level = is_variant_covered("V600E", specificity)
        assert covered is True
        assert level == "codon"

        covered, _ = is_variant_covered("V600K", specificity)
        assert covered is True

        covered, level = is_variant_covered("K601E", specificity)
        assert covered is False
        assert level is None


class TestTumorMatching:
    """Tests for tumor type matching in populate_locus_variant_match."""

    def test_tumor_match_nsclc(self):
        """Test tumor matching for NSCLC."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match
        from oncomind.models.evidence.fda import FDALabelEvidence

        label = FDALabelEvidence(
            drug="sotorasib",
            gene="KRAS",
            indications_and_usage="indicated for KRAS G12C-mutated non-small cell lung cancer",
        )

        populate_locus_variant_match([label], query_variant="G12C", query_tumor="NSCLC")

        assert label.biomarker_match is not None
        assert label.biomarker_match.matched is True
        assert label.biomarker_match.match_level == "variant"
        assert label.biomarker_match.tumor_matched is True
        assert label.biomarker_match.tumor_match_type == "exact"

    def test_tumor_match_different_cancer(self):
        """Test tumor mismatch when cancer type differs."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match
        from oncomind.models.evidence.fda import FDALabelEvidence

        label = FDALabelEvidence(
            drug="sotorasib",
            gene="KRAS",
            indications_and_usage="indicated for KRAS G12C-mutated non-small cell lung cancer",
        )

        populate_locus_variant_match([label], query_variant="G12C", query_tumor="colorectal cancer")

        assert label.biomarker_match is not None
        assert label.biomarker_match.matched is True  # Variant matches
        assert label.biomarker_match.tumor_matched is False  # But tumor doesn't

    def test_tumor_match_no_tumor_provided(self):
        """Test that tumor_matched is None when no tumor type is provided."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match
        from oncomind.models.evidence.fda import FDALabelEvidence

        label = FDALabelEvidence(
            drug="sotorasib",
            gene="KRAS",
            indications_and_usage="indicated for KRAS G12C-mutated non-small cell lung cancer",
        )

        populate_locus_variant_match([label], query_variant="G12C", query_tumor=None)

        assert label.biomarker_match is not None
        assert label.biomarker_match.matched is True
        assert label.biomarker_match.tumor_matched is None
        assert label.biomarker_match.tumor_match_type is None

    def test_tumor_match_pan_cancer_msi_h(self):
        """Test pan-cancer MSI-H approval matches any solid tumor."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match
        from oncomind.models.evidence.fda import FDALabelEvidence

        label = FDALabelEvidence(
            drug="pembrolizumab",
            gene="MSH2",
            indications_and_usage="[FDA APPROVED FOR MSI-H/dMMR] for MSI-H or dMMR cancer that has progressed",
        )

        populate_locus_variant_match([label], query_variant="E17K", query_tumor="endometrial cancer")

        assert label.biomarker_match is not None
        # MSI-H is phenotype-based, not gene-based, so variant match is False
        assert label.biomarker_match.tumor_matched is True
        assert label.biomarker_match.tumor_match_type == "pan_cancer"

    def test_tumor_match_breast_cancer(self):
        """Test tumor matching for breast cancer with AKT1."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match
        from oncomind.models.evidence.fda import FDALabelEvidence

        label = FDALabelEvidence(
            drug="capivasertib",
            gene="AKT1",
            indications_and_usage="TRUQAP is indicated for HR-positive, HER2-negative breast cancer with PIK3CA/AKT1/PTEN -alteration",
        )

        populate_locus_variant_match([label], query_variant="E17K", query_tumor="breast cancer")

        assert label.biomarker_match is not None
        assert label.biomarker_match.matched is True  # Gene-level approval covers E17K
        assert label.biomarker_match.match_level == "gene"
        assert label.biomarker_match.tumor_matched is True
        assert label.biomarker_match.tumor_match_type == "exact"

    def test_tumor_match_melanoma_braf(self):
        """Test tumor matching for melanoma with BRAF V600."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match
        from oncomind.models.evidence.fda import FDALabelEvidence

        label = FDALabelEvidence(
            drug="vemurafenib",
            gene="BRAF",
            indications_and_usage="ZELBORAF is indicated for unresectable or metastatic melanoma with BRAF V600 mutation",
        )

        populate_locus_variant_match([label], query_variant="V600E", query_tumor="melanoma")

        assert label.biomarker_match is not None
        assert label.biomarker_match.matched is True
        assert label.biomarker_match.match_level == "codon"
        assert label.biomarker_match.tumor_matched is True

    def test_tumor_match_melanoma_wrong_tumor(self):
        """Test tumor mismatch when BRAF melanoma drug queried for NSCLC."""
        from oncomind.insight_builder.fda_processor import populate_locus_variant_match
        from oncomind.models.evidence.fda import FDALabelEvidence

        label = FDALabelEvidence(
            drug="vemurafenib",
            gene="BRAF",
            indications_and_usage="ZELBORAF is indicated for unresectable or metastatic melanoma with BRAF V600 mutation",
        )

        populate_locus_variant_match([label], query_variant="V600E", query_tumor="NSCLC")

        assert label.biomarker_match is not None
        assert label.biomarker_match.matched is True  # Variant matches at codon level
        assert label.biomarker_match.tumor_matched is False  # But tumor doesn't match

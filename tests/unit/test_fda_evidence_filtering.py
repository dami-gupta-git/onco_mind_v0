"""Tests for FDA evidence filtering in Evidence model."""

import pytest

from oncomind.models.evidence.evidence import Evidence, VariantIdentifiers
from oncomind.models.evidence.fda_biomarker import (
    FDABiomarkerEvidence,
    BiomarkerRequirement,
    SpecificityLevel,
)


def create_evidence() -> Evidence:
    """Create an Evidence object with required fields."""
    return Evidence(
        identifiers=VariantIdentifiers(
            variant_id="TEST:V600E",
            gene="TEST",
            variant="V600E",
        )
    )


class TestGetFilteredFDAEvidence:
    """Tests for Evidence.get_filtered_fda_evidence()."""

    def _create_fda_evidence(
        self,
        drug_name: str,
        gene: str,
        requirement: BiomarkerRequirement = BiomarkerRequirement.REQUIRED_POSITIVE,
        specificity: SpecificityLevel = SpecificityLevel.VARIANT,
        specified_variants: list[str] | None = None,
        tumor_types: list[str] | None = None,
        combination_partners: list[str] | None = None,
        codon: str | None = None,
    ) -> FDABiomarkerEvidence:
        """Helper to create FDABiomarkerEvidence for testing."""
        return FDABiomarkerEvidence(
            drug_name=drug_name,
            gene=gene,
            requirement=requirement,
            specificity=specificity,
            specified_variants=specified_variants or [],
            tumor_types=tumor_types or [],
            combination_partners=combination_partners or [],
            codon=codon,
        )

    def test_filters_by_gene(self):
        """Evidence for different gene is filtered out."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
            self._create_fda_evidence(
                drug_name="Osimertinib",
                gene="EGFR",
                specified_variants=["L858R"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600E")

        assert len(filtered) == 1
        assert filtered[0].drug_name == "Vemurafenib"

    def test_filters_unknown_drug_names(self):
        """Evidence with drug_name 'UNKNOWN' is filtered out."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="UNKNOWN",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600E")

        assert len(filtered) == 1
        assert filtered[0].drug_name == "Vemurafenib"

    def test_filters_required_negative(self):
        """Evidence with REQUIRED_NEGATIVE is filtered out."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Imjudo",
                gene="EGFR",
                requirement=BiomarkerRequirement.REQUIRED_NEGATIVE,
                specificity=SpecificityLevel.GENE,
            ),
            self._create_fda_evidence(
                drug_name="Osimertinib",
                gene="EGFR",
                specified_variants=["L858R"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("EGFR", "L858R")

        assert len(filtered) == 1
        assert filtered[0].drug_name == "Osimertinib"

    def test_filters_by_tumor_type(self):
        """Evidence for non-matching tumor type is filtered out when tumor specified."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
                tumor_types=["Melanoma"],
            ),
            self._create_fda_evidence(
                drug_name="Dabrafenib",
                gene="BRAF",
                specified_variants=["V600E"],
                tumor_types=["Non-Small Cell Lung Cancer"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600E", "Melanoma")

        assert len(filtered) == 1
        assert filtered[0].drug_name == "Vemurafenib"

    def test_no_tumor_filter_returns_all(self):
        """When no tumor specified, all matching drugs are returned."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
                tumor_types=["Melanoma"],
            ),
            self._create_fda_evidence(
                drug_name="Dabrafenib",
                gene="BRAF",
                specified_variants=["V600E"],
                tumor_types=["Non-Small Cell Lung Cancer"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600E")

        assert len(filtered) == 2

    def test_variant_match_exact(self):
        """Exact variant match returns match_type 'exact'."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600E")

        assert len(filtered) == 1
        assert filtered[0].variant_match_result == "exact"

    def test_variant_match_gene_level(self):
        """Gene-level approval matches any variant."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Sotorasib",
                gene="KRAS",
                specificity=SpecificityLevel.GENE,
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("KRAS", "G12D")

        assert len(filtered) == 1
        assert filtered[0].variant_match_result == "gene"

    def test_variant_match_codon_level(self):
        """Codon-level approval matches variants at same codon."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Encorafenib",
                gene="BRAF",
                specificity=SpecificityLevel.CODON,
                codon="600",
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600K")

        assert len(filtered) == 1
        assert filtered[0].variant_match_result == "codon"

    def test_non_matching_variant_filtered_out(self):
        """Non-matching variant is filtered out."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "G469A")

        assert len(filtered) == 0

    def test_deduplicates_same_drug(self):
        """Duplicate drugs are deduplicated."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
                tumor_types=["Melanoma"],
            ),
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
                tumor_types=["Erdheim-Chester Disease"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600E")

        assert len(filtered) == 1
        assert filtered[0].drug_name == "Vemurafenib"

    def test_deduplicates_combination_drugs_different_order(self):
        """Combination drugs with different ordering are deduplicated."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="TRAMETINIB",
                gene="BRAF",
                specified_variants=["V600E"],
                combination_partners=["dabrafenib"],
            ),
            self._create_fda_evidence(
                drug_name="DABRAFENIB",
                gene="BRAF",
                specified_variants=["V600E"],
                combination_partners=["trametinib"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600E")

        # Should only get one entry since the combination is the same
        assert len(filtered) == 1

    def test_sets_filtered_count(self):
        """filtered_fda_biomarker_count is set after filtering."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
            self._create_fda_evidence(
                drug_name="Dabrafenib",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
        ]

        evidence.get_filtered_fda_evidence("BRAF", "V600E")

        assert evidence.filtered_fda_biomarker_count == 2

    def test_normalizes_three_letter_variant(self):
        """Three-letter variant notation (p.Leu858Arg) matches single-letter (L858R)."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Osimertinib",
                gene="EGFR",
                specified_variants=["L858R"],
            ),
        ]

        # Query with three-letter notation
        filtered = evidence.get_filtered_fda_evidence("EGFR", "p.Leu858Arg")

        assert len(filtered) == 1
        assert filtered[0].drug_name == "Osimertinib"
        assert filtered[0].variant_match_result == "exact"

    def test_normalizes_hgvs_variant(self):
        """HGVS variant notation (p.V600E) matches short form (V600E)."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("BRAF", "p.V600E")

        assert len(filtered) == 1
        assert filtered[0].variant_match_result == "exact"

    def test_case_insensitive_gene(self):
        """Gene matching is case-insensitive."""
        evidence = create_evidence()
        evidence.fda_biomarker_evidence = [
            self._create_fda_evidence(
                drug_name="Vemurafenib",
                gene="BRAF",
                specified_variants=["V600E"],
            ),
        ]

        filtered = evidence.get_filtered_fda_evidence("braf", "V600E")

        assert len(filtered) == 1

    def test_empty_evidence_returns_empty_list(self):
        """Empty fda_biomarker_evidence returns empty list."""
        evidence = create_evidence()

        filtered = evidence.get_filtered_fda_evidence("BRAF", "V600E")

        assert filtered == []
        assert evidence.filtered_fda_biomarker_count == 0


class TestToShortForm:
    """Tests for to_short_form variant normalization function."""

    def test_single_letter_passthrough(self):
        """Single-letter format is passed through."""
        from oncomind.utils.variant_normalization import to_short_form

        assert to_short_form("V600E") == "V600E"
        assert to_short_form("L858R") == "L858R"
        assert to_short_form("G12C") == "G12C"

    def test_three_letter_to_single(self):
        """Three-letter format is converted to single-letter."""
        from oncomind.utils.variant_normalization import to_short_form

        assert to_short_form("Val600Glu") == "V600E"
        assert to_short_form("Leu858Arg") == "L858R"
        assert to_short_form("Gly12Cys") == "G12C"

    def test_hgvs_prefix_stripped(self):
        """p. prefix is stripped."""
        from oncomind.utils.variant_normalization import to_short_form

        assert to_short_form("p.V600E") == "V600E"
        assert to_short_form("p.Leu858Arg") == "L858R"

    def test_non_missense_returns_none(self):
        """Non-missense variants return None."""
        from oncomind.utils.variant_normalization import to_short_form

        assert to_short_form("fusion") is None
        assert to_short_form("amplification") is None
        assert to_short_form("exon 19 deletion") is None

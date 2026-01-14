"""Tests for input parser functions."""

import pytest
from oncomind.normalization.input_parser import (
    ParsedVariant,
    parse_variant_input,
    parse_variant_row,
    parse_vcf_variant,
    _resolve_gene_alias,
    _extract_tumor_type,
)


class TestParseVariantInput:
    """Tests for parse_variant_input function."""

    def test_simple_gene_variant(self):
        """Test parsing simple gene + variant format."""
        result = parse_variant_input("BRAF V600E")

        assert result.gene == "BRAF"
        assert result.variant == "V600E"
        assert result.parse_confidence == 1.0
        assert result.tumor_type is None
        assert len(result.parse_warnings) == 0

    def test_colon_separated_format(self):
        """Test parsing colon-separated gene:variant format."""
        result = parse_variant_input("BRAF:V600E")

        assert result.gene == "BRAF"
        assert result.variant == "V600E"
        assert result.parse_confidence == 1.0

    def test_with_p_notation(self):
        """Test parsing with HGVS p. notation."""
        result = parse_variant_input("BRAF p.V600E")

        assert result.gene == "BRAF"
        assert result.variant == "V600E"
        assert result.parse_confidence == 1.0

    def test_three_letter_amino_acid(self):
        """Test parsing with three-letter amino acid codes."""
        result = parse_variant_input("BRAF p.Val600Glu")

        assert result.gene == "BRAF"
        assert result.variant == "Val600Glu"
        assert result.variant_normalized == "V600E"

    def test_nonsense_variant(self):
        """Test parsing nonsense (stop codon) variant."""
        result = parse_variant_input("TP53 R248*")

        assert result.gene == "TP53"
        assert result.variant == "R248*"
        assert result.variant_type == "nonsense"

    def test_complex_variant_delins(self):
        """Test parsing complex deletion-insertion variant."""
        result = parse_variant_input("EGFR L747_P753delinsS")

        assert result.gene == "EGFR"
        assert result.variant == "L747_P753delinsS"

    def test_tumor_type_extraction_with_in(self):
        """Test extraction of tumor type with 'in' keyword."""
        result = parse_variant_input("EGFR L858R in lung cancer")

        assert result.gene == "EGFR"
        assert result.variant == "L858R"
        assert result.tumor_type == "lung cancer"

    def test_tumor_type_extraction_with_with(self):
        """Test extraction of tumor type with 'with' keyword."""
        result = parse_variant_input("BRAF V600E with melanoma")

        assert result.gene == "BRAF"
        assert result.variant == "V600E"
        assert result.tumor_type == "melanoma"

    def test_tumor_type_extraction_with_for(self):
        """Test extraction of tumor type with 'for' keyword."""
        result = parse_variant_input("KRAS G12C for NSCLC")

        assert result.gene == "KRAS"
        assert result.variant == "G12C"
        assert result.tumor_type == "NSCLC"

    def test_explicit_tumor_type_overrides_extracted(self):
        """Test that explicit tumor_type parameter overrides extracted."""
        result = parse_variant_input("EGFR L858R in lung cancer", tumor_type="NSCLC")

        assert result.gene == "EGFR"
        assert result.variant == "L858R"
        assert result.tumor_type == "NSCLC"

    def test_case_insensitivity(self):
        """Test case insensitive parsing."""
        result = parse_variant_input("braf v600e")

        assert result.gene == "BRAF"
        assert result.variant == "v600e"

    def test_whitespace_handling(self):
        """Test handling of extra whitespace."""
        result = parse_variant_input("  BRAF   V600E  ")

        assert result.gene == "BRAF"
        assert result.variant == "V600E"

    def test_fallback_parsing(self):
        """Test fallback parsing for non-standard formats."""
        result = parse_variant_input("BRAF mutation_xyz")

        assert result.gene == "BRAF"
        assert result.variant == "mutation_xyz"
        assert result.parse_confidence == 0.7
        assert "fallback" in result.parse_warnings[0].lower()

    def test_raw_input_preserved(self):
        """Test that raw input is preserved."""
        input_str = "EGFR L858R in NSCLC"
        result = parse_variant_input(input_str)

        assert result.raw_input == input_str

    def test_invalid_input_raises_error(self):
        """Test that invalid input raises ValueError."""
        with pytest.raises(ValueError, match="Could not parse"):
            parse_variant_input("invalid")

    def test_empty_input_raises_error(self):
        """Test that empty input raises ValueError."""
        with pytest.raises(ValueError, match="Could not parse"):
            parse_variant_input("")

    def test_gene_alias_resolution(self):
        """Test that gene aliases are resolved to canonical names."""
        # NEU should resolve to ERBB2
        result = parse_variant_input("NEU V777L")

        assert result.gene == "ERBB2"
        assert result.variant == "V777L"

    def test_variant_type_classification(self):
        """Test that variant type is classified."""
        result = parse_variant_input("BRAF V600E")

        assert result.variant_type == "missense"

    def test_variant_normalization(self):
        """Test that variant is normalized."""
        result = parse_variant_input("BRAF V600E")

        assert result.variant_normalized is not None

    def test_to_dict_method(self):
        """Test conversion to dictionary."""
        result = parse_variant_input("BRAF V600E")
        d = result.to_dict()

        assert d["gene"] == "BRAF"
        assert d["variant"] == "V600E"
        assert "variant_normalized" in d
        assert "variant_type" in d
        assert "tumor_type" in d
        assert "raw_input" in d
        assert "parse_confidence" in d
        assert "parse_warnings" in d


class TestResolveGeneAlias:
    """Tests for _resolve_gene_alias helper function."""

    def test_canonical_gene_unchanged(self):
        """Test that canonical gene names are unchanged."""
        assert _resolve_gene_alias("BRAF") == "BRAF"
        assert _resolve_gene_alias("EGFR") == "EGFR"
        assert _resolve_gene_alias("KRAS") == "KRAS"

    def test_case_insensitivity(self):
        """Test case insensitive lookup."""
        assert _resolve_gene_alias("braf") == "BRAF"
        assert _resolve_gene_alias("Braf") == "BRAF"

    def test_unknown_gene_returned_uppercase(self):
        """Test that unknown genes are returned uppercase."""
        assert _resolve_gene_alias("unknowngene") == "UNKNOWNGENE"


class TestExtractTumorType:
    """Tests for _extract_tumor_type helper function."""

    def test_extracts_with_in_keyword(self):
        """Test extraction with 'in' keyword."""
        cleaned, tumor = _extract_tumor_type("BRAF V600E in melanoma")

        assert cleaned == "BRAF V600E"
        assert tumor == "melanoma"

    def test_extracts_with_with_keyword(self):
        """Test extraction with 'with' keyword."""
        cleaned, tumor = _extract_tumor_type("EGFR L858R with NSCLC")

        assert cleaned == "EGFR L858R"
        assert tumor == "NSCLC"

    def test_extracts_with_for_keyword(self):
        """Test extraction with 'for' keyword."""
        cleaned, tumor = _extract_tumor_type("KRAS G12C for lung cancer")

        assert cleaned == "KRAS G12C"
        assert tumor == "lung cancer"

    def test_no_tumor_type_returns_none(self):
        """Test that input without tumor type returns None."""
        cleaned, tumor = _extract_tumor_type("BRAF V600E")

        assert cleaned == "BRAF V600E"
        assert tumor is None


class TestParseVariantRow:
    """Tests for parse_variant_row function."""

    def test_basic_row_parsing(self):
        """Test parsing a basic row dictionary."""
        row = {"gene": "BRAF", "variant": "V600E", "tumor_type": "Melanoma"}
        result = parse_variant_row(row)

        assert result.gene == "BRAF"
        assert result.variant == "V600E"
        assert result.tumor_type == "Melanoma"
        assert result.parse_confidence == 1.0

    def test_custom_column_names(self):
        """Test parsing with custom column names."""
        row = {"gene_symbol": "EGFR", "mutation": "L858R"}
        result = parse_variant_row(
            row,
            gene_col="gene_symbol",
            variant_col="mutation",
            tumor_col=None,
        )

        assert result.gene == "EGFR"
        assert result.variant == "L858R"
        assert result.tumor_type is None

    def test_whitespace_in_values(self):
        """Test that whitespace in values is stripped."""
        row = {"gene": "  BRAF  ", "variant": " V600E "}
        result = parse_variant_row(row)

        assert result.gene == "BRAF"
        assert result.variant == "V600E"

    def test_missing_gene_raises_error(self):
        """Test that missing gene raises ValueError."""
        row = {"variant": "V600E"}

        with pytest.raises(ValueError, match="Missing gene"):
            parse_variant_row(row)

    def test_missing_variant_raises_error(self):
        """Test that missing variant raises ValueError."""
        row = {"gene": "BRAF"}

        with pytest.raises(ValueError, match="Missing gene or variant"):
            parse_variant_row(row)

    def test_gene_alias_resolved(self):
        """Test that gene aliases are resolved in row parsing."""
        row = {"gene": "NEU", "variant": "V777L"}
        result = parse_variant_row(row)

        assert result.gene == "ERBB2"


class TestParseVcfVariant:
    """Tests for parse_vcf_variant function."""

    def test_basic_snv(self):
        """Test parsing a basic SNV."""
        result = parse_vcf_variant("chr7", 140453136, "A", "T")

        assert result.variant == "chr7:140453136:A>T"
        assert result.variant_type == "snv"
        assert result.gene == "UNKNOWN"
        assert result.parse_confidence == 0.5

    def test_indel_detection(self):
        """Test that indels are detected."""
        result = parse_vcf_variant("chr7", 140453136, "AT", "A")

        assert result.variant_type == "indel"

    def test_gene_from_info(self):
        """Test extraction of gene from INFO field."""
        result = parse_vcf_variant(
            "chr7", 140453136, "A", "T",
            info={"GENE": "BRAF"}
        )

        assert result.gene == "BRAF"
        assert result.parse_confidence == 0.8

    def test_gene_from_info_lowercase(self):
        """Test extraction of gene from lowercase INFO key."""
        result = parse_vcf_variant(
            "chr7", 140453136, "A", "T",
            info={"gene": "EGFR"}
        )

        assert result.gene == "EGFR"

    def test_warning_about_basic_parsing(self):
        """Test that warning about basic parsing is present."""
        result = parse_vcf_variant("chr7", 140453136, "A", "T")

        assert any("basic" in w.lower() for w in result.parse_warnings)


class TestParsedVariantDataclass:
    """Tests for ParsedVariant dataclass."""

    def test_default_values(self):
        """Test default values in ParsedVariant."""
        pv = ParsedVariant(gene="BRAF", variant="V600E")

        assert pv.variant_normalized is None
        assert pv.variant_type is None
        assert pv.tumor_type is None
        assert pv.raw_input == ""
        assert pv.parse_confidence == 1.0
        assert pv.parse_warnings == []

    def test_to_dict_includes_all_fields(self):
        """Test that to_dict includes all fields."""
        pv = ParsedVariant(
            gene="BRAF",
            variant="V600E",
            variant_normalized="V600E",
            variant_type="missense",
            tumor_type="Melanoma",
            raw_input="BRAF V600E in Melanoma",
            parse_confidence=1.0,
            parse_warnings=["test warning"],
        )
        d = pv.to_dict()

        assert d["gene"] == "BRAF"
        assert d["variant"] == "V600E"
        assert d["variant_normalized"] == "V600E"
        assert d["variant_type"] == "missense"
        assert d["tumor_type"] == "Melanoma"
        assert d["raw_input"] == "BRAF V600E in Melanoma"
        assert d["parse_confidence"] == 1.0
        assert d["parse_warnings"] == ["test warning"]


class TestRealWorldVariants:
    """Tests with real-world variant examples."""

    @pytest.mark.parametrize("input_str,expected_gene,expected_variant", [
        ("BRAF V600E", "BRAF", "V600E"),
        ("EGFR L858R", "EGFR", "L858R"),
        ("KRAS G12C", "KRAS", "G12C"),
        ("KRAS G12D", "KRAS", "G12D"),
        ("TP53 R175H", "TP53", "R175H"),
        ("PIK3CA H1047R", "PIK3CA", "H1047R"),
        ("NRAS Q61K", "NRAS", "Q61K"),
        ("IDH1 R132H", "IDH1", "R132H"),
        ("KIT D816V", "KIT", "D816V"),
    ])
    def test_common_cancer_variants(self, input_str, expected_gene, expected_variant):
        """Test parsing common cancer variants."""
        result = parse_variant_input(input_str)

        assert result.gene == expected_gene
        assert result.variant == expected_variant
        assert result.variant_type == "missense"

    @pytest.mark.parametrize("input_str,expected_gene,expected_variant", [
        ("BRAF:V600E", "BRAF", "V600E"),
        ("EGFR:L858R", "EGFR", "L858R"),
        ("KRAS:G12C", "KRAS", "G12C"),
    ])
    def test_colon_format_variants(self, input_str, expected_gene, expected_variant):
        """Test parsing colon-separated format."""
        result = parse_variant_input(input_str)

        assert result.gene == expected_gene
        assert result.variant == expected_variant

    @pytest.mark.parametrize("input_str,expected_tumor", [
        ("EGFR L858R in NSCLC", "NSCLC"),
        ("BRAF V600E in melanoma", "melanoma"),
        ("KRAS G12C with lung cancer", "lung cancer"),
        ("PIK3CA H1047R for breast cancer", "breast cancer"),
    ])
    def test_variants_with_tumor_types(self, input_str, expected_tumor):
        """Test parsing variants with tumor types."""
        result = parse_variant_input(input_str)

        assert result.tumor_type == expected_tumor

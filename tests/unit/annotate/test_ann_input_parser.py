"""Unit tests for AnnInputParser."""

import pytest

from oncomind.annotator.ann_input_parser import AnnInputParser
from oncomind.annotator.ann_parsed_variant import AnnParsedVariant


class TestAnnInputParser:
    """Tests for AnnInputParser class."""

    def test_parser_creates_parsed_variant(self):
        """Test parser creates AnnParsedVariant from gene and variant."""
        parser = AnnInputParser(gene="BRAF", variant="V600E")
        result = parser.to_parsed_variant()

        assert isinstance(result, AnnParsedVariant)
        assert result.gene == "BRAF"
        assert result.variant == "V600E"

    def test_parser_uppercases_gene(self):
        """Test parser converts gene to uppercase."""
        parser = AnnInputParser(gene="braf", variant="V600E")
        result = parser.to_parsed_variant()

        assert result.gene == "BRAF"

    def test_parser_sets_tumor_type(self):
        """Test parser sets tumor type."""
        parser = AnnInputParser(gene="BRAF", variant="V600E", tumor_type="NSCLC")
        result = parser.to_parsed_variant()

        assert result.tumor_type == "NSCLC"

    def test_parser_populates_protein_field(self):
        """Test parser populates protein field with p. notation."""
        parser = AnnInputParser(gene="BRAF", variant="V600E")
        result = parser.to_parsed_variant()

        assert result.protein is not None
        assert result.protein.startswith("p.")

    def test_parser_sets_raw_input(self):
        """Test parser sets raw_input field."""
        parser = AnnInputParser(gene="BRAF", variant="V600E")
        result = parser.to_parsed_variant()

        assert result.raw_input == "BRAF V600E"

    def test_parser_sets_variant_type(self):
        """Test parser sets variant_type field."""
        parser = AnnInputParser(gene="BRAF", variant="V600E")
        result = parser.to_parsed_variant()

        assert result.variant_type is not None

    def test_parser_default_confidence(self):
        """Test parser sets default parse confidence."""
        parser = AnnInputParser(gene="BRAF", variant="V600E")
        result = parser.to_parsed_variant()

        assert result.parse_confidence == 1.0

    def test_parser_empty_warnings(self):
        """Test parser initializes empty warnings list."""
        parser = AnnInputParser(gene="BRAF", variant="V600E")
        result = parser.to_parsed_variant()

        assert result.parse_warnings == []

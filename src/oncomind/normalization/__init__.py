"""Normalization module for variant input parsing and standardization.

This module provides tools to:
- Parse free-text variant inputs (e.g., "BRAF V600E")
- Normalize variant notations to standard forms
- Interface with VEP for functional consequence annotation

Example usage:
    >>> from oncomind.normalization import parse_variant_input, normalize_variant
    >>> parsed = parse_variant_input("BRAF V600E")
    >>> print(parsed.gene, parsed.variant)
    BRAF V600E
"""

from oncomind.normalization.input_parser import (
    parse_variant_input,
    parse_variant_row,
    ParsedVariant,
)
from oncomind.normalization.hgvs_utils import (
    VariantNormalizer,
    normalize_variant,
    normalize_protein_change,
    classify_variant_type,
    is_missense_variant,
    get_protein_position,
    to_hgvs_protein,
    is_snp_or_small_indel,
)

__all__ = [
    # Input parsing
    "parse_variant_input",
    "parse_variant_row",
    "ParsedVariant",
    # Normalization
    "VariantNormalizer",
    "normalize_variant",
    "normalize_protein_change",
    "classify_variant_type",
    "is_missense_variant",
    "get_protein_position",
    "to_hgvs_protein",
    "is_snp_or_small_indel",
]

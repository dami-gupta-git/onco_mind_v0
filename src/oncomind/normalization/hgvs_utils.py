"""HGVS notation utilities for variant normalization.

This module re-exports the VariantNormalizer class and convenience functions
from the utils module, providing a clean interface for the normalization package.
"""

# Re-export from existing utils module to maintain compatibility
from oncomind.utils.variant_normalization import (
    VariantNormalizer,
    normalize_variant,
    is_missense_variant,
    get_protein_position,
    to_hgvs_protein,
    is_snp_or_small_indel,
)


def normalize_protein_change(variant: str) -> dict:
    """Normalize a protein change to multiple standard formats.

    Args:
        variant: Protein change in any format (V600E, Val600Glu, p.V600E, etc.)

    Returns:
        Dictionary with normalized representations:
        - short_form: One-letter code (V600E)
        - hgvs_protein: HGVS protein notation (p.V600E)
        - long_form: Three-letter code (Val600Glu)
        - position: Position number (600)
        - ref_aa: Reference amino acid one-letter (V)
        - alt_aa: Alternate amino acid one-letter (E)
        - is_missense: Boolean indicating if this is a missense variant

    Examples:
        >>> normalize_protein_change("V600E")
        {'short_form': 'V600E', 'hgvs_protein': 'p.V600E', ...}
    """
    return VariantNormalizer.normalize_protein_change(variant)


def classify_variant_type(variant: str) -> str:
    """Classify the type of variant.

    Args:
        variant: Variant string in any format

    Returns:
        Variant type: 'missense', 'nonsense', 'frameshift', 'deletion',
                     'insertion', 'duplication', 'fusion', 'amplification',
                     'splice', 'truncating', or 'unknown'

    Examples:
        >>> classify_variant_type("V600E")
        'missense'
        >>> classify_variant_type("R248*")
        'nonsense'
        >>> classify_variant_type("fusion")
        'fusion'
    """
    return VariantNormalizer.classify_variant_type(variant)


__all__ = [
    "VariantNormalizer",
    "normalize_variant",
    "normalize_protein_change",
    "classify_variant_type",
    "is_missense_variant",
    "get_protein_position",
    "to_hgvs_protein",
    "is_snp_or_small_indel",
]

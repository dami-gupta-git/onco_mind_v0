"""Parsed variant container for the annotator workflow.

This module provides the AnnParsedVariant dataclass for storing
parsed variant information from VCF files.
"""

from typing import Any
from dataclasses import dataclass, field


@dataclass
class AnnParsedVariant:
    """Container for parsed variant information.

    Attributes:
        gene: Gene symbol (uppercase)
        variant: Variant notation
        protein: Protein notation (e.g., p.V600E)
        variant_normalized: Normalized variant (if applicable)
        variant_type: Classified variant type
        tumor_type: Tumor type (if provided)
        raw_input: Original input string
        parse_confidence: Confidence in parsing (0-1)
        parse_warnings: Any warnings from parsing
    """
    gene: str
    variant: str
    protein: str | None = None
    variant_normalized: str | None = None
    variant_type: str | None = None
    tumor_type: str | None = None
    raw_input: str = ""
    parse_confidence: float = 1.0
    parse_warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "gene": self.gene,
            "variant": self.variant,
            "protein": self.protein,
            "variant_normalized": self.variant_normalized,
            "variant_type": self.variant_type,
            "tumor_type": self.tumor_type,
            "raw_input": self.raw_input,
            "parse_confidence": self.parse_confidence,
            "parse_warnings": self.parse_warnings,
        }


__all__ = [
    "AnnParsedVariant",
]

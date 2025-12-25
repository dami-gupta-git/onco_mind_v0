"""Input parsing for variant strings and tabular data.

This module provides functions to parse variant inputs from:
- Free-text strings (e.g., "BRAF V600E", "EGFR L858R in NSCLC")
- CSV/DataFrame rows
- VCF records (basic support)

Examples:
    >>> from oncomind.normalization import parse_variant_input
    >>> parsed = parse_variant_input("BRAF V600E")
    >>> print(parsed.gene, parsed.variant)
    BRAF V600E

    >>> parsed = parse_variant_input("EGFR L858R in lung cancer")
    >>> print(parsed.gene, parsed.variant, parsed.tumor_type)
    EGFR L858R lung cancer
"""

import re
from typing import Any
from dataclasses import dataclass, field

from oncomind.utils.variant_normalization import normalize_variant, classify_variant_type
from oncomind.constants import GENE_ALIASES


@dataclass
class ParsedVariant:
    """Container for parsed variant information.

    Attributes:
        gene: Gene symbol (uppercase)
        variant: Variant notation
        variant_normalized: Normalized variant (if applicable)
        variant_type: Classified variant type
        tumor_type: Tumor type (if provided)
        raw_input: Original input string
        parse_confidence: Confidence in parsing (0-1)
        parse_warnings: Any warnings from parsing
    """
    gene: str
    variant: str
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
            "variant_normalized": self.variant_normalized,
            "variant_type": self.variant_type,
            "tumor_type": self.tumor_type,
            "raw_input": self.raw_input,
            "parse_confidence": self.parse_confidence,
            "parse_warnings": self.parse_warnings,
        }


# Common variant patterns
VARIANT_PATTERNS = [
    # Gene + variant (BRAF V600E, EGFR L858R)
    re.compile(r'^([A-Z][A-Z0-9]+)\s+([A-Z]\d+[A-Z*])$', re.IGNORECASE),
    # Gene + variant with p. notation (BRAF p.V600E)
    re.compile(r'^([A-Z][A-Z0-9]+)\s+p\.([A-Z][a-z]{0,2}\d+[A-Z][a-z]{0,2})$', re.IGNORECASE),
    # Gene + complex variant (EGFR L747_P753delinsS)
    re.compile(r'^([A-Z][A-Z0-9]+)\s+([A-Z]\d+[_A-Z\d]*(?:del|ins|dup|fs)[A-Z0-9]*)$', re.IGNORECASE),
    # Gene + nonsense (TP53 R248*)
    re.compile(r'^([A-Z][A-Z0-9]+)\s+([A-Z]\d+\*)$', re.IGNORECASE),
    # Gene colon variant (BRAF:V600E)
    re.compile(r'^([A-Z][A-Z0-9]+):([A-Z0-9_*]+)$', re.IGNORECASE),
]

# Tumor type keywords and patterns
TUMOR_KEYWORDS = [
    "in", "with", "cancer", "carcinoma", "tumor", "tumour",
    "melanoma", "nsclc", "sclc", "lung", "breast", "colorectal",
    "gist", "aml", "cml", "all", "cll", "lymphoma", "leukemia",
    "glioma", "glioblastoma", "pancreatic", "ovarian", "prostate",
]

TUMOR_SEPARATOR_PATTERN = re.compile(
    r'\s+(?:in|with|for)\s+(.+)$',
    re.IGNORECASE
)


def _resolve_gene_alias(gene: str) -> str:
    """Resolve gene aliases to canonical name.

    GENE_ALIASES maps canonical_gene -> [aliases].
    This function does the reverse lookup: alias -> canonical_gene.
    """
    gene_upper = gene.upper()

    # First, check if the gene is already a canonical name
    if gene_upper in GENE_ALIASES:
        return gene_upper

    # Otherwise, search for the gene in the aliases lists
    for canonical, aliases in GENE_ALIASES.items():
        if gene_upper in aliases:
            return canonical

    # Not found - return as-is
    return gene_upper


def _extract_tumor_type(text: str) -> tuple[str, str | None]:
    """Extract tumor type from input string.

    Returns:
        Tuple of (cleaned_text, tumor_type)
    """
    # Try to match common separator patterns
    match = TUMOR_SEPARATOR_PATTERN.search(text)
    if match:
        tumor_type = match.group(1).strip()
        cleaned = text[:match.start()].strip()
        return cleaned, tumor_type

    return text, None


def parse_variant_input(
    variant_str: str,
    tumor_type: str | None = None,
) -> ParsedVariant:
    """Parse a free-text variant input string.

    Args:
        variant_str: Free-text variant input (e.g., "BRAF V600E", "EGFR L858R in NSCLC")
        tumor_type: Optional tumor type (overrides any extracted from string)

    Returns:
        ParsedVariant with extracted information

    Raises:
        ValueError: If the variant string cannot be parsed

    Examples:
        >>> parse_variant_input("BRAF V600E")
        ParsedVariant(gene='BRAF', variant='V600E', ...)

        >>> parse_variant_input("EGFR L858R in lung cancer")
        ParsedVariant(gene='EGFR', variant='L858R', tumor_type='lung cancer', ...)

        >>> parse_variant_input("BRAF:V600E")
        ParsedVariant(gene='BRAF', variant='V600E', ...)
    """
    warnings: list[str] = []
    confidence = 1.0

    # Clean input
    text = variant_str.strip()
    raw_input = text

    # Extract tumor type from string if not provided
    extracted_tumor = None
    if tumor_type is None:
        text, extracted_tumor = _extract_tumor_type(text)

    final_tumor = tumor_type or extracted_tumor

    # Try each pattern
    gene = None
    variant = None

    for pattern in VARIANT_PATTERNS:
        match = pattern.match(text)
        if match:
            gene, variant = match.groups()
            break

    if gene is None or variant is None:
        # Try more aggressive parsing - split on whitespace
        parts = text.split()
        if len(parts) >= 2:
            gene = parts[0]
            variant = parts[1]
            confidence = 0.7
            warnings.append("Parsed with fallback whitespace splitting")
        else:
            raise ValueError(f"Could not parse variant string: {variant_str}")

    # Normalize gene alias
    gene = _resolve_gene_alias(gene)

    # Normalize variant
    norm_result = normalize_variant(gene, variant)
    variant_normalized = norm_result.get("variant_normalized", variant)
    variant_type = norm_result.get("variant_type", classify_variant_type(variant))

    return ParsedVariant(
        gene=gene,
        variant=variant,
        variant_normalized=variant_normalized,
        variant_type=variant_type,
        tumor_type=final_tumor,
        raw_input=raw_input,
        parse_confidence=confidence,
        parse_warnings=warnings,
    )


def parse_variant_row(
    row: dict[str, Any],
    gene_col: str = "gene",
    variant_col: str = "variant",
    tumor_col: str | None = "tumor_type",
) -> ParsedVariant:
    """Parse a variant from a DataFrame row or dictionary.

    Args:
        row: Dictionary with variant information
        gene_col: Column name for gene
        variant_col: Column name for variant
        tumor_col: Column name for tumor type (optional)

    Returns:
        ParsedVariant with extracted information

    Examples:
        >>> row = {"gene": "BRAF", "variant": "V600E", "tumor_type": "Melanoma"}
        >>> parse_variant_row(row)
        ParsedVariant(gene='BRAF', variant='V600E', tumor_type='Melanoma', ...)
    """
    gene = row.get(gene_col, "")
    variant = row.get(variant_col, "")
    tumor_type = row.get(tumor_col) if tumor_col else None

    if not gene or not variant:
        raise ValueError(f"Missing gene or variant in row: {row}")

    gene = _resolve_gene_alias(str(gene).strip())
    variant = str(variant).strip()

    # Normalize
    norm_result = normalize_variant(gene, variant)
    variant_normalized = norm_result.get("variant_normalized", variant)
    variant_type = norm_result.get("variant_type", classify_variant_type(variant))

    return ParsedVariant(
        gene=gene,
        variant=variant,
        variant_normalized=variant_normalized,
        variant_type=variant_type,
        tumor_type=tumor_type,
        raw_input=f"{gene} {variant}",
        parse_confidence=1.0,
        parse_warnings=[],
    )


def parse_vcf_variant(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    info: dict[str, Any] | None = None,
) -> ParsedVariant:
    """Parse a variant from VCF-style coordinates.

    Note: This is a basic implementation. For full VCF support,
    use VEP annotation through the vep_client module.

    Args:
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternate allele
        info: Optional INFO field dictionary

    Returns:
        ParsedVariant with basic information (gene may be unknown)
    """
    # Basic variant representation
    variant_id = f"{chrom}:{pos}:{ref}>{alt}"

    # Check INFO for gene annotation
    gene = "UNKNOWN"
    if info:
        gene = info.get("GENE", info.get("Gene", info.get("gene", "UNKNOWN")))

    return ParsedVariant(
        gene=gene,
        variant=variant_id,
        variant_normalized=None,
        variant_type="snv" if len(ref) == 1 and len(alt) == 1 else "indel",
        tumor_type=None,
        raw_input=variant_id,
        parse_confidence=0.5 if gene == "UNKNOWN" else 0.8,
        parse_warnings=["VCF parsing is basic - use VEP for full annotation"],
    )


__all__ = [
    "ParsedVariant",
    "parse_variant_input",
    "parse_variant_row",
    "parse_vcf_variant",
]

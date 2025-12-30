"""Base class for all evidence items with common metadata fields."""

from typing import Literal
from pydantic import BaseModel, Field

from oncomind.config.constants import BROAD_VARIANTS


# Type aliases for the allowed values
VariantGeneLevel = Literal["variant", "gene"]
Scope = Literal["specific", "unspecified", "ambiguous"]


def is_ambiguous_variant(gene: str, variant: str) -> bool:
    """Check if a gene/variant pair is in the ambiguous variants set.

    Args:
        gene: Gene symbol (e.g., "KRAS")
        variant: Variant notation (e.g., "G12", "V600")

    Returns:
        True if this gene/variant pair is ambiguous
    """
    return (gene.upper(), variant.upper()) in {
        (g.upper(), v.upper()) for g, v in BROAD_VARIANTS
    }


Origin = Literal["kb", "trial", "inferred"]

# CancerSpecificity can be:
# - "cancer_specific": evidence matches the user's queried tumor type
# - "pan_cancer": evidence is tumor-agnostic or unknown tumor type
# - Any other string: the specific cancer type the evidence applies to (e.g., "ovarian cancer")
#   This is used when evidence doesn't match the queried tumor but we know what cancer it's for
CancerSpecificity = str  # Relaxed from Literal to allow specific cancer names


class EvidenceLevel(BaseModel):
    """Tracks the level, scope, and origin of evidence.

    Used for both variant/gene level and cancer type specificity.
    """

    level: VariantGeneLevel | CancerSpecificity | None = Field(
        default=None,
        description="Evidence level: 'variant'/'gene' for variant level, 'cancer_specific'/'pan_cancer' for cancer type"
    )
    scope: Scope | None = Field(
        default=None,
        description="Scope: 'specific' (exact match) or 'unspecified' (broad category) or ambiguous"
    )
    origin: Origin | None = Field(
        default=None,
        description="Origin: 'kb' (knowledge base), 'trial' (clinical trial), or 'inferred'"
    )


class EvidenceItemBase(BaseModel):
    """Base class for all evidence items.

    Provides common metadata fields for tracking evidence provenance:
    - variant_level: Tracks whether evidence is at variant or gene level
    - cancer_type_level: Tracks whether evidence is cancer-specific or pan-cancer
    """

    variant_level: EvidenceLevel | None = Field(
        default=None,
        description="Variant/gene level evidence metadata"
    )
    cancer_type_level: EvidenceLevel | None = Field(
        default=None,
        description="Cancer type specificity evidence metadata"
    )

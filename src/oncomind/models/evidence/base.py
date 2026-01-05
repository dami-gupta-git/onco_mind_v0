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
        description="Scope: 'specific' (exact match) or 'unspecified' or 'ambiguous'"
    )
    origin: Origin | None = Field(
        default=None,
        description="Origin: 'kb' (knowledge base), 'trial' (clinical trial), or 'inferred'"
    )


class EvidenceItemBase(BaseModel):
    """Base class for all evidence items.

    Provides common metadata fields for tracking evidence provenance:
    - locus_match: Tracks whether evidence is at variant or gene level
    - cancer_type_match: Tracks whether evidence is cancer-specific or pan-cancer
    """

    locus_match: EvidenceLevel | None = Field(
        default=None,
        description="Variant/gene level evidence metadata"
    )
    cancer_type_match: EvidenceLevel | None = Field(
        default=None,
        description="Cancer type specificity evidence metadata"
    )

    @property
    def is_tumor_match(self) -> bool | None:
        """Check if evidence matches the queried tumor type.

        Returns:
            True if cancer_specific, False if not cancer_specific, None if unknown.
        """
        if self.cancer_type_match and self.cancer_type_match.level:
            return self.cancer_type_match.level == "cancer_specific"
        return None

    @property
    def cancer_specificity(self) -> str | None:
        """Get the cancer specificity level.

        Returns:
            'cancer_specific' if matches queried tumor,
            'pan_cancer' if tumor-agnostic,
            or the specific cancer name if different tumor type,
            or None if unknown.
        """
        if self.cancer_type_match and self.cancer_type_match.level:
            return self.cancer_type_match.level
        return None

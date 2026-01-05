"""Base class for all evidence items with common metadata fields."""

from typing import Literal
from pydantic import BaseModel, Field

from oncomind.config.constants import BROAD_VARIANTS, PAN_CANCER_TERMS, PAN_CANCER_BIOMARKERS, TUMOR_TYPE_MAPPINGS


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


def is_pan_cancer_term(disease: str | None) -> bool:
    """Check if disease is a generic pan-cancer term or contains a tumor-agnostic biomarker.

    Pan-cancer diseases like "Cancer" or "Solid Tumor" apply broadly
    and should not be considered a specific tumor type match.

    Uses two matching strategies:
    - EXACT match for generic terms (e.g., "Cancer", "Solid Tumor") to avoid
      false positives on strings like "Lung Cancer"
    - SUBSTRING match for FDA tumor-agnostic biomarkers (e.g., "MSI-H", "NTRK")
      to catch FDA indications like "MSI-H Solid Tumors"

    Args:
        disease: Disease string from any knowledge base

    Returns:
        True if disease is a pan-cancer term or contains a tumor-agnostic biomarker
    """
    if not disease:
        return False

    disease_lower = disease.lower().strip()

    # Check for exact match on generic pan-cancer terms
    # (e.g., "Cancer", "Solid Tumor", "Malignant Neoplasm")
    if disease_lower in PAN_CANCER_TERMS:
        return True

    # Check for substring match on FDA tumor-agnostic biomarkers
    # (e.g., "MSI-H Solid Tumors", "NTRK fusion-positive tumors")
    return any(biomarker in disease_lower for biomarker in PAN_CANCER_BIOMARKERS)


def tumor_types_match(
    source_disease: str | None,
    queried_tumor: str | None,
    min_substring_length: int = 3,
) -> bool:
    """Check if two tumor types match via substring or TUMOR_TYPE_MAPPINGS.

    This is the centralized tumor matching logic used across all API clients.
    It handles:
    - Direct substring matches (bidirectional)
    - TUMOR_TYPE_MAPPINGS aliases (e.g., "NSCLC" matches "lung adenocarcinoma")
    - Minimum substring length to avoid false positives from short abbreviations

    Args:
        source_disease: Disease/tumor from the knowledge base (e.g., "Lung Adenocarcinoma")
        queried_tumor: User-provided tumor type (e.g., "NSCLC")
        min_substring_length: Minimum length for direct substring matching (default 3)

    Returns:
        True if tumor types match
    """
    if not source_disease or not queried_tumor:
        return False

    source_lower = source_disease.lower().strip()
    query_lower = queried_tumor.lower().strip()

    # Direct substring match (both directions)
    # Only if both strings are long enough to avoid false positives
    if len(source_lower) >= min_substring_length and len(query_lower) >= min_substring_length:
        if query_lower in source_lower or source_lower in query_lower:
            return True

    # Check via TUMOR_TYPE_MAPPINGS
    for abbrev, aliases in TUMOR_TYPE_MAPPINGS.items():
        # Check if queried tumor matches this mapping category
        query_matches = (
            query_lower == abbrev or
            any(query_lower in alias or alias in query_lower for alias in aliases)
        )
        # Check if source disease matches this mapping category
        source_matches = (
            source_lower == abbrev or
            any(source_lower in alias or alias in source_lower for alias in aliases)
        )
        # If both match the same category, they're related
        if query_matches and source_matches:
            return True

    return False


def compute_cancer_specificity(
    source_disease: str | None,
    queried_tumor: str | None,
) -> str:
    """Determine cancer specificity level for evidence.

    This is the centralized function for computing cancer_type_match.level.
    It determines whether evidence is:
    - "cancer_specific": Evidence matches the user's queried tumor type
    - "pan_cancer": Evidence is tumor-agnostic (e.g., "Solid Tumor", "Cancer")
    - specific cancer name: Evidence is for a different specific tumor type

    Args:
        source_disease: Disease/tumor from the knowledge base
        queried_tumor: User-provided tumor type to match against

    Returns:
        "cancer_specific" if matches queried tumor,
        "pan_cancer" if source is tumor-agnostic,
        or the source disease name if it's a different specific tumor
    """
    # Handle missing source disease
    if not source_disease:
        return "pan_cancer"

    # Check if source is a generic pan-cancer term
    if is_pan_cancer_term(source_disease):
        return "pan_cancer"

    # If no queried tumor, we can't determine specificity
    # Return the source disease as the specific cancer type
    if not queried_tumor:
        return source_disease

    # Check if there's a match
    if tumor_types_match(source_disease, queried_tumor):
        return "cancer_specific"

    # No match - return the source disease as the specific cancer type
    return source_disease


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

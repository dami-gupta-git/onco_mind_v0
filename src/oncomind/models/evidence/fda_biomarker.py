"""FDA biomarker evidence model for parsed FDA label indications.

This model represents structured biomarker-drug-tumor associations parsed
from FDA drug labels using FDALabelParser. It captures:
- Negation patterns (REQUIRED_POSITIVE vs REQUIRED_NEGATIVE)
- Variant specificity levels (variant, codon, gene, pathway)
- Tumor type context
- Combination therapy information
- Line of therapy
"""

from enum import Enum
from typing import Literal

from pydantic import Field

from oncomind.models.evidence.base import EvidenceItemBase


class BiomarkerRequirement(str, Enum):
    """Whether a biomarker is required present or absent for the indication."""
    REQUIRED_POSITIVE = "required_positive"  # Must HAVE the mutation
    REQUIRED_NEGATIVE = "required_negative"  # Must NOT have the mutation (wild-type)
    NOT_SPECIFIED = "not_specified"


class SpecificityLevel(str, Enum):
    """How specific the biomarker requirement is."""
    VARIANT = "variant"  # e.g., "EGFR L858R", "KRAS G12C"
    CODON = "codon"  # e.g., "EGFR exon 19 deletions"
    GENE = "gene"  # e.g., "EGFR mutation-positive"
    PATHWAY = "pathway"  # e.g., "HRR-deficient"


# Match type for variant matching results
MatchType = Literal[
    "exact",  # Exact variant match (e.g., L858R matches L858R)
    "codon",  # Same codon, different variant (e.g., V600K matches V600E approval)
    "gene",  # Gene-level approval (any mutation in gene)
    "excluded",  # This indication EXCLUDES the biomarker (REQUIRED_NEGATIVE)
    "same_codon_different_variant",  # Same codon but not covered (variant-specific approval)
    "different_variant",  # Different variant entirely
    "different_codon",  # Different codon
]


def _extract_codon(variant: str) -> str | None:
    """Extract codon number from variant notation. E.g., V600E -> 600"""
    import re
    match = re.search(r'[A-Z](\d+)[A-Z]?', variant.upper())
    return match.group(1) if match else None


def match_variant_to_indication(
    query_gene: str,
    query_variant: str,
    indication_gene: str,
    indication_requirement: BiomarkerRequirement,
    indication_specificity: SpecificityLevel,
    indication_specified_variants: list[str],
    indication_codon: str | None,
    normalize_variant: bool = False,
) -> dict:
    """Match a query variant against an FDA indication.

    This is the canonical matching logic used by both BiomarkerIndication (dataclass)
    and FDABiomarkerEvidence (Pydantic model).

    Args:
        query_gene: Gene symbol to match (e.g., "EGFR")
        query_variant: Variant notation to match (e.g., "L858R", "p.Leu858Arg")
        indication_gene: Gene from FDA indication
        indication_requirement: Whether biomarker must be present or absent
        indication_specificity: Specificity level (gene, codon, variant)
        indication_specified_variants: List of approved variants
        indication_codon: Codon number if specified
        normalize_variant: If True, normalize query_variant using to_short_form

    Returns:
        Dict with keys:
        - matches: bool - whether the variant is approved
        - match_type: MatchType - type of match
        - reason: str - explanation
    """
    query_gene = query_gene.upper()

    # Optionally normalize query variant to short form (e.g., "p.Leu858Arg" -> "L858R")
    if normalize_variant:
        from oncomind.utils.variant_normalization import to_short_form
        normalized = to_short_form(query_variant)
        query_variant = normalized.upper() if normalized else query_variant.upper()
    else:
        query_variant = query_variant.upper()

    if indication_gene.upper() != query_gene:
        return {"matches": False, "match_type": None, "reason": "Different gene"}

    # If this indication EXCLUDES the biomarker, it's a negative match
    if indication_requirement == BiomarkerRequirement.REQUIRED_NEGATIVE:
        return {
            "matches": False,
            "match_type": "excluded",
            "reason": f"This indication is for patients WITHOUT {indication_gene} mutations"
        }

    # Check specificity levels
    if indication_specificity == SpecificityLevel.GENE:
        return {
            "matches": True,
            "match_type": "gene",
            "reason": f"Gene-level approval: any {indication_gene} mutation"
        }

    if indication_specificity == SpecificityLevel.VARIANT:
        # Check for exact variant match
        if query_variant in [v.upper() for v in indication_specified_variants]:
            return {
                "matches": True,
                "match_type": "exact",
                "reason": f"Exact variant match: {query_variant}"
            }

        # Check for same codon
        query_codon = _extract_codon(query_variant)
        # Use indication_codon if set, otherwise extract from specified_variants
        codon = indication_codon
        if not codon and indication_specified_variants:
            for v in indication_specified_variants:
                codon = _extract_codon(v)
                if codon:
                    break
        if codon and query_codon == codon:
            return {
                "matches": False,
                "match_type": "same_codon_different_variant",
                "reason": f"Same codon ({codon}) but different variant. Approved: {indication_specified_variants}, Query: {query_variant}"
            }

        return {
            "matches": False,
            "match_type": "different_variant",
            "reason": f"Different variant. Approved: {indication_specified_variants}, Query: {query_variant}"
        }

    if indication_specificity == SpecificityLevel.CODON:
        query_codon = _extract_codon(query_variant)
        if indication_codon and query_codon == indication_codon:
            return {
                "matches": True,
                "match_type": "codon",
                "reason": f"Codon-level match: codon {indication_codon}"
            }
        return {
            "matches": False,
            "match_type": "different_codon",
            "reason": f"Different codon. Approved: {indication_codon}, Query: {query_codon}"
        }

    return {"matches": False, "match_type": None, "reason": "Unknown specificity"}


class FDABiomarkerEvidence(EvidenceItemBase):
    """FDA biomarker-drug indication parsed from drug label.

    This model represents a single biomarker-drug-tumor association extracted
    from an FDA drug label. This comes directly from parsing FDA label text
    sources, this comes directly from parsing FDA label text.

    Inherits from EvidenceItemBase which provides:
    - locus_variant_match: EvidenceLevel tracking variant/gene level match
    - cancer_type_match: EvidenceLevel tracking cancer specificity
    - locus_match: Computed field returning 'variant', 'codon', or 'gene'
    - tumor_match: Computed field returning True/False/None for cancer match
    """

    # Drug identification
    drug_name: str = Field(..., description="Generic drug name")
    brand_name: str | None = Field(default=None, description="Brand name")

    # FDA label identifiers for linking
    set_id: str | None = Field(default=None, description="FDA label set_id (UUID)")
    spl_version: str | None = Field(default=None, description="FDA label SPL version")

    # Biomarker details
    gene: str = Field(..., description="Gene symbol (e.g., EGFR, BRAF)")
    requirement: BiomarkerRequirement = Field(
        ...,
        description="Whether biomarker must be present or absent"
    )
    specificity: SpecificityLevel = Field(
        ...,
        description="Specificity level of the biomarker requirement"
    )

    # Variant details (populated when specificity is VARIANT or CODON)
    specified_variants: list[str] = Field(
        default_factory=list,
        description="Specific variants covered (e.g., ['L858R', 'exon 19 del'])"
    )
    codon: str | None = Field(
        default=None,
        description="Codon number if relevant (e.g., '600' for V600E)"
    )

    # Tumor context
    tumor_types: list[str] = Field(
        default_factory=list,
        description="Tumor types covered by this indication"
    )
    tumor_stage: str | None = Field(
        default=None,
        description="Tumor stage/setting (e.g., 'metastatic', 'locally advanced')"
    )

    # Combination therapy
    combination_partners: list[str] = Field(
        default_factory=list,
        description="Partner drugs for combination therapy"
    )
    is_monotherapy: bool = Field(
        default=True,
        description="True if drug is used as monotherapy"
    )

    # Line of therapy
    line_of_therapy: str | None = Field(
        default=None,
        description="Line of therapy (e.g., 'first-line', 'post-progression')"
    )

    # Source text
    indication_text: str = Field(
        default="",
        description="Original indication text from FDA label (truncated)"
    )
    section: str = Field(
        default="indications_and_usage",
        description="FDA label section this was parsed from"
    )

    # Match result (populated when matching against a query variant)
    variant_match_result: MatchType | None = Field(
        default=None,
        description="Result of matching query variant against this indication"
    )
    variant_match_reason: str | None = Field(
        default=None,
        description="Explanation of the match result"
    )

    def matches_variant(self, query_gene: str, query_variant: str) -> dict:
        """Check if a user's variant query matches this indication.

        Delegates to match_variant_to_indication() with variant normalization enabled.

        Args:
            query_gene: Gene symbol to match (e.g., "EGFR")
            query_variant: Variant notation to match (e.g., "L858R", "p.Leu858Arg")

        Returns:
            Dict with keys:
            - matches: bool - whether the variant is approved
            - match_type: MatchType - type of match
            - reason: str - explanation
        """
        return match_variant_to_indication(
            query_gene=query_gene,
            query_variant=query_variant,
            indication_gene=self.gene,
            indication_requirement=self.requirement,
            indication_specificity=self.specificity,
            indication_specified_variants=self.specified_variants,
            indication_codon=self.codon,
            normalize_variant=True,  # Pydantic model normalizes variants
        )

    @staticmethod
    def _extract_codon(variant: str) -> str | None:
        """Extract codon number from variant notation. E.g., V600E -> 600"""
        import re
        match = re.search(r'[A-Z](\d+)[A-Z]?', variant.upper())
        return match.group(1) if match else None

    @classmethod
    def from_biomarker_indication(
        cls,
        indication: "BiomarkerIndication",
        set_id: str | None = None,
        spl_version: str | None = None,
    ) -> "FDABiomarkerEvidence":
        """Create FDABiomarkerEvidence from a parsed BiomarkerIndication.

        This is the factory method to convert from the dataclass used by
        FDALabelParser to this Pydantic model.

        Args:
            indication: BiomarkerIndication from fda_label_parser.py
            set_id: FDA label set_id for linking
            spl_version: FDA label version

        Returns:
            FDABiomarkerEvidence instance
        """
        # Import here to avoid circular imports
        from oncomind.api.fda_label_parser import (
            BiomarkerIndication as ParserIndication,
            BiomarkerRequirement as ParserRequirement,
            SpecificityLevel as ParserSpecificity,
        )

        # Map parser enums to model enums
        requirement_map = {
            ParserRequirement.REQUIRED_POSITIVE: BiomarkerRequirement.REQUIRED_POSITIVE,
            ParserRequirement.REQUIRED_NEGATIVE: BiomarkerRequirement.REQUIRED_NEGATIVE,
            ParserRequirement.NOT_SPECIFIED: BiomarkerRequirement.NOT_SPECIFIED,
        }

        specificity_map = {
            ParserSpecificity.VARIANT: SpecificityLevel.VARIANT,
            ParserSpecificity.CODON: SpecificityLevel.CODON,
            ParserSpecificity.GENE: SpecificityLevel.GENE,
            ParserSpecificity.PATHWAY: SpecificityLevel.PATHWAY,
        }

        return cls(
            drug_name=indication.drug_name,
            brand_name=indication.brand_name,
            set_id=set_id,
            spl_version=spl_version,
            gene=indication.gene,
            requirement=requirement_map[indication.requirement],
            specificity=specificity_map[indication.specificity],
            specified_variants=indication.specified_variants,
            codon=indication.codon,
            tumor_types=indication.tumor_types,
            tumor_stage=indication.tumor_stage,
            combination_partners=indication.combination_partners,
            is_monotherapy=indication.is_monotherapy,
            line_of_therapy=indication.line_of_therapy,
            indication_text=indication.indication_text,
            section=indication.section,
        )

    @property
    def is_exclusion(self) -> bool:
        """Check if this indication excludes the biomarker (wild-type requirement)."""
        return self.requirement == BiomarkerRequirement.REQUIRED_NEGATIVE

    @property
    def fda_label_url(self) -> str | None:
        """Generate FDA DailyMed URL for this drug label."""
        if self.set_id:
            return f"https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={self.set_id}"
        return None

    def get_display_drug_name(self) -> str:
        """Get drug name for display, including combination partners."""
        if self.combination_partners:
            partners = " + ".join(self.combination_partners)
            return f"{self.drug_name} + {partners}"
        return self.drug_name

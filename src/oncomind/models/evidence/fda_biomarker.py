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


class FDABiomarkerEvidence(EvidenceItemBase):
    """FDA biomarker-drug indication parsed from drug label.

    This model represents a single biomarker-drug-tumor association extracted
    from an FDA drug label. Unlike FDAApproval which is derived from multiple
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

        This mirrors the matches_variant method from BiomarkerIndication
        but operates on the Pydantic model.

        Args:
            query_gene: Gene symbol to match (e.g., "EGFR")
            query_variant: Variant notation to match (e.g., "L858R")

        Returns:
            Dict with keys:
            - matches: bool - whether the variant is approved
            - match_type: MatchType - type of match
            - reason: str - explanation
        """
        import re

        query_gene = query_gene.upper()
        query_variant = query_variant.upper()

        if self.gene.upper() != query_gene:
            return {"matches": False, "match_type": None, "reason": "Different gene"}

        # If this indication EXCLUDES the biomarker, it's a negative match
        if self.requirement == BiomarkerRequirement.REQUIRED_NEGATIVE:
            return {
                "matches": False,
                "match_type": "excluded",
                "reason": f"This indication is for patients WITHOUT {self.gene} mutations"
            }

        # Check specificity levels
        if self.specificity == SpecificityLevel.GENE:
            return {
                "matches": True,
                "match_type": "gene",
                "reason": f"Gene-level approval: any {self.gene} mutation"
            }

        if self.specificity == SpecificityLevel.VARIANT:
            # Check for exact variant match
            if query_variant in [v.upper() for v in self.specified_variants]:
                return {
                    "matches": True,
                    "match_type": "exact",
                    "reason": f"Exact variant match: {query_variant}"
                }

            # Check for same codon
            query_codon = self._extract_codon(query_variant)
            if self.codon and query_codon == self.codon:
                return {
                    "matches": False,
                    "match_type": "same_codon_different_variant",
                    "reason": f"Same codon ({self.codon}) but different variant. Approved: {self.specified_variants}, Query: {query_variant}"
                }

            return {
                "matches": False,
                "match_type": "different_variant",
                "reason": f"Different variant. Approved: {self.specified_variants}, Query: {query_variant}"
            }

        if self.specificity == SpecificityLevel.CODON:
            query_codon = self._extract_codon(query_variant)
            if self.codon and query_codon == self.codon:
                return {
                    "matches": True,
                    "match_type": "codon",
                    "reason": f"Codon-level match: codon {self.codon}"
                }
            return {
                "matches": False,
                "match_type": "different_codon",
                "reason": f"Different codon. Approved: {self.codon}, Query: {query_codon}"
            }

        return {"matches": False, "match_type": None, "reason": "Unknown specificity"}

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

from pydantic import Field, computed_field

from oncomind.models.evidence.base import EvidenceItemBase


class VICCEvidence(EvidenceItemBase):
    """Evidence from VICC MetaKB (harmonized multi-KB interpretations)."""

    description: str | None = None
    gene: str | None = None
    variant: str | None = None
    disease: str | None = None
    drugs: list[str] = Field(default_factory=list)
    evidence_level: str | None = None
    response_type: str | None = None
    source: str | None = None
    publication_url: str | list[str] | None = None
    oncogenic: str | None = None
    is_sensitivity: bool = False
    is_resistance: bool = False
    oncokb_level: str | None = None
    # Additional fields for molecular profile context
    molecular_profile: str | None = None
    molecular_profile_score: float | None = None
    # Additional match tracking fields (not in base class)
    matched_profile: str | None = Field(
        default=None,
        description="The molecular profile that was actually matched"
    )

    @computed_field
    @property
    def match_level(self) -> str:
        """Get the locus match level: 'variant', 'codon', or 'gene'.

        Uses variant_level.level from EvidenceItemBase for consistency
        with ClinicalTrialEvidence, FDAApproval, CIViC, and CGI models.
        """
        if self.variant_level and self.variant_level.level:
            return self.variant_level.level
        return "gene"

    @property
    def is_tumor_match(self) -> bool | None:
        """Check if evidence matches the queried tumor type.

        Uses cancer_type_level from EvidenceItemBase for consistency
        with other therapeutic evidence models.

        Returns:
            True if cancer_specific, False if not, None if unknown.
        """
        if self.cancer_type_level and self.cancer_type_level.level:
            return self.cancer_type_level.level == "cancer_specific"
        return None

    @property
    def cancer_specificity(self) -> str | None:
        """Get the cancer specificity level.

        Returns:
            'cancer_specific' if matches queried tumor,
            'pan_cancer' if tumor-agnostic,
            or the specific cancer name if different tumor type.
        """
        if self.cancer_type_level and self.cancer_type_level.level:
            return self.cancer_type_level.level
        return None

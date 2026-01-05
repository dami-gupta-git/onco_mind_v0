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

        Uses locus_match.level from EvidenceItemBase for consistency
        with ClinicalTrialEvidence, FDAApproval, CIViC, and CGI models.
        """
        if self.locus_match and self.locus_match.level:
            return self.locus_match.level
        return "gene"

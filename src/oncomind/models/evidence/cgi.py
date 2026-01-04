from pydantic import Field, computed_field

from oncomind.models.evidence.base import EvidenceItemBase


class CGIBiomarkerEvidence(EvidenceItemBase):
    """Evidence from Cancer Genome Interpreter biomarkers database."""

    gene: str | None = None
    alteration: str | None = None
    drug: str | None = None
    drug_status: str | None = None
    association: str | None = None
    evidence_level: str | None = None
    source: str | None = None
    tumor_type: str | None = None
    fda_approved: bool = False
    # Additional match tracking fields (not in base class)
    matched_alteration: str | None = Field(
        default=None,
        description="The alteration that was actually matched"
    )

    @computed_field
    @property
    def match_level(self) -> str:
        """Get the locus match level: 'variant', 'codon', or 'gene'.

        Uses variant_level.level from EvidenceItemBase for consistency
        with ClinicalTrialEvidence, FDAApproval, and CIViC models.
        """
        if self.variant_level and self.variant_level.level:
            return self.variant_level.level
        return "gene"

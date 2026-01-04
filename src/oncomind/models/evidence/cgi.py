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

    @property
    def is_tumor_match(self) -> bool | None:
        """Check if biomarker matches the queried tumor type.

        Uses cancer_type_level from EvidenceItemBase for consistency
        with ClinicalTrialEvidence, FDAApproval, and CIViC models.

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

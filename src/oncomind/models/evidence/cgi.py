from pydantic import Field

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
    # Match specificity tracking
    match_level: str | None = Field(
        default=None,
        description="Level of match: 'variant' (exact), 'codon' (same position), 'gene' (gene-only)"
    )
    matched_alteration: str | None = Field(
        default=None,
        description="The alteration that was actually matched"
    )

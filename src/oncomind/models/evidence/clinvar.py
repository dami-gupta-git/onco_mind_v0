from pydantic import Field

from oncomind.models.evidence.base import EvidenceItemBase


class ClinVarEvidence(EvidenceItemBase):
    """Evidence from ClinVar."""

    clinical_significance: str | None = None
    review_status: str | None = None
    conditions: list[str] = Field(default_factory=list)
    last_evaluated: str | None = None
    variation_id: str | None = None
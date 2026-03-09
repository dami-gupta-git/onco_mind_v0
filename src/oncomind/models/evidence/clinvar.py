from pydantic import Field

from oncomind.models.evidence.base import EvidenceItemBase


class ClinVarEvidence(EvidenceItemBase):
    """Evidence from ClinVar."""

    clinical_significance: str | None = None
    review_status: str | None = None
    conditions: list[str] = Field(default_factory=list)
    last_evaluated: str | None = None
    variation_id: str | None = None

    def get_url(self) -> str | None:
        """Get the ClinVar variation URL."""
        if self.variation_id:
            return (
                f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{self.variation_id}/"
            )
        return None

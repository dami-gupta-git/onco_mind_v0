from pydantic import Field, computed_field

from oncomind.models.evidence.base import EvidenceItemBase


class CIViCEvidence(EvidenceItemBase):
    """Evidence from CIViC (Clinical Interpretations of Variants in Cancer).

    Each evidence item has a unique EID (Evidence Item ID) in the format "EID{number}".
    Example: EID5586
    """

    evidence_id: int | None = Field(default=None, description="CIViC evidence item ID (numeric)")
    evidence_type: str | None = None
    evidence_level: str | None = None
    evidence_direction: str | None = None
    clinical_significance: str | None = None
    disease: str | None = None
    drugs: list[str] = Field(default_factory=list)
    description: str | None = None
    source: str | None = None
    rating: int | None = None
    # Additional fields for provenance
    pmid: str | None = None
    source_url: str | None = None
    trust_rating: int | None = None  # 1-5 star rating

    @computed_field
    @property
    def eid(self) -> str | None:
        """Formatted Evidence Item ID (e.g., 'EID5586')."""
        if self.evidence_id is not None:
            return f"EID{self.evidence_id}"
        return None

    @computed_field
    @property
    def civic_url(self) -> str | None:
        """URL to the CIViC evidence item page."""
        if self.evidence_id is not None:
            return f"https://civicdb.org/evidence/{self.evidence_id}/summary"
        return None





class CIViCAssertionEvidence(EvidenceItemBase):
    """Evidence from CIViC Assertions (curated AMP/ASCO/CAP classifications).

    Each assertion has a unique AID (Assertion ID) in the format "AID{number}".
    Example: AID20
    """

    assertion_id: int | None = Field(default=None, description="CIViC assertion ID (numeric)")
    name: str | None = None
    amp_level: str | None = None
    amp_tier: str | None = None
    amp_level_letter: str | None = None
    assertion_type: str | None = None
    significance: str | None = None
    status: str | None = None
    molecular_profile: str | None = None
    disease: str | None = None
    therapies: list[str] = Field(default_factory=list)
    fda_companion_test: bool | None = None
    nccn_guideline: str | None = None
    description: str | None = None
    is_sensitivity: bool = False
    is_resistance: bool = False

    @computed_field
    @property
    def aid(self) -> str | None:
        """Formatted Assertion ID (e.g., 'AID20')."""
        if self.assertion_id is not None:
            return f"AID{self.assertion_id}"
        return None

    @computed_field
    @property
    def civic_url(self) -> str | None:
        """URL to the CIViC assertion page."""
        if self.assertion_id is not None:
            return f"https://civicdb.org/assertions/{self.assertion_id}/summary"
        return None


from pydantic import Field, computed_field

from oncomind.models.evidence.base import EvidenceItemBase


class CIViCEvidence(EvidenceItemBase):
    """Evidence from CIViC (Clinical Interpretations of Variants in Cancer).

    Each evidence item has a unique EID (Evidence Item ID) in the format "EID{number}".
    Example: EID5586
    """

    evidence_id: int | None = Field(
        default=None, description="CIViC evidence item ID (numeric)"
    )
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
    # Additional match tracking fields (not in base class)
    matched_profile: str | None = Field(
        default=None,
        description="The molecular profile that was actually matched (e.g., 'EGFR L858R')",
    )

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

    @computed_field
    @property
    def is_sensitivity(self) -> bool:
        """Check if this evidence SUPPORTS sensitivity/response.

        CIViC has two fields:
        - clinical_significance: What the evidence is ABOUT (e.g., "Sensitivity/Response")
        - evidence_direction: Whether the study SUPPORTS or DOES_NOT_SUPPORT that significance

        True sensitivity requires:
        - clinical_significance contains "sensitiv" or "response"
        - AND evidence_direction is "SUPPORTS" (or missing, for backward compatibility)
        """
        if not self.clinical_significance:
            return False
        significance = self.clinical_significance.upper()
        direction = (self.evidence_direction or "").upper()

        if "SENSITIV" in significance or "RESPONSE" in significance:
            # If direction is missing/empty, assume SUPPORTS for backward compatibility
            if not direction or direction == "SUPPORTS":
                return True
            else:
                return False
        return False

    @computed_field
    @property
    def is_resistance(self) -> bool:
        """Check if this evidence SUPPORTS resistance.

        CIViC has two fields:
        - clinical_significance: What the evidence is ABOUT (e.g., "Resistance")
        - evidence_direction: Whether the study SUPPORTS or DOES_NOT_SUPPORT that significance

        True resistance requires:
        - clinical_significance contains "resist"
        - AND evidence_direction is "SUPPORTS" (or missing, for backward compatibility)
        """
        if not self.clinical_significance:
            return False
        sig_upper = self.clinical_significance.upper()
        direction = (self.evidence_direction or "").upper()

        is_resist_topic = "RESIST" in sig_upper
        # If direction is missing/empty, assume SUPPORTS for backward compatibility
        if is_resist_topic and (not direction or direction == "SUPPORTS"):
            return True

        return False


class CIViCAssertionEvidence(EvidenceItemBase):
    """Evidence from CIViC Assertions (curated AMP/ASCO/CAP classifications).

    Each assertion has a unique AID (Assertion ID) in the format "AID{number}".
    Example: AID20
    """

    assertion_id: int | None = Field(
        default=None, description="CIViC assertion ID (numeric)"
    )
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
    # Additional match tracking fields (not in base class)
    matched_profile: str | None = Field(
        default=None,
        description="The molecular profile that was actually matched (e.g., 'EGFR L858R')",
    )

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

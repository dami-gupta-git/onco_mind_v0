from pydantic import BaseModel, Field


class RecommendedTherapy(BaseModel):
    """Recommended therapy based on variant."""

    drug_name: str = Field(..., description="Name of the therapeutic agent")
    evidence_level: str | None = Field(None, description="Level of supporting evidence")
    approval_status: str | None = Field(None, description="FDA approval status for this indication")
    clinical_context: str | None = Field(
        None, description="Clinical context (e.g., first-line, resistant)"
    )

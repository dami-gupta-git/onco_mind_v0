"""Variant insight and annotation models."""

from pydantic import BaseModel, Field

from oncomind.models.annotations import VariantAnnotations


class RecommendedTherapy(BaseModel):
    """Recommended therapy based on variant."""

    drug_name: str = Field(..., description="Name of the therapeutic agent")
    evidence_level: str | None = Field(None, description="Level of supporting evidence")
    approval_status: str | None = Field(None, description="FDA approval status for this indication")
    clinical_context: str | None = Field(
        None, description="Clinical context (e.g., first-line, resistant)"
    )


class LLMInsight(VariantAnnotations):
    """LLM-generated narrative insight for a variant.

    This contains the LLM's synthesized interpretation of the evidence,
    including clinical summary, recommended therapies, and rationale.
    """

    gene: str
    variant: str
    tumor_type: str | None
    llm_summary: str = Field(..., description="LLM-generated narrative summary of the variant")
    recommended_therapies: list[RecommendedTherapy] = Field(default_factory=list)
    rationale: str = Field(..., description="Detailed rationale for clinical interpretation")
    evidence_strength: str | None = Field(
        None, description="Overall strength of evidence (Strong/Moderate/Weak)"
    )
    clinical_trials_available: bool = Field(
        default=False, description="Whether relevant clinical trials exist"
    )
    references: list[str] = Field(
        default_factory=list, description="Key references supporting the insight"
    )

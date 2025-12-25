"""Variant insight and annotation models."""

from pydantic import BaseModel, Field

from oncomind.models.recommended_therapies import RecommendedTherapy


class LLMInsight(BaseModel):
    """LLM-generated narrative insight for a variant.

    This contains the LLM's synthesized interpretation of the evidence,
    including clinical summary, recommended therapies, and rationale.

    This is embedded within the parent Insight object, which contains
    all the structured evidence and annotations.
    """

    llm_summary: str = Field(..., description="LLM-generated narrative summary of the variant")
    recommended_therapies: list[RecommendedTherapy] = Field(default_factory=list)
    rationale: str = Field(..., description="Detailed rationale for clinical interpretation")
    clinical_trials_available: bool = Field(
        default=False, description="Whether relevant clinical trials exist"
    )
    references: list[str] = Field(
        default_factory=list, description="Key references supporting the insight"
    )

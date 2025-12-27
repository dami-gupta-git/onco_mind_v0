"""Variant insight and annotation models."""

from pydantic import BaseModel, Field

from oncomind.models.therapeutic_evidence import TherapeuticEvidence


class LLMInsight(BaseModel):
    """LLM-generated narrative insight for a variant.

    This contains the LLM's synthesized interpretation of the evidence,
    including clinical summary, recommended therapies, and rationale.

    This is embedded within the parent Insight object, which contains
    all the structured evidence and annotations.
    """

    llm_summary: str = Field(..., description="LLM-generated narrative summary of the variant")
    rationale: str = Field("", description="Reasoning behind the summary")
    clinical_trials_available: bool = Field(False, description="Whether trials exist")

    # Raw component fields (UI layer handles formatting)
    functional_summary: str | None = Field(None, description="Functional impact of the variant")
    biological_context: str | None = Field(None, description="Biological context and mechanism")
    therapeutic_landscape: dict | None = Field(
        None,
        description="Therapeutic landscape: fda_approved, clinical_evidence, preclinical, resistance_mechanisms"
    )

    # Changed from recommended_therapies to therapeutic_evidence
    therapeutic_evidence: list[TherapeuticEvidence] = Field(
        default_factory=list,
        description="Therapeutic evidence at all levels"
    )

    references: list[str] = Field(default_factory=list, description="Key references")

    # New research-focused fields
    evidence_quality: str | None = Field(
        None,
        description="Overall evidence quality assessment (comprehensive/moderate/limited/minimal)"
    )
    knowledge_gaps: list[str] = Field(
        default_factory=list,
        description="Identified gaps in knowledge"
    )
    well_characterized: list[str] = Field(
        default_factory=list,
        description="Aspects with strong, consistent evidence"
    )
    conflicting_evidence: list[str] = Field(
        default_factory=list,
        description="Areas where sources disagree or suggest different interpretations"
    )
    research_implications: str | None = Field(
        None,
        description="Implications for future research"
    )
    evidence_tags: list[str] = Field(
        default_factory=list,
        description="Labels indicating evidence types (e.g., 'direct clinical data', 'preclinical only')"
    )
    research_hypotheses: list[str] = Field(
        default_factory=list,
        description="Testable research hypotheses generated from evidence gaps"
    )

    # Backwards compatibility
    @property
    def recommended_therapies(self) -> list[TherapeuticEvidence]:
        """Backwards-compatible alias."""
        return self.therapeutic_evidence

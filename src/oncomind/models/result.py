"""Result - Combined evidence and LLM insight container.

This is the output when LLM synthesis is enabled, combining:
- Evidence: Structured data from databases and APIs
- LLMInsight: LLM-generated narrative and recommendations
"""

from pydantic import BaseModel, Field

from oncomind.models.evidence import Evidence
from oncomind.models.llm_insight import LLMInsight


class Result(BaseModel):
    """Combined evidence and LLM insight for a variant.

    This is returned by Conductor when LLM mode is enabled.
    It wraps the structured Evidence with optional LLM narrative.

    Example:
        >>> result = await conductor.run("BRAF V600E", tumor_type="Melanoma")
        >>> print(result.evidence.identifiers.gene)
        BRAF
        >>> if result.llm:
        ...     print(result.llm.llm_summary)
    """

    evidence: Evidence = Field(..., description="Structured evidence from databases and APIs")
    llm: LLMInsight | None = Field(
        default=None, description="LLM-generated narrative insight"
    )

    @property
    def identifiers(self):
        """Shortcut to evidence.identifiers."""
        return self.evidence.identifiers

    @property
    def clinical(self):
        """Shortcut to evidence.clinical."""
        return self.evidence.clinical

    @property
    def kb(self):
        """Shortcut to evidence.kb."""
        return self.evidence.kb

    @property
    def functional(self):
        """Shortcut to evidence.functional."""
        return self.evidence.functional

    @property
    def literature(self):
        """Shortcut to evidence.literature."""
        return self.evidence.literature

    def has_evidence(self) -> bool:
        """Check if any evidence was found."""
        return self.evidence.has_evidence()

    def get_summary(self) -> str:
        """Generate a summary of the variant."""
        return self.evidence.get_summary()

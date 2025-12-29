"""Evidence gaps detection for research context.

Identifies what's unknown or understudied about a variant,
helping researchers identify opportunities for investigation.
"""

from pydantic import BaseModel, Field
from enum import Enum


class GapSeverity(str, Enum):
    """How significant is this evidence gap?"""
    CRITICAL = "critical"      # No data at all in key area
    SIGNIFICANT = "significant"  # Limited data, needs more research
    MINOR = "minor"            # Some data exists but could be deeper


class GapCategory(str, Enum):
    """Category of evidence gap."""
    FUNCTIONAL = "functional"           # Mechanism unknown
    CLINICAL = "clinical"               # No clinical trials/outcomes
    TUMOR_TYPE = "tumor_type"           # Not studied in this tumor type
    DRUG_RESPONSE = "drug_response"     # No drug sensitivity data
    RESISTANCE = "resistance"           # Resistance mechanisms unknown
    PRECLINICAL = "preclinical"         # No cell line/model data
    PREVALENCE = "prevalence"           # Frequency unknown
    PROGNOSTIC = "prognostic"           # Survival impact unknown
    DISCORDANT = "discordant"           # Conflicting evidence between sources
    VALIDATION = "validation"           # Strong oncogenic signal but lacks therapeutic validation


class CharacterizedAspect(BaseModel):
    """A well-characterized aspect of the variant with its basis."""

    aspect: str = Field(..., description="What is well characterized (e.g., 'clinical actionability')")
    basis: str = Field(..., description="Why we think so (e.g., 'FDA-approved therapies exist')")
    category: GapCategory | None = Field(None, description="Category this aspect belongs to (for grouping)")
    matches_on: str | None = Field(None, description="Match level breakdown (e.g., '2 variant, 1 gene')")


class EvidenceGap(BaseModel):
    """A specific gap in the evidence for a variant."""

    category: GapCategory
    severity: GapSeverity
    description: str = Field(..., description="Human-readable description of the gap")

    # What would fill this gap?
    suggested_studies: list[str] = Field(
        default_factory=list,
        description="Types of studies that would address this gap"
    )

    # Is this gap addressable with existing tools/data?
    addressable_with: list[str] = Field(
        default_factory=list,
        description="Databases or methods that could fill this gap"
    )


class EvidenceGaps(BaseModel):
    """Collection of evidence gaps for a variant."""

    gaps: list[EvidenceGap] = Field(default_factory=list)

    # Overall assessment
    overall_evidence_quality: str = Field(
        "unknown",
        description="comprehensive | moderate | limited | minimal | none"
    )

    # What's well-characterized vs not
    well_characterized: list[str] = Field(
        default_factory=list,
        description="Aspects with strong evidence (legacy: simple strings)"
    )

    # Structured version with basis/reasoning
    well_characterized_detailed: list[CharacterizedAspect] = Field(
        default_factory=list,
        description="Aspects with strong evidence, including basis for each"
    )

    poorly_characterized: list[str] = Field(
        default_factory=list,
        description="Aspects needing more research"
    )

    # Research priority
    research_priority: str = Field(
        "unknown",
        description="high | medium | low - based on gap severity and variant importance"
    )

    def has_critical_gaps(self) -> bool:
        """Check if there are any critical gaps."""
        return any(g.severity == GapSeverity.CRITICAL for g in self.gaps)

    def get_gaps_by_category(self, category: GapCategory) -> list[EvidenceGap]:
        """Get gaps filtered by category."""
        return [g for g in self.gaps if g.category == category]

    def get_gaps_by_severity(self, severity: GapSeverity) -> list[EvidenceGap]:
        """Get gaps filtered by severity."""
        return [g for g in self.gaps if g.severity == severity]

    def to_summary(self) -> str:
        """Generate human-readable summary of evidence gaps."""
        if not self.gaps:
            return "No significant evidence gaps identified."

        lines = [f"Evidence Quality: {self.overall_evidence_quality.upper()}"]
        lines.append(f"Research Priority: {self.research_priority.upper()}")

        if self.well_characterized:
            lines.append(f"\nWell characterized: {', '.join(self.well_characterized)}")

        if self.poorly_characterized:
            lines.append(f"Needs research: {', '.join(self.poorly_characterized)}")

        critical = self.get_gaps_by_severity(GapSeverity.CRITICAL)
        if critical:
            lines.append(f"\nCritical gaps ({len(critical)}):")
            for g in critical:
                lines.append(f"  • {g.description}")

        significant = self.get_gaps_by_severity(GapSeverity.SIGNIFICANT)
        if significant:
            lines.append(f"\nSignificant gaps ({len(significant)}):")
            for g in significant:
                lines.append(f"  • {g.description}")

        return "\n".join(lines)

    def to_dict_for_llm(self) -> dict:
        """Convert to dict optimized for LLM prompt."""
        return {
            "overall_quality": self.overall_evidence_quality,
            "research_priority": self.research_priority,
            "well_characterized": self.well_characterized,
            "knowledge_gaps": self.poorly_characterized,  # Alias for prompt compatibility
            "conflicting_evidence": [
                g.description for g in self.get_gaps_by_category(GapCategory.DISCORDANT)
            ],
            "critical_gaps": [
                {"description": g.description, "suggested_studies": g.suggested_studies}
                for g in self.get_gaps_by_severity(GapSeverity.CRITICAL)
            ],
            "significant_gaps": [
                {"description": g.description, "suggested_studies": g.suggested_studies}
                for g in self.get_gaps_by_severity(GapSeverity.SIGNIFICANT)
            ],
        }

    def top_gaps(self, n: int = 3) -> list[EvidenceGap]:
        """Get the top N most important gaps, sorted by severity.

        Args:
            n: Maximum number of gaps to return (default 3)

        Returns:
            List of up to N gaps, prioritized by severity (CRITICAL > SIGNIFICANT > MINOR)
        """
        # Define severity order
        severity_order = {
            GapSeverity.CRITICAL: 0,
            GapSeverity.SIGNIFICANT: 1,
            GapSeverity.MINOR: 2,
        }

        # Sort by severity (critical first)
        sorted_gaps = sorted(
            self.gaps,
            key=lambda g: severity_order.get(g.severity, 99)
        )

        return sorted_gaps[:n]

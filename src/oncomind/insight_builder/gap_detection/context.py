"""Gap detection context for accumulating results.

Provides GapDetectionContext dataclass for tracking gaps, well-characterized
aspects, and flags during gap detection analysis.
"""

from dataclasses import dataclass, field

from oncomind.models.extracted.evidence_gaps import (
    EvidenceGap, GapCategory, GapSeverity, CharacterizedAspect
)


@dataclass
class GapDetectionContext:
    """Context accumulated during gap detection.

    Tracks gaps found, well-characterized aspects, and flags for use
    across multiple detection functions.
    """
    gene: str
    variant: str
    tumor_type: str | None
    is_cancer_gene: bool
    has_pathogenic_signal: bool

    # Accumulated results
    gaps: list[EvidenceGap] = field(default_factory=list)
    well_characterized: list[str] = field(default_factory=list)
    well_characterized_detailed: list[CharacterizedAspect] = field(default_factory=list)
    poorly_characterized: list[str] = field(default_factory=list)

    # Flags set during detection (used by later checks)
    has_clinical: bool = False
    has_drug_data: bool = False

    def add_well_characterized(
        self,
        aspect: str,
        basis: str,
        category: GapCategory | None = None,
        matches_on: str | None = None,
        tumor_match: str | None = None,
    ) -> None:
        """Add a well-characterized aspect with its basis and category."""
        # Use title case for aspect
        aspect_title = aspect.title()
        self.well_characterized.append(aspect_title)
        self.well_characterized_detailed.append(
            CharacterizedAspect(
                aspect=aspect_title,
                basis=basis,
                category=category,
                matches_on=matches_on,
                tumor_match=tumor_match,
            )
        )

    def add_gap(
        self,
        category: GapCategory,
        severity: GapSeverity,
        description: str,
        suggested_studies: list[str] | None = None,
        addressable_with: list[str] | None = None,
    ) -> None:
        """Add an evidence gap."""
        self.gaps.append(EvidenceGap(
            category=category,
            severity=severity,
            description=description,
            suggested_studies=suggested_studies or [],
            addressable_with=addressable_with or [],
        ))

    def add_poorly_characterized(self, aspect: str) -> None:
        """Add a poorly-characterized aspect."""
        self.poorly_characterized.append(aspect)

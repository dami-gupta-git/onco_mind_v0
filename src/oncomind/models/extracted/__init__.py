"""Extracted models - unified/derived from evidence items.

These models are created by extracting and unifying data from raw evidence items.
They are not evidence sources themselves.
"""

from oncomind.models.extracted.therapeutic_data import TherapeuticData
from oncomind.models.extracted.literature_knowledge import (
    LiteratureKnowledge,
    LitDrugResistance,
    LitDrugSensitivity,
)
from oncomind.models.extracted.evidence_gaps import (
    EvidenceGaps,
    EvidenceGap,
    GapCategory,
    GapSeverity,
    CharacterizedAspect,
)

__all__ = [
    # Therapeutic data
    "TherapeuticData",
    # Literature knowledge
    "LiteratureKnowledge",
    "LitDrugResistance",
    "LitDrugSensitivity",
    # Evidence gaps
    "EvidenceGaps",
    "EvidenceGap",
    "GapCategory",
    "GapSeverity",
    "CharacterizedAspect",
]

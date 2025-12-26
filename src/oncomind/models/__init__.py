"""Data models for OncoMind."""

from oncomind.models.llm_insight import LLMInsight
from oncomind.models.therapeutic_evidence import TherapeuticEvidence, RecommendedTherapy
from oncomind.models.result import Result
from oncomind.models.evidence import Evidence
from oncomind.models.evidence.civic import CIViCEvidence
from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.cosmic import COSMICEvidence
from oncomind.models.variant import VariantInput

__all__ = [
    "Evidence",
    "Result",
    "VariantInput",
    "CIViCEvidence",
    "ClinVarEvidence",
    "COSMICEvidence",
    "TherapeuticEvidence",
    "RecommendedTherapy",  # Backwards compatibility alias
    "LLMInsight",
]

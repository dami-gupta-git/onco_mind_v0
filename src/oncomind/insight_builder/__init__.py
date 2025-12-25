"""Insight building module.

This module provides:
- EvidenceAggregator: Aggregates evidence from multiple data sources → Evidence
- Conductor: Orchestrates evidence aggregation + LLM insight generation → Result

Example (evidence only):
    >>> from oncomind.insight_builder import build_evidence
    >>> evidence = await build_evidence("BRAF V600E", tumor_type="Melanoma")
    >>> print(evidence.clinical.get_approved_drugs())
    ['Dabrafenib', 'Vemurafenib']

Example (with LLM insight):
    >>> from oncomind.insight_builder import conduct, ConductorConfig
    >>> config = ConductorConfig(enable_llm=True)
    >>> result = await conduct("BRAF V600E", tumor_type="Melanoma", config=config)
    >>> print(result.llm.llm_summary)
"""

from oncomind.insight_builder.aggregator import (
    EvidenceAggregator,
    EvidenceAggregatorConfig,
    build_evidence,
    build_evidences,
)
from oncomind.insight_builder.conductor import (
    Conductor,
    ConductorConfig,
    conduct,
    conduct_batch,
)


__all__ = [
    # Aggregator (evidence only)
    "EvidenceAggregator",
    "EvidenceAggregatorConfig",
    "build_evidence",
    "build_evidences",
    # Conductor (evidence + LLM)
    "Conductor",
    "ConductorConfig",
    "conduct",
    "conduct_batch",
]

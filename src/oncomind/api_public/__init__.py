"""Public API for OncoMind variant insight.

This module provides the high-level public API for variant insight:

    >>> from oncomind import get_insight, get_insights
    >>> panel = await get_insight("BRAF V600E", tumor_type="Melanoma")
    >>> panels = await get_insights(["BRAF V600E", "EGFR L858R"])

The API returns strongly-typed EvidencePanel objects with all aggregated
evidence organized into sections:
- identifiers: Variant IDs, HGVS notation
- kb: Knowledgebase evidence (CIViC, ClinVar, COSMIC, etc.)
- functional: Computational predictions (AlphaMissense, CADD, etc.)
- clinical: FDA approvals, clinical trials, gene context
- literature: PubMed articles, extracted knowledge
- meta: Processing metadata
"""

from oncomind.api_public.insight import (
    get_insight,
    get_insights,
    InsightConfig,
)

__all__ = [
    "get_insight",
    "get_insights",
    "InsightConfig",
]

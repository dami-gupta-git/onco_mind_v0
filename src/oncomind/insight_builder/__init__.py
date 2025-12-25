"""Insight building module.

This module provides the InsightBuilder class that orchestrates
API calls to multiple data sources and aggregates results into
an Insight object.

Example:
    >>> from oncomind.insight_builder import build_insight
    >>> insight = await build_insight("BRAF V600E", tumor_type="Melanoma")
    >>> print(insight.clinical.get_approved_drugs())
    ['Dabrafenib', 'Vemurafenib']

For batch processing:
    >>> from oncomind.insight_builder import build_insights
    >>> insights = await build_insights(["BRAF V600E", "EGFR L858R"])
"""

from oncomind.insight_builder.builder import (
    InsightBuilder,
    InsightBuilderConfig,
    build_insight,
    build_insights,
)


__all__ = [
    "InsightBuilder",
    "InsightBuilderConfig",
    "build_insight",
    "build_insights",
]

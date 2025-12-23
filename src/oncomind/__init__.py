"""OncoMind - AI-powered cancer variant insight and evidence synthesis.

Public API:
    >>> from oncomind import get_insight, get_insights
    >>> panel = await get_insight("BRAF V600E", tumor_type="Melanoma")
    >>> panels = await get_insights(["BRAF V600E", "EGFR L858R"])

For synchronous usage:
    >>> from oncomind import get_insight_sync
    >>> panel = get_insight_sync("BRAF V600E")
"""

__version__ = "0.1.0"

# Public API
from oncomind.api_public.annotate import (
    get_insight,
    get_insights,
    get_insight_sync,
    get_insights_sync,
    AnnotationConfig,
)

# Core models
from oncomind.models.evidence.evidence_panel import EvidencePanel
from oncomind.normalization import ParsedVariant, parse_variant_input

__all__ = [
    # Version
    "__version__",
    # Public API
    "get_insight",
    "get_insights",
    "get_insight_sync",
    "get_insights_sync",
    "AnnotationConfig",
    # Core models
    "EvidencePanel",
    "ParsedVariant",
    "parse_variant_input",
]

"""OncoMind - LLM-powered cancer variant annotation and evidence synthesis.

Public API:
    >>> from oncomind import process_variant, process_variants
    >>> panel = await process_variant("BRAF V600E", tumor_type="Melanoma")
    >>> panels = await process_variants(["BRAF V600E", "EGFR L858R"])

For synchronous usage:
    >>> from oncomind import process_variant_sync
    >>> panel = process_variant_sync("BRAF V600E")
"""

__version__ = "0.1.0"

# Public API
from oncomind.api_public.annotate import (
    process_variant,
    process_variants,
    process_variant_sync,
    process_variants_sync,
    AnnotationConfig,
)

# Core models
from oncomind.models.evidence.evidence_panel import EvidencePanel
from oncomind.normalization import ParsedVariant, parse_variant_input

__all__ = [
    # Version
    "__version__",
    # Public API
    "process_variant",
    "process_variants",
    "process_variant_sync",
    "process_variants_sync",
    "AnnotationConfig",
    # Core models
    "EvidencePanel",
    "ParsedVariant",
    "parse_variant_input",
]

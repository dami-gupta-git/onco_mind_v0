"""Public API for OncoMind variant annotation.

This module provides the high-level public API for variant annotation:

    >>> from oncomind import process_variant, process_variants
    >>> panel = await process_variant("BRAF V600E", tumor_type="Melanoma")
    >>> panels = await process_variants(["BRAF V600E", "EGFR L858R"])

The API returns strongly-typed EvidencePanel objects with all aggregated
evidence organized into sections:
- identifiers: Variant IDs, HGVS notation
- kb: Knowledgebase evidence (CIViC, ClinVar, COSMIC, etc.)
- functional: Computational predictions (AlphaMissense, CADD, etc.)
- clinical: FDA approvals, clinical trials, gene context
- literature: PubMed articles, extracted knowledge
- meta: Processing metadata
"""

from oncomind.api_public.annotate import (
    process_variant,
    process_variants,
    AnnotationConfig,
)

__all__ = [
    "process_variant",
    "process_variants",
    "AnnotationConfig",
]

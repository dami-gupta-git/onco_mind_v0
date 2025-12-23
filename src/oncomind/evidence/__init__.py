"""Evidence aggregation module.

This module provides the EvidenceBuilder class that orchestrates
API calls to multiple data sources and aggregates results into
an EvidencePanel.

Example:
    >>> from oncomind.evidence import build_evidence_panel
    >>> panel = await build_evidence_panel("BRAF V600E", tumor_type="Melanoma")
    >>> print(panel.clinical.get_approved_drugs())
    ['Dabrafenib', 'Vemurafenib']

For batch processing:
    >>> from oncomind.evidence import build_evidence_panels
    >>> panels = await build_evidence_panels(["BRAF V600E", "EGFR L858R"])
"""

from oncomind.evidence.builder import (
    EvidenceBuilder,
    EvidenceBuilderConfig,
    build_evidence_panel,
    build_evidence_panels,
)

__all__ = [
    "EvidenceBuilder",
    "EvidenceBuilderConfig",
    "build_evidence_panel",
    "build_evidence_panels",
]

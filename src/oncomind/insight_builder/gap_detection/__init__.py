"""Gap detection subpackage for evidence analysis.

This package analyzes evidence to identify gaps and areas needing further research.
"""

from .detector import detect_evidence_gaps
from .context import GapDetectionContext

__all__ = [
    "detect_evidence_gaps",
    "GapDetectionContext",
]

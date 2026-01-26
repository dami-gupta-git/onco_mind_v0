"""Evidence gap detection from aggregated evidence.

Analyzes what's missing from the evidence to guide research priorities.

This module is a thin wrapper that re-exports from the gap_detection subpackage.
The actual implementation has been refactored into:
  - gap_detection/helpers.py: Match counting and utility functions
  - gap_detection/checks.py: Individual gap detection check functions
  - gap_detection/discordance.py: Discordant evidence detection
  - gap_detection/context.py: GapDetectionContext dataclass
  - gap_detection/detector.py: Main detect_evidence_gaps entry point
"""

from .gap_detection import detect_evidence_gaps, GapDetectionContext

__all__ = ["detect_evidence_gaps", "GapDetectionContext"]

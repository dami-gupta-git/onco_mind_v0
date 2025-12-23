"""Experimental features module.

This module contains experimental features that are not part of the
stable public API. These features may change or be removed without notice.

Currently includes:
- tiering: Experimental AMP/ASCO/CAP tier computation

IMPORTANT: Do not rely on these features for production use.
They are provided for research and testing purposes only.
"""

from oncomind.experimental.tiering import (
    compute_experimental_tier,
    TierResult,
)

__all__ = [
    "compute_experimental_tier",
    "TierResult",
]

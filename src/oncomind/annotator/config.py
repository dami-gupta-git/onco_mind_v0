"""Annotation configuration."""

from dataclasses import dataclass


@dataclass
class AnnotatorConfig:
    """Configuration for the Annotator.

    Controls annotation sources and processing options.
    """

    # Future: annotation source toggles
    # enable_vep: bool = True
    # enable_clinvar: bool = True
    # enable_cosmic: bool = True

    # Logging
    log_level: str = "INFO"


__all__ = [
    "AnnotatorConfig",
]

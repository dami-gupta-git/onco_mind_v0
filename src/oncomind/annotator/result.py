"""Annotation result data model."""

from dataclasses import dataclass, field
from typing import Any


@dataclass
class AnnotationResult:
    """Result of annotating a VCF file."""

    variants: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "variants": self.variants,
        }


__all__ = [
    "AnnotationResult",
]

"""Embeddings module for numeric feature extraction.

This module provides tools to convert Insight objects into
numeric feature vectors suitable for:
- Similarity search
- Clustering
- ML model input
- Visualization

This is a placeholder for future implementation.

Example (future):
    >>> from oncomind.embeddings import extract_features, create_embedding
    >>> from oncomind import Conductor
    >>> async with Conductor() as conductor:
    ...     result = await conductor.run("BRAF V600E")
    >>> features = extract_features(result.evidence)  # dict of numeric features
    >>> embedding = create_embedding(result.evidence)  # dense vector
"""

from oncomind.embeddings.features import (
    extract_features,
    FeatureConfig,
)

__all__ = [
    "extract_features",
    "FeatureConfig",
]

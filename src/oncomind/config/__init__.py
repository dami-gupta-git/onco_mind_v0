"""Configuration module for OncoMind.

Constants are available via: from oncomind.config.constants import ...
Variant classes are available via: from oncomind.config import VariantClassConfig
"""

from oncomind.config.variant_classes import VariantClassConfig, load_variant_classes

__all__ = ["VariantClassConfig", "load_variant_classes"]

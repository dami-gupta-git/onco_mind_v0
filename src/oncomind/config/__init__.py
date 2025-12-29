"""Configuration module for OncoMind.

Constants are available via: from oncomind.config.constants import ...
Variant classes are available via: from oncomind.config import VariantClassConfig
Debug/logging: from oncomind.config.debug import get_logger, set_log_level
"""

from oncomind.config.variant_classes import VariantClassConfig, load_variant_classes
from oncomind.config.debug import (
    get_logger,
    set_log_level,
    get_log_level,
    is_debug,
    debug,
    info,
    warn,
    error,
)

__all__ = [
    "VariantClassConfig",
    "load_variant_classes",
    "get_logger",
    "set_log_level",
    "get_log_level",
    "is_debug",
    "debug",
    "info",
    "warn",
    "error",
]

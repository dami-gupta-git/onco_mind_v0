"""Debug and logging configuration for OncoMind.

Provides application-wide logging with configurable log levels.
Log level can be set via:
1. Environment variable: ONCOMIND_LOG_LEVEL=DEBUG|INFO|WARN|ERROR
2. Programmatically: set_log_level("DEBUG")
3. CLI flag: --debug (sets DEBUG level)

Default level is INFO.
"""

import logging
import os
import sys
from enum import Enum
from typing import Literal

# Log level type
LogLevel = Literal["DEBUG", "INFO", "WARN", "WARNING", "ERROR"]


class LogLevels(str, Enum):
    """Available log levels."""
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARN = "WARN"
    WARNING = "WARNING"  # Alias for WARN
    ERROR = "ERROR"


# Map string levels to logging module levels
_LEVEL_MAP = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARN": logging.WARNING,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}

# Default log level
DEFAULT_LOG_LEVEL = "INFO"

# Global logger instance
_logger: logging.Logger | None = None
_current_level: str = DEFAULT_LOG_LEVEL


def _get_level_from_env() -> str:
    """Get log level from environment variable."""
    env_level = os.environ.get("ONCOMIND_LOG_LEVEL", "").upper()
    if env_level in _LEVEL_MAP:
        return env_level
    return DEFAULT_LOG_LEVEL


def _create_logger(level: str = DEFAULT_LOG_LEVEL) -> logging.Logger:
    """Create and configure the OncoMind logger."""
    logger = logging.getLogger("oncomind")

    # Clear existing handlers to avoid duplicates
    logger.handlers.clear()

    # Set level
    logger.setLevel(_LEVEL_MAP.get(level.upper(), logging.INFO))

    # Console handler with colored output
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(_LEVEL_MAP.get(level.upper(), logging.INFO))

    # Format: [LEVEL] module: message
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S"
    )
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Prevent propagation to root logger
    logger.propagate = False

    return logger


def get_logger(name: str | None = None) -> logging.Logger:
    """Get a logger instance for the given module.

    Args:
        name: Optional module name (e.g., "oncomind.api.civic").
              If None, returns the root oncomind logger.

    Returns:
        Logger instance configured with the current log level.

    Example:
        from oncomind.config.debug import get_logger
        logger = get_logger(__name__)
        logger.debug("Fetching data from CIViC...")
        logger.info("Found 5 assertions")
        logger.warning("Rate limit approaching")
        logger.error("API request failed")
    """
    global _logger, _current_level

    # Initialize root logger if needed
    if _logger is None:
        _current_level = _get_level_from_env()
        _logger = _create_logger(_current_level)

    # Return root logger or child logger
    if name is None or name == "oncomind":
        return _logger

    # Create child logger that inherits from root
    if name.startswith("oncomind."):
        return logging.getLogger(name)
    else:
        return logging.getLogger(f"oncomind.{name}")


def set_log_level(level: LogLevel) -> None:
    """Set the log level for all OncoMind loggers.

    Args:
        level: One of "DEBUG", "INFO", "WARN", "ERROR"

    Example:
        from oncomind.config.debug import set_log_level
        set_log_level("DEBUG")  # Enable verbose output
    """
    global _logger, _current_level

    level_upper = level.upper()
    if level_upper not in _LEVEL_MAP:
        raise ValueError(f"Invalid log level: {level}. Must be one of: DEBUG, INFO, WARN, ERROR")

    _current_level = level_upper

    # Update root logger
    if _logger is not None:
        _logger.setLevel(_LEVEL_MAP[level_upper])
        for handler in _logger.handlers:
            handler.setLevel(_LEVEL_MAP[level_upper])
    else:
        # Initialize with the new level
        _logger = _create_logger(level_upper)

    # Also set environment variable for child processes
    os.environ["ONCOMIND_LOG_LEVEL"] = level_upper


def get_log_level() -> str:
    """Get the current log level.

    Returns:
        Current log level as a string (DEBUG, INFO, WARN, ERROR)
    """
    global _current_level
    return _current_level


def is_debug() -> bool:
    """Check if debug logging is enabled.

    Returns:
        True if current log level is DEBUG
    """
    return get_log_level() == "DEBUG"


def reset_logger() -> None:
    """Reset the logger (mainly for testing)."""
    global _logger, _current_level
    if _logger is not None:
        _logger.handlers.clear()
    _logger = None
    _current_level = DEFAULT_LOG_LEVEL


# Convenience function for quick debug prints during development
def debug(msg: str, *args, **kwargs) -> None:
    """Quick debug logging without getting a logger instance.

    Example:
        from oncomind.config.debug import debug
        debug("Processing variant: %s %s", gene, variant)
    """
    get_logger().debug(msg, *args, **kwargs)


def info(msg: str, *args, **kwargs) -> None:
    """Quick info logging without getting a logger instance."""
    get_logger().info(msg, *args, **kwargs)


def warn(msg: str, *args, **kwargs) -> None:
    """Quick warning logging without getting a logger instance."""
    get_logger().warning(msg, *args, **kwargs)


def error(msg: str, *args, **kwargs) -> None:
    """Quick error logging without getting a logger instance."""
    get_logger().error(msg, *args, **kwargs)

"""Tests for debug/logging configuration."""

import logging
import os
import pytest
from unittest.mock import patch

from oncomind.config.debug import (
    get_logger,
    set_log_level,
    get_log_level,
    is_debug,
    reset_logger,
    debug,
    info,
    warn,
    error,
    DEFAULT_LOG_LEVEL,
    _LEVEL_MAP,
)


class TestLogLevel:
    """Tests for log level management."""

    def setup_method(self):
        """Reset logger state before each test."""
        reset_logger()
        # Clear environment variable if set
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def teardown_method(self):
        """Clean up after each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def test_default_log_level(self):
        """Test that default log level is INFO."""
        assert DEFAULT_LOG_LEVEL == "INFO"
        logger = get_logger()
        assert get_log_level() == "INFO"

    def test_set_log_level_debug(self):
        """Test setting log level to DEBUG."""
        set_log_level("DEBUG")
        assert get_log_level() == "DEBUG"
        assert is_debug() is True

    def test_set_log_level_info(self):
        """Test setting log level to INFO."""
        set_log_level("DEBUG")  # First set to DEBUG
        set_log_level("INFO")   # Then back to INFO
        assert get_log_level() == "INFO"
        assert is_debug() is False

    def test_set_log_level_warn(self):
        """Test setting log level to WARN."""
        set_log_level("WARN")
        assert get_log_level() == "WARN"
        assert is_debug() is False

    def test_set_log_level_warning_alias(self):
        """Test that WARNING is accepted as alias for WARN."""
        set_log_level("WARNING")
        assert get_log_level() == "WARNING"

    def test_set_log_level_error(self):
        """Test setting log level to ERROR."""
        set_log_level("ERROR")
        assert get_log_level() == "ERROR"
        assert is_debug() is False

    def test_set_log_level_case_insensitive(self):
        """Test that log level setting is case insensitive."""
        set_log_level("debug")
        assert get_log_level() == "DEBUG"

        set_log_level("Info")
        assert get_log_level() == "INFO"

        set_log_level("wArN")
        assert get_log_level() == "WARN"

    def test_set_invalid_log_level_raises(self):
        """Test that invalid log level raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            set_log_level("INVALID")
        assert "Invalid log level" in str(exc_info.value)

    def test_environment_variable_sets_level(self):
        """Test that ONCOMIND_LOG_LEVEL env var sets the level."""
        os.environ["ONCOMIND_LOG_LEVEL"] = "DEBUG"
        reset_logger()  # Reset to pick up env var
        logger = get_logger()  # Initialize logger
        assert get_log_level() == "DEBUG"

    def test_environment_variable_ignored_if_invalid(self):
        """Test that invalid env var value falls back to default."""
        os.environ["ONCOMIND_LOG_LEVEL"] = "INVALID"
        reset_logger()
        logger = get_logger()
        assert get_log_level() == "INFO"

    def test_set_log_level_updates_environment(self):
        """Test that set_log_level updates the environment variable."""
        set_log_level("DEBUG")
        assert os.environ.get("ONCOMIND_LOG_LEVEL") == "DEBUG"


class TestGetLogger:
    """Tests for get_logger function."""

    def setup_method(self):
        """Reset logger state before each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def teardown_method(self):
        """Clean up after each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def test_get_root_logger(self):
        """Test getting the root oncomind logger."""
        logger = get_logger()
        assert logger.name == "oncomind"

    def test_get_root_logger_with_explicit_name(self):
        """Test getting root logger with explicit 'oncomind' name."""
        logger = get_logger("oncomind")
        assert logger.name == "oncomind"

    def test_get_root_logger_with_none(self):
        """Test getting root logger with None name."""
        logger = get_logger(None)
        assert logger.name == "oncomind"

    def test_get_child_logger(self):
        """Test getting a child logger."""
        logger = get_logger("oncomind.api.civic")
        assert logger.name == "oncomind.api.civic"

    def test_get_child_logger_auto_prefix(self):
        """Test that non-prefixed names get 'oncomind.' prefix."""
        logger = get_logger("mymodule")
        assert logger.name == "oncomind.mymodule"

    def test_logger_has_handler(self):
        """Test that the root logger has at least one handler."""
        logger = get_logger()
        assert len(logger.handlers) > 0

    def test_logger_respects_level(self):
        """Test that logger respects the set level."""
        set_log_level("DEBUG")
        logger = get_logger()
        assert logger.level == logging.DEBUG

        set_log_level("ERROR")
        assert logger.level == logging.ERROR

    def test_multiple_get_logger_returns_same_root(self):
        """Test that multiple calls return the same root logger instance."""
        logger1 = get_logger()
        logger2 = get_logger()
        assert logger1 is logger2

    def test_logger_does_not_propagate(self):
        """Test that the root logger does not propagate to parent."""
        logger = get_logger()
        assert logger.propagate is False


class TestConvenienceFunctions:
    """Tests for convenience logging functions."""

    def setup_method(self):
        """Reset logger state before each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def teardown_method(self):
        """Clean up after each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def test_debug_function(self):
        """Test the debug() convenience function."""
        set_log_level("DEBUG")
        with patch.object(get_logger(), 'debug') as mock_debug:
            debug("Test message")
            mock_debug.assert_called_once_with("Test message")

    def test_info_function(self):
        """Test the info() convenience function."""
        with patch.object(get_logger(), 'info') as mock_info:
            info("Test info message")
            mock_info.assert_called_once_with("Test info message")

    def test_warn_function(self):
        """Test the warn() convenience function."""
        with patch.object(get_logger(), 'warning') as mock_warn:
            warn("Test warning message")
            mock_warn.assert_called_once_with("Test warning message")

    def test_error_function(self):
        """Test the error() convenience function."""
        with patch.object(get_logger(), 'error') as mock_error:
            error("Test error message")
            mock_error.assert_called_once_with("Test error message")

    def test_debug_with_args(self):
        """Test debug function with formatting args."""
        set_log_level("DEBUG")
        with patch.object(get_logger(), 'debug') as mock_debug:
            debug("Processing %s %s", "BRAF", "V600E")
            mock_debug.assert_called_once_with("Processing %s %s", "BRAF", "V600E")


class TestIsDebug:
    """Tests for is_debug function."""

    def setup_method(self):
        """Reset logger state before each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def teardown_method(self):
        """Clean up after each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def test_is_debug_true_when_debug(self):
        """Test is_debug returns True when level is DEBUG."""
        set_log_level("DEBUG")
        assert is_debug() is True

    def test_is_debug_false_when_info(self):
        """Test is_debug returns False when level is INFO."""
        set_log_level("INFO")
        assert is_debug() is False

    def test_is_debug_false_when_warn(self):
        """Test is_debug returns False when level is WARN."""
        set_log_level("WARN")
        assert is_debug() is False

    def test_is_debug_false_when_error(self):
        """Test is_debug returns False when level is ERROR."""
        set_log_level("ERROR")
        assert is_debug() is False


class TestResetLogger:
    """Tests for reset_logger function."""

    def setup_method(self):
        """Reset logger state before each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def teardown_method(self):
        """Clean up after each test."""
        reset_logger()
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]

    def test_reset_clears_handlers(self):
        """Test that reset_logger clears handlers."""
        logger = get_logger()
        assert len(logger.handlers) > 0
        reset_logger()
        # After reset, getting logger again should reinitialize
        logger2 = get_logger()
        assert len(logger2.handlers) > 0

    def test_reset_restores_default_level(self):
        """Test that reset_logger restores default log level."""
        set_log_level("DEBUG")
        assert get_log_level() == "DEBUG"
        reset_logger()
        # Clear env var since set_log_level updates it
        if "ONCOMIND_LOG_LEVEL" in os.environ:
            del os.environ["ONCOMIND_LOG_LEVEL"]
        # After reset, level should be back to default
        logger = get_logger()
        assert get_log_level() == DEFAULT_LOG_LEVEL


class TestLevelMap:
    """Tests for level mapping."""

    def test_all_levels_mapped(self):
        """Test that all expected levels are in the mapping."""
        expected_levels = ["DEBUG", "INFO", "WARN", "WARNING", "ERROR"]
        for level in expected_levels:
            assert level in _LEVEL_MAP

    def test_debug_maps_to_logging_debug(self):
        """Test DEBUG maps to logging.DEBUG."""
        assert _LEVEL_MAP["DEBUG"] == logging.DEBUG

    def test_info_maps_to_logging_info(self):
        """Test INFO maps to logging.INFO."""
        assert _LEVEL_MAP["INFO"] == logging.INFO

    def test_warn_maps_to_logging_warning(self):
        """Test WARN maps to logging.WARNING."""
        assert _LEVEL_MAP["WARN"] == logging.WARNING

    def test_warning_maps_to_logging_warning(self):
        """Test WARNING maps to logging.WARNING."""
        assert _LEVEL_MAP["WARNING"] == logging.WARNING

    def test_error_maps_to_logging_error(self):
        """Test ERROR maps to logging.ERROR."""
        assert _LEVEL_MAP["ERROR"] == logging.ERROR

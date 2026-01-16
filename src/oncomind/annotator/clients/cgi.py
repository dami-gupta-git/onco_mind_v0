"""Annotator CGI (Cancer Genome Interpreter) client."""

from typing import Any

from oncomind.config.debug import get_logger

logger = get_logger(__name__)


class AnnotatorCGIClient:
    """CGI client for the annotator workflow."""

    def __init__(self, timeout: float = 30.0) -> None:
        """Initialize the client.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout

    async def __aenter__(self) -> "AnnotatorCGIClient":
        """Async context manager entry."""
        return self

    async def __aexit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Async context manager exit."""
        pass

    async def close(self) -> None:
        """Close the client."""
        pass

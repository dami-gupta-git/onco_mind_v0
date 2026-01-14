"""Integration tests for Streamlit backend functions.

Tests validate the backend functions used by the Streamlit UI.
"""

import pytest
import sys
from pathlib import Path

# Add streamlit directory to path for imports
streamlit_dir = Path(__file__).parent.parent.parent.parent / "streamlit"
sys.path.insert(0, str(streamlit_dir))

from backend import get_variant_insight, batch_get_variant_insights


class TestGetVariantInsight:
    """Tests for get_variant_insight function."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_basic_annotation(self):
        """Basic annotation should work."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["gene"] == "BRAF"
        assert result["variant"]["variant"] == "V600E"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_annotation_with_tumor_type(self):
        """Annotation with tumor type should work."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["tumor_type"] == "Melanoma"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fast_mode_annotation(self):
        """Fast mode (no literature) should work."""
        result = await get_variant_insight(
            gene="EGFR",
            variant="L858R",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["gene"] == "EGFR"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_annotation_returns_expected_keys(self):
        """Annotation should return expected structure."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "variant" in result
        assert "insight" in result
        assert "identifiers" in result
        assert "hgvs" in result
        assert "clinvar" in result
        assert "annotations" in result
        assert "recommended_therapies" in result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_annotation_insight_section(self):
        """Annotation insight section should have summary."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )

        assert "insight" in result
        assert "summary" in result["insight"]
        assert isinstance(result["insight"]["summary"], str)


class TestBatchGetVariantInsights:
    """Tests for batch_get_variant_insights function."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_basic(self):
        """Basic batch processing should work."""
        variants = [
            {"gene": "BRAF", "variant": "V600E"},
            {"gene": "EGFR", "variant": "L858R"},
        ]

        results = await batch_get_variant_insights(
            variants,
            enable_llm=False,
            enable_literature=False,
        )

        assert len(results) == 2
        assert results[0]["variant"]["gene"] == "BRAF"
        assert results[1]["variant"]["gene"] == "EGFR"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_with_tumor_types(self):
        """Batch with tumor types should work."""
        variants = [
            {"gene": "BRAF", "variant": "V600E", "tumor_type": "Melanoma"},
            {"gene": "EGFR", "variant": "L858R", "tumor_type": "NSCLC"},
        ]

        results = await batch_get_variant_insights(
            variants,
            enable_llm=False,
            enable_literature=False,
        )

        assert len(results) == 2
        assert results[0]["variant"]["tumor_type"] == "Melanoma"
        assert results[1]["variant"]["tumor_type"] == "NSCLC"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_fast_mode(self):
        """Batch fast mode should work."""
        variants = [{"gene": "KRAS", "variant": "G12C"}]

        results = await batch_get_variant_insights(
            variants,
            enable_llm=False,
            enable_literature=False,
        )

        assert len(results) == 1
        assert "error" not in results[0]

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_with_progress_callback(self):
        """Batch with progress callback should work."""
        variants = [
            {"gene": "BRAF", "variant": "V600E"},
            {"gene": "TP53", "variant": "R248W"},
        ]

        progress_calls = []

        def progress_callback(current, total):
            progress_calls.append((current, total))

        results = await batch_get_variant_insights(
            variants,
            enable_llm=False,
            enable_literature=False,
            progress_callback=progress_callback,
        )

        assert len(results) == 2
        assert len(progress_calls) == 2
        assert progress_calls[-1] == (2, 2)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_handles_errors_gracefully(self):
        """Batch should handle individual errors gracefully."""
        variants = [
            {"gene": "BRAF", "variant": "V600E"},
            {"gene": "INVALID_GENE_SYMBOL_12345", "variant": "X999Y"},
        ]

        results = await batch_get_variant_insights(
            variants,
            enable_llm=False,
            enable_literature=False,
        )

        # Should return results for both, even if one has an error
        assert len(results) == 2
        # First should succeed
        assert results[0]["variant"]["gene"] == "BRAF"


class TestBackendOptions:
    """Tests for backend option combinations."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_llm_disabled_literature_enabled(self):
        """LLM disabled with literature enabled should work."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=True,
        )

        assert "error" not in result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_both_disabled(self):
        """Both LLM and literature disabled should be fastest."""
        import time

        start = time.time()
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            enable_llm=False,
            enable_literature=False,
        )
        elapsed = time.time() - start

        assert "error" not in result
        # Should complete in under 30 seconds
        assert elapsed < 30, f"Fastest mode took {elapsed:.1f}s"

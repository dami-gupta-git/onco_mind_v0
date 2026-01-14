"""Integration tests for the public API (get_insight, get_insights).

Tests validate fast mode, LLM options, and batch processing work correctly.
"""

import pytest
import json
import tempfile
from pathlib import Path

from oncomind import get_insight, get_insights, InsightConfig, Result


class TestGetInsightFastMode:
    """Tests for get_insight with fast mode (no literature)."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fast_mode_skips_literature(self):
        """Fast mode should skip literature search."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        panel = await get_insight("BRAF V600E", config=config)

        assert isinstance(panel, Result)
        assert panel.identifiers.gene == "BRAF"
        assert panel.identifiers.variant == "V600E"
        # Literature should be empty in fast mode
        assert len(panel.evidence.pubmed_articles) == 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fast_mode_still_fetches_databases(self):
        """Fast mode should still fetch database evidence."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        panel = await get_insight("BRAF V600E", tumor_type="Melanoma", config=config)

        # Should still have database evidence
        assert panel.identifiers.gene == "BRAF"
        # BRAF V600E is well-characterized, should have evidence (FDA approvals, KB data)
        assert panel.has_evidence()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fast_mode_performance(self):
        """Fast mode should complete quickly."""
        import time

        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        start = time.time()
        panel = await get_insight("EGFR L858R", config=config)
        elapsed = time.time() - start

        assert isinstance(panel, Result)
        # Fast mode should complete in under 30 seconds
        assert elapsed < 30, f"Fast mode took {elapsed:.1f}s, expected < 30s"


class TestGetInsightFullMode:
    """Tests for get_insight with full mode (with literature)."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_full_mode_includes_literature(self):
        """Full mode should include literature search."""
        config = InsightConfig(
            enable_literature=True,
            enable_llm=False,
        )

        panel = await get_insight("BRAF V600E", tumor_type="Melanoma", config=config)

        assert isinstance(panel, Result)
        assert panel.identifiers.gene == "BRAF"
        # Full mode should attempt literature search
        # (may or may not find articles depending on API availability)
        # Literature is accessed via panel.evidence.pubmed_articles
        assert hasattr(panel.evidence, 'pubmed_articles')


class TestGetInsightLLMOptions:
    """Tests for get_insight with LLM options."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_llm_disabled_by_default(self):
        """LLM should be disabled by default."""
        config = InsightConfig()

        assert config.enable_llm is False

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_config_llm_model_setting(self):
        """LLM model should be configurable."""
        config = InsightConfig(
            enable_llm=True,
            llm_model="gpt-4o-mini",
            llm_temperature=0.2,
        )

        assert config.enable_llm is True
        assert config.llm_model == "gpt-4o-mini"
        assert config.llm_temperature == 0.2


class TestGetInsights:
    """Tests for batch get_insights."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_with_list_of_strings(self):
        """Batch should work with list of variant strings."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        panels = await get_insights(
            ["BRAF V600E", "EGFR L858R"],
            config=config,
        )

        assert len(panels) == 2
        assert all(isinstance(p, Result) for p in panels)
        assert panels[0].identifiers.gene == "BRAF"
        assert panels[1].identifiers.gene == "EGFR"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_with_progress_callback(self):
        """Batch should call progress callback."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        progress_calls = []

        def progress_callback(current: int, total: int):
            progress_calls.append((current, total))

        panels = await get_insights(
            ["BRAF V600E", "KRAS G12C"],
            config=config,
            progress_callback=progress_callback,
        )

        assert len(panels) == 2
        assert len(progress_calls) >= 2
        # Should have called with (1, 2), (2, 2)
        assert progress_calls[-1] == (2, 2)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_fast_mode(self):
        """Batch fast mode should skip literature for all variants."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        panels = await get_insights(
            ["BRAF V600E", "TP53 R248W"],
            config=config,
        )

        assert len(panels) == 2
        # All panels should have empty literature in fast mode
        for panel in panels:
            assert len(panel.evidence.pubmed_articles) == 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_batch_with_tumor_type(self):
        """Batch should apply tumor type to all variants."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        panels = await get_insights(
            ["BRAF V600E", "EGFR L858R"],
            tumor_type="NSCLC",
            config=config,
        )

        assert len(panels) == 2
        # Tumor type should be applied
        for panel in panels:
            assert panel.context.tumor_type == "NSCLC"


class TestInsightConfig:
    """Tests for InsightConfig."""

    def test_default_config(self):
        """Default config should have sensible defaults."""
        config = InsightConfig()

        assert config.enable_vicc is True
        assert config.enable_civic_assertions is True
        assert config.enable_clinical_trials is True
        assert config.enable_literature is True
        assert config.enable_llm is False
        assert config.llm_model == "gpt-4o-mini"
        assert config.llm_temperature == 0.1

    def test_fast_mode_config(self):
        """Fast mode config should disable literature."""
        config = InsightConfig(
            enable_literature=False,
        )

        assert config.enable_literature is False
        assert config.enable_vicc is True  # Other sources still enabled

    def test_llm_config(self):
        """LLM config should be customizable."""
        config = InsightConfig(
            enable_llm=True,
            llm_model="gpt-4o",
            llm_temperature=0.5,
        )

        assert config.enable_llm is True
        assert config.llm_model == "gpt-4o"
        assert config.llm_temperature == 0.5


class TestInsightOutput:
    """Tests for Insight output structure."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_panel_has_required_sections(self):
        """Insight should have all required sections."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        panel = await get_insight("BRAF V600E", config=config)

        # Check all sections exist via Result's evidence field
        assert hasattr(panel, 'evidence')
        assert hasattr(panel.evidence, 'identifiers')
        assert hasattr(panel.evidence, 'functional')
        assert hasattr(panel.evidence, 'pubmed_articles')
        # Result also has property shortcuts
        assert hasattr(panel, 'identifiers')
        assert hasattr(panel, 'functional')

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_panel_serialization(self):
        """Insight should serialize to JSON."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        panel = await get_insight("BRAF V600E", config=config)

        # Should serialize without error
        json_data = panel.model_dump(mode="json")
        assert isinstance(json_data, dict)

        # Should be valid JSON
        json_str = json.dumps(json_data)
        assert len(json_str) > 0

        # Should round-trip
        parsed = json.loads(json_str)
        # In new Result model, identifiers are nested under evidence
        assert parsed['evidence']['identifiers']['gene'] == "BRAF"
        assert parsed['evidence']['identifiers']['variant'] == "V600E"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_panel_summary(self):
        """Result should provide a summary."""
        config = InsightConfig(
            enable_literature=False,
            enable_llm=False,
        )

        panel = await get_insight("BRAF V600E", tumor_type="Melanoma", config=config)

        # Result now has get_summary() method and has_evidence()
        summary = panel.get_summary()
        assert isinstance(summary, str)
        assert len(summary) > 0

        # Should have evidence for well-characterized variant
        assert panel.has_evidence()

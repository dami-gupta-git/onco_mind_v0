"""Tests for cross-source analysis functionality.

Tests for:
1. create_cross_source_prompt - prompt structure and content
2. get_cross_source_analysis - LLM service method
"""

import json
import pytest
from unittest.mock import AsyncMock, patch

from oncomind.llm.service import LLMService
from oncomind.llm.prompts import create_cross_source_prompt


class TestCreateCrossSourcePrompt:
    """Tests for create_cross_source_prompt function."""

    def test_returns_list_of_messages(self):
        """Should return a list of message dicts."""
        messages = create_cross_source_prompt(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            gene_context="EGFR is a receptor tyrosine kinase",
            cross_source_synthesis="CGI: erlotinib responsive\nCIViC: gefitinib sensitive",
        )

        assert isinstance(messages, list)
        assert len(messages) == 2
        assert all("role" in m and "content" in m for m in messages)

    def test_message_roles(self):
        """Should have system and user roles in correct order."""
        messages = create_cross_source_prompt(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            gene_context="BRAF is a kinase in MAPK pathway",
            cross_source_synthesis="FDA approved: vemurafenib",
        )

        assert messages[0]["role"] == "system"
        assert messages[1]["role"] == "user"

    def test_user_prompt_contains_gene_variant(self):
        """User prompt should contain the gene and variant."""
        messages = create_cross_source_prompt(
            gene="KRAS",
            variant="G12C",
            tumor_type="NSCLC",
            gene_context="KRAS is a GTPase",
            cross_source_synthesis="Sotorasib is FDA approved",
        )

        user_content = messages[1]["content"]
        assert "KRAS" in user_content
        assert "G12C" in user_content
        assert "NSCLC" in user_content

    def test_user_prompt_contains_cross_source_synthesis(self):
        """User prompt should contain the cross-source synthesis data."""
        synthesis = "CGI: Vemurafenib FDA-approved\nCIViC: Dabrafenib sensitive\nLiterature: Encorafenib effective"

        messages = create_cross_source_prompt(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            gene_context="BRAF kinase",
            cross_source_synthesis=synthesis,
        )

        user_content = messages[1]["content"]
        assert "Vemurafenib" in user_content or "FDA-approved" in user_content
        assert synthesis in user_content or all(
            term in user_content for term in ["CGI", "CIViC", "Literature"]
        )

    def test_handles_none_tumor_type(self):
        """Should handle None tumor type gracefully."""
        messages = create_cross_source_prompt(
            gene="PIK3CA",
            variant="H1047R",
            tumor_type=None,
            gene_context="PIK3CA is catalytic subunit of PI3K",
            cross_source_synthesis="Alpelisib is FDA approved",
        )

        user_content = messages[1]["content"]
        # Should use "Pan-cancer" or similar fallback
        assert "PIK3CA" in user_content
        assert "H1047R" in user_content
        # The function uses "Pan-cancer" as fallback
        assert "Pan-cancer" in user_content or "pan-cancer" in user_content.lower()

    def test_includes_sensitivity_resistance_summaries(self):
        """Should include sensitivity and resistance summaries when provided."""
        messages = create_cross_source_prompt(
            gene="KIT",
            variant="D816V",
            tumor_type="Mastocytosis",
            gene_context="KIT receptor tyrosine kinase",
            cross_source_synthesis="Avapritinib FDA approved",
            sensitivity_summary="Sensitive to avapritinib",
            resistance_summary="Resistant to imatinib",
        )

        user_content = messages[1]["content"]
        assert "avapritinib" in user_content.lower()
        assert "imatinib" in user_content.lower()

    def test_system_prompt_has_expected_structure(self):
        """System prompt should contain expected analysis guidance."""
        messages = create_cross_source_prompt(
            gene="ALK",
            variant="EML4-ALK",
            tumor_type="NSCLC",
            gene_context="ALK fusion oncogene",
            cross_source_synthesis="Crizotinib, alectinib FDA approved",
        )

        system_content = messages[0]["content"].lower()
        # Should contain guidance about cross-source analysis
        assert any(term in system_content for term in [
            "cross-source", "evidence", "synthesis", "analyze", "drug"
        ])


class TestGetCrossSourceAnalysis:
    """Tests for LLMService.get_cross_source_analysis method."""

    @pytest.mark.asyncio
    async def test_returns_none_when_no_synthesis_data(self):
        """Should return None when no cross-source synthesis is provided."""
        service = LLMService(model="gpt-4o-mini")

        result = await service.get_cross_source_analysis(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            gene_context="EGFR kinase",
            cross_source_synthesis="",  # Empty
        )

        assert result is None

    @pytest.mark.asyncio
    async def test_returns_none_for_placeholder_text(self):
        """Should return None for placeholder text."""
        service = LLMService(model="gpt-4o-mini")

        result = await service.get_cross_source_analysis(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            gene_context="EGFR kinase",
            cross_source_synthesis="No cross-source synthesis available.",
        )

        assert result is None

    @pytest.mark.asyncio
    async def test_calls_llm_with_correct_messages(self):
        """Should call LLM with properly formatted messages."""
        captured_messages = None

        async def mock_acompletion(**kwargs):
            nonlocal captured_messages
            captured_messages = kwargs["messages"]

            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps({
                        "strongest_evidence": [],
                        "conflicting_signals": [],
                        "emerging_targets": [],
                        "key_gaps": [],
                        "summary": "Test summary"
                    })),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            await service.get_cross_source_analysis(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                gene_context="BRAF kinase",
                cross_source_synthesis="CGI: Vemurafenib responsive",
            )

        assert captured_messages is not None
        assert len(captured_messages) == 2
        assert captured_messages[0]["role"] == "system"
        assert captured_messages[1]["role"] == "user"
        assert "BRAF" in captured_messages[1]["content"]
        assert "V600E" in captured_messages[1]["content"]

    @pytest.mark.asyncio
    async def test_returns_structured_result(self):
        """Should return properly structured result dict."""
        llm_response = {
            "strongest_evidence": [
                {"drug": "vemurafenib", "evidence_level": "FDA", "sources": ["CGI", "CIViC"]}
            ],
            "conflicting_signals": [],
            "emerging_targets": [
                {"drug": "encorafenib", "rationale": "Newer BRAF inhibitor"}
            ],
            "key_gaps": ["Resistance mechanisms"],
            "summary": "V600E is well-characterized actionable target"
        }

        async def mock_acompletion(**kwargs):
            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps(llm_response)),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            result = await service.get_cross_source_analysis(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                gene_context="BRAF kinase in MAPK pathway",
                cross_source_synthesis="CGI: Vemurafenib FDA approved",
            )

        assert result is not None
        assert "strongest_evidence" in result
        assert "conflicting_signals" in result
        assert "emerging_targets" in result
        assert len(result["strongest_evidence"]) == 1
        assert result["strongest_evidence"][0]["drug"] == "vemurafenib"

    @pytest.mark.asyncio
    async def test_returns_none_on_llm_error(self):
        """Should return None when LLM call fails."""
        async def mock_acompletion(**kwargs):
            raise Exception("LLM API error")

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            result = await service.get_cross_source_analysis(
                gene="EGFR",
                variant="L858R",
                tumor_type="NSCLC",
                gene_context="EGFR kinase",
                cross_source_synthesis="Erlotinib sensitive",
            )

        assert result is None

    @pytest.mark.asyncio
    async def test_handles_json_in_markdown(self):
        """Should handle JSON wrapped in markdown code blocks."""
        llm_response = {
            "strongest_evidence": [],
            "conflicting_signals": [],
            "emerging_targets": [],
            "key_gaps": [],
            "summary": "Test"
        }
        markdown_wrapped = f"```json\n{json.dumps(llm_response)}\n```"

        async def mock_acompletion(**kwargs):
            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=markdown_wrapped),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            result = await service.get_cross_source_analysis(
                gene="KRAS",
                variant="G12C",
                tumor_type="NSCLC",
                gene_context="KRAS GTPase",
                cross_source_synthesis="Sotorasib FDA approved",
            )

        assert result is not None
        assert "summary" in result

    @pytest.mark.asyncio
    async def test_includes_sensitivity_resistance_in_prompt(self):
        """Should include sensitivity and resistance summaries in the prompt."""
        captured_messages = None

        async def mock_acompletion(**kwargs):
            nonlocal captured_messages
            captured_messages = kwargs["messages"]

            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps({
                        "strongest_evidence": [],
                        "conflicting_signals": [],
                        "emerging_targets": [],
                        "key_gaps": [],
                        "summary": "Test"
                    })),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            await service.get_cross_source_analysis(
                gene="KIT",
                variant="D816V",
                tumor_type="GIST",
                gene_context="KIT receptor",
                cross_source_synthesis="Avapritinib responsive",
                sensitivity_summary="Avapritinib shows efficacy",
                resistance_summary="Imatinib resistance observed",
            )

        assert captured_messages is not None
        user_content = captured_messages[1]["content"]
        assert "Avapritinib" in user_content
        assert "Imatinib" in user_content or "imatinib" in user_content.lower()

"""Tests for LLM service."""

import json
import pytest
from unittest.mock import AsyncMock, patch

from oncomind.llm.service import LLMService


class TestLLMService:
    """Tests for LLMService."""

    @pytest.mark.asyncio
    async def test_get_llm_insight(self):
        """Test variant insight generation."""
        service = LLMService()

        # Mock the acompletion call
        with patch("oncomind.llm.service.acompletion", new_callable=AsyncMock) as mock_call:
            # Create mock response object with new format (raw component fields)
            response_json = {
                "functional_summary": "BRAF V600E is a well-characterized oncogenic mutation.",
                "biological_context": "Constitutive MAPK pathway activation.",
                "research_implications": "Multiple FDA-approved BRAF inhibitors exist for melanoma.",
                "therapeutic_landscape": {
                    "fda_approved": ["Vemurafenib", "Dabrafenib"],
                    "clinical_evidence": [],
                    "preclinical": [],
                    "resistance_mechanisms": []
                },
                "evidence_assessment": {
                    "overall_quality": "comprehensive",
                    "well_characterized": ["therapeutic response"],
                    "knowledge_gaps": [],
                    "conflicting_evidence": []
                },
                "key_references": ["PMID:12345"]
            }
            mock_response = AsyncMock()
            mock_response.choices = [AsyncMock()]
            mock_response.choices[0].message.content = json.dumps(response_json)
            mock_call.return_value = mock_response

            insight = await service.get_llm_insight(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                evidence_summary="BRAF V600E with FDA approved therapies.",
                has_clinical_trials=True,
            )

            assert "BRAF V600E" in insight.llm_summary
            assert insight.clinical_trials_available is True
            assert insight.functional_summary == "BRAF V600E is a well-characterized oncogenic mutation."
            assert insight.therapeutic_landscape["fda_approved"] == ["Vemurafenib", "Dabrafenib"]

    @pytest.mark.asyncio
    async def test_get_llm_insight_with_markdown(self):
        """Test insight generation with markdown-wrapped JSON."""
        service = LLMService()

        response_json = {
            "functional_summary": "Test summary for the variant.",
            "biological_context": "",
            "research_implications": "",
            "therapeutic_landscape": {},
            "evidence_assessment": {},
            "key_references": []
        }

        markdown_response = f"```json\n{json.dumps(response_json)}\n```"

        with patch("oncomind.llm.service.acompletion", new_callable=AsyncMock) as mock_call:
            mock_response = AsyncMock()
            mock_response.choices = [AsyncMock()]
            mock_response.choices[0].message.content = markdown_response
            mock_call.return_value = mock_response

            insight = await service.get_llm_insight(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                evidence_summary="Test evidence summary.",
            )

            assert insight.llm_summary == "Test summary for the variant."
            assert insight.functional_summary == "Test summary for the variant."

    @pytest.mark.asyncio
    async def test_llm_service_with_custom_temperature(self):
        """Test LLM service with custom temperature parameter."""
        custom_temp = 0.5
        service = LLMService(model="gpt-4o-mini", temperature=custom_temp)

        assert service.temperature == custom_temp
        assert service.model == "gpt-4o-mini"

        response_json = {
            "functional_summary": "Test summary for the variant.",
            "biological_context": "",
            "research_implications": "",
            "therapeutic_landscape": {},
            "evidence_assessment": {},
            "key_references": []
        }

        with patch("oncomind.llm.service.acompletion", new_callable=AsyncMock) as mock_call:
            mock_response = AsyncMock()
            mock_response.choices = [AsyncMock()]
            mock_response.choices[0].message.content = json.dumps(response_json)
            mock_call.return_value = mock_response

            await service.get_llm_insight(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                evidence_summary="Test evidence summary.",
            )

            # Verify temperature was passed to acompletion
            mock_call.assert_called_once()
            call_kwargs = mock_call.call_args[1]
            assert call_kwargs["temperature"] == custom_temp
            assert call_kwargs["model"] == "gpt-4o-mini"

    @pytest.mark.asyncio
    async def test_llm_failure_fallback(self):
        """Test that LLM failure still returns valid insight."""
        service = LLMService()

        with patch("oncomind.llm.service.acompletion", new_callable=AsyncMock) as mock_call:
            mock_call.side_effect = Exception("LLM API error")

            insight = await service.get_llm_insight(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                evidence_summary="Test evidence summary.",
            )

            # Should still get a valid insight with fallback values
            assert "BRAF V600E" in insight.llm_summary
            assert "failed" in insight.rationale.lower() or "error" in insight.rationale.lower()

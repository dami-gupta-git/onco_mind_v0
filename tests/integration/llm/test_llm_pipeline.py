"""Integration tests for the LLM pipeline.

Tests that:
1. A single LLM call is made (single-stage research dossier pipeline)
2. The prompt contains expected structure and keywords
3. The response is parsed correctly into LLMInsight
"""

import pytest
from unittest.mock import AsyncMock, patch
from oncomind.llm.service import LLMService


class TestSingleStagePipeline:
    """Tests for the single-stage research dossier LLM pipeline."""

    @pytest.mark.asyncio
    async def test_single_llm_call(self):
        """Should make exactly one LLM call (single-stage pipeline)."""
        synthesis_response = {
            "functional_impact": "BRAF V600E is an activating mutation.",
            "tumor_biology": "Found in 50% of melanomas.",
            "therapeutic_landscape_prose": "Vemurafenib is FDA-approved.",
            "therapeutic_landscape": {
                "fda_approved_exact": ["vemurafenib (melanoma)"],
                "fda_approved_codon": [],
                "fda_approved_gene": [],
                "clinical_evidence": [],
                "preclinical": [],
                "resistance_mechanisms": [],
            },
            "evidence_quality": "Comprehensive evidence.",
            "open_questions": ["CNS penetration?"],
            "research_program": [
                {
                    "aim_title": "Aim 1: CNS activity",
                    "rationale": "Test CNS penetration.",
                    "approaches": ["Use BrainMets model"],
                }
            ],
            "key_references": ["PMID:12345"],
            "evidence_tags": ["direct clinical data"],
        }

        call_count = 0
        captured_messages = []

        async def mock_acompletion(**kwargs):
            nonlocal call_count
            call_count += 1
            captured_messages.append(kwargs["messages"])

            import json

            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps(synthesis_response)),
                    finish_reason="stop",
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            result = await service.get_llm_insight(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                evidence_summary="FDA-approved, CIViC evidence",
                biological_context="50% melanoma prevalence",
                evidence_assessment={
                    "overall_quality": "comprehensive",
                    "well_characterized": ["FDA therapies"],
                    "knowledge_gaps": ["CNS resistance"],
                    "conflicting_evidence": [],
                },
                literature_summary="10 papers",
            )

        # Should have made exactly 1 call (single-stage)
        assert call_count == 1, f"Expected 1 LLM call, got {call_count}"

        # Result should be populated
        assert result.functional_summary == "BRAF V600E is an activating mutation."
        assert result.biological_context == "Found in 50% of melanomas."
        assert result.therapeutic_summary == "Vemurafenib is FDA-approved."
        assert len(result.research_hypotheses) >= 1

    @pytest.mark.asyncio
    async def test_llm_call_structure(self):
        """The LLM call should use correct message structure (system + user)."""
        captured_messages = []

        async def mock_acompletion(**kwargs):
            captured_messages.append(kwargs["messages"])

            import json

            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(
                        content=json.dumps(
                            {
                                "functional_impact": "test",
                                "tumor_biology": "",
                                "therapeutic_landscape_prose": "",
                                "therapeutic_landscape": {},
                                "evidence_quality": "",
                                "open_questions": [],
                                "research_program": [],
                                "key_references": [],
                                "evidence_tags": [],
                            }
                        )
                    ),
                    finish_reason="stop",
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            await service.get_llm_insight(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                evidence_summary="test",
                evidence_assessment={"overall_quality": "limited", "knowledge_gaps": []},
            )

        assert len(captured_messages) == 1
        messages = captured_messages[0]
        assert len(messages) == 2
        assert messages[0]["role"] == "system"
        assert messages[1]["role"] == "user"

        # User prompt contains gene/variant/tumor
        user_content = messages[1]["content"]
        assert "BRAF" in user_content
        assert "V600E" in user_content
        assert "Melanoma" in user_content

    @pytest.mark.asyncio
    async def test_result_on_llm_failure(self):
        """Should return a fallback LLMInsight when the LLM call fails."""

        async def mock_acompletion(**kwargs):
            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content="not valid json"),
                    finish_reason="stop",
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            result = await service.get_llm_insight(
                gene="TP53",
                variant="R248W",
                tumor_type="Breast",
                evidence_summary="test",
                evidence_assessment={"overall_quality": "limited", "knowledge_gaps": []},
            )

        # Should return a non-None LLMInsight (fallback)
        assert result is not None
        assert result.llm_summary is not None


class TestPromptContent:
    """Tests for prompt content structure and keywords."""

    @pytest.mark.asyncio
    async def test_synthesis_prompt_contains_calibration_rules(self):
        """Synthesis system prompt should contain calibration/constraint rules."""
        captured_system_prompt = None

        async def mock_acompletion(**kwargs):
            nonlocal captured_system_prompt
            captured_system_prompt = kwargs["messages"][0]["content"]

            import json

            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(
                        content=json.dumps(
                            {
                                "functional_impact": "",
                                "tumor_biology": "",
                                "therapeutic_landscape_prose": "",
                                "therapeutic_landscape": {},
                                "evidence_quality": "",
                                "open_questions": [],
                                "research_program": [],
                                "key_references": [],
                                "evidence_tags": [],
                            }
                        )
                    ),
                    finish_reason="stop",
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService()
            await service.get_llm_insight(
                gene="EGFR",
                variant="L858R",
                tumor_type="NSCLC",
                evidence_summary="test",
                evidence_assessment={
                    "overall_quality": "limited",
                    "knowledge_gaps": [],
                },
            )

        assert captured_system_prompt is not None
        # Prompt should contain key constraint language
        prompt_lower = captured_system_prompt.lower()
        assert "only" in prompt_lower  # "Use ONLY facts from..."
        assert "not" in prompt_lower  # "Do NOT use outside knowledge..."
        assert "research" in prompt_lower

    @pytest.mark.asyncio
    async def test_synthesis_prompt_contains_match_specificity(self):
        """Synthesis system prompt should contain match specificity guidance."""
        captured_system_prompt = None

        async def mock_acompletion(**kwargs):
            nonlocal captured_system_prompt
            captured_system_prompt = kwargs["messages"][0]["content"]

            import json

            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(
                        content=json.dumps(
                            {
                                "functional_impact": "",
                                "tumor_biology": "",
                                "therapeutic_landscape_prose": "",
                                "therapeutic_landscape": {},
                                "evidence_quality": "",
                                "open_questions": [],
                                "research_program": [],
                                "key_references": [],
                                "evidence_tags": [],
                            }
                        )
                    ),
                    finish_reason="stop",
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService()
            await service.get_llm_insight(
                gene="GNAQ",
                variant="Q209L",
                tumor_type="Uveal Melanoma",
                evidence_summary="test",
                evidence_assessment={
                    "overall_quality": "moderate",
                    "knowledge_gaps": [],
                },
            )

        assert captured_system_prompt is not None
        # Check match specificity guidance is present
        assert "variant-level" in captured_system_prompt.lower()
        assert "gene-level" in captured_system_prompt.lower()
        assert "codon-level" in captured_system_prompt.lower()

    @pytest.mark.asyncio
    async def test_synthesis_prompt_contains_research_program_spec(self):
        """Synthesis system prompt should describe the research program (aims) section."""
        captured_system_prompt = None

        async def mock_acompletion(**kwargs):
            nonlocal captured_system_prompt
            captured_system_prompt = kwargs["messages"][0]["content"]

            import json

            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(
                        content=json.dumps(
                            {
                                "functional_impact": "",
                                "tumor_biology": "",
                                "therapeutic_landscape_prose": "",
                                "therapeutic_landscape": {},
                                "evidence_quality": "",
                                "open_questions": [],
                                "research_program": [],
                                "key_references": [],
                                "evidence_tags": [],
                            }
                        )
                    ),
                    finish_reason="stop",
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService()
            await service.get_llm_insight(
                gene="KRAS",
                variant="G12C",
                tumor_type="NSCLC",
                evidence_summary="test",
                evidence_assessment={
                    "overall_quality": "moderate",
                    "knowledge_gaps": ["some gaps"],
                    "well_characterized": [],
                    "conflicting_evidence": [],
                },
            )

        assert captured_system_prompt is not None
        prompt_lower = captured_system_prompt.lower()
        # Should describe the aims/research program section
        assert "aim" in prompt_lower
        assert "research" in prompt_lower

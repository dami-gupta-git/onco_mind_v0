"""Integration tests for the two-stage LLM pipeline.

Tests that:
1. Two LLM calls are made (synthesis + hypothesis) when knowledge gaps exist
2. The prompts contain expected structure and keywords
3. Stage 2 is skipped when no knowledge gaps exist
"""

import pytest
from unittest.mock import AsyncMock, patch, call
from oncomind.llm.service import LLMService


class TestTwoStagePipeline:
    """Tests for the two-stage LLM pipeline."""

    @pytest.mark.asyncio
    async def test_two_llm_calls_when_gaps_exist(self):
        """Should make two LLM calls when knowledge gaps exist."""
        # Mock LLM responses
        synthesis_response = {
            "functional_summary": "BRAF V600E is an activating mutation.",
            "biological_context": "Found in 50% of melanomas.",
            "therapeutic_landscape": {
                "fda_approved": ["vemurafenib (variant-level, melanoma)"],
                "clinical_evidence": [],
                "preclinical": [],
                "resistance_mechanisms": [],
            },
            "evidence_assessment": {
                "overall_quality": "comprehensive",
                "well_characterized": ["FDA-approved therapies"],
                "knowledge_gaps": ["resistance patterns in CNS metastases"],
                "conflicting_evidence": [],
            },
            "key_references": ["PMID:12345"],
            "evidence_tags": ["direct clinical data"],
        }

        hypothesis_response = {
            "research_hypotheses": [
                "[Direct Clinical Data] Given CNS penetration data, test intracranial response rates.",
            ],
            "research_implications": "CNS-specific studies needed.",
        }

        call_count = 0
        captured_messages = []

        async def mock_acompletion(**kwargs):
            nonlocal call_count
            call_count += 1
            captured_messages.append(kwargs["messages"])

            # Return appropriate response based on call number
            if call_count == 1:
                content = synthesis_response
            else:
                content = hypothesis_response

            import json
            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps(content)),
                    finish_reason="stop"
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
                generate_hypotheses=True,
            )

        # Should have made exactly 2 calls
        assert call_count == 2, f"Expected 2 LLM calls, got {call_count}"

        # Verify synthesis prompt (call 1) structure
        synthesis_messages = captured_messages[0]
        assert len(synthesis_messages) == 2
        assert synthesis_messages[0]["role"] == "system"
        assert synthesis_messages[1]["role"] == "user"

        # Check synthesis system prompt keywords
        synthesis_system = synthesis_messages[0]["content"]
        assert "cancer genomics researcher" in synthesis_system.lower()
        assert "synthesis" in synthesis_system.lower() or "functional" in synthesis_system.lower()
        assert "calibration" in synthesis_system.lower()

        # Check synthesis user prompt keywords
        synthesis_user = synthesis_messages[1]["content"]
        assert "BRAF" in synthesis_user
        assert "V600E" in synthesis_user
        assert "Melanoma" in synthesis_user

        # Verify hypothesis prompt (call 2) structure
        hypothesis_messages = captured_messages[1]
        assert len(hypothesis_messages) == 2
        assert hypothesis_messages[0]["role"] == "system"
        assert hypothesis_messages[1]["role"] == "user"

        # Check hypothesis system prompt keywords
        hypothesis_system = hypothesis_messages[0]["content"]
        assert "hypothes" in hypothesis_system.lower()
        assert "research" in hypothesis_system.lower()
        assert "testable" in hypothesis_system.lower()

        # Check hypothesis user prompt keywords
        hypothesis_user = hypothesis_messages[1]["content"]
        assert "BRAF" in hypothesis_user
        assert "V600E" in hypothesis_user
        assert "gap" in hypothesis_user.lower()

        # Verify result contains hypothesis
        assert len(result.research_hypotheses) >= 1

    @pytest.mark.asyncio
    async def test_single_call_when_no_gaps(self):
        """Should make only one LLM call when no knowledge gaps exist."""
        synthesis_response = {
            "functional_summary": "BRAF V600E is well-characterized.",
            "biological_context": "Extensively studied.",
            "therapeutic_landscape": {
                "fda_approved": ["vemurafenib"],
                "clinical_evidence": [],
                "preclinical": [],
                "resistance_mechanisms": [],
            },
            "evidence_assessment": {
                "overall_quality": "comprehensive",
                "well_characterized": ["everything"],
                "knowledge_gaps": [],
                "conflicting_evidence": [],
            },
            "key_references": [],
            "evidence_tags": [],
        }

        call_count = 0

        async def mock_acompletion(**kwargs):
            nonlocal call_count
            call_count += 1

            import json
            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps(synthesis_response)),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            result = await service.get_llm_insight(
                gene="BRAF",
                variant="V600E",
                tumor_type="Melanoma",
                evidence_summary="FDA-approved",
                evidence_assessment={
                    "overall_quality": "comprehensive",
                    "well_characterized": ["everything"],
                    "knowledge_gaps": [],  # No gaps!
                    "conflicting_evidence": [],
                },
                generate_hypotheses=True,
            )

        # Should have made only 1 call (no stage 2 without gaps)
        assert call_count == 1, f"Expected 1 LLM call when no gaps, got {call_count}"
        assert result.research_hypotheses == []

    @pytest.mark.asyncio
    async def test_single_call_when_hypotheses_disabled(self):
        """Should make only one LLM call when generate_hypotheses=False."""
        synthesis_response = {
            "functional_summary": "Test summary",
            "biological_context": "",
            "therapeutic_landscape": {},
            "evidence_assessment": {"overall_quality": "limited"},
            "key_references": [],
            "evidence_tags": [],
        }

        call_count = 0

        async def mock_acompletion(**kwargs):
            nonlocal call_count
            call_count += 1

            import json
            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps(synthesis_response)),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService(model="gpt-4o-mini")
            await service.get_llm_insight(
                gene="TP53",
                variant="R248W",
                tumor_type="Breast",
                evidence_summary="Some evidence",
                evidence_assessment={
                    "overall_quality": "limited",
                    "well_characterized": [],
                    "knowledge_gaps": ["many gaps exist"],
                    "conflicting_evidence": [],
                },
                generate_hypotheses=False,  # Disabled!
            )

        # Should have made only 1 call (hypothesis disabled)
        assert call_count == 1, f"Expected 1 LLM call when hypotheses disabled, got {call_count}"


class TestPromptContent:
    """Tests for prompt content structure and keywords."""

    @pytest.mark.asyncio
    async def test_synthesis_prompt_contains_calibration_rules(self):
        """Synthesis system prompt should contain calibration rules."""
        captured_system_prompt = None

        async def mock_acompletion(**kwargs):
            nonlocal captured_system_prompt
            captured_system_prompt = kwargs["messages"][0]["content"]

            import json
            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps({
                        "functional_summary": "", "biological_context": "",
                        "therapeutic_landscape": {}, "evidence_assessment": {},
                        "key_references": [], "evidence_tags": []
                    })),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService()
            await service.get_llm_insight(
                gene="EGFR", variant="L858R", tumor_type="NSCLC",
                evidence_summary="test",
                evidence_assessment={"overall_quality": "limited", "knowledge_gaps": []},
            )

        # Check calibration rules are present
        assert "calibration" in captured_system_prompt.lower()
        assert "limited" in captured_system_prompt.lower() or "minimal" in captured_system_prompt.lower()
        assert "generic" in captured_system_prompt.lower()

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
                    message=AsyncMock(content=json.dumps({
                        "functional_summary": "", "biological_context": "",
                        "therapeutic_landscape": {}, "evidence_assessment": {},
                        "key_references": [], "evidence_tags": []
                    })),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService()
            await service.get_llm_insight(
                gene="GNAQ", variant="Q209L", tumor_type="Uveal Melanoma",
                evidence_summary="test",
                evidence_assessment={"overall_quality": "moderate", "knowledge_gaps": []},
            )

        # Check match specificity guidance
        assert "match specificity" in captured_system_prompt.lower()
        assert "variant-level" in captured_system_prompt.lower()
        assert "gene-level" in captured_system_prompt.lower()
        assert "codon-level" in captured_system_prompt.lower()

    @pytest.mark.asyncio
    async def test_hypothesis_prompt_contains_evidence_tags(self):
        """Hypothesis system prompt should require evidence basis tags."""
        captured_hypothesis_system = None

        call_count = 0
        async def mock_acompletion(**kwargs):
            nonlocal captured_hypothesis_system, call_count
            call_count += 1

            if call_count == 2:  # Hypothesis call
                captured_hypothesis_system = kwargs["messages"][0]["content"]

            import json
            if call_count == 1:
                response = {
                    "functional_summary": "test", "biological_context": "",
                    "therapeutic_landscape": {}, "evidence_assessment": {},
                    "key_references": [], "evidence_tags": []
                }
            else:
                response = {"research_hypotheses": [], "research_implications": ""}

            mock_response = AsyncMock()
            mock_response.choices = [
                AsyncMock(
                    message=AsyncMock(content=json.dumps(response)),
                    finish_reason="stop"
                )
            ]
            return mock_response

        with patch("oncomind.llm.service.acompletion", side_effect=mock_acompletion):
            service = LLMService()
            await service.get_llm_insight(
                gene="KRAS", variant="G12C", tumor_type="NSCLC",
                evidence_summary="test",
                evidence_assessment={
                    "overall_quality": "moderate",
                    "knowledge_gaps": ["some gaps"],
                    "well_characterized": [],
                    "conflicting_evidence": [],
                },
                generate_hypotheses=True,
            )

        assert captured_hypothesis_system is not None
        # Check evidence basis tags are required
        hypothesis_lower = captured_hypothesis_system.lower()
        assert "direct clinical data" in hypothesis_lower
        assert "preclinical data" in hypothesis_lower
        assert "pan-cancer extrapolation" in hypothesis_lower

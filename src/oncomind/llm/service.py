"""LLM service for variant insight generation.

Single-stage pipeline:
1. Synthesis - integrates evidence into a five-section research dossier covering
   functional impact, tumor biology, therapeutic landscape, evidence quality,
   and an emerging research program with concrete aims.
"""

import json
import re
import time
from dataclasses import dataclass

from litellm import acompletion

from oncomind.config.constants import (
    LLM_DEFAULT_MODEL,
    LLM_FAST_MODEL,
    LLM_MAX_PAPERS_FOR_EXTRACTION,
    LLM_MAX_TOKENS_KNOWLEDGE_EXTRACTION,
    LLM_MAX_TOKENS_PAPER_SCORING,
    LLM_MAX_TOKENS_SYNTHESIS,
    LLM_PAPER_ABSTRACT_TRUNCATION,
    LLM_PAPER_CONTENT_TRUNCATION,
    LLM_PAPER_RELEVANCE_THRESHOLD,
    LLM_TIMEOUT_CLAUDE,
    LLM_TIMEOUT_DEFAULT,
)
from oncomind.config.debug import get_logger
from oncomind.llm.prompts import (
    create_cross_source_prompt,
    create_synthesis_prompt,
)
from oncomind.llm.llm_insight import LLMInsight

logger = get_logger(__name__)


@dataclass
class LLMInsightInput:
    """Input data for LLM insight generation."""

    gene: str
    variant: str
    tumor_type: str | None
    evidence_summary: str
    biological_context: str = ""
    evidence_assessment: dict | None = None
    literature_summary: str = ""
    has_clinical_trials: bool = False
    data_availability: dict | None = None
    resistance_summary: str = ""
    sensitivity_summary: str = ""
    locus_match_summary: dict | None = None


def _empty_paper_scoring_result(key_finding: str = "", confidence: float = 0.0) -> dict:
    """Return empty/error result for paper scoring."""
    return {
        "relevance_score": 0.0,
        "is_relevant": False,
        "signal_type": "unclear",
        "drugs_mentioned": [],
        "key_finding": key_finding,
        "confidence": confidence,
    }


class LLMService:
    """LLM service for generating variant annotation narratives.

    The LLM synthesizes evidence from multiple databases into a research dossier
    covering functional impact, tumor biology, therapeutic landscape, evidence
    quality, and an emerging research program with concrete aims.
    """

    def __init__(self, model: str = LLM_DEFAULT_MODEL, temperature: float = 0.0):
        self.model = model
        self.temperature = temperature
        logger.debug(
            f"LLMService initialized with model={model}, temperature={temperature}"
        )

    async def _call_llm(
        self, messages: list[dict], max_tokens: int = LLM_MAX_TOKENS_SYNTHESIS
    ) -> dict | None:
        """Make LLM API call and parse JSON response.

        Args:
            messages: List of message dicts with role/content
            max_tokens: Maximum tokens for response

        Returns:
            Parsed JSON dict or None on error
        """
        timeout = (
            LLM_TIMEOUT_CLAUDE
            if "claude" in self.model.lower()
            else LLM_TIMEOUT_DEFAULT
        )

        completion_kwargs = {
            "model": self.model,
            "messages": messages,
            "temperature": self.temperature,
            "max_tokens": max_tokens,
            "timeout": timeout,
        }

        # Use JSON mode for OpenAI models
        if "gpt" in self.model.lower():
            completion_kwargs["response_format"] = {"type": "json_object"}

        input_chars = sum(len(m.get("content", "")) for m in messages)
        logger.debug(f"LLM request: model={self.model}, payload={input_chars} chars")

        try:
            t0 = time.time()
            response = await acompletion(**completion_kwargs)
            llm_time = time.time() - t0
            logger.info(f"LLM call completed in {llm_time:.2f}s")

            raw_content = response.choices[0].message.content.strip()
            finish_reason = response.choices[0].finish_reason

            if finish_reason == "length":
                logger.warning("LLM response truncated (finish_reason=length)")

            logger.debug(
                f"LLM raw response (finish_reason={finish_reason}):\n{raw_content[:500]}..."
            )

            # Parse JSON response
            return self._parse_json_response(raw_content)

        except Exception as e:
            logger.error(f"LLM call failed: {e}")
            return None

    def _parse_json_response(self, content: str) -> dict | None:
        """Parse JSON from LLM response, handling common formatting issues."""
        # Strip markdown code blocks if present
        if "```" in content:
            parts = content.split("```")
            for part in parts:
                stripped = part.strip()
                if stripped.lower().startswith("json"):
                    stripped = stripped[4:].lstrip()
                if stripped.startswith("{"):
                    content = stripped
                    break

        # Find JSON object if not at start
        if not content.strip().startswith("{"):
            start_idx = content.find("{")
            end_idx = content.rfind("}")
            if start_idx != -1 and end_idx != -1 and end_idx > start_idx:
                content = content[start_idx : end_idx + 1]

        content = content.strip()

        try:
            return json.loads(content)
        except json.JSONDecodeError as e:
            logger.warning(f"JSON parse failed: {e}. Attempting repair...")

            # Common fixes
            repaired = content
            repaired = repaired.replace("\\'", "'")
            repaired = re.sub(
                r'(?<!\\)\n(?=(?:[^"]*"[^"]*")*[^"]*"[^"]*$)', "\\n", repaired
            )
            repaired = re.sub(r",\s*([}\]])", r"\1", repaired)

            try:
                return json.loads(repaired)
            except json.JSONDecodeError:
                logger.error("JSON repair failed")
                return None

    async def get_llm_insight(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
        evidence_summary: str,
        biological_context: str = "",
        evidence_assessment: dict | None = None,
        literature_summary: str = "",
        has_clinical_trials: bool = False,
        data_availability: dict | None = None,
        resistance_summary: str = "",
        sensitivity_summary: str = "",
        locus_match_summary: dict | None = None,
        tumor_match_summary: dict | None = None,
    ) -> LLMInsight:
        """Generate variant insight using single-stage research dossier pipeline.

        Args:
            gene: Gene symbol (e.g., BRAF)
            variant: Variant notation (e.g., V600E)
            tumor_type: Tumor type context
            evidence_summary: Compact text summary of evidence
            biological_context: cBioPortal prevalence/co-mutation data
            evidence_assessment: Dict with keys: overall_quality, research_priority,
                well_characterized, knowledge_gaps, significant_gaps, conflicting_evidence
            literature_summary: PubMed literature findings
            has_clinical_trials: Whether clinical trials are available
            data_availability: Dict with boolean flags for data presence
            resistance_summary: Concise summary of resistance evidence
            sensitivity_summary: Concise summary of sensitivity evidence
            locus_match_summary: Dict with locus match specificity info
            tumor_match_summary: Dict with tumor match specificity info

        Returns:
            LLMInsight with five-section research dossier
        """
        # Default empty dicts if not provided
        if evidence_assessment is None:
            evidence_assessment = {
                "overall_quality": "unknown",
                "research_priority": "unknown",
                "well_characterized": [],
                "knowledge_gaps": [],
                "significant_gaps": [],
                "conflicting_evidence": [],
            }

        if data_availability is None:
            data_availability = {
                "has_tumor_specific_cbioportal": False,
                "has_civic_assertions": False,
                "has_fda_biomarker_evidence": False,
                "has_vicc_evidence": False,
            }

        # =====================================================================
        # SYNTHESIS: research dossier (single stage)
        # =====================================================================
        logger.info(f"Synthesizing research dossier for {gene} {variant}")

        synthesis_messages = create_synthesis_prompt(
            gene=gene,
            variant=variant,
            tumor_type=tumor_type,
            biological_context=biological_context,
            evidence_summary=evidence_summary,
            evidence_assessment=evidence_assessment,
            literature_summary=literature_summary,
            data_availability=data_availability,
            resistance_summary=resistance_summary,
            sensitivity_summary=sensitivity_summary,
            locus_match_summary=locus_match_summary,
            tumor_match_summary=tumor_match_summary,
        )

        print("\n" + "=" * 80)
        print("LLM CALL: RESEARCH DOSSIER SYNTHESIS")
        # print("=" * 80)
        # for msg in synthesis_messages:
        #     print(f"\n--- {msg['role'].upper()} ---")
        #     print(msg["content"])
        # print("=" * 80 + "\n")

        synthesis_data = await self._call_llm(
            synthesis_messages, max_tokens=LLM_MAX_TOKENS_SYNTHESIS
        )

        if synthesis_data is None:
            logger.error("Research dossier synthesis failed")
            return LLMInsight(
                llm_summary=f"Evidence summary for {gene} {variant}. See database annotations below.",
                rationale="LLM synthesis failed",
                clinical_trials_available=has_clinical_trials,
                therapeutic_evidence=[],
                references=[],
            )

        # =====================================================================
        # EXTRACT DOSSIER COMPONENTS
        # =====================================================================
        functional_impact = synthesis_data.get("functional_impact", "")
        tumor_biology = synthesis_data.get("tumor_biology", "")
        therapeutic_landscape_prose = synthesis_data.get("therapeutic_landscape_prose", "")
        therapeutic = synthesis_data.get("therapeutic_landscape", {})
        evidence_quality_text = synthesis_data.get("evidence_quality", "")
        open_questions = synthesis_data.get("open_questions", [])
        research_program = synthesis_data.get("research_program", [])
        evidence_tags = synthesis_data.get("evidence_tags", [])
        key_references = synthesis_data.get("key_references", [])

        # Build plain text summary from sections 1–3
        summary_parts = [p for p in [functional_impact, tumor_biology, therapeutic_landscape_prose] if p]
        llm_summary = " ".join(summary_parts) if summary_parts else "No summary available"

        return LLMInsight(
            llm_summary=llm_summary,
            rationale="",
            clinical_trials_available=has_clinical_trials,
            therapeutic_evidence=[],
            references=key_references,
            functional_summary=functional_impact or None,
            biological_context=tumor_biology or None,
            therapeutic_summary=therapeutic_landscape_prose or None,
            therapeutic_landscape=therapeutic or None,
            evidence_quality=evidence_assessment.get("overall_quality"),
            knowledge_gaps=open_questions,
            well_characterized=evidence_assessment.get("well_characterized", []),
            conflicting_evidence=evidence_assessment.get("conflicting_evidence", []),
            research_implications=evidence_quality_text or None,
            evidence_tags=evidence_tags,
            research_hypotheses=[aim.get("aim_title", "") for aim in research_program if aim.get("aim_title")],
            research_program=research_program,
        )

    async def get_cross_source_analysis(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
        gene_context: str,
        cross_source_synthesis: str,
        sensitivity_summary: str = "",
        resistance_summary: str = "",
    ) -> dict | None:
        """Generate cross-source synthesis analysis.

        This is a separate LLM call focused specifically on analyzing
        drug evidence across CGI, CIViC, VICC, and Literature sources.

        Args:
            gene: Gene symbol
            variant: Variant notation
            tumor_type: Tumor type context
            gene_context: Gene role, pathway, function information
            cross_source_synthesis: Pre-computed cross-source drug grouping
            sensitivity_summary: Summary of sensitivity signals
            resistance_summary: Summary of resistance signals

        Returns:
            Dict with keys: strongest_evidence, conflicting_signals,
            emerging_targets, key_gaps, summary. Or None on error.
        """
        if (
            not cross_source_synthesis
            or cross_source_synthesis == "No cross-source synthesis available."
        ):
            logger.debug("Skipping cross-source analysis - no synthesis data")
            return None

        logger.info(f"Generating cross-source analysis for {gene} {variant}")

        messages = create_cross_source_prompt(
            gene=gene,
            variant=variant,
            tumor_type=tumor_type,
            gene_context=gene_context,
            cross_source_synthesis=cross_source_synthesis,
            sensitivity_summary=sensitivity_summary,
            resistance_summary=resistance_summary,
        )

        result = await self._call_llm(messages, max_tokens=LLM_MAX_TOKENS_SYNTHESIS)

        if result is None:
            logger.warning("Cross-source analysis LLM call failed")
            return None

        logger.info(
            f"Cross-source analysis complete: {len(result.get('strongest_evidence', []))} strong, "
            f"{len(result.get('conflicting_signals', []))} conflicts, "
            f"{len(result.get('emerging_targets', []))} emerging"
        )

        return result

    async def score_paper_relevance(
        self,
        title: str,
        abstract: str | None,
        tldr: str | None,
        gene: str,
        variant: str,
        tumor_type: str | None,
    ) -> dict:
        """Score a paper's relevance to a specific gene/variant/tumor context.

        Uses LLM to determine if a paper is actually relevant to the specific
        clinical context, extracting key findings about resistance/sensitivity.

        Args:
            title: Paper title
            abstract: Paper abstract (may be None)
            tldr: AI-generated summary from Semantic Scholar (may be None)
            gene: Gene symbol (e.g., "KIT")
            variant: Variant notation (e.g., "D816V")
            tumor_type: Tumor type (e.g., "GIST", "Gastrointestinal Stromal Tumor")

        Returns:
            dict with keys:
                - relevance_score: float 0-1 (1 = directly about this variant in this tumor)
                - is_relevant: bool (True if score >= 0.6)
                - signal_type: "resistance", "sensitivity", "mixed", "prognostic", or "unclear"
                - drugs_mentioned: list of drug names affected by this variant
                - key_finding: one sentence summary of the paper's relevance
                - confidence: float 0-1 for the extraction confidence
        """
        # Use the best available text
        text_content = tldr or abstract or ""
        if not text_content:
            return _empty_paper_scoring_result("No abstract or summary available")

        tumor_context = tumor_type or "cancer (unspecified)"

        system_prompt = """You are an expert oncology literature analyst. Your task is to evaluate whether a scientific paper is relevant to understanding a specific gene variant in a specific tumor type.

Be INCLUSIVE for clinically relevant papers:
- Papers about the SAME EXON or SAME CODON are highly relevant (e.g., exon 17 papers for D816V)
- Papers about drugs targeting this mutation class are relevant (e.g., avapritinib for KIT mutations in GIST)
- Papers about resistance mechanisms in this tumor type are relevant
- Papers about related variants in the SAME gene and SAME tumor are relevant

Be STRICT only about tumor type:
- A paper about KIT D816V in mastocytosis is NOT relevant if we're asking about GIST
- A paper about a completely different gene is NOT relevant

CRITICAL - Distinguish PREDICTIVE vs PROGNOSTIC signals:
- PREDICTIVE (resistance/sensitivity): Paper shows variant PREDICTS response or resistance to a SPECIFIC drug
- PROGNOSTIC: Paper shows variant is associated with OUTCOMES (survival, recurrence) but NOT specific drug response

Return valid JSON only, no markdown."""

        user_prompt = f"""Evaluate this paper's relevance to {gene} {variant} in {tumor_context}:

TITLE: {title}

CONTENT: {text_content[:LLM_PAPER_CONTENT_TRUNCATION]}

Return JSON with these exact fields:
{{
    "relevance_score": <float 0-1>,
    "signal_type": "<resistance|sensitivity|mixed|prognostic|unclear>",
    "drugs_mentioned": [<list of specific drug names mentioned>],
    "key_finding": "<one sentence summary>",
    "confidence": <float 0-1>
}}"""

        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ]

        print("\n" + "=" * 80)
        print("LLM CALL: PAPER RELEVANCE SCORING")
        # print("=" * 80)
        # for msg in messages:
        #     print(f"\n--- {msg['role'].upper()} ---")
        #     print(msg["content"])
        # print("=" * 80 + "\n")

        try:
            # Use fast/cheap model for paper scoring (high volume operation)
            response = await acompletion(
                model=LLM_FAST_MODEL,
                messages=messages,
                temperature=0.0,
                max_tokens=LLM_MAX_TOKENS_PAPER_SCORING,
                response_format={"type": "json_object"},
            )

            content = response.choices[0].message.content.strip()
            data = self._parse_json_response(content)

            if data is None:
                return _empty_paper_scoring_result("Failed to parse LLM response")

            relevance_score = float(data.get("relevance_score", 0.0))
            relevance_score = max(0.0, min(1.0, relevance_score))

            return {
                "relevance_score": relevance_score,
                "is_relevant": relevance_score >= LLM_PAPER_RELEVANCE_THRESHOLD,
                "signal_type": data.get("signal_type", "unclear"),
                "drugs_mentioned": data.get("drugs_mentioned", []),
                "key_finding": data.get("key_finding", ""),
                "confidence": float(data.get("confidence", 0.5)),
            }

        except Exception as e:
            logger.error(f"Paper relevance scoring error: {e}")
            return _empty_paper_scoring_result(f"Error during analysis: {str(e)[:100]}")

    async def extract_variant_knowledge(
        self,
        gene: str,
        variant: str,
        tumor_type: str,
        paper_contents: list[dict],
    ) -> dict:
        """Extract structured knowledge about a variant from multiple papers.

        Uses LLM to synthesize information from paper abstracts/TLDRs into
        structured knowledge about the variant's clinical significance.

        Args:
            gene: Gene symbol
            variant: Variant notation
            tumor_type: Tumor type
            paper_contents: List of dicts with keys: title, abstract, tldr, pmid, url

        Returns:
            dict with structured knowledge about therapeutic implications
        """
        if not paper_contents:
            return None

        # Format paper contents for the prompt
        papers_text = []
        for i, paper in enumerate(paper_contents[:LLM_MAX_PAPERS_FOR_EXTRACTION], 1):
            content = paper.get("tldr") or paper.get("abstract") or ""
            pmid = paper.get("pmid", "Unknown")
            title = paper.get("title", "Untitled")
            papers_text.append(f"""
Paper {i} (PMID: {pmid}):
Title: {title}
Content: {content[:LLM_PAPER_ABSTRACT_TRUNCATION]}
""")

        papers_combined = "\n".join(papers_text)

        system_prompt = """You are an expert oncology researcher synthesizing knowledge from scientific literature.

Your task is to extract structured, clinically actionable information about a specific gene variant from research papers.

CRITICAL DISTINCTION - PREDICTIVE vs PROGNOSTIC:
- PREDICTIVE: Variant predicts response to a SPECIFIC TARGETED THERAPY
- PROGNOSTIC: Variant associated with outcomes but NOT specific drug response

CRITICAL DISTINCTION - MATCH LEVEL:
For each drug response finding, determine the MATCH LEVEL:
- "variant": The paper discusses THIS EXACT VARIANT (e.g., "EGFR L858R shows sensitivity to...")
- "codon": The paper discusses the same codon position but different amino acid (e.g., paper discusses V600K when queried about V600E)
- "gene": The paper discusses the gene generally (e.g., "EGFR mutations respond to...") without specifying this exact variant

Default to "gene" if the paper doesn't clearly specify the exact variant.

Be PRECISE and EVIDENCE-BASED:
- Only report findings directly supported by the papers
- Distinguish between preclinical and clinical evidence
- Note the strength of evidence

Return valid JSON only, no markdown."""

        user_prompt = f"""Extract structured knowledge about {gene} {variant} in {tumor_type} from these papers:

{papers_combined}

Return JSON with these exact fields:
{{
    "mutation_type": "<primary|secondary|both|unknown>",
    "resistant_to": [
        {{"drug": "<drug name>", "evidence": "<preclinical|clinical|FDA-labeled>", "mechanism": "<mechanism>", "locus_match": "<variant|codon|gene>"}}
    ],
    "sensitive_to": [
        {{"drug": "<drug name>", "evidence": "<preclinical|clinical|FDA-labeled>", "locus_match": "<variant|codon|gene>"}}
    ],
    "clinical_significance": "<2-3 sentence summary>",
    "evidence_level": "<FDA-approved|Phase 3|Phase 2|Preclinical|Case reports|None>",
    "references": ["<PMID1>", "<PMID2>"],
    "key_findings": ["<finding 1>", "<finding 2>"],
    "confidence": <0.0-1.0>
}}"""

        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ]

        print("\n" + "=" * 80)
        print("LLM CALL: VARIANT KNOWLEDGE EXTRACTION")
        print("=" * 80)
        # for msg in messages:
        #     print(f"\n--- {msg['role'].upper()} ---")
        #     print(msg["content"])
        # print("=" * 80 + "\n")

        try:
            # Use fast/cheap model for knowledge extraction (high volume operation)
            response = await acompletion(
                model=LLM_FAST_MODEL,
                messages=messages,
                temperature=0.0,
                max_tokens=LLM_MAX_TOKENS_KNOWLEDGE_EXTRACTION,
                response_format={"type": "json_object"},
            )

            content = response.choices[0].message.content.strip()
            data = self._parse_json_response(content)

            if data is None:
                logger.error("Failed to parse variant knowledge extraction response")
                return None

            return {
                "mutation_type": data.get("mutation_type", "unknown"),
                "resistant_to": data.get("resistant_to", []),
                "sensitive_to": data.get("sensitive_to", []),
                "clinical_significance": data.get("clinical_significance", ""),
                "evidence_level": data.get("evidence_level", "None"),
                "references": data.get("references", []),
                "key_findings": data.get("key_findings", []),
                "confidence": float(data.get("confidence", 0.5)),
            }

        except Exception as e:
            logger.error(f"Variant knowledge extraction error: {e}")
            return None

    async def test_llm_baseline_knowledge(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None = None,
    ) -> dict:
        """Test what the LLM knows from training data alone (no evidence provided).

        This is a diagnostic method to understand what the LLM's baseline knowledge
        is for a given variant, without any database evidence. Useful for detecting
        when LLM outputs come from training data vs provided evidence.

        Args:
            gene: Gene symbol (e.g., "PIK3CA")
            variant: Variant notation (e.g., "H1047R")
            tumor_type: Optional tumor type context

        Returns:
            dict with LLM's baseline knowledge about the variant
        """
        tumor_context = f" in {tumor_type}" if tumor_type else ""

        system_prompt = """You are a cancer genomics expert. Answer based ONLY on your training knowledge.
Do NOT say you need more information. Tell me everything you know from your training data."""

        user_prompt = f"""What do you know about {gene} {variant}{tumor_context}?

Return JSON with:
{{
    "variant_function": "What does this variant do molecularly?",
    "oncogenic_mechanism": "How does it drive cancer?",
    "clinical_significance": "What is known clinically?",
    "fda_approved_drugs": ["List any FDA-approved drugs for this variant"],
    "resistance_mechanisms": ["Known resistance mechanisms"],
    "sensitivity_mechanisms": ["Known drug sensitivities"],
    "prevalence": "How common is this variant?",
    "key_facts": ["Other important facts"],
    "confidence": "How confident are you in this information? (high/medium/low)"
}}"""

        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ]

        # print("\n" + "=" * 80)
        # print("LLM BASELINE KNOWLEDGE TEST (NO EVIDENCE PROVIDED)")
        # print(f"Gene: {gene}, Variant: {variant}, Tumor: {tumor_type or 'None'}")
        # print("=" * 80 + "\n")

        try:
            response = await acompletion(
                model=self.model,
                messages=messages,
                temperature=0.0,
                max_tokens=2000,
            )

            content = response.choices[0].message.content.strip()
            # print("RAW LLM RESPONSE:")
            # print(content)
            # print("=" * 80 + "\n")

            data = self._parse_json_response(content)
            return (
                data if data else {"error": "Failed to parse response", "raw": content}
            )

        except Exception as e:
            logger.error(f"Baseline knowledge test error: {e}")
            return {"error": str(e)}

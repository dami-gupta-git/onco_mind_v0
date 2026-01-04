"""LLM service for variant insight generation.

Two-stage pipeline:
1. Synthesis - integrates evidence into structured summary
2. Hypothesis - generates testable research questions (optional, based on synthesis)
"""

import json
import re
import time
from litellm import acompletion
from oncomind.llm.prompts import create_synthesis_prompt, create_hypothesis_prompt
from oncomind.models.llm_insight import LLMInsight

from oncomind.config.debug import get_logger

logger = get_logger(__name__)


class LLMService:
    """LLM service for generating variant annotation narratives.

    The LLM synthesizes evidence from multiple databases into a clear,
    human-readable summary of the variant's clinical significance.

    Uses a two-stage pipeline:
    1. Synthesis: evidence → structured summary (functional, biological, therapeutic)
    2. Hypothesis: synthesis + gaps → research hypotheses
    """

    def __init__(self, model: str = "gpt-4o-mini", temperature: float = 0.0):
        self.model = model
        self.temperature = temperature
        logger.debug(f"LLMService initialized with model={model}, temperature={temperature}")

    async def _call_llm(self, messages: list[dict], max_tokens: int = 1500) -> dict | None:
        """Make LLM API call and parse JSON response.

        Args:
            messages: List of message dicts with role/content
            max_tokens: Maximum tokens for response

        Returns:
            Parsed JSON dict or None on error
        """
        timeout = 120 if "claude" in self.model.lower() else 60

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

            logger.debug(f"LLM raw response (finish_reason={finish_reason}):\n{raw_content[:500]}...")

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
                content = content[start_idx:end_idx + 1]

        content = content.strip()

        try:
            return json.loads(content)
        except json.JSONDecodeError as e:
            logger.warning(f"JSON parse failed: {e}. Attempting repair...")

            # Common fixes
            repaired = content
            repaired = repaired.replace("\\'", "'")
            repaired = re.sub(r'(?<!\\)\n(?=(?:[^"]*"[^"]*")*[^"]*"[^"]*$)', '\\n', repaired)
            repaired = re.sub(r',\s*([}\]])', r'\1', repaired)

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
        match_level_summary: dict | None = None,
        generate_hypotheses: bool = True,
    ) -> LLMInsight:
        """Generate variant insight using two-stage LLM pipeline.

        Stage 1: Synthesis - integrates evidence into structured summary
        Stage 2: Hypothesis - generates research questions (if generate_hypotheses=True)

        Args:
            gene: Gene symbol (e.g., BRAF)
            variant: Variant notation (e.g., V600E)
            tumor_type: Tumor type context
            evidence_summary: Compact text summary of evidence
            biological_context: cBioPortal prevalence/co-mutation data
            evidence_assessment: Dict with keys: overall_quality, well_characterized,
                knowledge_gaps, conflicting_evidence
            literature_summary: PubMed literature findings
            has_clinical_trials: Whether clinical trials are available
            data_availability: Dict with boolean flags for data presence
            resistance_summary: Concise summary of resistance evidence
            sensitivity_summary: Concise summary of sensitivity evidence
            match_level_summary: Dict with match specificity info
            generate_hypotheses: Whether to run stage 2 (default True)

        Returns:
            LLMInsight with synthesized narrative and research hypotheses
        """
        # Default empty dicts if not provided
        if evidence_assessment is None:
            evidence_assessment = {
                "overall_quality": "unknown",
                "well_characterized": [],
                "knowledge_gaps": [],
                "conflicting_evidence": [],
            }

        if data_availability is None:
            data_availability = {
                "has_tumor_specific_cbioportal": False,
                "has_civic_assertions": False,
                "has_fda_approvals": False,
                "has_vicc_evidence": False,
            }

        # =====================================================================
        # STAGE 1: SYNTHESIS
        # =====================================================================
        logger.info(f"Stage 1: Synthesizing evidence for {gene} {variant}")

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
            match_level_summary=match_level_summary,
        )

        synthesis_data = await self._call_llm(synthesis_messages, max_tokens=1500)

        if synthesis_data is None:
            logger.error("Stage 1 synthesis failed")
            return LLMInsight(
                llm_summary=f"Evidence summary for {gene} {variant}. See database annotations below.",
                rationale="LLM synthesis failed",
                clinical_trials_available=has_clinical_trials,
                therapeutic_evidence=[],
                references=[],
            )

        # Extract synthesis components
        functional_summary = synthesis_data.get("functional_summary", "")
        biological_context_text = synthesis_data.get("biological_context", "")
        therapeutic = synthesis_data.get("therapeutic_landscape", {})
        evidence_assess = synthesis_data.get("evidence_assessment", {})
        evidence_tags = synthesis_data.get("evidence_tags", [])
        key_references = synthesis_data.get("key_references", [])

        # =====================================================================
        # STAGE 2: HYPOTHESIS GENERATION (optional)
        # =====================================================================
        research_hypotheses = []
        research_implications = ""

        if generate_hypotheses and evidence_assessment.get("knowledge_gaps"):
            logger.info(f"Stage 2: Generating hypotheses for {gene} {variant}")

            # Build therapeutic signals summary for hypothesis context
            therapeutic_signals = f"Sensitivity: {sensitivity_summary}\nResistance: {resistance_summary}"

            hypothesis_messages = create_hypothesis_prompt(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                synthesis_result=synthesis_data,
                evidence_assessment=evidence_assessment,
                therapeutic_signals=therapeutic_signals,
            )

            hypothesis_data = await self._call_llm(hypothesis_messages, max_tokens=800)

            if hypothesis_data:
                research_hypotheses = hypothesis_data.get("research_hypotheses", [])
                research_implications = hypothesis_data.get("research_implications", "")
                logger.info(f"Generated {len(research_hypotheses)} hypotheses")
            else:
                logger.warning("Stage 2 hypothesis generation failed")
        else:
            logger.debug("Skipping hypothesis generation (no knowledge gaps or disabled)")

        # =====================================================================
        # BUILD FINAL INSIGHT
        # =====================================================================

        # Build plain text summary
        summary_parts = []
        if functional_summary:
            summary_parts.append(functional_summary)
        if biological_context_text:
            summary_parts.append(biological_context_text)
        if research_implications:
            summary_parts.append(research_implications)

        llm_summary = " ".join(summary_parts) if summary_parts else "No summary available"

        return LLMInsight(
            llm_summary=llm_summary,
            rationale=research_implications,
            clinical_trials_available=has_clinical_trials,
            therapeutic_evidence=[],
            references=key_references,
            # Raw component data for UI formatting
            functional_summary=functional_summary or None,
            biological_context=biological_context_text or None,
            therapeutic_landscape=therapeutic or None,
            # Research assessment fields
            evidence_quality=evidence_assess.get("overall_quality"),
            knowledge_gaps=evidence_assess.get("knowledge_gaps", []),
            well_characterized=evidence_assess.get("well_characterized", []),
            conflicting_evidence=evidence_assess.get("conflicting_evidence", []),
            research_implications=research_implications,
            evidence_tags=evidence_tags,
            research_hypotheses=research_hypotheses,
        )

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
            return {
                "relevance_score": 0.0,
                "is_relevant": False,
                "signal_type": "unclear",
                "drugs_mentioned": [],
                "key_finding": "No abstract or summary available",
                "confidence": 0.0,
            }

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

CONTENT: {text_content[:1500]}

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

        try:
            response = await acompletion(
                model="gpt-4o-mini",
                messages=messages,
                temperature=0.0,
                max_tokens=500,
                response_format={"type": "json_object"},
            )

            content = response.choices[0].message.content.strip()

            if content.startswith("```"):
                parts = content.split("```")
                content = parts[1] if len(parts) > 1 else parts[0]
                if content.lower().startswith("json"):
                    content = content[4:].lstrip()

            data = json.loads(content)

            relevance_score = float(data.get("relevance_score", 0.0))
            relevance_score = max(0.0, min(1.0, relevance_score))

            return {
                "relevance_score": relevance_score,
                "is_relevant": relevance_score >= 0.6,
                "signal_type": data.get("signal_type", "unclear"),
                "drugs_mentioned": data.get("drugs_mentioned", []),
                "key_finding": data.get("key_finding", ""),
                "confidence": float(data.get("confidence", 0.5)),
            }

        except Exception as e:
            logger.error(f"Paper relevance scoring error: {e}")
            return {
                "relevance_score": None,
                "is_relevant": False,
                "signal_type": "unclear",
                "drugs_mentioned": [],
                "key_finding": f"Error during analysis: {str(e)[:100]}",
                "confidence": 0.0,
            }

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
        for i, paper in enumerate(paper_contents[:5], 1):
            content = paper.get("tldr") or paper.get("abstract") or ""
            pmid = paper.get("pmid", "Unknown")
            title = paper.get("title", "Untitled")
            papers_text.append(f"""
Paper {i} (PMID: {pmid}):
Title: {title}
Content: {content[:1000]}
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
        {{"drug": "<drug name>", "evidence": "<preclinical|clinical|FDA-labeled>", "mechanism": "<mechanism>", "match_level": "<variant|codon|gene>"}}
    ],
    "sensitive_to": [
        {{"drug": "<drug name>", "evidence": "<preclinical|clinical|FDA-labeled>", "match_level": "<variant|codon|gene>"}}
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

        try:
            response = await acompletion(
                model="gpt-4o-mini",
                messages=messages,
                temperature=0.0,
                max_tokens=1500,
                response_format={"type": "json_object"},
            )

            content = response.choices[0].message.content.strip()

            if content.startswith("```"):
                parts = content.split("```")
                content = parts[1] if len(parts) > 1 else parts[0]
                if content.lower().startswith("json"):
                    content = content[4:].lstrip()

            data = json.loads(content)

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

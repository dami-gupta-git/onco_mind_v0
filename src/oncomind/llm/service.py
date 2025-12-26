"""LLM service for variant insight generation."""

import json
from litellm import acompletion
from oncomind.llm.prompts import create_research_prompt
from oncomind.models import RecommendedTherapy
from oncomind.models.llm_insight import LLMInsight
from oncomind.models.gene_context import get_oncogene_mutation_class

from oncomind.utils.logging_config import get_logger


class LLMService:
    """LLM service for generating variant annotation narratives.

    The LLM synthesizes evidence from multiple databases into a clear,
    human-readable summary of the variant's clinical significance.
    """

    def __init__(self, model: str = "gpt-4o-mini", temperature: float = 0.0, enable_logging: bool = False):
        self.model = model
        self.temperature = temperature
        self.enable_logging = enable_logging
        self.logger = get_logger(enable_console_logging=enable_logging) if enable_logging else None

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
    ) -> LLMInsight:
        """Generate variant insight by synthesizing evidence with LLM.

        Args:
            gene: Gene symbol (e.g., BRAF)
            variant: Variant notation (e.g., V600E)
            tumor_type: Tumor type context
            evidence_summary: Compact text summary of evidence
            biological_context: cBioPortal prevalence/co-mutation data
            evidence_assessment: Dict with keys: overall_quality, well_characterized,
                knowledge_gaps, conflicting_evidence, critical_gaps, significant_gaps
            literature_summary: PubMed literature findings (for full mode)
            has_clinical_trials: Whether clinical trials are available

        Returns:
            LLMInsight with LLM-generated research-focused narrative
        """
        # Default empty evidence assessment if not provided
        if evidence_assessment is None:
            evidence_assessment = {
                "overall_quality": "unknown",
                "well_characterized": [],
                "knowledge_gaps": [],
                "conflicting_evidence": [],
                "critical_gaps": [],
                "significant_gaps": [],
            }

        # Create research-focused prompt
        messages = create_research_prompt(
            gene=gene,
            variant=variant,
            tumor_type=tumor_type,
            biological_context=biological_context,
            evidence_summary=evidence_summary,
            evidence_assessment=evidence_assessment,
            literature_summary=literature_summary,
        )

        # Call LLM for narrative generation
        completion_kwargs = {
            "model": self.model,
            "messages": messages,
            "temperature": self.temperature,
            "max_tokens": 1000,
        }

        # Use JSON mode for OpenAI models
        if "gpt" in self.model.lower():
            completion_kwargs["response_format"] = {"type": "json_object"}

        try:
            response = await acompletion(**completion_kwargs)
            raw_content = response.choices[0].message.content.strip()

            # Parse JSON response
            content = raw_content
            if content.startswith("```"):
                parts = content.split("```")
                content = parts[1] if len(parts) > 1 else parts[0]
                if content.lower().startswith("json"):
                    content = content[4:].lstrip()

            data = json.loads(content)

            # Build research-focused LLM summary from response components
            functional_summary = data.get("functional_summary", "")
            biological_context_text = data.get("biological_context", "")
            research_implications = data.get("research_implications", "")

            # Combine into a comprehensive summary
            summary_parts = []
            if functional_summary:
                summary_parts.append(f"**Functional Impact:** {functional_summary}")
            if biological_context_text:
                summary_parts.append(f"**Biological Context:** {biological_context_text}")

            # Add therapeutic landscape
            therapeutic = data.get("therapeutic_landscape", {})
            if therapeutic:
                therapy_parts = []
                if therapeutic.get("fda_approved"):
                    therapy_parts.append(f"FDA-approved: {', '.join(therapeutic['fda_approved'])}")
                if therapeutic.get("clinical_evidence"):
                    therapy_parts.append(f"Clinical evidence: {', '.join(therapeutic['clinical_evidence'])}")
                if therapeutic.get("preclinical"):
                    therapy_parts.append(f"Preclinical: {', '.join(therapeutic['preclinical'])}")
                if therapeutic.get("resistance_mechanisms"):
                    therapy_parts.append(f"Resistance: {', '.join(therapeutic['resistance_mechanisms'])}")
                if therapy_parts:
                    summary_parts.append(f"**Therapeutic Landscape:** {'; '.join(therapy_parts)}")

            if research_implications:
                summary_parts.append(f"**Research Implications:** {research_implications}")

            llm_summary = "\n\n".join(summary_parts) if summary_parts else "No summary available"

            # Extract evidence assessment
            evidence_assessment = data.get("evidence_assessment", {})
            evidence_quality = evidence_assessment.get("overall_quality")
            knowledge_gaps = evidence_assessment.get("knowledge_gaps", [])
            well_characterized = evidence_assessment.get("well_characterized", [])
            conflicting_evidence = evidence_assessment.get("conflicting_evidence", [])

            #print(f"{evidence_assessment=},{evidence_quality=}, {knowledge_gaps=},{well_characterized=},{conflicting_evidence=}")

            # Extract evidence tags for transparency
            evidence_tags = data.get("evidence_tags", [])

            # Build insight with research-focused narrative
            return LLMInsight(
                llm_summary=llm_summary,
                rationale=research_implications,
                clinical_trials_available=has_clinical_trials,
                therapeutic_evidence=[],  # Therapeutic evidence comes from Evidence.get_therapeutic_evidence()
                references=data.get("key_references", []),
                evidence_quality=evidence_quality,
                knowledge_gaps=knowledge_gaps,
                well_characterized=well_characterized,
                conflicting_evidence=conflicting_evidence,
                research_implications=research_implications,
                evidence_tags=evidence_tags,
            )

        except Exception as e:
            # On LLM failure, return insight with basic error message
            return LLMInsight(
                llm_summary=f"Evidence summary for {gene} {variant}. See database annotations below.",
                rationale=f"LLM narrative generation failed: {str(e)}",
                clinical_trials_available=has_clinical_trials,
                therapeutic_evidence=[],
                references=[],
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
            print(f"Paper relevance scoring error: {e}")
            return {
                "relevance_score": 0.5,
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
            return {
                "mutation_type": "unknown",
                "resistant_to": [],
                "sensitive_to": [],
                "clinical_significance": "No literature available for analysis",
                "evidence_level": "None",
                "references": [],
                "confidence": 0.0,
            }

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
        {{"drug": "<drug name>", "evidence": "<preclinical|clinical|FDA-labeled>", "mechanism": "<mechanism>"}}
    ],
    "sensitive_to": [
        {{"drug": "<drug name>", "evidence": "<preclinical|clinical|FDA-labeled>"}}
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
            print(f"Variant knowledge extraction error: {e}")
            return {
                "mutation_type": "unknown",
                "resistant_to": [],
                "sensitive_to": [],
                "clinical_significance": f"Error during extraction: {str(e)[:100]}",
                "evidence_level": "None",
                "references": [],
                "key_findings": [],
                "confidence": 0.0,
            }

"""LLM service for variant insight generation."""

import json
from litellm import acompletion
from oncomind.llm.prompts import create_annotation_prompt
from oncomind.models import EvidenceForLLM
from oncomind.models.insight import LLMInsight, RecommendedTherapy
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

    async def get_variant_insight(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
        evidence: EvidenceForLLM,
    ) -> LLMInsight:
        """Generate variant insight by synthesizing evidence with LLM.

        Args:
            gene: Gene symbol (e.g., BRAF)
            variant: Variant notation (e.g., V600E)
            tumor_type: Tumor type context
            evidence: Aggregated evidence from databases

        Returns:
            LLMInsight with annotations and LLM-generated summary
        """
        # Get evidence summary for context
        evidence_summary = evidence.summary_compact(tumor_type=tumor_type)

        # Check for oncogene mutation class therapy notes
        therapy_note = None
        mutation_class = get_oncogene_mutation_class(gene, variant)
        if mutation_class:
            # Check for tumor-specific therapy note first (most relevant)
            tumor_specific = mutation_class.get("tumor_specific", {})
            tumor_note = None
            if tumor_type:
                tumor_lower = tumor_type.lower()
                for tumor_key, note in tumor_specific.items():
                    if tumor_key in tumor_lower or tumor_lower in tumor_key:
                        tumor_note = note
                        break

            # Build therapy note from mutation class info
            notes = []
            if tumor_note:
                notes.append(tumor_note)
            else:
                if mutation_class.get("note"):
                    notes.append(mutation_class["note"])
                if mutation_class.get("mechanism"):
                    notes.append(f"Mechanism: {mutation_class['mechanism']}")
                if mutation_class.get("drugs"):
                    drugs_str = ", ".join(mutation_class["drugs"][:3])
                    notes.append(f"Therapeutic options: {drugs_str}")
            if notes:
                therapy_note = " | ".join(notes)

        # Create prompt for LLM
        messages = create_annotation_prompt(
            gene=gene,
            variant=variant,
            tumor_type=tumor_type,
            evidence_summary=evidence_summary,
            therapy_note=therapy_note,
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

            # Extract recommended therapies from LLM response
            therapies = []
            for therapy_data in data.get("recommended_therapies", []):
                therapies.append(RecommendedTherapy(
                    drug_name=therapy_data.get("drug_name", ""),
                    evidence_level=therapy_data.get("evidence_level"),
                    approval_status=therapy_data.get("approval_status"),
                    clinical_context=therapy_data.get("clinical_context"),
                ))

            # Determine evidence strength from evidence
            evidence_strength = self._compute_evidence_strength(evidence, tumor_type)

            # Build insight with LLM narrative
            insight = LLMInsight(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                summary=data.get("summary", "No summary available"),
                rationale=data.get("rationale", ""),
                evidence_strength=evidence_strength,
                clinical_trials_available=bool(evidence.clinical_trials),
                recommended_therapies=therapies,
                references=data.get("references", []),
                **evidence.model_dump(include={
                    'cosmic_id', 'ncbi_gene_id', 'dbsnp_id', 'clinvar_id',
                    'clinvar_clinical_significance', 'clinvar_accession',
                    'hgvs_genomic', 'hgvs_protein', 'hgvs_transcript',
                    'snpeff_effect', 'polyphen2_prediction', 'cadd_score', 'gnomad_exome_af',
                    'alphamissense_score', 'alphamissense_prediction',
                    'transcript_id', 'transcript_consequence'
                })
            )

            return insight

        except Exception as e:
            # On LLM failure, return insight with basic evidence summary
            return LLMInsight(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                summary=f"Evidence summary for {gene} {variant}. See database annotations below.",
                rationale=f"LLM narrative generation failed: {str(e)}",
                evidence_strength=self._compute_evidence_strength(evidence, tumor_type),
                clinical_trials_available=bool(evidence.clinical_trials),
                recommended_therapies=[],
                references=[],
                **evidence.model_dump(include={
                    'cosmic_id', 'ncbi_gene_id', 'dbsnp_id', 'clinvar_id',
                    'clinvar_clinical_significance', 'clinvar_accession',
                    'hgvs_genomic', 'hgvs_protein', 'hgvs_transcript',
                    'snpeff_effect', 'polyphen2_prediction', 'cadd_score', 'gnomad_exome_af',
                    'alphamissense_score', 'alphamissense_prediction',
                    'transcript_id', 'transcript_consequence'
                })
            )

    def _compute_evidence_strength(self, evidence: EvidenceForLLM, tumor_type: str | None) -> str:
        """Compute evidence strength from aggregated evidence.

        Returns "Strong", "Moderate", or "Weak" based on:
        - FDA approvals
        - CIViC evidence levels
        - ClinVar clinical significance
        """
        # Check for FDA approval
        if evidence.fda_approvals:
            for approval in evidence.fda_approvals:
                if tumor_type:
                    parsed = approval.parse_indication_for_tumor(tumor_type)
                    if parsed.get('tumor_match'):
                        return "Strong"
                else:
                    return "Strong"

        # Check CIViC evidence levels
        has_level_a = any(ev.evidence_level == 'A' for ev in evidence.civic)
        has_level_b = any(ev.evidence_level == 'B' for ev in evidence.civic)

        if has_level_a:
            return "Strong"
        if has_level_b:
            return "Moderate"

        # Check ClinVar pathogenicity
        if evidence.clinvar_clinical_significance:
            sig = evidence.clinvar_clinical_significance.lower()
            if 'pathogenic' in sig and 'likely' not in sig:
                return "Moderate"

        # Check CGI biomarkers
        if any(b.fda_approved for b in evidence.cgi_biomarkers):
            return "Strong"

        return "Weak"

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

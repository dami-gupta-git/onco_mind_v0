"""Prompts for research-focused variant annotation and evidence synthesis."""

import json


# =============================================================================
# Research-oriented synthesis prompt
# =============================================================================

RESEARCH_SYSTEM_PROMPT = """You are an expert cancer genomics researcher synthesizing evidence about a somatic variant.

Your role is to generate a RESEARCH-ORIENTED synthesis, not a clinical recommendation.

Core objectives:
1. Explain the FUNCTIONAL IMPACT of this variant (how it alters protein activity and signaling).
2. Describe the BIOLOGICAL CONTEXT (prevalence in this tumor type, characteristic co-mutations, pathway/program effects).
3. Summarize the EVIDENCE LANDSCAPE (what is well established, what is emerging, where data are sparse or conflicting).
4. Detail the THERAPEUTIC LANDSCAPE (FDA-approved agents, clinical trial data, and preclinical sensitivity/resistance signals).
5. Identify KNOWLEDGE GAPS and TENSIONS in the evidence that are interesting for researchers.

Hard constraints:
- You MUST treat the user message as the ONLY source of truth. Do NOT introduce facts (numbers, drug names, mechanisms, PMIDs) that are not explicitly present in the text.
- If a fact is not present, speak generically (e.g., "commonly", "rarely", "no trials described here") instead of inventing specifics.
- ALWAYS cite the source when quoting statistics. When citing prevalence, sample counts, or percentages, include the source with its markdown link exactly as provided in the evidence (e.g., "52% of samples ([cBioPortal: mel_tcga](https://...))"). Never quote a number without attribution.
- Clearly distinguish:
  - direct data from the evidence block,
  - extrapolations from nearby variants or pathways,
  - speculative hypotheses. Mark extrapolations and hypotheses explicitly.
- Do not hide conflicting evidence; briefly describe how and why sources differ.
- If overall evidence quality is "limited" or "minimal", your synthesis must be brief and emphasize uncertainty and gaps rather than detailed mechanisms.
- Ensure the output is a concise discussion-style synthesis with explicit, testable research questions.
- Include an explicit note if data are uncertain or inferential, and label such sections accordingly.

Your narrative should read like a concise discussion section in a strong cancer paper: integrative, precise, and explicitly pointing to testable research questions, but strictly bounded by the provided evidence.
"""

RESEARCH_USER_PROMPT = """Synthesize research knowledge about this variant using ONLY the evidence below.

Gene: {gene}
Variant: {variant}
Tumor Type: {tumor_type}

## BIOLOGICAL CONTEXT
{biological_context}

## DATABASE EVIDENCE
{evidence_summary}

## LITERATURE FINDINGS
{literature_summary}

## EVIDENCE GAPS
{evidence_gaps}

## EVIDENCE ASSESSMENT
Overall quality: {overall_quality}
Well-characterized aspects: {well_characterized_text}
Known gaps: {known_gaps_text}
Conflicting evidence: {conflicting_evidence_text}

Task:
1. Extract and integrate the key functional, biological, and therapeutic signals.
2. Explicitly note where evidence is strong vs limited vs conflicting, in a way that is consistent with the evidence assessment above.
3. Propose ONE clear, hypothesis-generating research direction that ties together at least two distinct evidence elements, and clearly label this final hypothesis as exploratory.

Additional guidelines:
- Explicitly label extrapolations and hypotheses.
- Prefer hypotheses that are testable and tied to at least two evidence elements.
- Maintain a strict JSON output schema as shown below.

CRITICAL: For biological_context, you MUST include the source citation INLINE with each statistic. Look for "Source: [cBioPortal: ...](...)" in the BIOLOGICAL CONTEXT section above and copy that markdown link directly into your text next to the numbers.

Respond ONLY with valid JSON:
{{
  "functional_summary": "2–3 sentences on what this variant does to protein function and downstream signaling, using only provided evidence.",
  "biological_context": "2–3 sentences with INLINE source citations. Example: 'BRAF mutations occur in 52% of cases ([cBioPortal: mel_tcga](https://...)), with V600E in 35%.' Copy the markdown link exactly from the Source line in the evidence.",
  "therapeutic_landscape": {{
    "fda_approved": ["Drug names with FDA approval in this variant/tumor context, or empty if none are mentioned in the evidence blocks."],
    "clinical_evidence": ["Drugs with Phase 2/3 (or similar) data described in the evidence blocks."],
    "preclinical": ["Agents with preclinical sensitivity/resistance data described in the evidence blocks."],
    "resistance_mechanisms": ["Resistance mechanisms explicitly mentioned in the evidence blocks (e.g., secondary mutations, bypass pathways, histologic transformation)."]
  }},
  "evidence_assessment": {{
    "overall_quality": "{overall_quality}",
    "well_characterized": {well_characterized_json},
    "knowledge_gaps": {known_gaps_json},
    "conflicting_evidence": {conflicting_evidence_json}
  }},
  "research_implications": "2–3 sentences describing what researchers should investigate next, including ONE concrete, testable hypothesis or study design that explicitly references at least two evidence elements (for example: co-mutation pattern + drug response, or functional score + pathway membership). Clearly mark which parts are hypothesis.",
  "key_references": ["PMIDs, trial IDs, or database/source names appearing in the evidence text."],
  "evidence_tags": [
    "Short labels describing which parts of your synthesis are: 'direct clinical data', 'preclinical only', 'nearby-variant extrapolation', 'pathway-level hypothesis', etc."
  ]
}}
Use empty arrays or null where appropriate.
"""


def create_research_prompt(
    gene: str,
    variant: str,
    tumor_type: str | None,
    biological_context: str,
    evidence_summary: str,
    evidence_assessment: dict,
    literature_summary: str = "",
) -> list[dict]:
    """
    Create prompt for research-focused variant synthesis.

    Args:
        gene: Gene symbol
        variant: Variant notation
        tumor_type: Tumor type context
        biological_context: cBioPortal prevalence/co-mutation data
        evidence_summary: Aggregated database evidence
        evidence_assessment: Dict with keys overall_quality, well_characterized, knowledge_gaps, conflicting_evidence
        literature_summary: PubMed/Semantic Scholar findings (optional)

    Returns:
        Messages list for LLM API call.
    """
    tumor_display = tumor_type or "Pan-cancer"

    overall_quality = evidence_assessment.get("overall_quality", "minimal")
    well_char = evidence_assessment.get("well_characterized", []) or []
    gaps = evidence_assessment.get("knowledge_gaps", []) or []
    conflicts = evidence_assessment.get("conflicting_evidence", []) or []

    # Build the JSON-safe strings for the template
    user_content = RESEARCH_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        biological_context=biological_context or "No cBioPortal or co-mutation data available.",
        evidence_summary=(evidence_summary or "No database evidence available.").strip()[:4000],
        literature_summary=literature_summary or "No literature search performed.",
        evidence_gaps="\n".join(gaps) or "No explicit gaps analysis available.",
        overall_quality=overall_quality,
        well_characterized_text="; ".join(well_char) or "None specified.",
        known_gaps_text="; ".join(gaps) or "None specified.",
        conflicting_evidence_text="; ".join(conflicts) or "None specified.",
        well_characterized_json=json.dumps(well_char),
        known_gaps_json=json.dumps(gaps),
        conflicting_evidence_json=json.dumps(conflicts),
    )

    return [
        {"role": "system", "content": RESEARCH_SYSTEM_PROMPT},
        {"role": "user", "content": user_content},
    ]


# =============================================================================
# Legacy clinical-annotation prompts (kept for backwards compatibility)
# =============================================================================

ANNOTATION_SYSTEM_PROMPT = """You are an expert molecular tumor board pathologist writing a concise clinical summary for a cancer variant.

Your task is to synthesize the provided evidence into a clear, actionable summary that:
1. States the clinical significance of this variant based on the evidence provided.
2. Notes any therapeutic implications (approved therapies, contraindications, or relevant trials).
3. Identifies key database annotations and their clinical relevance.
4. Is suitable for inclusion in a clinical-style report.

Use ONLY the supplied evidence. Keep the tone focused and evidence-based. Do not speculate beyond the data.
"""

ANNOTATION_USER_PROMPT = """Write a clinical annotation summary for this variant:

Gene: {gene}
Variant: {variant}
Tumor Type: {tumor_type}
{therapy_note_section}
Evidence Summary:
{evidence_summary}

Respond ONLY with valid JSON:
{{
  "summary": "3–5 sentence clinical-style summary covering significance and therapeutic implications, using only the provided evidence.",
  "rationale": "Brief explanation of how the evidence supports this summary.",
  "recommended_therapies": [
    {{
      "drug_name": "Drug name",
      "evidence_level": "FDA-approved | Clinical trial | Preclinical",
      "approval_status": "Approved | Investigational | Off-label",
      "clinical_context": "e.g., first-line, second-line, post-resistance, biomarker-selected"
    }}
  ],
  "references": ["Key PMIDs, guidelines, or database sources mentioned in the evidence block."]
}}
"""

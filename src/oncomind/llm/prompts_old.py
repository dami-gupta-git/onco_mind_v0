"""Prompts for research-focused variant annotation and evidence synthesis."""

RESEARCH_SYSTEM_PROMPT = """You are an expert cancer genomics researcher synthesizing evidence about a somatic variant.

Your task is to provide a comprehensive RESEARCH-ORIENTED summary that:

1. Explains the FUNCTIONAL IMPACT of this variant (what does it do to the protein?)
2. Describes the BIOLOGICAL CONTEXT (prevalence, co-mutations, pathway effects)
3. Summarizes the EVIDENCE LANDSCAPE (what's known, what's uncertain, what's controversial)
4. Notes ALL THERAPEUTIC EVIDENCE (FDA-approved, clinical trials, AND preclinical)
5. Identifies KNOWLEDGE GAPS (what's missing, where is more research needed)

This is for researchers, not clinical decision-making. Be comprehensive and nuanced.
Distinguish between established findings and emerging hypotheses.
Note conflicting evidence rather than hiding it.
Include preclinical and mechanistic data, not just approved therapies."""


RESEARCH_USER_PROMPT = """Synthesize research knowledge about this variant:

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

Respond with JSON:
{{
  "functional_summary": "2-3 sentences on what this variant does to protein function",
  "biological_context": "2-3 sentences on prevalence, co-mutations, pathway effects",
  "therapeutic_landscape": {{
    "fda_approved": ["Drug names with FDA approval"],
    "clinical_evidence": ["Drugs with Phase 2/3 data"],
    "preclinical": ["Drugs with preclinical sensitivity/resistance data"],
    "resistance_mechanisms": ["Known resistance mechanisms"]
  }},
  "evidence_assessment": {{
    "overall_quality": "comprehensive | moderate | limited | minimal",
    "well_characterized": ["Aspects with strong evidence"],
    "knowledge_gaps": ["What's unknown or understudied"],
    "conflicting_evidence": ["Where sources disagree"]
  }},
  "research_implications": "What should researchers investigate next?",
  "key_references": ["PMIDs or database sources"]
}}
"""


def create_research_prompt(
    gene: str,
    variant: str,
    tumor_type: str | None,
    biological_context: str,
    evidence_summary: str,
    evidence_gaps: str,
    literature_summary: str = "",
) -> list[dict]:
    """Create prompt for research-focused variant synthesis.

    Args:
        gene: Gene symbol
        variant: Variant notation
        tumor_type: Tumor type context
        biological_context: cBioPortal prevalence/co-mutation data
        evidence_summary: Database evidence
        evidence_gaps: Evidence gaps summary
        literature_summary: PubMed literature findings (for full mode)

    Returns:
        Messages list for LLM API call
    """
    tumor_display = tumor_type or "Pan-cancer"

    user_content = RESEARCH_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        biological_context=biological_context or "No cBioPortal data available",
        evidence_summary=evidence_summary.strip()[:4000],
        literature_summary=literature_summary or "No literature search performed",
        evidence_gaps=evidence_gaps or "No gaps analysis available",
    )

    return [
        {"role": "system", "content": RESEARCH_SYSTEM_PROMPT},
        {"role": "user", "content": user_content}
    ]


# === KEEP OLD PROMPTS for backwards compatibility ===

ANNOTATION_SYSTEM_PROMPT = """You are an expert molecular tumor board pathologist writing a concise clinical summary for a cancer variant.

Your task is to synthesize evidence from multiple databases into a clear, actionable summary that:
1. States the clinical significance of this variant based on the evidence provided
2. Notes any therapeutic implications (approved therapies, contraindications, or trials)
3. Identifies key database annotations and their clinical relevance
4. Is suitable for inclusion in a clinical report

Keep it focused and evidence-based. Prioritize actionable information."""

ANNOTATION_USER_PROMPT = """Write a clinical annotation summary for this variant:

Gene: {gene}
Variant: {variant}
Tumor Type: {tumor_type}
{therapy_note_section}
Evidence Summary:
{evidence_summary}

Respond with JSON:
{{
  "summary": "3-5 sentence clinical summary covering significance and therapeutic implications",
  "rationale": "Brief explanation of the evidence supporting this summary",
  "recommended_therapies": [
    {{
      "drug_name": "Drug name",
      "evidence_level": "FDA-approved/Clinical trial/Preclinical",
      "approval_status": "Approved/Investigational/Off-label",
      "clinical_context": "First-line/Second-line/Resistance etc."
    }}
  ],
  "references": ["Key reference PMIDs or sources"]
}}
"""


# def create_annotation_prompt(
#     gene: str,
#     variant: str,
#     tumor_type: str | None,
#     evidence_summary: str,
#     therapy_note: str | None = None,
# ) -> list[dict]:
#     """Create prompt for clinical annotation (backwards compatible)."""
#     tumor_display = tumor_type if tumor_type else "Unspecified"
#
#     therapy_note_section = ""
#     if therapy_note:
#         therapy_note_section = f"\nTherapeutic Context: {therapy_note}\n"
#
#     user_content = ANNOTATION_USER_PROMPT.format(
#         gene=gene,
#         variant=variant,
#         tumor_type=tumor_display,
#         therapy_note_section=therapy_note_section,
#         evidence_summary=evidence_summary.strip()[:3000],
#     )
#
#     return [
#         {"role": "system", "content": ANNOTATION_SYSTEM_PROMPT},
#         {"role": "user", "content": user_content}
#     ]

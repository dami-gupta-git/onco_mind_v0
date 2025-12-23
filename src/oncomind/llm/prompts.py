# oncomind/prompts.py
"""
Prompts for variant annotation and evidence synthesis.
"""

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


def create_annotation_prompt(
    gene: str,
    variant: str,
    tumor_type: str | None,
    evidence_summary: str,
    therapy_note: str | None = None,
) -> list[dict]:
    """
    Create a prompt for the LLM to synthesize evidence into a variant annotation.

    Args:
        gene: Gene symbol
        variant: Variant notation
        tumor_type: Patient's tumor type
        evidence_summary: Formatted evidence from databases
        therapy_note: Optional note about therapeutic implications (e.g., mutation class info)

    Returns:
        Messages list for LLM API call
    """
    tumor_display = tumor_type if tumor_type else "Unspecified"

    # Format therapy note section if provided
    therapy_note_section = ""
    if therapy_note:
        therapy_note_section = f"\nTherapeutic Context: {therapy_note}\n"

    user_content = ANNOTATION_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        therapy_note_section=therapy_note_section,
        evidence_summary=evidence_summary.strip()[:3000],  # Limit context size
    )

    return [
        {"role": "system", "content": ANNOTATION_SYSTEM_PROMPT},
        {"role": "user", "content": user_content}
    ]

"""Prompts for research-focused variant annotation and evidence synthesis.

Two-stage LLM pipeline:
1. Synthesis (SYNTHESIS_*) - integrates evidence into functional/biological/therapeutic summary
2. Hypothesis (HYPOTHESIS_*) - generates testable research questions from synthesis + gaps
"""

import json


# =============================================================================
# STAGE 1: EVIDENCE SYNTHESIS
# =============================================================================

SYNTHESIS_SYSTEM_PROMPT = """You are an expert cancer genomics researcher synthesizing evidence about a somatic variant.

Generate a RESEARCH-ORIENTED synthesis (not clinical recommendations) covering:
1. FUNCTIONAL IMPACT - how the variant alters protein activity
2. BIOLOGICAL CONTEXT - prevalence, co-mutations, pathway effects
3. THERAPEUTIC LANDSCAPE - FDA-approved, clinical, preclinical, resistance
4. EVIDENCE QUALITY - what's established vs sparse vs conflicting

=== CALIBRATION RULES (STRICT) ===

When overall_quality is "limited" or "minimal":
- functional_summary MUST be GENERIC (gene function only, not variant-specific effects)
- Do NOT assign oncogene/tumor-suppressor roles unless GENE ROLE section states it
- Do NOT predict drug response unless FDA/CIViC/VICC evidence exists for THIS variant
- Keep sections brief (2-3 sentences), focus on what's UNKNOWN

When has_tumor_specific_cbioportal_data is FALSE:
- State: "Pan-cancer data; no {tumor_type}-specific prevalence available"
- Do NOT extrapolate pan-cancer co-mutations to specific tumor context

=== MATCH SPECIFICITY (CRITICAL) ===

Always indicate evidence match level:
- VARIANT-LEVEL: "osimertinib (this variant, {tumor_type})"
- CODON-LEVEL: "MEK inhibitors (studied for Q209P, not Q209L)" - same codon, different amino acid
- GENE-LEVEL: "erlotinib (gene-level)"
- OTHER CANCER: "sotorasib (approved for NSCLC, not {tumor_type})"

CODON-LEVEL EVIDENCE WARNING:
When evidence comes from OTHER variants at the same codon (e.g., Q209P data applied to Q209L):
- State explicitly: "Evidence from [other variant] at same codon; [queried variant]-specific data limited"
- Note that different amino acid substitutions can have distinct signaling and drug-response profiles
- Flag this as a knowledge gap requiring variant-specific validation

If THERAPEUTIC SIGNALS says "FDA-approved for OTHER cancers (NOT {tumor_type})", report it as approved for that other cancer, NOT {tumor_type}.

=== CONFLICTING EVIDENCE ===

Distinguish expected biology from true conflicts:
- EXPECTED: T790M resistant to erlotinib but sensitive to osimertinib (sequential therapy)
- TRUE CONFLICT: Same drug, same setting, contradictory outcomes → flag as gap

=== HARD CONSTRAINTS ===

- Use ONLY evidence from the user message. Do NOT invent facts.
- ALWAYS cite sources for statistics with markdown links as provided
- Include resistance signals in therapeutic_landscape.resistance_mechanisms
- Include sensitivity signals in therapeutic_landscape.clinical_evidence or preclinical
"""

SYNTHESIS_USER_PROMPT = """Synthesize evidence for this variant. Use ONLY the data below.

Gene: {gene}
Variant: {variant}
Tumor Type: {tumor_type}

## DATA FLAGS
has_tumor_specific_cbioportal_data: {has_tumor_specific_cbioportal}
has_civic_assertions: {has_civic_assertions}
has_fda_approvals: {has_fda_approvals}
has_vicc_evidence: {has_vicc_evidence}

## BIOLOGICAL CONTEXT
{biological_context}

## THERAPEUTIC SIGNALS
Sensitivity: {sensitivity_summary}
Resistance: {resistance_summary}

## DATABASE EVIDENCE
{evidence_summary}

## LITERATURE
{literature_summary}

## EVIDENCE ASSESSMENT
Overall quality: {overall_quality}
Well-characterized: {well_characterized_text}
Gaps: {known_gaps_text}
Conflicts: {conflicting_evidence_text}

Respond with valid JSON only:
{{
  "functional_summary": "Gene function. If quality is limited/minimal: generic only. If moderate/comprehensive: variant-specific with citations.",
  "biological_context": "2-3 sentences. Start with 'As per [STUDY] - '. If no tumor-specific data, state 'Pan-cancer data shown.'",
  "therapeutic_landscape": {{
    "fda_approved": ["drug (match level, cancer type)"],
    "clinical_evidence": ["drug - source"],
    "preclinical": ["drug - source"],
    "resistance_mechanisms": ["drug - mechanism"]
  }},
  "evidence_assessment": {{
    "overall_quality": "{overall_quality}",
    "well_characterized": {well_characterized_json},
    "knowledge_gaps": {known_gaps_json},
    "conflicting_evidence": {conflicting_evidence_json}
  }},
  "key_references": ["PMIDs, databases from evidence"],
  "evidence_tags": ["direct clinical data | preclinical only | pan-cancer extrapolation | limited evidence"]
}}
"""


# =============================================================================
# STAGE 2: HYPOTHESIS GENERATION
# =============================================================================

HYPOTHESIS_SYSTEM_PROMPT = """You are a cancer genomics researcher generating testable research hypotheses.

You will receive:
1. A synthesis of variant evidence (from stage 1)
2. Evidence gaps that need investigation

Generate 2-3 SPECIFIC, TESTABLE research hypotheses that:
- Address the identified knowledge gaps
- Build on existing evidence (not speculation)
- Are experimentally tractable
- Focus on research questions (NOT clinical recommendations)

=== HYPOTHESIS REQUIREMENTS ===

Each hypothesis MUST:
1. Start with an EVIDENCE BASIS TAG:
   - [Direct Clinical Data] - builds on FDA/CIViC/Phase 2-3 trials for THIS variant
   - [Preclinical Data] - builds on DepMap/cell line data
   - [Pan-Cancer Extrapolation] - extrapolates from other tumor types
   - [Nearby-Variant Inference] - extrapolates from other variants in same gene
   - [Pathway-Level Inference] - infers from pathway biology

2. Be SPECIFIC and TESTABLE (not vague)
3. Connect a GAP to existing EVIDENCE
4. Suggest concrete experimental approach

=== EXAMPLES ===

Good: "[Preclinical Data] Given DepMap shows BRAF V600E dependency in melanoma but no drug sensitivity data exists for this variant, systematic testing of BRAF inhibitors in isogenic models could establish therapeutic vulnerability."

Good: "[Pan-Cancer Extrapolation] EGFR L858R shows osimertinib sensitivity in NSCLC; testing cross-histology response in breast cancer models would determine tissue-specific effects."

Bad: "More research is needed" (vague)
Bad: "Patients should receive this drug" (clinical recommendation)
Bad: "Test JAK inhibitors" (no evidence basis tag)
"""

HYPOTHESIS_USER_PROMPT = """Generate research hypotheses based on this synthesis and gaps.

Gene: {gene}
Variant: {variant}
Tumor Type: {tumor_type}

## SYNTHESIS (from stage 1)
Functional: {functional_summary}
Biological: {biological_context}
Therapeutic: {therapeutic_landscape}
Evidence quality: {overall_quality}

## KNOWLEDGE GAPS TO ADDRESS
{knowledge_gaps}

## AVAILABLE EVIDENCE TO BUILD ON
Well-characterized: {well_characterized}
Therapeutic signals: {therapeutic_signals}

Generate 2-3 hypotheses. Respond with valid JSON only:
{{
  "research_hypotheses": [
    "[Evidence Tag] Specific testable hypothesis connecting gap to evidence...",
    "[Evidence Tag] Another hypothesis..."
  ],
  "research_implications": "2-3 sentence summary of key research directions."
}}
"""


# =============================================================================
# PROMPT BUILDERS
# =============================================================================

def create_synthesis_prompt(
    gene: str,
    variant: str,
    tumor_type: str | None,
    biological_context: str,
    evidence_summary: str,
    evidence_assessment: dict,
    literature_summary: str = "",
    data_availability: dict | None = None,
    resistance_summary: str = "",
    sensitivity_summary: str = "",
    match_level_summary: dict | None = None,
) -> list[dict]:
    """Create prompt for stage 1: evidence synthesis."""
    tumor_display = tumor_type or "Pan-cancer"

    overall_quality = evidence_assessment.get("overall_quality", "minimal")
    well_char = evidence_assessment.get("well_characterized", []) or []
    gaps = evidence_assessment.get("knowledge_gaps", []) or []
    conflicts = evidence_assessment.get("conflicting_evidence", []) or []

    if data_availability is None:
        data_availability = {}

    _ = match_level_summary  # Reserved for future use

    user_content = SYNTHESIS_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        has_tumor_specific_cbioportal=str(data_availability.get("has_tumor_specific_cbioportal", False)).upper(),
        has_civic_assertions=str(data_availability.get("has_civic_assertions", False)).upper(),
        has_fda_approvals=str(data_availability.get("has_fda_approvals", False)).upper(),
        has_vicc_evidence=str(data_availability.get("has_vicc_evidence", False)).upper(),
        biological_context=biological_context or "No cBioPortal data available.",
        resistance_summary=resistance_summary or "No resistance signals.",
        sensitivity_summary=sensitivity_summary or "No sensitivity signals.",
        evidence_summary=(evidence_summary or "No database evidence.").strip()[:4000],
        literature_summary=literature_summary or "No literature search performed.",
        overall_quality=overall_quality,
        well_characterized_text="; ".join(well_char) or "None.",
        known_gaps_text="; ".join(gaps) or "None.",
        conflicting_evidence_text="; ".join(conflicts) or "None.",
        well_characterized_json=json.dumps(well_char),
        known_gaps_json=json.dumps(gaps),
        conflicting_evidence_json=json.dumps(conflicts),
    )

    return [
        {"role": "system", "content": SYNTHESIS_SYSTEM_PROMPT},
        {"role": "user", "content": user_content},
    ]


def create_hypothesis_prompt(
    gene: str,
    variant: str,
    tumor_type: str | None,
    synthesis_result: dict,
    evidence_assessment: dict,
    therapeutic_signals: str = "",
) -> list[dict]:
    """Create prompt for stage 2: hypothesis generation.

    Args:
        gene: Gene symbol
        variant: Variant notation
        tumor_type: Tumor type context
        synthesis_result: Output from stage 1 synthesis
        evidence_assessment: Dict with knowledge_gaps, well_characterized
        therapeutic_signals: Summary of sensitivity/resistance signals
    """
    tumor_display = tumor_type or "Pan-cancer"

    gaps = evidence_assessment.get("knowledge_gaps", []) or []
    well_char = evidence_assessment.get("well_characterized", []) or []

    # Format therapeutic landscape for context
    therapeutic = synthesis_result.get("therapeutic_landscape", {})
    therapeutic_str = json.dumps(therapeutic, indent=2) if therapeutic else "None"

    user_content = HYPOTHESIS_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        functional_summary=synthesis_result.get("functional_summary", "Not available"),
        biological_context=synthesis_result.get("biological_context", "Not available"),
        therapeutic_landscape=therapeutic_str,
        overall_quality=synthesis_result.get("evidence_assessment", {}).get("overall_quality", "unknown"),
        knowledge_gaps="\n".join(f"- {g}" for g in gaps) or "None identified.",
        well_characterized="; ".join(well_char) or "None.",
        therapeutic_signals=therapeutic_signals or "No therapeutic signals.",
    )

    return [
        {"role": "system", "content": HYPOTHESIS_SYSTEM_PROMPT},
        {"role": "user", "content": user_content},
    ]


# =============================================================================
# BACKWARDS COMPATIBILITY
# =============================================================================

def create_research_prompt(
    gene: str,
    variant: str,
    tumor_type: str | None,
    biological_context: str,
    evidence_summary: str,
    evidence_assessment: dict,
    literature_summary: str = "",
    data_availability: dict | None = None,
    resistance_summary: str = "",
    sensitivity_summary: str = "",
    match_level_summary: dict | None = None,
) -> list[dict]:
    """Create prompt for research-focused variant synthesis.

    DEPRECATED: Use create_synthesis_prompt + create_hypothesis_prompt instead.
    This function is kept for backwards compatibility.
    """
    return create_synthesis_prompt(
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

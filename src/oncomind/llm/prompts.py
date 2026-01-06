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

Always indicate both LOCUS match and TUMOR match level for therapeutic evidence:

LOCUS MATCH (variant specificity):
- VARIANT-LEVEL: Evidence directly studied for THIS specific variant
- CODON-LEVEL: Evidence from OTHER variants at same codon (e.g., Q209P data applied to Q209L)
- GENE-LEVEL: Evidence from gene-level studies (any mutation in the gene)

TUMOR MATCH (cancer specificity):
- CANCER-SPECIFIC: Evidence from the patient's tumor type ({tumor_type})
- PAN-CANCER: Tumor-agnostic evidence (e.g., MSI-H, Solid Tumor)
- OTHER CANCER: Evidence from a different specific cancer type

Format therapeutic entries as: "drug (locus-level, approved for VARIANT, tumor-context)"
IMPORTANT: When evidence says "approved for [VARIANT]", ALWAYS include the variant info.
NEVER say "approved for any [GENE] mutation" - most targeted therapies are approved for SPECIFIC variants only.
Examples:
- "osimertinib (variant-level, approved for T790M, NSCLC)" - best evidence, variant-specific
- "sotorasib (gene-level, approved for G12C, NSCLC)" - G12C-specific drug used for other KRAS variant
- "adagrasib (gene-level, approved for G12C, not pancreatic)" - different tumor AND different variant
- "MEK inhibitors (codon-level from Q209P, melanoma)" - same codon, different AA
- "gefitinib (codon-level, approved for exon 19 del/L858R, NSCLC)" - approved for specific sensitizing mutations, NOT for T790M or other variants
- "erlotinib (codon-level, approved for exon 19 del/L858R, NSCLC)" - same as gefitinib, specific variants only

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
- NEVER describe HOW a variant works mechanistically. Only state THAT it is oncogenic/pathogenic if evidence says so. Delete any phrases about membrane localization, pathway activation, signaling, or protein function mechanisms.
- NEVER say "approved for any [GENE] mutation" - most targeted therapies are approved for SPECIFIC variants only. Read the indication text in the evidence and specify the exact approved variants.

=== CRITICAL: NO HALLUCINATION ===

- NEVER mention drugs, clinical trials, or treatments NOT explicitly listed in the DATABASE EVIDENCE section
- If no FDA approvals are listed → say "No FDA-approved therapies"
- If no clinical trials are listed → say "No clinical trials found"
- Do NOT invent drug names, trial identifiers (NCT numbers), or phase information
- Do NOT use your training knowledge about drugs or trials - ONLY use what's provided below
- If evidence is sparse, say "Limited evidence available" - do NOT fill gaps with training data
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

## LOCUS MATCH SUMMARY
{locus_match_text}

## TUMOR MATCH SUMMARY
{tumor_match_text}

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
    "fda_approved": ["drug (locus-level, approved for VARIANT if gene/codon-level, tumor)"],
    "clinical_evidence": ["drug (locus-level, tumor) - source"],
    "preclinical": ["drug (locus-level) - source"],
    "resistance_mechanisms": ["drug - mechanism (locus-level)"]
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
    locus_match_summary: dict | None = None,
    tumor_match_summary: dict | None = None,
) -> list[dict]:
    """Create prompt for stage 1: evidence synthesis."""
    tumor_display = tumor_type or "Pan-cancer"

    overall_quality = evidence_assessment.get("overall_quality", "minimal")
    well_char = evidence_assessment.get("well_characterized", []) or []
    gaps = evidence_assessment.get("knowledge_gaps", []) or []
    conflicts = evidence_assessment.get("conflicting_evidence", []) or []

    if data_availability is None:
        data_availability = {}

    # Build locus match text for LLM context
    if locus_match_summary:
        locus_match_text = locus_match_summary.get("summary_text", "No locus match data available.")
        # Add detail about what's variant-specific vs gene-level
        if locus_match_summary.get("is_all_gene_level"):
            locus_match_text += " CAUTION: No variant-specific evidence found - use gene-level inferences carefully."
        elif not locus_match_summary.get("has_variant_specific"):
            locus_match_text += " WARNING: Limited variant-specific data."
    else:
        locus_match_text = "No locus match data available."

    # Build tumor match text for LLM context
    if tumor_match_summary:
        tumor_match_text = tumor_match_summary.get("summary_text", "No tumor match data available.")
        # Add warnings about evidence from other tumor types
        if tumor_match_summary.get("is_all_other"):
            other_tumors = ", ".join(tumor_match_summary.get("other_tumor_types", [])[:3])
            tumor_match_text += f" CAUTION: All evidence is from other tumor types ({other_tumors}) - extrapolation required."
        elif not tumor_match_summary.get("has_tumor_specific"):
            tumor_match_text += " WARNING: Limited tumor-specific evidence."
    else:
        tumor_match_text = "No tumor match data available."

    user_content = SYNTHESIS_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        has_tumor_specific_cbioportal=str(data_availability.get("has_tumor_specific_cbioportal", False)).upper(),
        has_civic_assertions=str(data_availability.get("has_civic_assertions", False)).upper(),
        has_fda_approvals=str(data_availability.get("has_fda_approvals", False)).upper(),
        has_vicc_evidence=str(data_availability.get("has_vicc_evidence", False)).upper(),
        locus_match_text=locus_match_text,
        tumor_match_text=tumor_match_text,
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
    locus_match_summary: dict | None = None,
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
        locus_match_summary=locus_match_summary,
    )

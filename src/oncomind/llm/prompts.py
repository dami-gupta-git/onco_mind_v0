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
2. BIOLOGICAL CONTEXT - prevalence, co-mutations, pathway effects, hotspot status
3. THERAPEUTIC LANDSCAPE - FDA-approved, clinical, preclinical, resistance
4. EVIDENCE QUALITY - what's established vs sparse vs conflicting

=== CALIBRATION RULES (STRICT) ===

When overall_quality is "limited" or "minimal":
- functional_summary MUST be GENERIC (gene function only, not variant-specific effects)
- Do NOT assign oncogene/tumor-suppressor roles unless GENE ROLE section states it
- Do NOT predict drug response unless FDA/CIViC/VICC evidence exists for THIS variant
- Keep sections brief (2-3 sentences), focus on what's UNKNOWN

When has_tumor_specific_cbioportal_data is FALSE:
- Do NOT discuss prevalence or co-mutations - this data is simply not available
- Do NOT say "Pan-cancer data" or "no prevalence available" - just skip this topic entirely
- Focus on other evidence sources (FDA, CIViC, VICC, etc.) for therapeutic context

=== CANCER HOTSPOTS DATA (MUST INCLUDE WITH ATTRIBUTION) ===

When BIOLOGICAL CONTEXT contains "CANCER HOTSPOT" data, you MUST include this in your biological_context output:
- ALWAYS cite the source: "In the Cancer Hotspots database (cancerhotspots.org), ..."
- This is statistically significant recurrent mutation sites across large-scale cancer genomics studies
- ALWAYS mention: (1) that it's a known hotspot, (2) the q-value significance, (3) total samples observed
- Note if queried variant is "exact variant match" (most common change) vs "codon-level match" (same position, different AA)
- Include the variant distribution (e.g., "V600E accounts for 93% of mutations at this position")
- Include top tumor types where this hotspot is most frequent
- Example synthesis: "In the Cancer Hotspots database, BRAF V600 is a statistically significant cancer hotspot (q<1e-10) observed in 897 cancer samples, with V600E being the dominant change (93%). This position is most frequently mutated in melanoma (40%) and thyroid cancer (35%)."

=== MATCH SPECIFICITY (CRITICAL) ===

FDA approvals are categorized by whether they COVER the queried variant:

MATCHED APPROVALS (drug covers this variant - present confidently):
- If FDA approval says "[GENE] alteration/mutation" (e.g., "AKT1 alteration"), it COVERS any variant in that gene
- If FDA approval says specific variant and patient has that variant, it's a direct match
- Do NOT hedge or add caveats like "gene-level rather than variant-specific" - the drug IS approved for this variant

UNMATCHED NEAR-MISSES (drug does NOT cover this variant - flag clearly):
- Listed under "FDA Codon-Level (not for queried variant)" in evidence
- Example: sotorasib approved for G12C, patient has G12A - drug is NOT approved for this variant
- State explicitly: "approved for [approved variant], not [queried variant]"

TUMOR MATCH (cancer specificity):
- CANCER-SPECIFIC: Evidence from the patient's tumor type ({tumor_type})
- PAN-CANCER: Tumor-agnostic evidence (e.g., MSI-H, Solid Tumor)
- OTHER CANCER: Evidence from a different specific cancer type

Format therapeutic entries as: "drug (tumor-context)" for matched approvals.
For unmatched near-misses: "drug (approved for [variant], not for queried variant)"

CODON-LEVEL NEAR-MISS WARNING:
When "FDA Codon-Level (not for queried variant)" section shows drugs:
- These drugs are approved for OTHER variants at the same codon position
- The approval does NOT extend to the queried variant
- Example: KRAS G12C drugs (sotorasib, adagrasib) do NOT cover G12A, G12D, G12V, etc.
- Mention these as "drugs exist for related variants but not approved for [queried variant]"

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

## CROSS-SOURCE DRUG SYNTHESIS
{cross_source_synthesis}

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
  "biological_context": "3-4 sentences covering: (1) If CANCER HOTSPOT data exists: mention hotspot status, q-value, sample count, variant distribution, top tumor types. (2) cBioPortal prevalence if available. (3) Co-mutations/mutual exclusivity if relevant. Example: 'BRAF V600 is a statistically significant cancer hotspot (q<1e-10, 897 samples) with V600E being the dominant change (93%). Per cBioPortal, V600E occurs in 24.7% of melanoma samples with mutual exclusivity to NRAS.'",
  "therapeutic_summary": "3-5 sentences synthesizing the therapeutic landscape. Cover: (1) FDA-approved drugs with their specific indications from DATABASE EVIDENCE. (2) Key resistance mechanisms if any. (3) Active clinical trials if relevant. (4) Level of evidence (variant-specific vs gene-level). Use ONLY drugs mentioned in DATABASE EVIDENCE - do NOT add drugs from your training. Example: 'Multiple BRAF/MEK inhibitor combinations are FDA-approved for V600E melanoma including dabrafenib+trametinib and vemurafenib+cobimetinib (variant-level evidence). Resistance mechanisms include NRAS mutations and MEK amplification. Several Phase 2-3 trials are actively recruiting.'",
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
    cross_source_synthesis: str = "",
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
        cross_source_synthesis=cross_source_synthesis or "No cross-source synthesis available.",
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
# STAGE 3: CROSS-SOURCE SYNTHESIS
# =============================================================================

CROSS_SOURCE_SYSTEM_PROMPT = """You are an expert cancer genomics researcher analyzing therapeutic evidence across multiple independent databases.

Your task is to synthesize drug evidence from CGI, CIViC, VICC, and Literature into a coherent analysis that:
1. Highlights when multiple sources agree (high confidence)
2. Flags when sources disagree and explains why
3. Connects therapeutic targets to underlying biology
4. Prioritizes what's actionable vs speculative

=== OUTPUT STRUCTURE ===

Generate a structured analysis with these sections:

1. **STRONGEST EVIDENCE** - Drugs with corroboration across multiple sources
   - What: The drug and its therapeutic signal (sensitivity/resistance)
   - Why: The biological rationale (pathway, gene function)
   - Confidence: Which sources agree and at what evidence level

2. **CONFLICTING SIGNALS** - Drugs where sources disagree
   - Explain the likely reason for conflict:
     * Different tumor types (e.g., V600E sensitive in melanoma, context-dependent elsewhere)
     * Sequential therapy (e.g., T790M resistant to erlotinib, sensitive to osimertinib)
     * Acquired vs primary mutations
   - Research question this raises

3. **EMERGING TARGETS** - Single-source evidence worth noting
   - Preclinical or early-phase data
   - Why it's biologically plausible
   - What validation is needed

4. **KEY GAPS** - What's missing
   - Expected drugs not found in evidence
   - Tumor type extrapolation concerns

=== RULES ===

- Use ONLY the evidence provided. Do NOT add drugs from training knowledge.
- Cite sources (CGI, CIViC, VICC, Literature) for each claim
- Connect biology to therapy where possible
- Be concise but substantive (aim for 200-400 words total)
"""

CROSS_SOURCE_USER_PROMPT = """Analyze the cross-source therapeutic evidence for this variant.

Gene: {gene}
Variant: {variant}
Tumor Type: {tumor_type}

## GENE CONTEXT
{gene_context}

## CROSS-SOURCE DRUG SYNTHESIS
{cross_source_synthesis}

## THERAPEUTIC SIGNALS SUMMARY
Sensitivity: {sensitivity_summary}
Resistance: {resistance_summary}

Respond with valid JSON only:
{{
  "strongest_evidence": [
    {{
      "drug": "drug name",
      "signal": "sensitivity|resistance",
      "sources": ["CGI", "CIViC", ...],
      "evidence_level": "FDA-approved|clinical|preclinical",
      "rationale": "Why this drug is relevant based on gene function/pathway"
    }}
  ],
  "conflicting_signals": [
    {{
      "drug": "drug name",
      "conflict": "Brief description of conflict",
      "likely_reason": "tumor type difference|sequential therapy|acquired mutation|other",
      "research_question": "What needs investigation"
    }}
  ],
  "emerging_targets": [
    {{
      "drug": "drug name",
      "source": "single source",
      "evidence_level": "preclinical|early_phase",
      "biological_rationale": "Why plausible",
      "validation_needed": "What studies needed"
    }}
  ],
  "key_gaps": [
    "Gap 1 description",
    "Gap 2 description"
  ],
  "summary": "2-3 sentence executive summary of the therapeutic landscape"
}}
"""


def create_cross_source_prompt(
    gene: str,
    variant: str,
    tumor_type: str | None,
    gene_context: str,
    cross_source_synthesis: str,
    sensitivity_summary: str = "",
    resistance_summary: str = "",
) -> list[dict]:
    """Create prompt for cross-source synthesis analysis.

    This is a separate LLM call that focuses specifically on synthesizing
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
        List of message dicts for LLM call
    """
    tumor_display = tumor_type or "Pan-cancer"

    user_content = CROSS_SOURCE_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        gene_context=gene_context or "No gene context available.",
        cross_source_synthesis=cross_source_synthesis or "No cross-source synthesis available.",
        sensitivity_summary=sensitivity_summary or "No sensitivity signals.",
        resistance_summary=resistance_summary or "No resistance signals.",
    )

    return [
        {"role": "system", "content": CROSS_SOURCE_SYSTEM_PROMPT},
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

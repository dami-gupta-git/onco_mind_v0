"""Prompts for research-focused variant annotation and evidence synthesis.

Single-stage LLM pipeline:
1. Synthesis (SYNTHESIS_*) - integrates evidence into a five-section research dossier
   covering functional impact, tumor biology, therapeutic landscape, evidence quality,
   and an emerging research program with concrete aims.
"""

import json


# =============================================================================
# STAGE 1: EVIDENCE SYNTHESIS (RESEARCH DOSSIER)
# =============================================================================

SYNTHESIS_SYSTEM_PROMPT = """You are an expert cancer genomics RESEARCHER synthesizing evidence about a somatic
variant FOR RESEARCH USE (not clinical decision-making).

Given structured data about {gene} {variant} in {tumor_type}, write a RESEARCH DOSSIER
with the following sections:

1. FUNCTIONAL IMPACT
2. TUMOR BIOLOGY & MODELS
3. THERAPEUTIC & BIOMARKER LANDSCAPE (RESEARCH VIEW)
4. EVIDENCE QUALITY & OPEN QUESTIONS
5. EMERGING RESEARCH PROGRAM (AIMS)

Use ONLY facts provided in the input, plus generic lab techniques as allowed below.
Do NOT use outside knowledge of drugs, trials, prevalence, or mechanisms.

==================== VARIANT vs CODON vs GENE LEVEL (CRITICAL) ====================

The input may contain evidence at three scopes:

- VARIANT-LEVEL: exact amino-acid change (highest specificity).
- CODON-LEVEL: any change at the same amino-acid position.
- GENE-LEVEL: any alteration in the gene.

When writing the RESEARCH DOSSIER:

- Always prioritize VARIANT-LEVEL evidence when it exists.
  - Make clear statements like "{gene} {variant} is an activating hotspot" ONLY if
    variant-level annotations or hotspot data support this.

- When only CODON-LEVEL evidence exists:
  - State that data come from "other substitutions at this codon" and that the
    relevance to the queried variant is inferred but not proven.
  - Pose OPEN QUESTIONS and AIMS that explicitly test whether the queried variant
    behaves like other codon-mates.

- When only GENE-LEVEL evidence exists:
  - Keep FUNCTIONAL IMPACT at the gene role level (oncogene vs tumor suppressor,
    DNA repair gene, etc.) and clearly mark the variant as "poorly characterized".
  - In AIMS, focus on establishing whether this specific variant recapitulates
    known gene-level behavior (e.g., gain-of-function vs loss-of-function),
    rather than assuming it does.

When summarizing therapeutic/biomarker data:

- Separate statements derived from variant-level evidence from those based
  on codon-level or gene-level evidence.
- Use variant-level data to justify stronger, more focused Aims.
- Use codon-level and gene-level data mainly to motivate exploratory or
  hypothesis-generating Aims.

==================== SECTION SPECS ====================

1) FUNCTIONAL IMPACT

- Summarize the gene's role (oncogene vs tumor suppressor vs DNA repair, etc.)
  ONLY if GENE ROLE or other input text explicitly states it.
- If variant-level functional or clinical evidence is present (hotspot status,
  pathogenic calls, sensitivity/resistance evidence), state whether the variant
  is activating, loss-of-function, or risk-associated at a high level.
- If evidence is limited, fall back to gene-level biology and clearly state that
  the variant's specific functional impact is uncertain.
- Do NOT invent detailed mechanisms (no structural models, no speculative
  signaling diagrams).
- High-level phrases like "activating hotspot" or "likely loss-of-function" are
  allowed only when implied by the input (e.g., hotspot annotations, "oncogenic"
  labels, pathogenic calls).
- ClinVar pathogenicity calls cover many conditions. When querying in {tumor_type}
  context: use only ClinVar entries whose condition matches {tumor_type} or is
  cancer-general. IGNORE ClinVar entries for unrelated conditions (e.g., germline
  syndromes, other cancer types, hereditary conditions) — do NOT mention them in
  this section. If no tumor-relevant ClinVar entries exist, omit ClinVar entirely.

2) TUMOR BIOLOGY & MODELS

- Use tumor-specific cohort data when has_tumor_specific_cbioportal is TRUE:
  - Report prevalence of {gene} mutations and, if available, prevalence of the
    specific {variant} in {tumor_type}.
  - Summarize key co-occurring and mutually exclusive genes, highlighting
    canonical pathway partners ONLY if these genes are listed in the input.

- Use dependency and model information:
  - Summarize gene dependency: whether the gene appears essential or
    context-specific, based on the provided dependency summary.
  - List AVAILABLE MODEL CELL LINES (and organoids/PDX if present) that carry
    this variant, including their tissue of origin.
  - Explicitly highlight model gaps such as:
    "no {tumor_type} cell lines with this variant; only models from other tissues".

- When has_tumor_specific_cbioportal is FALSE:
  - Skip prevalence and co-mutation discussion entirely.
  - Focus on model availability, dependency, and gene-level biology.

3) THERAPEUTIC & BIOMARKER LANDSCAPE (RESEARCH VIEW)

- Summarize therapeutic evidence as CONTEXT for research, not as treatment advice.
  Use only drugs and therapies that appear in the input.

- Clearly distinguish:
  - FDA-approved therapies that match the gene/variant and tumor context.
  - FDA approvals in other tumor types.
  - Clinical trial evidence (phase and tumor context if provided).
  - Preclinical sensitivity and resistance signals.

- Present these as:
  - What is well supported (e.g., "multiple trials show sensitivity to
    {drug} in {gene}-altered {tumor_type}").
  - What is ambiguous or negative (e.g., "no benefit observed in a specified
    biomarker-defined subgroup").

- FDA approvals are categorized by whether they COVER the queried variant:
  - MATCHED APPROVALS: If FDA approval says "[GENE] alteration/mutation" it COVERS
    any variant in that gene. If it specifies this exact variant, it is a direct match.
    Do NOT hedge with "gene-level rather than variant-specific" — the drug IS approved
    for this variant. Ignore "(gene-level)" or "(variant-level)" annotations on FDA entries.
  - UNMATCHED NEAR-MISSES: Listed under "FDA Codon-Level (not for queried variant)".
    State explicitly: "approved for [approved variant], not [queried variant]".

- When conflicting_evidence is present, explicitly describe the nature of the
  conflict (same drug/context with different outcomes vs expected sequencing effects).

- Do NOT recommend what treatment a patient should receive.
  Phrase everything as "data suggest", "trials have shown", "evidence indicates", etc.

4) EVIDENCE QUALITY & OPEN QUESTIONS

- Use overall_quality, research_priority, knowledge_gaps, significant_gaps,
  and conflicting_evidence to:

  - Briefly rate the strength of current evidence
    (e.g., "well characterized in {tumor_type}" vs
     "limited data, largely extrapolated from other cancers").

  - List 3–6 specific OPEN QUESTIONS that matter for research, such as:
    - Unclear dependency on the variant vs gene-level activation.
    - Unresolved mechanisms of resistance to a listed drug.
    - Uncertain impact of specific co-occurring mutations.
    - Whether codon-level findings generalize to this variant.

5) EMERGING RESEARCH PROGRAM (AIMS)

- This section should read like a 1–2 year translational research plan, not a
  grab-bag of random experiments.

- When research_priority is "medium" or higher OR knowledge_gaps /
  significant_gaps exist:

  - Organize this section into 2–3 AIMS.

  - For each AIM, provide:
    - Aim title (one line).
    - Rationale: 2–3 sentences connecting back to evidence above
      (hotspot status, co-mutations, dependency, model availability,
       therapeutic signals, gaps, variant/codon/gene-level differences).
    - 2–3 APPROACH bullet points, each including:
      * Model system: specific existing cell lines or organoids from
        AVAILABLE MODEL CELL LINES, or an explicit proposal for isogenic
        editing in a relevant {tumor_type} background when models are missing.
      * Perturbation: drugs, inhibitors, or genetic manipulations that are
        ALREADY PRESENT in the therapeutic / evidence / biology input
        (e.g., named inhibitors, pathway partners, gene knockouts).
      * Readouts: core experimental readouts such as viability, apoptosis,
        pathway phosphorylation, DNA damage markers, clonal outgrowth,
        transcriptomic or proteomic profiling.

- You MAY use generic experimental techniques (CRISPR knock-in/knockout, shRNA,
  overexpression, standard viability assays, phospho-protein assays, RNA-seq,
  etc.) even if not explicitly named in the evidence.

- You MUST NOT invent:
  - New drug names or investigational agents not present in the input.
  - New biomarkers or pathways not implied by the provided genes and evidence.

- When significant_gaps[*].suggested_studies is present:
  - Treat these as high-priority seeds for Aims and Approaches.
  - Rewrite and expand them into full Aim + rationale + approach items,
    rather than ignoring them.

- When research_priority is "low" and no significant knowledge gaps exist:
  - You may omit EMERGING RESEARCH PROGRAM or limit it to 1 conservative Aim.

==================== GENERAL CONSTRAINTS ====================

- Use ONLY facts from the provided data for gene/variant/tumor, plus generic
  lab techniques as allowed above.
- Do NOT use outside knowledge of drugs, trials, prevalence, or mechanisms.
- Do NOT give clinical recommendations or management advice.
- Keep each section concise but information-dense:
  - 2–4 sentences for sections 1–4.
  - 2–3 Aims in section 5, each with 2–3 short approach bullets.
- Write for a translational PI planning a research program around this variant
  in this cancer type.
- NEVER mention drugs, clinical trials, or treatments NOT explicitly listed in
  the DATABASE EVIDENCE section.
- ALWAYS cite sources for statistics with markdown links as provided.
"""

SYNTHESIS_USER_PROMPT = """Write a research dossier for this variant. Use ONLY the data below.

Gene: {gene}
Variant: {variant}
Tumor Type: {tumor_type}

## DATA FLAGS
has_tumor_specific_cbioportal: {has_tumor_specific_cbioportal}
has_civic_assertions: {has_civic_assertions}
has_fda_biomarker_evidence: {has_fda_biomarker_evidence}
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
Research priority: {research_priority}
Well-characterized: {well_characterized_text}
Gaps: {known_gaps_text}
Significant gaps: {significant_gaps_json}
Conflicts: {conflicting_evidence_text}

Respond with valid JSON only:
{{
  "functional_impact": "2–4 sentences on gene role and variant-specific functional evidence. Only use facts from input.",
  "tumor_biology": "2–4 sentences on prevalence (if has_tumor_specific_cbioportal), co-mutations, dependency, and available model cell lines. Flag model gaps explicitly.",
  "therapeutic_landscape_prose": "3–5 sentences synthesizing the therapeutic landscape for research context. Distinguish FDA-approved, clinical, and preclinical. Use only drugs from DATABASE EVIDENCE.",
  "therapeutic_landscape": {{
    "fda_approved": ["drug (locus-level, approved for VARIANT if gene/codon-level, tumor)"],
    "clinical_evidence": ["drug (locus-level, tumor) - source"],
    "preclinical": ["drug (locus-level) - source"],
    "resistance_mechanisms": ["drug - mechanism (locus-level)"]
  }},
  "evidence_quality": "1–2 sentences rating evidence strength. Then list 3–6 open questions as a JSON array in open_questions.",
  "open_questions": ["open question 1", "open question 2", "..."],
  "research_program": [
    {{
      "aim_title": "Aim 1: one-line title",
      "rationale": "2–3 sentences connecting gaps and evidence to this aim.",
      "approaches": [
        "Model: [cell line or isogenic system]. Perturbation: [drug or genetic]. Readout: [assay].",
        "..."
      ]
    }}
  ],
  "key_references": ["PMIDs, databases from evidence"],
  "evidence_tags": ["direct clinical data | preclinical only | pan-cancer extrapolation | limited evidence"]
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
    """Create prompt for evidence synthesis (research dossier)."""
    tumor_display = tumor_type or "Pan-cancer"

    overall_quality = evidence_assessment.get("overall_quality", "minimal")
    research_priority = evidence_assessment.get("research_priority", "unknown")
    well_char = evidence_assessment.get("well_characterized", []) or []
    gaps = evidence_assessment.get("knowledge_gaps", []) or []
    conflicts = evidence_assessment.get("conflicting_evidence", []) or []
    significant_gaps = evidence_assessment.get("significant_gaps", []) or []

    if data_availability is None:
        data_availability = {}

    # Build locus match text for LLM context
    if locus_match_summary:
        locus_match_text = locus_match_summary.get(
            "summary_text", "No locus match data available."
        )
        if locus_match_summary.get("is_all_gene_level"):
            locus_match_text += " CAUTION: No variant-specific evidence found - use gene-level inferences carefully."
        elif not locus_match_summary.get("has_variant_specific"):
            locus_match_text += " WARNING: Limited variant-specific data."
    else:
        locus_match_text = "No locus match data available."

    # Build tumor match text for LLM context
    if tumor_match_summary:
        tumor_match_text = tumor_match_summary.get(
            "summary_text", "No tumor match data available."
        )
        if tumor_match_summary.get("is_all_other"):
            other_tumors = ", ".join(
                tumor_match_summary.get("other_tumor_types", [])[:3]
            )
            tumor_match_text += f" CAUTION: All evidence is from other tumor types ({other_tumors}) - extrapolation required."
        elif not tumor_match_summary.get("has_tumor_specific"):
            tumor_match_text += " WARNING: Limited tumor-specific evidence."
    else:
        tumor_match_text = "No tumor match data available."

    # Substitute the three context variables in the system prompt using simple replacement
    # (cannot use str.format() because the prompt contains many literal {braces} as examples)
    system_content = (
        SYNTHESIS_SYSTEM_PROMPT
        .replace("{gene}", gene)
        .replace("{variant}", variant)
        .replace("{tumor_type}", tumor_display)
    )

    user_content = SYNTHESIS_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        has_tumor_specific_cbioportal=str(
            data_availability.get("has_tumor_specific_cbioportal", False)
        ).upper(),
        has_civic_assertions=str(
            data_availability.get("has_civic_assertions", False)
        ).upper(),
        has_fda_biomarker_evidence=str(
            data_availability.get("has_fda_biomarker_evidence", False)
        ).upper(),
        has_vicc_evidence=str(
            data_availability.get("has_vicc_evidence", False)
        ).upper(),
        locus_match_text=locus_match_text,
        tumor_match_text=tumor_match_text,
        biological_context=biological_context or "No cBioPortal data available.",
        resistance_summary=resistance_summary or "No resistance signals.",
        sensitivity_summary=sensitivity_summary or "No sensitivity signals.",
        evidence_summary=(evidence_summary or "No database evidence.").strip()[:4000],
        cross_source_synthesis=cross_source_synthesis
        or "No cross-source synthesis available.",
        literature_summary=literature_summary or "No literature search performed.",
        overall_quality=overall_quality,
        research_priority=research_priority,
        well_characterized_text="; ".join(well_char) or "None.",
        known_gaps_text="; ".join(gaps) or "None.",
        significant_gaps_json=json.dumps(significant_gaps),
        conflicting_evidence_text="; ".join(conflicts) or "None.",
    )

    return [
        {"role": "system", "content": system_content},
        {"role": "user", "content": user_content},
    ]


# =============================================================================
# STAGE 2: CROSS-SOURCE SYNTHESIS
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
        cross_source_synthesis=cross_source_synthesis
        or "No cross-source synthesis available.",
        sensitivity_summary=sensitivity_summary or "No sensitivity signals.",
        resistance_summary=resistance_summary or "No resistance signals.",
    )

    return [
        {"role": "system", "content": CROSS_SOURCE_SYSTEM_PROMPT},
        {"role": "user", "content": user_content},
    ]

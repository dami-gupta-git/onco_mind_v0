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
6. Generate TESTABLE RESEARCH HYPOTHESES based on evidence gaps (see hypothesis generation rules below).

=== CRITICAL: CALIBRATE CONFIDENCE TO EVIDENCE QUALITY ===

EVIDENCE QUALITY RULES (STRICT - MUST FOLLOW):

When overall_quality is "limited" or "minimal":
- Your functional_summary MUST remain GENERIC (e.g., "JAK1 is involved in cytokine signaling" NOT "JAK1 V657F causes loss of function")
- You MUST NOT assign oncogene/tumor-suppressor roles unless the GENE ROLE section explicitly states it
- You MUST NOT make drug-response predictions unless FDA approvals, CIViC assertions, or VICC evidence explicitly exist for THIS variant
- Keep synthesis brief (2-3 sentences per section), focused on what's UNKNOWN
- Your evidence_tags MUST NOT include "direct clinical data" unless CIViC/FDA/VICC evidence exists for THIS specific variant

When has_tumor_specific_cbioportal_data is FALSE (check the BIOLOGICAL CONTEXT section):
- You MUST NOT describe prevalence as being "in [tumor type]" - state it's pan-cancer data
- State explicitly: "No [tumor_type]-specific data; pan-cancer statistics shown"
- Do not extrapolate pan-cancer co-mutations to the specific tumor context

EVIDENCE TAG VALIDATION (STRICT):
- Use "direct clinical data" ONLY if: CIViC assertions, FDA approvals, or Phase 2/3 trial data exist for THIS variant in THIS tumor
- Use "preclinical only" if: only DepMap/cell line/in vitro data exists
- Use "pan-cancer extrapolation" if: cBioPortal data is pan-cancer, not tumor-specific
- Use "pathway-level inference" if: inferring from general pathway knowledge (NOT direct evidence)
- Use "nearby-variant extrapolation" if: making inferences from other variants in same gene
- Do NOT use "direct clinical data" while simultaneously listing "clinical evidence" as a knowledge gap - this is contradictory

=== HARD CONSTRAINTS ===

- You MUST treat the user message as the ONLY source of truth. Do NOT introduce facts (numbers, drug names, mechanisms, PMIDs) that are not explicitly present in the text.
- If a fact is not present, speak generically (e.g., "commonly", "rarely", "no trials described here") instead of inventing specifics.
- ALWAYS cite the source when quoting statistics. When citing prevalence, sample counts, or percentages, include the source with its markdown link exactly as provided in the evidence (e.g., "52% of samples ([cBioPortal: mel_tcga](https://...))"). Never quote a number without attribution.
- Clearly distinguish:
  - direct data from the evidence block,
  - extrapolations from nearby variants or pathways,
  - speculative hypotheses. Mark extrapolations and hypotheses explicitly.
- Do not hide conflicting evidence; briefly describe how and why sources differ.
- Ensure the output is a concise discussion-style synthesis with explicit, testable research questions.
- Include an explicit note if data are uncertain or inferential, and label such sections accordingly.

Your narrative should read like a concise discussion section in a strong cancer paper: integrative, precise, and explicitly pointing to testable research questions, but strictly bounded by the provided evidence.

=== THERAPEUTIC SIGNALS HANDLING ===

When the THERAPEUTIC SIGNALS section contains Sensitivity or Resistance data:
- ALWAYS include sensitivity signals in therapeutic_landscape.clinical_evidence or therapeutic_landscape.preclinical as appropriate
- ALWAYS include resistance signals in therapeutic_landscape.resistance_mechanisms
- Generate at least one research hypothesis that addresses resistance mechanisms, bypass pathways, or sensitivity validation when resistance or sensitivity evidence exists

=== RESEARCH HYPOTHESIS GENERATION RULES ===

Generate 2-3 testable research hypotheses based on the evidence gaps. Focus on:
- BIOLOGICAL MECHANISM: What functional experiments could clarify the variant's effect?
- PRECLINICAL TESTING: What cell line models or drug sensitivity assays are warranted?
- CO-MUTATION EFFECTS: How might co-occurring alterations modify the phenotype?
- CROSS-HISTOLOGY: Could findings from other tumor types apply here?

HYPOTHESIS FORMAT:
- Each hypothesis must be specific and testable (not vague directions)
- Emphasize novelty and tractability
- Connect at least two evidence elements (e.g., gap + existing data point)
- Avoid clinical treatment recommendations - focus on research questions
- CRITICAL: Each hypothesis MUST explicitly state its EVIDENCE BASIS using one of these tags:
  * [Direct Clinical Data] - hypothesis builds on FDA approvals, CIViC assertions, or Phase 2/3 trials for THIS variant in THIS tumor
  * [Preclinical Data] - hypothesis builds on DepMap, cell line, or in vitro data
  * [Pan-Cancer Extrapolation] - hypothesis extrapolates from data in other tumor types
  * [Nearby-Variant Inference] - hypothesis extrapolates from other variants in the same gene/domain
  * [Pathway-Level Inference] - hypothesis infers from general pathway biology, not direct evidence

EXAMPLES:
- Good: "[Preclinical Data] Given the lack of functional data for JAK1 V657F despite its recurrence in T-ALL, isogenic knock-in models could determine whether this variant causes gain- or loss-of-function signaling."
- Good: "[Nearby-Variant Inference] The absence of preclinical drug sensitivity data for this variant, combined with its structural similarity to JAK2 V617F, suggests testing JAK inhibitor panels in cell lines harboring this mutation."
- Good: "[Pan-Cancer Extrapolation] While EGFR L858R shows sensitivity to osimertinib in NSCLC, testing this response in breast cancer models would determine cross-histology applicability."
- Bad: "Patients with this variant should be treated with JAK inhibitors" (clinical recommendation, not hypothesis)
- Bad: "More research is needed" (vague, not testable)
- Bad: "JAK inhibitors should be tested" (missing evidence basis tag)
"""

RESEARCH_USER_PROMPT = """Synthesize research knowledge about this variant using ONLY the evidence below.

Gene: {gene}
Variant: {variant}
Tumor Type: {tumor_type}

## DATA AVAILABILITY FLAGS
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
3. Generate 2-3 specific, TESTABLE research hypotheses based on the evidence gaps. Each hypothesis should:
   - Be concrete and experimentally addressable
   - Connect multiple evidence elements (gap + existing data point)
   - Focus on biological mechanism, preclinical testing, or co-mutation effects
   - Avoid clinical recommendations
   - MUST include an evidence basis tag at the start: [Direct Clinical Data], [Preclinical Data], [Pan-Cancer Extrapolation], [Nearby-Variant Inference], or [Pathway-Level Inference]

CALIBRATION RULES (MUST FOLLOW):
- If overall_quality is "limited" or "minimal": Keep functional_summary GENERIC (describe the gene's general function, not variant-specific effects). Focus on gaps.
- If has_tumor_specific_cbioportal_data is FALSE: State "pan-cancer data shown; no {{tumor_type}}-specific prevalence available."
- If has_civic_assertions, has_fda_approvals, AND has_vicc_evidence are ALL FALSE: Do NOT use "direct clinical data" in evidence_tags.
- Do NOT assign oncogene/tumor-suppressor roles unless GENE ROLE section explicitly states it.
- CRITICAL: If THERAPEUTIC SIGNALS section contains resistance data, you MUST include those resistance mechanisms in therapeutic_landscape.resistance_mechanisms.
- CRITICAL: If THERAPEUTIC SIGNALS section contains sensitivity data, you MUST include those drugs in therapeutic_landscape.clinical_evidence or therapeutic_landscape.preclinical.

CRITICAL: For biological_context, include the source citation ONCE at the end of all cBioPortal statistics. Look for "cite as: ([cBioPortal: ...](...)" in the BIOLOGICAL CONTEXT section and copy that markdown link.

Respond ONLY with valid JSON:
{{
  "functional_summary": "IF overall_quality is 'limited'/'minimal': Generic gene function only (e.g., 'JAK1 encodes a tyrosine kinase involved in cytokine receptor signaling'). IF 'moderate'/'comprehensive': Variant-specific effects with citations.",
  "biological_context": "2–3 sentences. IF has_tumor_specific_cbioportal_data is FALSE, state 'Pan-cancer data; no {{tumor_type}}-specific data available.' Cite source ONCE at end.",
  "therapeutic_landscape": {{
    "fda_approved": ["ONLY drugs explicitly listed in FDA APPROVALS section. Empty array if none."],
    "clinical_evidence": ["ONLY drugs with CIViC/VICC clinical evidence OR from THERAPEUTIC SIGNALS Sensitivity. Empty array if none."],
    "preclinical": ["Drugs from DepMap or preclinical sources only."],
    "resistance_mechanisms": ["MUST include signals from THERAPEUTIC SIGNALS Resistance section. Also include any mechanisms from evidence. Do not infer beyond what's stated."]
  }},
  "evidence_assessment": {{
    "overall_quality": "{overall_quality}",
    "well_characterized": {well_characterized_json},
    "knowledge_gaps": {known_gaps_json},
    "conflicting_evidence": {conflicting_evidence_json}
  }},
  "research_implications": "2–3 sentences. For limited/minimal quality: focus on basic characterization needs. For moderate/comprehensive: specific hypotheses.",
  "research_hypotheses": [
    "[Evidence Basis Tag] REQUIRED: 2-3 specific, testable hypotheses based on evidence gaps.",
    "Each MUST start with one of: [Direct Clinical Data], [Preclinical Data], [Pan-Cancer Extrapolation], [Nearby-Variant Inference], or [Pathway-Level Inference]",
    "Example: '[Pan-Cancer Extrapolation] Given KRAS G12C sensitivity to sotorasib in NSCLC, testing this response in breast cancer models would determine cross-histology applicability.'"
  ],
  "key_references": ["PMIDs, trial IDs, or database/source names appearing in the evidence text."],
  "evidence_tags": [
    "MUST be consistent with DATA AVAILABILITY FLAGS above.",
    "Use 'direct clinical data' ONLY if has_civic_assertions OR has_fda_approvals is TRUE.",
    "Use 'preclinical only' if only DepMap data exists.",
    "Use 'pan-cancer extrapolation' if has_tumor_specific_cbioportal_data is FALSE.",
    "Use 'limited evidence' if overall_quality is 'limited' or 'minimal'."
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
    data_availability: dict | None = None,
    resistance_summary: str = "",
    sensitivity_summary: str = "",
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
        data_availability: Dict with boolean flags for data presence (tumor_specific_cbioportal, civic_assertions, etc.)
        resistance_summary: Concise summary of resistance signals
        sensitivity_summary: Concise summary of sensitivity signals

    Returns:
        Messages list for LLM API call.
    """
    tumor_display = tumor_type or "Pan-cancer"

    overall_quality = evidence_assessment.get("overall_quality", "minimal")
    well_char = evidence_assessment.get("well_characterized", []) or []
    gaps = evidence_assessment.get("knowledge_gaps", []) or []
    conflicts = evidence_assessment.get("conflicting_evidence", []) or []

    # Default data availability flags
    if data_availability is None:
        data_availability = {}

    has_tumor_specific_cbioportal = data_availability.get("has_tumor_specific_cbioportal", False)
    has_civic_assertions = data_availability.get("has_civic_assertions", False)
    has_fda_approvals = data_availability.get("has_fda_approvals", False)
    has_vicc_evidence = data_availability.get("has_vicc_evidence", False)

    # Build the JSON-safe strings for the template
    user_content = RESEARCH_USER_PROMPT.format(
        gene=gene,
        variant=variant,
        tumor_type=tumor_display,
        has_tumor_specific_cbioportal=str(has_tumor_specific_cbioportal).upper(),
        has_civic_assertions=str(has_civic_assertions).upper(),
        has_fda_approvals=str(has_fda_approvals).upper(),
        has_vicc_evidence=str(has_vicc_evidence).upper(),
        biological_context=biological_context or "No cBioPortal or co-mutation data available.",
        resistance_summary=resistance_summary or "No resistance signals identified.",
        sensitivity_summary=sensitivity_summary or "No sensitivity signals identified.",
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

# OncoMind LLM Integration

This document describes how OncoMind uses Large Language Models (LLMs) for evidence synthesis and research analysis.

---

## Table of Contents

1. [Overview](#overview)
2. [LLM Research Synthesis](#llm-research-synthesis)
3. [Cross-Source Drug Analysis](#cross-source-drug-analysis)
4. [Paper Relevance Scoring](#paper-relevance-scoring)
5. [Variant Knowledge Extraction](#variant-knowledge-extraction)

---

## Overview

OncoMind uses LLMs for tasks where they excel: synthesis, narrative generation, and extracting structured knowledge from unstructured text. All LLM calls are optional and run through a unified service layer.

### Multi-Provider Support

OncoMind uses [litellm](https://github.com/BerriAI/litellm) for multi-provider LLM support:

```python
# Supported models
"claude-sonnet-4-20250514"  # Default, recommended
"claude-3-5-haiku-20241022" # Faster, lower cost
"gpt-4o-mini"               # OpenAI
"gpt-4o"                    # OpenAI
"gpt-4-turbo"               # OpenAI
```

### Parallel Execution Model

LLM Research Synthesis and Cross-Source Drug Analysis run **in parallel** using `asyncio.gather()`:

```python
# Both LLM calls run concurrently
llm_insight_task = self._get_llm_insight(evidence, gene, variant, tumor_type)
cross_source_task = self._get_cross_source_analysis(evidence, gene, variant, tumor_type)

llm_insight, cross_source_analysis = await asyncio.gather(
    llm_insight_task,
    cross_source_task,
    return_exceptions=True
)
```

This parallel execution reduces total latency compared to sequential calls.

---

## LLM Research Synthesis

**Purpose:** Generate research-focused narrative synthesis calibrated to evidence quality.

### Research-Oriented Prompt Design

The LLM receives structured context with explicit data availability flags:

```
## DATA AVAILABILITY FLAGS
has_tumor_specific_cbioportal_data: TRUE/FALSE
has_civic_assertions: TRUE/FALSE
has_fda_approvals: TRUE/FALSE
has_vicc_evidence: TRUE/FALSE

## BIOLOGICAL CONTEXT
{cBioPortal prevalence, co-mutations, DepMap essentiality}

## THERAPEUTIC SIGNALS
Sensitivity: {synthesized from Evidence.get_sensitivity_summary()}
Resistance: {synthesized from Evidence.get_resistance_summary()}

## DATABASE EVIDENCE
{compact summary from Evidence.get_evidence_summary_for_llm()}

## LITERATURE FINDINGS
{from Evidence.get_literature_summary_for_llm()}

## EVIDENCE GAPS
{from EvidenceGaps.to_dict_for_llm()}
```

### Evidence Quality Calibration

The LLM is instructed to calibrate confidence based on evidence quality:

| Quality | LLM Behavior |
|---------|--------------|
| `limited` / `minimal` | Generic gene function only, no variant-specific claims |
| No tumor-specific cBioPortal | Must state "pan-cancer data; no {tumor}-specific data" |
| No CIViC/FDA/VICC | Must NOT use "direct clinical data" evidence tag |

### Research Hypothesis Generation

The LLM generates 2-3 testable hypotheses, each with an evidence basis tag:

| Tag | Meaning |
|-----|---------|
| `[Direct Clinical Data]` | Builds on FDA/CIViC/Phase 2-3 for THIS variant |
| `[Preclinical Data]` | Builds on DepMap/cell line/in vitro data |
| `[Pan-Cancer Extrapolation]` | Extrapolates from other tumor types |
| `[Nearby-Variant Inference]` | Extrapolates from other variants in same gene |
| `[Pathway-Level Inference]` | Infers from general pathway biology |

Example hypotheses:
```
[Preclinical Data] Given the lack of functional data for JAK1 V657F despite
its recurrence in T-ALL, isogenic knock-in models could determine whether
this variant causes gain- or loss-of-function signaling.

[Pan-Cancer Extrapolation] While EGFR L858R shows sensitivity to osimertinib
in NSCLC, testing this response in breast cancer models would determine
cross-histology applicability.
```

### LLMInsight Model

The LLM returns structured insight via the `LLMInsight` model:

```python
class LLMInsight(BaseModel):
    llm_summary: str                    # Combined narrative
    rationale: str                      # Research implications
    clinical_trials_available: bool
    therapeutic_evidence: list          # Empty (comes from Evidence)
    references: list[str]               # PMIDs, trial IDs, sources

    # Raw component data for UI formatting
    functional_summary: str | None
    biological_context: str | None
    therapeutic_landscape: dict | None  # fda_approved, clinical_evidence, preclinical, resistance_mechanisms

    # Research assessment fields
    evidence_quality: str | None        # "comprehensive" | "moderate" | "limited" | "minimal"
    knowledge_gaps: list[str]
    well_characterized: list[str]
    conflicting_evidence: list[str]
    research_implications: str
    evidence_tags: list[str]            # Transparency tags
    research_hypotheses: list[str]      # With evidence basis tags
```

---

## Cross-Source Drug Analysis

**Purpose:** Synthesize therapeutic evidence across multiple sources to identify corroboration, conflicts, and emerging signals.

### Overview

Cross-Source Drug Analysis is a separate LLM analysis that runs **in parallel** with the main LLM synthesis. It groups drug evidence across CGI, CIViC, VICC, and Literature sources to provide a unified view of the therapeutic landscape.

### Input: Drug Evidence Grouping

Before calling the LLM, evidence is grouped by drug across all sources:

```python
def _build_cross_source_drug_summary(evidence: Evidence) -> str:
    """
    Groups drugs from:
    - CGI biomarkers (FDA guidelines, NCCN guidelines, clinical trials, preclinical)
    - CIViC assertions and evidence items
    - VICC MetaKB entries
    - Literature mentions (from PubMed articles)

    Returns formatted summary like:

    === VEMURAFENIB ===
    Sources: CGI, CIViC, VICC, Literature

    CGI:
      - FDA guidelines: Responsive in Melanoma (variant-level match)
      - NCCN guidelines: Responsive in NSCLC (gene-level match)

    CIViC:
      - Tier I Level A: Sensitivity/Response in Melanoma (BRAF V600E)
      - Evidence: 3 items supporting sensitivity

    VICC:
      - oncokb: Level 1, Responsive in Melanoma
      - civic: Level A, Sensitive in Melanoma

    Literature:
      - 5 papers mention this drug
      - Signal types: sensitivity (3), resistance (2)
    """
```

### Output Structure

The LLM returns structured analysis:

```python
cross_source_analysis = {
    "strongest_evidence": [
        {
            "drug": "vemurafenib",
            "sources": ["CGI", "CIViC", "VICC", "Literature"],
            "corroboration_summary": "FDA-approved with consistent sensitivity signals across 4 sources",
            "biological_rationale": "BRAF V600E causes constitutive kinase activation; vemurafenib directly inhibits mutant BRAF"
        }
    ],
    "conflicting_signals": [
        {
            "drug": "dabrafenib",
            "conflict_type": "response_vs_resistance",
            "sources_supporting": ["CGI", "CIViC"],
            "sources_conflicting": ["Literature"],
            "likely_explanation": "Literature reports acquired resistance after initial response; not primary resistance"
        }
    ],
    "emerging_targets": [
        {
            "drug": "PLX8394",
            "source": "Literature",
            "evidence_type": "preclinical",
            "potential": "Paradox-breaker RAF inhibitor showing activity in BRAF V600E cell lines"
        }
    ],
    "key_gaps": [
        "No combination therapy data for BRAF + MEK inhibitors in this tumor type",
        "Resistance mechanisms not characterized for third-generation inhibitors"
    ]
}
```

### Analysis Categories

| Category | Description | When Flagged |
|----------|-------------|--------------|
| **Strongest Evidence** | Drugs with corroboration across ≥2 independent sources | Multiple sources agree on response type |
| **Conflicting Signals** | Drugs where sources disagree | One source says sensitive, another says resistant |
| **Emerging Targets** | Single-source preclinical or early evidence | Only Literature or CGI preclinical mentions |
| **Key Gaps** | Expected drugs not found or extrapolation concerns | Known drug class missing, or only gene-level match |

### Conflict Detection Logic

The LLM is instructed to identify likely explanations for conflicts:

| Conflict Pattern | Likely Explanation |
|------------------|-------------------|
| CIViC sensitive + Literature resistant | Acquired resistance after initial response |
| Gene-level sensitive + Variant-level resistant | Variant-specific resistance mechanism |
| Different tumor types | Tissue-specific response differences |
| Different evidence levels | Early preclinical vs mature clinical data |

### UI Display

In Streamlit, cross-source analysis appears in the "LLM Analysis" tab:

```
### Cross-Source Drug Analysis

#### Strongest Evidence
- **Vemurafenib**: Corroborated across CGI, CIViC, VICC, Literature
  - FDA-approved with consistent sensitivity signals
  - Biological rationale: Direct inhibition of V600E-mutant BRAF kinase

#### Conflicting Signals
- **Dabrafenib**: CGI/CIViC show sensitivity, Literature reports resistance
  - Likely explanation: Acquired resistance after initial response

#### Emerging Targets
- **PLX8394** (Literature, preclinical): Paradox-breaker RAF inhibitor

#### Key Gaps
- No combination therapy data for this tumor type
```

### Result Model Integration

Cross-source analysis is stored in the `Result` model:

```python
class Result(BaseModel):
    evidence: Evidence
    llm: LLMInsight | None
    cross_source_analysis: dict | None  # New field
```

---

## Paper Relevance Scoring

The LLM service scores individual papers for relevance to a query variant:

```python
async def score_paper_relevance(
    title: str,
    abstract: str | None,
    tldr: str | None,
    gene: str,
    variant: str,
    tumor_type: str | None,
) -> dict:
    """
    Returns:
        relevance_score: float 0-1
        is_relevant: bool (True if >= 0.6)
        signal_type: "resistance" | "sensitivity" | "mixed" | "prognostic" | "unclear"
        drugs_mentioned: list[str]
        key_finding: str
        confidence: float 0-1
    """
```

---

## Variant Knowledge Extraction

Extracts structured knowledge from multiple papers:

```python
async def extract_variant_knowledge(
    gene: str,
    variant: str,
    tumor_type: str,
    paper_contents: list[dict],
) -> dict:
    """
    Returns:
        mutation_type: "primary" | "secondary" | "both" | "unknown"
        resistant_to: list[{drug, evidence, mechanism}]
        sensitive_to: list[{drug, evidence}]
        clinical_significance: str
        evidence_level: str
        key_findings: list[str]
        confidence: float
    """
```

---

## Code References

- **LLM service**: [service.py](src/oncomind/llm/service.py) - All LLM calls go through here
- **Prompts**: [prompts.py](src/oncomind/llm/prompts.py) - Prompt templates
- **Drug grouping**: [insight_builder.py](src/oncomind/insight_builder/insight_builder.py) - `_build_cross_source_drug_summary()`
- **Parallel execution**: [insight_builder.py](src/oncomind/insight_builder/insight_builder.py) - `_run_llm_layer()`
- **LLMInsight model**: [llm_insight.py](src/oncomind/llm/llm_insight.py)
- **Result model**: [result.py](src/oncomind/models/result.py)

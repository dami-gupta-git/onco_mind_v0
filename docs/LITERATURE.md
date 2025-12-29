# Literature Synthesis Tools & Strategy

This document evaluates external tools for literature synthesis and outlines OncoMind's approach.

---

## Tool Evaluation

| Tool | What It Is | Status |
|------|-----------|--------|
| **litellm** | LLM API router (OpenAI, Anthropic, Cohere, etc.) | **✅ Implemented** - powers multi-provider LLM support |
| **hint-lab/pubmed-agent** | Agentic PubMed search with tool use | No - adds complexity, we have the APIs already |
| **llm-literature-review-tool** | Generic literature review automation | No - overkill for focused oncology use case |
| **LitSense** | NIH semantic search for PubMed | Maybe - supplementary data source |
| **PubMedSummarizer** | Basic LLM summarization of abstracts | No - trivial to build ourselves |

### litellm ✅ IMPLEMENTED

**What it does:** Unified interface to 100+ LLM providers. One API call, swap models via config.

**Current implementation in OncoMind:**
- Integrated in `src/oncomind/llm/service.py`
- Supports Claude (Sonnet 4, Haiku), GPT-4o, GPT-4o-mini, GPT-4-turbo
- Handles provider differences via unified `acompletion()` interface
- Configured via `ONCOMIND_MODEL` environment variable

**Supported models:**
```python
SUPPORTED_MODELS = [
    "claude-sonnet-4-20250514",    # Default, recommended
    "claude-3-5-haiku-20241022",   # Faster, lower cost
    "gpt-4o-mini",                 # OpenAI
    "gpt-4o",                      # OpenAI
    "gpt-4-turbo",                 # OpenAI
]
```

**Usage:**
```python
from litellm import acompletion

# Same code, different models
response = await acompletion(model="claude-sonnet-4-20250514", messages=messages)
```

**Link:** https://github.com/BerriAI/litellm

---

### hint-lab/pubmed-agent

**What it does:** Agentic PubMed querying with LLM-driven search refinement.

**Why we're skipping it:**
- We already have PubMed and Semantic Scholar APIs integrated
- Adds agent framework dependency (LangChain or similar)
- Our value is in cancer-specific extraction, not generic search
- Can build targeted search logic ourselves with less overhead

**Link:** https://github.com/Hint-lab/pubmed-agent

---

### llm-literature-review-tool

**What it does:** Automated literature review pipeline - search, filter, summarize, synthesize.

**Why we're skipping it:**
- Generic tool, not cancer-focused
- We need domain-specific extraction (resistance signals, drug sensitivity)
- Overkill for our focused use case
- Would require significant customization anyway

---

### LitSense

**What it does:** NIH's AI-powered semantic search for PubMed. Finds related sentences across papers.

**Potential use:** Supplementary search when Semantic Scholar + PubMed don't find enough.

**Considerations:**
- API access unclear (may be web-only)
- Could help find obscure case reports
- Lower priority than core functionality

**Link:** https://www.ncbi.nlm.nih.gov/research/litsense/

---

### PubMedSummarizer (Generic)

**What it does:** LLM wrapper to summarize PubMed abstracts.

**Why we're skipping it:**
- Trivial to implement ourselves
- We need structured extraction, not generic summaries
- Our `LiteratureSynthesis` model requires specific fields (key papers, drug signals, etc.)

---

## OncoMind's Literature Strategy

### Current Implementation ✅

1. **Search:** Semantic Scholar API + PubMed E-utilities
2. **Filter:** Cancer-relevant papers, recency weighting
3. **Score:** LLM-powered relevance scoring with signal type extraction
4. **Extract:** Structured knowledge extraction into `LiteratureKnowledge` model
5. **Synthesize:** Summary with key findings, therapeutic signals, and evidence gaps

### Implemented Models

```python
class DrugResistance(BaseModel):
    drug: str
    evidence: str  # "in vitro" | "preclinical" | "clinical" | "FDA-labeled"
    mechanism: str | None
    is_predictive: bool  # True = affects drug selection, False = just prognostic

class DrugSensitivity(BaseModel):
    drug: str
    evidence: str  # "in vitro" | "preclinical" | "clinical" | "FDA-labeled"
    ic50_nM: str | None

class LiteratureKnowledge(BaseModel):
    mutation_type: str  # "primary" | "secondary" | "both" | "unknown"
    is_prognostic_only: bool  # True if variant only prognostic, not predictive
    resistant_to: list[DrugResistance]
    sensitive_to: list[DrugSensitivity]
    clinical_significance: str
    evidence_level: str  # "FDA-approved" | "Phase 3" | ... | "None"
    references: list[str]  # PMIDs
    key_findings: list[str]
    confidence: float  # 0-1
```

### Paper Relevance Scoring

```python
# LLM service returns structured scoring
{
    "relevance_score": 0.85,  # 0-1
    "is_relevant": True,      # >= 0.6 threshold
    "signal_type": "resistance",  # resistance | sensitivity | mixed | prognostic | unclear
    "drugs_mentioned": ["osimertinib", "gefitinib"],
    "key_finding": "T790M confers resistance to first-gen EGFR TKIs",
    "confidence": 0.9
}
```

### Search Strategy

**Primary sources:**
- Semantic Scholar (better relevance ranking, citation data)
- PubMed (comprehensive, authoritative)

**Query construction:**
```
"{gene} {variant}" AND ({tumor_type} OR cancer OR oncology)
```

**Filtering:**
- Exclude reviews, meta-analyses (want primary data)
- Weight by recency (recent papers more relevant for resistance)
- Boost clinical studies over preclinical
- Prioritize papers with drug mentions

### Extraction Prompts

**Resistance/Sensitivity extraction:**
```
Given this abstract about {gene} {variant} in {tumor_type}:

{abstract}

Extract:
1. Any drugs mentioned with sensitivity or resistance findings
2. The evidence type (clinical trial, case report, preclinical)
3. Any resistance mechanisms described
4. Confidence level based on study design

Return structured JSON.
```

**Key finding extraction:**
```
Summarize the single most important finding from this paper
regarding {gene} {variant} treatment implications.
One sentence, be specific about outcomes if mentioned.
```

---

## Future Considerations

### Ensemble Literature Extraction (v0.6)

Run extraction across multiple LLMs and compare outputs:
- Claude Sonnet 4 (primary)
- GPT-4o-mini (secondary)
- Claude 3.5 Haiku (fast/cheap)

Compare outputs, flag disagreements, report confidence.

**Implementation:** Already have litellm integrated. Need ensemble voting logic.

### Pre-fetched Literature Cache

For top 50-100 cancer genes, pre-fetch papers weekly:
- Reduces API latency for common queries
- SQLite or Redis backend
- TTL: 7-30 days depending on gene activity

### Citation Network Analysis

Future feature: trace citation chains to find seminal papers and emerging work.
- Which papers does everyone cite? (foundational)
- Which papers cite our key papers? (emerging)

---

## API Rate Limits

| Source | Rate Limit | Strategy |
|--------|-----------|----------|
| Semantic Scholar | 100 req/5 min | Batch queries, cache results |
| PubMed E-utilities | 3 req/sec (with API key) | Sequential with delays |
| LLM APIs | Varies by provider | LitLLM handles fallbacks |

---

## Decision Log

| Decision | Rationale |
|----------|-----------|
| Build own extraction, not use generic tools | Cancer-specific logic is our IP |
| Use litellm for multi-provider support | ✅ Implemented - simplifies provider orchestration |
| Skip agentic frameworks for now | Adds complexity, not needed for literature search |
| Structured output over summaries | Enables downstream analysis and conflict detection |
| Predictive vs prognostic distinction | Critical for clinical relevance - prognostic markers don't guide drug selection |
| Signal type extraction | Enables automated resistance/sensitivity classification |

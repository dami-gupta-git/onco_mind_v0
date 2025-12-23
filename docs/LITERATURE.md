# Literature Synthesis Tools & Strategy

This document evaluates external tools for literature synthesis and outlines OncoMind's approach.

---

## Tool Evaluation

| Tool | What It Is | Recommendation |
|------|-----------|----------------|
| **LitLLM** | LLM API router (OpenAI, Anthropic, Cohere, etc.) | **Yes, for v0.6** - simplifies ensemble LLM orchestration |
| **hint-lab/pubmed-agent** | Agentic PubMed search with tool use | No - adds complexity, we have the APIs already |
| **llm-literature-review-tool** | Generic literature review automation | No - overkill for focused oncology use case |
| **LitSense** | NIH semantic search for PubMed | Maybe - supplementary data source |
| **PubMedSummarizer** | Basic LLM summarization of abstracts | No - trivial to build ourselves |

### LitLLM

**What it does:** Unified interface to 100+ LLM providers. One API call, swap models via config.

**Why it's useful for OncoMind:**
- v0.6 requires running extraction across multiple models (GPT-4o-mini, Claude Haiku, etc.)
- LitLLM handles the provider differences, rate limits, fallbacks
- Makes ensemble voting straightforward

**Example:**
```python
from litellm import completion

# Same code, different models
response1 = completion(model="gpt-4o-mini", messages=messages)
response2 = completion(model="claude-3-haiku-20240307", messages=messages)
response3 = completion(model="gemini/gemini-1.5-flash", messages=messages)
```

**When to add:** v0.6 (Ensemble LLM & Uncertainty)

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

### Current Implementation (v0.1)

1. **Search:** Semantic Scholar API + PubMed E-utilities
2. **Filter:** Cancer-relevant papers, recency weighting
3. **Extract:** LLM extracts resistance/sensitivity signals from abstracts
4. **Synthesize:** Summary with key findings and evidence gaps

### Target Implementation (v0.2)

Structured `LiteratureSynthesis` output:

```python
class KeyPaper(BaseModel):
    pmid: str
    title: str
    year: int
    relevance_score: float  # 0-1
    key_finding: str

class DrugSignal(BaseModel):
    drug: str
    signal_type: Literal["sensitivity", "resistance", "mixed"]
    evidence_level: str  # "clinical", "preclinical", "case_report"
    mechanism: Optional[str]
    supporting_pmids: list[str]

class LiteratureSynthesis(BaseModel):
    key_papers: list[KeyPaper]
    resistance_signals: list[DrugSignal]
    sensitivity_signals: list[DrugSignal]
    summary: str
    evidence_gaps: list[str]
    confidence: float  # 0-1, based on paper count and agreement
    cross_paper_agreement: float  # 0-1
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

Run extraction across multiple LLMs:
- GPT-4o-mini
- Claude 3 Haiku
- Gemini 1.5 Flash

Compare outputs, flag disagreements, report confidence.

**Implementation:** Use LitLLM for provider abstraction.

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
| Use LitLLM for ensemble (v0.6) | Simplifies multi-provider orchestration |
| Skip agentic frameworks for now | Adds complexity, not needed for literature search |
| Structured output over summaries | Enables downstream analysis and conflict detection |

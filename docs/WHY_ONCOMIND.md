# Why OncoMind?

The honest case for what this tool does and doesn't do.

---

## The Problem We're Solving

Variant interpretation tools exist. CIViC, OncoKB, CancerVar, PCGR, OpenCRAVAT — they're excellent at annotating well-characterized mutations.

But for less common variants, researchers still spend hours:
- Searching PubMed for case reports and functional studies
- Digging through ClinicalTrials.gov for relevant trials
- Reconciling conflicting information across databases
- Synthesizing what it all means for a specific tumor type

**OncoMind automates that workflow.** Aggregate, search, synthesize, cite — in seconds instead of hours.

---

## What We're Not

Let's be clear about what OncoMind isn't:

| Not This | Why |
|----------|-----|
| A database | We query databases, we don't replace them |
| A classifier | We don't assign tiers or pathogenicity — we surface evidence |
| A clinical decision tool | Research use only, not for patient care |
| Better than CIViC/OncoKB | We use them as sources, not compete with them |
| A complete solution | We focus on SNPs/indels now; fusions coming |

---

## What Makes OncoMind Different

### 1. We Tell You What We Don't Know

Every tool tells you what evidence exists. None tell you what's missing.

```
Evidence Gaps for KRAS G12D in Pancreatic Cancer:
- No FDA-approved targeted therapies (G12C inhibitors don't work here)
- Limited resistance data (fewer treatment options studied)
- No active Phase 3 trials for this specific variant
- Conflicting prognostic data across studies
```

**Why it matters:** Transparency about uncertainty builds trust. Knowing where to focus manual effort is actionable.

### 2. We Surface Conflicts, Not Hide Them

When CIViC says "sensitive" and a 2024 paper says "resistance emerging," most tools pick one. We show both.

```
Conflict Detected:
- CIViC (2019): Sensitive to Drug X (Level B evidence)
- PMID:38291034 (2024): Resistance reported in 15% of cases after 12 months
- Resolution: Initial sensitivity with acquired resistance risk
```

**Why it matters:** Oncologists are trained to be skeptical. A tool that hides disagreement loses trust.

### 3. We Output Context, Not Data Dumps

Legacy tools produce 50-page PDFs or 10,000-row CSVs. An LLM chokes on that. A human skims it.

OncoMind outputs dense, grounded context blocks:

```
EGFR T790M in NSCLC:
- Acquired resistance mutation (50-60% of EGFR TKI failures)
- Resistant to: erlotinib, gefitinib, afatinib (1st/2nd gen TKIs)
- Sensitive to: osimertinib (FDA-approved, AURA3 trial)
- Next resistance: C797S (10-15%), MET amp (5-10%), histologic transformation
- Sources: FDA label, CIViC:assertion:42, PMID:25923549
```

**Why it matters:** This is the "fuel" for downstream AI systems. We're not building the car; we're providing high-octane context.

### 4. No Source, No Claim

Every assertion links to a PMID, FDA label, or database entry. If we can't cite it, we don't say it.

```python
class AttributedClaim(BaseModel):
    claim: str           # "Resistant to osimertinib"
    source: str          # "PMID:29151359"
    source_url: str      # "https://pubmed.ncbi.nlm.nih.gov/29151359/"
    evidence_level: str  # "Clinical (Phase 2)"
```

**Why it matters:** "Click through to verify" beats "trust the AI."

### 5. We Focus on What Happens Next

Most annotators tell you what the variant *is*. We tell you what happens *next*.

```
Resistance Escape Routes for BRAF V600E on Vemurafenib:
├── NRAS mutations (20%) → Consider immunotherapy
├── MEK1/2 mutations (5-10%) → MEK inhibitor may help
├── BRAF amplification (10%) → Dose escalation or combo
└── PI3K/AKT activation (10-15%) → Clinical trial

Median time to resistance: 6-12 months
Next-line options: Add MEK inhibitor, switch to immunotherapy, clinical trial
```

**Why it matters:** Clinicians ask "what do I watch for?" and "what do I do when it fails?" Nobody answers this cleanly.

---

## The Trust Stack

We don't say "trust the AI." We provide a verification framework:

| Layer | What It Proves |
|-------|----------------|
| **Source Attribution** | Every claim has a receipt |
| **Conflict Detection** | We're not hiding disagreement |
| **Evidence Gaps** | We're honest about what we don't know |
| **Ensemble Agreement** | Multiple models, surface disagreement |
| **Benchmarks** | Quantitative accuracy on gold standard |

---

## Who Should Use This

### Primary Users

#### 1. Grad Students & Postdocs in Cancer Genomics Labs

**The Problem:**
Your advisor says "look into this variant." You spend days searching PubMed, checking CIViC, cross-referencing ClinVar, digging through cBioPortal... just to understand the current state of knowledge.

**OncoMind Solution:**
- Instant evidence summary across 14+ sources
- Explicit gap analysis: "Here's what's NOT known"
- Research implications: "Here's what you could study"
- Suggested experiments based on gap severity

**Value:** Saves weeks of literature review. Helps scope thesis projects with clear, defensible knowledge gaps.

---

#### 2. PIs Writing Grant Proposals

**The Problem:**
R01 proposals need to justify "why this research question matters." That means showing the gap — not just what's known, but what's missing and why it matters.

**OncoMind Solution:**
- Evidence gap detection with severity ratings (critical/significant/minor)
- Gap categories: functional, clinical, resistance mechanisms, preclinical data
- Suggested studies mapped to each gap
- Source attribution for every claim (PMIDs, database entries)

**Value:** Write stronger proposals with clearly articulated gaps. "We propose to address the CRITICAL gap in resistance mechanism understanding for [variant] — no published studies characterize acquired resistance despite 3 FDA-approved therapies."

---

#### 3. Pharma/Biotech Researchers in Early Discovery

**The Problem:**
"Should we pursue this variant as a drug target? What do we actually know? What's the white space?"

**OncoMind Solution:**
- Gap analysis shows research opportunities
- DepMap integration: gene essentiality (is it a good target?), drug sensitivity (what works?), model systems (can we test it?)
- Resistance mechanism extraction from literature
- Preclinical vs clinical evidence stratification

**Value:** Prioritize research portfolio based on evidence landscape. Identify white space before competitors. Make go/no-go decisions with full context.

---

#### 4. Bioinformaticians Supporting Multiple PIs

**The Problem:**
Constant ad-hoc variant lookups for different projects. Each PI wants "just a quick look" at their pet variant. You become the bottleneck.

**OncoMind Solution:**
- Standardized, comprehensive reports in seconds
- Batch processing for multiple variants
- JSON/CSV output for downstream analysis
- Consistent format across all queries

**Value:** Don't become the bottleneck. Deliver comprehensive variant briefs quickly. Spend your time on analysis, not data wrangling.

---

### Secondary Users

- **Clinical trial designers** — Identifying patient populations with unmet therapeutic needs. OncoMind surfaces variants where clinical evidence is lacking but biological rationale exists.
- **Science journalists & medical writers** — Understanding what's known vs. unknown about a variant. OncoMind provides source-attributed facts, not speculation.
- **Biotech BD teams** — Competitive intelligence on the research landscape. Which variants are over-studied? Where are the opportunities?

---

### Good Fit vs Not a Fit

**Good fit:**
- Researchers exploring variant evidence quickly
- Data scientists building oncology AI applications
- Developers who need LLM-ready cancer variant context
- Anyone doing the "2-hour PubMed rabbit hole" repeatedly

**Not a fit:**
- Clinical labs needing certified tools
- Anyone expecting a diagnostic or treatment recommendation
- Users who need fusion/amplification support (coming soon)

---

## The Competitive Landscape

| Tool | Strength | OncoMind's Angle |
|------|----------|------------------|
| **CIViC** | Gold-standard curation | We query it, add literature synthesis |
| **OncoKB** | Actionability levels | We surface it via VICC, add conflict detection |
| **VICC MetaKB** | Aggregates everything | We use it, add LLM intelligence layer |
| **PCGR** | Full reports, AMP/ASCO tiers | We output context blocks, not PDFs |
| **CancerVar** | 13M pre-computed variants | We do real-time literature synthesis |
| **OpenCRAVAT** | 300+ annotation modules | We focus on intelligence, not extensibility |

**We don't compete on annotation.** The databases do that well. We compete on:
- Intelligence (synthesis, not just aggregation)
- Trust (gaps, conflicts, attribution)
- Format (LLM-native context, not reports)

---

## The Honest Limitations

1. **SNPs and indels only** — Fusions, amplifications, CNVs coming in v0.3
2. **LLM extraction is imperfect** — We're building validation benchmarks
3. **Not for clinical use** — Research tool, not a medical device
4. **Rate limits on literature APIs** — Semantic Scholar and PubMed have limits
5. **Dependent on upstream databases** — We're only as good as CIViC, OncoKB, etc.

---

## The Vision

> "OncoMind: Grounded context for cancer variant reasoning. What happens next, with receipts."

We want to be the **context layer** for AI-assisted oncology. Not the database, not the classifier, not the report generator — the intelligence layer that turns raw evidence into LLM-ready, fact-checked, source-attributed knowledge.

For well-known variants, we give you what the databases have.
For everything else, we do the literature search you'd do manually — and show our work.

---

## Try It

```python
from oncomind import get_insight

result = await get_insight("EGFR T790M", tumor_type="NSCLC")

# What do the databases say?
print(result.kb.civic_assertions)
print(result.clinical.fda_approvals)

# What does the literature say?
print(result.literature.literature_knowledge)

# LLM-generated clinical narrative (when enabled)
if result.llm:
    print(result.llm.llm_summary)
```

See [README.md](../README.md) for installation and full API documentation.

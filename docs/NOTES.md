# OncoMind: Design Notes & Context

Internal documentation capturing design decisions, competitive analysis, and implementation details discussed during planning.

---

## The Core Problem

Variant interpretation tools (CIViC, OncoKB, CancerVar, OpenCRAVAT) are excellent for well-characterized mutations. But for less common variants, researchers spend hours manually searching PubMed, ClinicalTrials.gov, and Google Scholar for case reports, functional studies, and resistance mechanisms.

**The gap:** Raw data (APIs) → Ready-to-use intelligence (what a developer/clinician wants)

**OncoMind's position:** Automate the literature search and synthesis step, with full source attribution.

---

## Competitive Landscape

### Standard Aggregators
- **VICC MetaKB**: Gold standard for harmonizing CIViC, OncoKB, CGI
- **MyVariant.info**: High-performance API bundling ClinVar, COSMIC, gnomAD, dbSNP
- **Genome Nexus**: cBioPortal's backend, 12+ resources, cancer-specific

### Mature Annotation Pipelines
- **PCGR**: Extends VEP with AMP/ASCO/CAP interpretation, HTML reports
- **CancerVar**: ML-based Tier I-IV classification, 13M pre-computed variants
- **OpenCRAVAT**: 300+ annotation modules, extensible plugin system

### Where OncoMind Fits

| Tool | Strength | OncoMind's Angle |
|------|----------|------------------|
| VICC MetaKB | Aggregation | We use it as a source |
| PCGR | Full reports | We output LLM-ready context, not PDFs |
| CancerVar | Pre-computed tiers | We do real-time literature synthesis |
| OpenCRAVAT | Extensibility | We focus on intelligence layer |

**OncoMind doesn't compete on annotation.** The differentiator is LLM-powered synthesis + trust features.

---

## Key Differentiators

### 1. LLM-Native Output (Knowledge Headers)

Legacy tools output 50-page PDFs or 10K-row CSVs. An LLM will choke on that.

OncoMind outputs dense, grounded context blocks:

```
BRAF V600E in melanoma:
- Oncogenic driver (constitutively activates MAPK)
- FDA-approved: vemurafenib, dabrafenib + trametinib
- Resistance: typically 6-12 months via NRAS, MEK1/2, or BRAF amp
- When resistant: consider immunotherapy or clinical trial
- Sources: CIViC:12, OncoKB:Level1, PMID:22735384
```

This is the "high-octane fuel" for downstream AI systems.

### 2. Anti-Hallucination by Design

| Mechanism | Implementation |
|-----------|----------------|
| No source → no claim | If FDA API returns nothing, LLM is forbidden from mentioning "FDA-approved" |
| Conflict surfacing | If CIViC says sensitive and CGI says resistant, show both |
| Forced attribution | Every field requires `source_url` or it's rejected |
| Judge LLM check | Secondary model verifies claims against raw data |

### 3. Resistance Escape Routes

Most tools tell you what the variant *is*. OncoMind tells you what happens *next*.

```python
result = await oncomind.get_escape_routes(
    target="BRAF V600E",
    therapy="vemurafenib",
    tumor="melanoma"
)

# Returns:
# - bypass_pathways: ["PI3K/AKT", "NRAS reactivation"]
# - gatekeeper_mutations: ["BRAF amp", "MEK1 mutations"]
# - frequencies: {"PI3K activation": "10-15%", "NRAS": "20%"}
# - countermeasures: ["Add MEK inhibitor", "Consider immunotherapy"]
```

### 4. Uncertainty Quantification

Ensemble LLM extraction with agreement scores:

```python
class EnsembleResult(BaseModel):
    value: str                    # "resistant"
    agreement: float              # 0.67 (2/3 models agree)
    model_votes: dict[str, str]   # {"gpt-4o": "resistant", "claude": "resistant", "gemini": "sensitive"}
    confidence: str               # "moderate"
```

If models disagree, that's a data point. Surface it, don't hide it.

---

## Data Models

### Attributed Claims

Every claim must have a receipt:

```python
class AttributedClaim(BaseModel):
    claim: str
    source_type: Literal["FDA", "CIViC", "OncoKB", "PMID", "clinical_trial"]
    source_id: str           # "PMID:22735384" or "CIViC:assertion:12"
    source_url: str          # Clickable link
    retrieval_date: date
```

### Resistance Mechanisms

```python
class ResistanceMechanism(BaseModel):
    initial_therapy: str                    # "osimertinib"
    initial_target: str                     # "EGFR L858R"
    resistance_variant: str                 # "EGFR C797S"
    resistance_type: Literal["gatekeeper", "bypass", "downstream", "phenotypic"]
    frequency: str | None                   # "10-15% of progressors"
    median_time_to_emergence: str | None    # "12-18 months"
    next_line_options: list[str]            # ["combo EGFR inhibitor", "chemotherapy"]
    source_pmids: list[str]
```

### Resistance Classification

```python
class ResistanceProfile(BaseModel):
    variant: str
    typically_acquired: bool      # T790M = usually acquired after TKI
    typically_intrinsic: bool     # KRAS G12C = usually present at diagnosis
    emergence_context: str        # "After EGFR TKI therapy"
    median_time_to_emergence: str | None
    source_pmids: list[str]
```

### V2P (Variant to Phenotype)

```python
class PhenotypeInference(BaseModel):
    msi_status: Literal["MSI-H", "MSS", "unknown"]
    hrd_status: Literal["HRD", "HRP", "unknown"]
    hypermutated: bool
    actionable_phenotypes: list[str]  # ["immunotherapy_responsive", "PARP_sensitive"]
    inferred_from: list[str]          # ["MLH1 p.R100*", "BRCA2 loss"]
    confidence: float
```

Rule-based phenotypes (MSI-H from MMR genes, HRD from BRCA1/2) are low-hanging fruit.

### Synthetic Lethality

```python
class SyntheticLethality(BaseModel):
    context: str                # "BRCA1/2 loss"
    vulnerability: str          # "PARP dependency"
    drug_class: str             # "PARP inhibitors"
    drugs: list[str]            # ["olaparib", "rucaparib", "niraparib"]
    evidence_level: str         # "FDA-approved"
    tumor_types: list[str]      # ["ovarian", "breast", "prostate", "pancreatic"]
```

Well-established pairs:
- BRCA loss → PARP inhibitors
- MSI-H → checkpoint inhibitors
- ATM loss → PARP/ATR inhibitors (emerging)
- MTAP loss → MAT2A inhibitors (emerging)

### Trial Eligibility

```python
class TrialEligibility(BaseModel):
    nct_id: str
    required_variants: list[str]      # ["EGFR exon 19 del", "EGFR L858R"]
    excluded_variants: list[str]      # ["EGFR T790M"]
    variant_logic: str                # "any_of" or "all_of"
    other_criteria: list[str]         # ["No prior EGFR TKI", "ECOG 0-1"]
```

### Knowledge Header

```python
class KnowledgeHeader(BaseModel):
    """LLM-optimized context block for a variant."""

    variant_id: str
    tumor_context: str
    headline: str                     # One-sentence summary
    mechanism: str                    # "Constitutive MAPK pathway activation"
    sensitivity: list[str]
    resistance_risk: str
    resistance_mechanisms: list[str]
    next_line_options: list[str]
    phenotype: str | None
    evidence_level: str
    conflicts: list[str] | None
    sources: list[str]

    def to_prompt_block(self) -> str:
        """Format for injection into any LLM prompt."""
        ...
```

---

## Validation Framework

### Gold Standard Benchmark

50 variants with known correct answers:

| Variant | Tumor | Expected | Source |
|---------|-------|----------|--------|
| EGFR T790M | NSCLC | Resistant to erlotinib/gefitinib, sensitive to osimertinib | CIViC, FDA |
| EGFR C797S | NSCLC | Resistant to osimertinib | CIViC |
| BRAF V600E | Melanoma | Sensitive to vemurafenib/dabrafenib | FDA label |
| KRAS G12C | NSCLC | Sensitive to sotorasib/adagrasib | FDA label |
| ALK fusion | NSCLC | Sensitive to crizotinib/alectinib/lorlatinib | FDA label |

### Metrics

| Metric | Target | Why |
|--------|--------|-----|
| Accuracy vs. gold standard | >85% | Core extraction works |
| Citation coverage | 100% | Every claim has a receipt |
| Cohen's Kappa vs. expert | >0.75 | Clinical-grade agreement |
| Hallucination rate | <5% | LLM isn't making things up |

### benchmark.py

```python
"""
Validation suite for OncoMind.
Run: python benchmark.py --dataset data/gold_standard.csv
"""

async def run_benchmark(dataset_path: str) -> BenchmarkResult:
    """
    Loads ground truth, runs OncoMind, compares results.

    Outputs:
    - Accuracy on sensitivity/resistance calls
    - Citation coverage (% claims with sources)
    - Cohen's Kappa vs. expert annotations
    - Confusion matrix
    """
    ...
```

### Faithfulness Checking

```python
class FaithfulnessCheck(BaseModel):
    claim: str                    # "Resistant to osimertinib"
    supported: bool               # False
    supporting_source: str | None # None
    verdict: Literal["grounded", "hallucinated", "uncertain"]

async def check_faithfulness(
    llm_output: KnowledgeHeader,
    raw_evidence: EvidencePanel
) -> list[FaithfulnessCheck]:
    """Use judge LLM to verify each claim against sources."""
    ...
```

---

## Ensemble LLM Implementation

### When to Use Ensemble

| Use Case | Ensemble Worth It? |
|----------|-------------------|
| Resistance/sensitivity extraction | Yes — high stakes |
| Paper relevance scoring | Maybe — disagreement flags edge cases |
| Narrative generation | No — no "correct" answer |
| Trial eligibility parsing | Yes — errors matter |

### Implementation

```python
async def extract_with_confidence(paper, variant):
    results = await asyncio.gather(
        extract_with_model("gpt-4o-mini", paper, variant),
        extract_with_model("claude-haiku", paper, variant),
    )

    if results[0].signal == results[1].signal:
        return {"signal": results[0].signal, "confidence": "high"}
    else:
        return {"signal": "uncertain", "confidence": "low", "disagreement": results}
```

### Cost/Latency Tradeoff

| Setup | Latency | Cost per variant |
|-------|---------|------------------|
| Single GPT-4o-mini | ~1s | ~$0.01 |
| GPT-4o + Claude | ~3s (parallel) | ~$0.10 |
| 3x prompts same model | ~3s | ~$0.03 |

Use ensemble for high-stakes extractions only.

---

## Structural Variants (v0.3)

### The Challenge

- SNPs: `BRAF V600E` → straightforward protein notation
- Fusions: `EML4-ALK` vs `ALK-EML4` vs `ALK fusion` — many representations
- Amplifications: `MET amp` vs `MET copy number gain` vs `MET CN ≥6`

### Pragmatic Approach

1. Accept simple input formats: `ALK fusion`, `EML4-ALK`, `MET amp`
2. Map to canonical representations in each database's format
3. Don't try to solve general SV normalization

### Database Coverage

| Source | Fusions | Amps/CNVs |
|--------|---------|-----------|
| CIViC | Yes (well-curated) | Limited |
| VICC MetaKB | Yes | Yes |
| OncoKB | Yes (excellent) | Yes |
| CGI | Yes | Yes |
| MyVariant.info | No | No |

MyVariant.info's functional scores don't apply to fusions anyway.

### NSCLC Priority List

- ALK fusions (EML4-ALK most common)
- ROS1 fusions (CD74-ROS1, etc.)
- RET fusions (KIF5B-RET, CCDC6-RET)
- NTRK1/2/3 fusions
- MET exon 14 skipping
- MET amplification
- ERBB2 (HER2) amplification/mutations

---

## Pathway Context

### What It Means

```python
result = await oncomind.get_pathway_context("EGFR", "L858R", tumor="NSCLC")

result.pathway              # "MAPK/ERK signaling"
result.position             # "receptor"
result.common_bypass_nodes  # ["MET", "HER2", "KRAS"]
result.known_escape_mechanisms # [
                            #   {"trigger": "EGFR TKI", "bypass": "MET amp", "frequency": "5-10%"},
                            # ]
result.co_targeting_options # ["EGFR + MET inhibitor"]
```

### Data Sources

| Data | Source | Availability |
|------|--------|--------------|
| Gene → Pathway | KEGG, Reactome | Structured, free |
| Bypass mechanisms | Literature | Needs LLM extraction |
| Frequency estimates | Clinical studies | Scattered |
| Co-targeting strategies | Trial data | Unstructured |

### Reality Check

The insight isn't "EGFR is in MAPK" — any database tells you that.

The value is: "When you block EGFR, here's what happens next, with frequencies and citations."

That's curated knowledge from literature, not pathway databases.

---

## Market Reality

### Who Would Use OncoMind?

| User | Likelihood | Why |
|------|------------|-----|
| AI/ML developers building oncology apps | High | Don't want to learn domain APIs |
| Researchers prototyping | High | Want to test ideas quickly |
| Startups without bioinformatics staff | Medium | Need quick integration |
| Hackathon participants | High | Need fast results |
| Clinical labs | Low | Need validated, certified tools |
| Large pharma | Low | Have internal infrastructure |

### The Honest Pitch

> "OncoMind saves you 2 weeks of API integration so you can focus on your actual problem. For well-known variants, we give you what the databases have. For everything else, we do the literature search you'd do manually — and show you the sources."

---

## Technical Notes

### Streamlit Literature Tab

Already have the components:
- `semantic_scholar.py` — literature search
- `pubmed.py` — fallback search
- `service.py:208-331` — paper relevance scoring
- `service.py:333-460` — knowledge extraction

Just need UI wiring:

```python
with st.expander("Literature & Trials"):
    papers = await semantic_scholar.search_variant_literature(gene, variant, tumor)
    scored = await llm_service.score_paper_relevance(papers, variant)

    for paper in sorted(scored, key=lambda x: x.relevance_score, reverse=True):
        st.markdown(f"**{paper.title}** ({paper.relevance_score:.0%})")
        st.markdown(f"[PubMed]({paper.url})")
```

### Embeddings (Future)

Possible use cases:
- Variant similarity search ("find variants with similar evidence")
- Literature corpus embedding (semantic search beyond keywords)
- Cross-variant reasoning ("this VUS is similar to known pathogenic variants")

**Don't build until there's a concrete use case.** Structured queries often beat vector search for structured data.

---

## Future Differentiators (Ideas)

### 1. Evidence Gap Identification

**The idea:** "Here's what we DON'T know about this variant."

No one does this well. Current tools tell you what evidence exists. They don't tell you what's missing.

```python
class EvidenceGaps(BaseModel):
    variant: str
    tumor_type: str

    # What's missing
    no_clinical_trials: bool
    no_resistance_data: bool
    no_functional_studies: bool
    limited_tumor_specific_data: bool
    conflicting_evidence: bool

    # Specific gaps
    gaps: list[str]  # ["No Phase 3 data", "Resistance mechanisms unknown", "No data in this tumor type"]

    # What would help
    suggested_searches: list[str]  # ["Check for related variants in same codon", "Look for pathway-level evidence"]
```

**Value:** Transparency about uncertainty. Researchers know where to focus their manual effort.

**Implementation:** Compare expected evidence (based on gene importance, variant type) against actual evidence retrieved.

### 2. Multi-Variant Reasoning

**The idea:** "I have these 3 variants together, what does that mean?"

Current tools annotate variants independently. Real tumors have multiple variants that interact.

```python
async def analyze_variant_profile(
    variants: list[str],  # ["EGFR L858R", "TP53 R248W", "MET amp"]
    tumor_type: str
) -> ProfileAnalysis:
    ...

class ProfileAnalysis(BaseModel):
    variants: list[str]
    tumor_type: str

    # Individual annotations
    individual_evidence: dict[str, EvidencePanel]

    # Combined analysis
    interactions: list[VariantInteraction]
    prognostic_impact: str  # "TP53 co-mutation associated with worse prognosis in EGFR+ NSCLC"
    therapeutic_implications: list[str]  # ["Consider osimertinib", "TP53 may reduce TKI duration"]
    resistance_risk: str  # "MET amp suggests potential primary resistance"

    # Literature specific to combination
    combination_papers: list[PubMedEvidence]

    sources: list[str]
```

**Concrete examples:**
- EGFR L858R + TP53 → worse prognosis, shorter TKI response
- EGFR L858R + MET amp → potential primary resistance, consider combo
- BRAF V600E + PIK3CA → may reduce BRAF inhibitor efficacy

**Implementation:**
1. Annotate each variant individually
2. Search literature for co-occurrence patterns
3. LLM synthesis of combined implications
4. Surface known interaction patterns (can curate top combinations)

### 3. Clinical Question Answering

**The idea:** Not a report, but an answer to a specific question.

Instead of: "Here's everything about BRAF V600E"
Answer: "Should I try pembrolizumab given BRAF V600E and melanoma?"

```python
async def answer_clinical_question(
    question: str,  # "Should I try drug X given variant Y and tumor Z?"
    context: EvidencePanel | None = None
) -> ClinicalAnswer:
    ...

class ClinicalAnswer(BaseModel):
    question: str

    # Direct answer
    answer: str  # "Pembrolizumab is not first-line for BRAF V600E melanoma. BRAF/MEK inhibitors are preferred."
    confidence: Literal["high", "moderate", "low"]

    # Supporting evidence
    supporting_evidence: list[AttributedClaim]
    contradicting_evidence: list[AttributedClaim]

    # Caveats
    caveats: list[str]  # ["Response to immunotherapy may be considered after BRAF inhibitor progression"]

    # What we don't know
    evidence_gaps: list[str]

    # Disclaimer
    disclaimer: str = "This is not medical advice. Discuss with treating physician."
```

**Question types:**
- "Should I try [drug] given [variant] and [tumor]?"
- "What's the expected response duration for [variant] on [drug]?"
- "What should I monitor for resistance?"
- "Is this patient eligible for [trial NCT#]?"

**Implementation:**
1. Parse question to extract drug, variant, tumor
2. Retrieve relevant evidence
3. LLM generates answer grounded in evidence
4. Require citations for every claim
5. Surface contradictions and gaps

**Risk:** This is close to clinical decision support. Be very clear about disclaimers and limitations.

### 4. Dynamic Evidence Synthesis

**The idea:** Not just "here's what we know" — reasoning across sources in real-time.

Static annotation: "CIViC says X, OncoKB says Y, here's a list."

Dynamic synthesis: "CIViC and OncoKB agree on sensitivity. However, three recent papers (2023-2024) report emerging resistance. The FDA label doesn't reflect this yet. Weight: moderate confidence in sensitivity, with caveats."

```python
class SynthesizedEvidence(BaseModel):
    variant: str
    tumor_type: str
    query_context: str  # What we're trying to answer

    # Synthesis result
    consensus: str | None  # None if no consensus
    confidence: Literal["high", "moderate", "low", "conflicting"]

    # Evidence weighting
    weighted_sources: list[WeightedSource]

    # Conflict analysis
    conflicts: list[EvidenceConflict]
    conflict_resolution: str | None  # "Recent literature suggests resistance emerging despite FDA approval"

    # Temporal analysis
    evidence_trend: str | None  # "Older sources say sensitive, newer sources report resistance"
    most_recent_evidence: str

    # Final synthesis
    synthesis_narrative: str  # Human-readable paragraph weighing all evidence
    key_caveats: list[str]

class WeightedSource(BaseModel):
    source: str
    claim: str
    evidence_level: str
    publication_date: date | None
    weight: float  # Based on recency, evidence level, specificity
    rationale: str  # Why this weight

class EvidenceConflict(BaseModel):
    source_a: str
    claim_a: str
    source_b: str
    claim_b: str
    conflict_type: Literal["direct_contradiction", "different_context", "different_evidence_level", "temporal"]
    resolution_hint: str | None
```

**Key insight:** Evidence has a temporal dimension. A 2018 FDA approval doesn't account for 2024 resistance reports. The synthesis should reflect this.

**Implementation:**
1. Retrieve all evidence with timestamps
2. Detect conflicts across sources
3. Weight by recency, evidence level, tumor specificity
4. LLM generates synthesis narrative
5. Surface conflicts explicitly, don't hide them

### 5. Contextual Question-Answering

**The idea:** "What does this variant mean for a patient who already failed X therapy?"

This is not static lookup. It's contextual reasoning with patient history.

```python
async def contextual_query(
    variant: str,
    tumor_type: str,
    clinical_context: ClinicalContext
) -> ContextualAnswer:
    ...

class ClinicalContext(BaseModel):
    prior_therapies: list[str]  # ["erlotinib", "osimertinib"]
    therapy_outcomes: dict[str, str]  # {"erlotinib": "progressed after 8 months"}
    current_status: str  # "progressing on osimertinib"
    other_variants: list[str]  # Co-occurring variants
    patient_factors: dict[str, str] | None  # {"ECOG": "1", "brain_mets": "no"}

class ContextualAnswer(BaseModel):
    question_interpreted: str  # "What options after osimertinib failure in EGFR L858R NSCLC?"

    # Context-aware answer
    answer: str
    confidence: Literal["high", "moderate", "low"]

    # Why this answer given the context
    reasoning: list[ReasoningStep]

    # Options ranked by context
    therapeutic_options: list[TherapeuticOption]

    # Trials specifically for this context
    relevant_trials: list[ClinicalTrialEvidence]  # Trials for post-osimertinib

    # What resistance mechanism to look for
    suggested_testing: list[str]  # ["Check for C797S", "Consider MET FISH"]

    sources: list[str]

class TherapeuticOption(BaseModel):
    option: str  # "Platinum-based chemotherapy"
    rationale: str  # "Standard after TKI failure"
    evidence_level: str
    context_fit: str  # "Appropriate given prior TKI exposure"
    caveats: list[str]
```

**Example query flow:**
- Input: "EGFR L858R, failed erlotinib (8mo), now failed osimertinib (14mo)"
- Reasoning:
  1. T790M likely emerged after erlotinib (explains osimertinib use)
  2. Osimertinib failure suggests C797S or bypass (MET, HER2)
  3. Check resistance mechanism to guide next therapy
  4. If C797S: consider combo or chemo
  5. If MET amp: consider tepotinib/capmatinib combo trials
- Output: Ranked options with rationale specific to this treatment history

**This is the holy grail of precision oncology tooling.** Moving from "here's everything about this variant" to "here's what matters for THIS patient's situation."

---

## Open Questions

1. **Structural variant normalization** — How much effort? Pragmatic approach vs. proper HGVS?
2. **LLM synthesis quality** — Need to benchmark against known variants
3. **Handling contradictions** — Paper A says sensitive, Paper B says resistant. Current approach: surface both.
4. **Pathway data source** — Curate manually for top variants, or extract from literature?
5. **Multi-variant interactions** — Curate known combinations, or rely on literature search?
6. **Clinical Q&A liability** — How to frame to avoid being perceived as medical advice?

---

## References

- AMP/ASCO/CAP 2017 Guidelines for somatic variant classification
- VICC Meta-Knowledgebase: https://search.cancervariants.org/
- CIViC: https://civicdb.org/
- OncoKB: https://www.oncokb.org/

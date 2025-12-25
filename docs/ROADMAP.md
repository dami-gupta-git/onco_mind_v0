# OncoMind Roadmap

This document outlines the development roadmap for OncoMind, organized by release phase.

## Vision

OncoMind aims to be the **grounded context layer** for AI-assisted oncology reasoning. Not a database, not a classifier — infrastructure that provides LLM-ready, fact-checked, source-attributed evidence about cancer variants.

### Core Principles

1. **No source, no claim** — Every assertion must link to a PMID, FDA label, or database entry
2. **Surface uncertainty** — Flag conflicts and disagreements, don't hide them
3. **What happens next** — Focus on resistance mechanisms and clinical trajectory, not just variant classification
4. **LLM-native output** — Designed for injection into downstream AI systems

---

## Current State (v0.1)

**Status:** Working

- SNP and small indel annotation
- Multi-source evidence aggregation (CIViC, ClinVar, COSMIC, VICC, CGI, FDA, ClinicalTrials.gov)
- Functional predictions (AlphaMissense, CADD, PolyPhen2, SIFT, gnomAD)
- Literature search (Semantic Scholar, PubMed)
- Basic LLM synthesis for literature analysis
- Pydantic `Result` output model (wraps `Evidence` + optional `LLMInsight`)
- CLI and Python API

---

## v0.2: Trust Foundation

**Goal:** Establish the trust layer that makes OncoMind credible.

### Source Attribution
- [ ] Consistent attribution across all evidence types
- [ ] Every extracted field includes `source`, `source_id`, `source_url`
- [ ] `AttributedClaim` model for LLM-generated content
- [ ] Retrieval timestamps on all evidence

### Conflict Detection
- [ ] Cross-source comparison for therapeutic predictions
- [ ] `result.evidence.meta.conflicts` field populated automatically
- [ ] Flag discrepancies between CIViC, VICC/OncoKB, CGI
- [ ] Surface evidence level disagreements (e.g., Tier I vs Tier II)

### Structured Literature Synthesis
- [ ] `LiteratureSynthesis` model with key papers, resistance/sensitivity signals
- [ ] `KeyPaper` model: PMID, title, year, relevance score, key finding
- [ ] `DrugSignal` model: drug, signal type, evidence level, mechanism, supporting PMIDs
- [ ] Summary paragraph with evidence gaps
- [ ] Confidence and cross-paper agreement metrics

### Evidence Gap Identification
- [ ] `result.evidence.meta.evidence_gaps` field populated automatically
- [ ] Detect missing evidence: no trials, no resistance data, no functional studies
- [ ] Detect limited tumor-specific data
- [ ] Surface conflicting evidence as a gap type
- [ ] Suggested next steps (where to focus manual effort)

### Deliverables
- [ ] Updated Pydantic models with attribution
- [ ] Conflict detection logic in `EvidenceBuilder`
- [ ] Evidence gap detection in `EvidenceBuilder`
- [ ] `LiteratureSynthesis` output in `result.literature`
- [ ] Unit tests for conflict and gap scenarios

---

## v0.3: Structural Variants

**Goal:** Support the variant types required for real clinical utility.

### Fusion Support
- [ ] Input parsing for fusion notation (`EML4-ALK`, `ALK fusion`)
- [ ] CIViC fusion queries
- [ ] VICC MetaKB fusion queries
- [ ] CGI fusion biomarkers

### Amplification Support
- [ ] Input parsing for amplification notation (`MET amp`, `ERBB2 amplification`)
- [ ] Database queries for copy number events

### Target Coverage (NSCLC focus)
- [ ] ALK fusions
- [ ] ROS1 fusions
- [ ] RET fusions
- [ ] NTRK1/2/3 fusions
- [ ] MET exon 14 skipping
- [ ] MET amplification
- [ ] ERBB2 (HER2) amplification

### Deliverables
- [ ] Extended `parse_variant_input()` for structural variants
- [ ] Updated API clients for fusion/amp queries
- [ ] Integration tests with known fusions

---

## v0.4: Validation & Benchmarking

**Goal:** Prove the LLM extraction works with quantitative metrics.

### Gold Standard Dataset
- [ ] Curate 50 variants with known correct answers
- [ ] Include sensitivity, resistance, evidence levels
- [ ] Cover well-characterized and less common variants
- [ ] Document sources for ground truth

### Benchmark Script
- [ ] `benchmark.py` in repo root
- [ ] Loads gold standard, runs OncoMind, compares results
- [ ] Outputs confusion matrix
- [ ] Reports accuracy, precision, recall

### Metrics
- [ ] Accuracy vs. gold standard (target: >85%)
- [ ] Citation coverage (target: 100%)
- [ ] Cohen's Kappa vs. expert annotations (target: >0.75)

### Faithfulness Checking
- [ ] Judge LLM to verify claims against raw evidence
- [ ] Hallucination rate tracking
- [ ] Flag unsupported claims

### Deliverables
- [ ] `data/gold_standard.csv` with curated variants
- [ ] `benchmark.py` with reproducible validation
- [ ] Metrics documented in README

---

## v0.5: Visualizations

**Goal:** Make evidence tangible and explorable.

### Evidence Visualizations
- [ ] **Evidence gap radar**: Spider chart showing evidence completeness across categories
- [ ] **Conflict heatmap**: Which sources agree/disagree for a variant
- [ ] **Resistance trajectory**: Timeline showing typical resistance evolution (e.g., L858R → T790M → C797S)
- [ ] **Evidence network**: Graph showing variant → drug → outcome relationships

### Embedding Visualizations
- [ ] **Variant landscape map**: t-SNE/UMAP of variant embeddings, colored by drug response
- [ ] Cluster similar variants by evidence profile
- [ ] Interactive tooltips with variant details

### Streamlit Integration
- [ ] Visualization components in Streamlit app
- [ ] Interactive exploration of evidence
- [ ] Export as static images for reports

### Deliverables
- [ ] Visualization module (`oncomind.viz`)
- [ ] Streamlit visualization tab
- [ ] Example outputs in documentation

---

## v0.6: Ensemble LLM & Uncertainty

**Goal:** Quantify uncertainty in LLM extractions.

### Multi-Model Extraction
- [ ] Run extraction across 2-3 models (GPT-4o-mini, Claude Haiku, etc.)
- [ ] Compare outputs field-by-field
- [ ] Require identical Pydantic schema for comparison

### Uncertainty Quantification
- [ ] Agreement score (0-1) for each extraction
- [ ] `confidence` field on LLM-generated content
- [ ] Flag low-agreement extractions for review

### Consensus Logic
- [ ] Majority voting for categorical fields
- [ ] Surface disagreements rather than forcing consensus
- [ ] `EnsembleResult` model with `model_votes` field

### Deliverables
- [ ] Ensemble extraction in `LLMService`
- [ ] Agreement metrics in output
- [ ] Cost/latency analysis documented

---

## v0.7: Clinical Intelligence

**Goal:** Move from "what the variant is" to "what happens next."

### Trial Eligibility Parsing
- [ ] LLM extraction of inclusion/exclusion criteria
- [ ] Structured eligibility output (`required_variants`, `excluded_variants`)
- [ ] Match variants against parsed criteria
- [ ] Flag variant-specific vs. gene-level trials

### Resistance Escape Routes
- [ ] Curated resistance mechanisms for top variants
- [ ] `ResistanceMechanism` model with drugs affected, bypass pathways
- [ ] Frequency estimates where available
- [ ] Next-line therapy options

### Resistance Classification
- [ ] `typically_acquired` vs `typically_intrinsic` flags
- [ ] Emergence context (e.g., "after EGFR TKI therapy")
- [ ] Gatekeeper mutation identification

### V2P (Variant to Phenotype)
- [ ] MSI-H inference from MMR gene mutations
- [ ] HRD inference from HR gene panel
- [ ] Hypermutated phenotype from POLE/POLD1
- [ ] `PhenotypeInference` model

### Deliverables
- [ ] Trial eligibility extraction
- [ ] Resistance data model and initial curation
- [ ] Phenotype inference logic
- [ ] Integration tests

---

## v0.8: LLM Infrastructure

**Goal:** Position OncoMind as the context layer for AI oncology applications.

### Knowledge Header Export
- [ ] `result.to_knowledge_header()` method
- [ ] Dense, grounded context block format
- [ ] Optimized for token efficiency
- [ ] All claims attributed inline

### Prompt-Ready Output
- [ ] `result.to_prompt_context()` for system prompt injection
- [ ] Structured JSON for agent tool use
- [ ] Markdown format for chat interfaces

### Anti-Hallucination Guarantees
- [ ] Constraint-based generation (no FDA mention if no FDA data)
- [ ] Required fields validation
- [ ] Source coverage assertion

### Deliverables
- [ ] Export methods on `Result`
- [ ] Documentation for downstream integration
- [ ] Example agent integration

---

## Future Considerations

These are tracked but not scheduled:

### Gene Summaries
- Gene-level context beyond variant annotation
- Background, function, pathway, hotspot variants
- Variant frequency by tumor type (from COSMIC/cBioPortal)
- Approved therapies, resistance patterns, common co-mutations
- Options: pre-curated for top 50 genes, or dynamic generation with LLM

### Clinical Question Answering
- Answer specific questions, not generate reports
- "Should I try drug X given variant Y and tumor Z?"
- "What's the expected response duration?"
- "What should I monitor for resistance?"
- Grounded in evidence with citations, not open-ended generation
- Surface contradictions and evidence gaps in answer

### Trial Matching Engine
- Reasoning across eligibility criteria, not just keyword matching
- LLM-parsed eligibility: required variants, excluded variants, prior therapy requirements
- Match patient profile against parsed criteria
- Output: match confidence, rationale, potential issues
- Scope: variant/therapy matching (avoid full PHI-based patient matching)

### Multi-Variant Reasoning
- Input: multiple variants together (e.g., EGFR L858R + TP53 + MET amp)
- Co-occurrence impact on prognosis/resistance
- Literature search for variant combinations
- Combined therapeutic implications

### Dynamic Evidence Synthesis
- Real-time reasoning across sources with temporal weighting
- Recent papers vs old FDA labels — recency matters
- Conflict analysis and resolution hints
- Synthesis narrative with explicit caveats

### Contextual Question-Answering
- "What does this variant mean for a patient who already failed X therapy?"
- Reasoning with treatment history context
- Ranked therapeutic options with rationale specific to prior therapy
- Suggested resistance testing based on clinical trajectory

### Embeddings & Vector Search
- Embed `Result` for similarity search
- Use case: "Find variants with similar evidence profiles"
- Depends on: concrete use case emerging

### Pathway Context
- Gene → pathway membership (KEGG, Reactome)
- Bypass pathway identification
- Cross-pathway resistance tracking

### Pre-fetched Literature Cache
- Daily/weekly cache of papers for top cancer genes
- Reduce API latency for common queries
- SQLite or Redis backend

---

## Contributing

See individual phase issues for specific tasks. PRs welcome for any roadmap item.

## Timeline

No fixed dates. This is a portfolio/research project with incremental releases as features are completed and validated.

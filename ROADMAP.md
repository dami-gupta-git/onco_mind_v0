# OncoMind Roadmap

## Vision

OncoMind is evolving from a variant annotation tool into a **research intelligence platform** for cancer genomics. The core insight: for well-characterized variants, the answers are in databases. For everything else, the interesting question is *what don't we know yet?*

Our differentiator is **evidence gap detection** — explicitly identifying what's missing, where sources conflict, and what research questions remain open.

---

## Current State (v0.1)

### Evidence Aggregation
- CIViC, VICC MetaKB, ClinVar, COSMIC, CGI, FDA labels
- ClinicalTrials.gov integration
- cBioPortal co-mutation patterns
- AlphaMissense, CADD, PolyPhen2, gnomAD

### Research Intelligence
- Evidence gap detection (no functional data, limited clinical evidence, etc.)
- Evidence quality assessment (comprehensive/moderate/limited/minimal)
- Source conflict detection
- LLM-powered hypothesis generation
- Cited sources for all claims

### Literature
- Semantic Scholar + PubMed search
- LLM relevance scoring
- Resistance/sensitivity signal extraction

---

## Near-Term (Q1 2025)

### Reactome Integration
**Goal:** Add pathway context to variant analysis

- Query pathways containing the mutated gene
- Show upstream/downstream pathway members
- Identify pathway-level therapeutic opportunities
- Surface when variants affect pathway bottlenecks

**Research value:** "This variant affects a gene in the PI3K-AKT pathway, which has 3 FDA-approved inhibitors targeting other nodes."

### DepMap/CCLE Integration
**Goal:** Add cell line dependency and drug sensitivity data

- Gene dependency scores across cancer cell lines
- Drug sensitivity (PRISM) for cell lines with the mutation
- Identify cell line models for functional studies
- Co-essential gene patterns

**Research value:** "Cell lines with this mutation show 2.3x higher sensitivity to MEK inhibitors (PRISM data, n=47 lines)."

### Enhanced Hypothesis Generation
- Cross-reference co-mutation patterns with drug sensitivity
- Suggest synthetic lethality opportunities based on DepMap
- Generate testable hypotheses with specific experimental designs

---

## Mid-Term (Q2-Q3 2025)

### UniProt Integration
- Protein domain annotations
- Post-translational modification sites
- Structural features affected by variants
- Protein-protein interaction context

### bioRxiv Integration
- Preprint search for cutting-edge findings
- Flag preprints not yet peer-reviewed
- Earlier access to emerging resistance mechanisms

### Structural Variant Support
- Fusions (ALK, ROS1, RET, NTRK)
- Gene amplifications (ERBB2, MET, MYC)
- Copy number variants
- Rearrangements

### NIH RePORTER Integration
- Active grants studying the gene/variant
- Research funding landscape
- Identify collaboration opportunities

---

## Long-Term (2025+)

### Research Tracking
- Watch list for variants of interest
- Alerts when new evidence appears (papers, trials, FDA approvals)
- Longitudinal evidence quality tracking

### Batch Analysis Features
- Cohort-level analysis
- Cross-variant pattern detection
- Research portfolio prioritization

### API & Integrations
- RESTful API for programmatic access
- LIMS/EHR integration hooks
- Jupyter notebook widgets

---

## Design Principles

1. **Gaps over facts** — Prioritize surfacing what's unknown
2. **Sources required** — Every claim must link to evidence
3. **Conflicts visible** — Never hide disagreements between sources
4. **Research-first** — Optimize for research workflows, not clinical decisions
5. **Hypothesis-generating** — Output should suggest next experiments

---

## Contributing

We welcome contributions! Priority areas:
- New data source integrations (see Near-Term goals)
- Evidence gap detection improvements
- Hypothesis generation algorithms
- Structural variant support

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

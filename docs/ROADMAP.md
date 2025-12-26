# OncoMind Roadmap

## Vision

OncoMind is evolving from a variant annotation tool into a **research intelligence platform** for cancer genomics. The core insight: for well-characterized variants, the answers are in databases. For everything else, the interesting question is *what don't we know yet?*

Our differentiator is **evidence gap detection** — explicitly identifying what's missing, where sources conflict, and what research questions remain open.

---

## Current State (v0.1) ✅

### Evidence Aggregation
- CIViC, VICC MetaKB, ClinVar, COSMIC, CGI, FDA labels
- ClinicalTrials.gov integration
- cBioPortal co-mutation patterns (prevalence, co-occurring, mutually exclusive)
- DepMap/CCLE integration (gene essentiality, PRISM drug sensitivity, cell line models)
- AlphaMissense, CADD, PolyPhen2, SIFT, gnomAD (via MyVariant.info)
- OncoTree tumor type normalization
- VEP fallback for variant normalization

### Research Intelligence
- Evidence gap detection (`gap_detector.py`)
  - Gap categories: functional, clinical, tumor_type, drug_response, resistance, preclinical, prevalence, prognostic
  - Severity levels: critical, significant, minor
  - Suggested studies and addressable data sources per gap
- Evidence quality assessment (comprehensive/moderate/limited/minimal)
- Research priority scoring (high/medium/low based on gene importance + gaps)
- Source conflict detection (in LLM output)
- LLM-powered hypothesis generation with evidence tags
- Cited sources for all claims (enforced in prompts)

### Therapeutic Evidence Model
- `TherapeuticEvidence` model with:
  - All evidence levels (FDA → preclinical → computational)
  - Response type (sensitivity/resistance)
  - Mechanism of action
  - Quantitative data (IC50, response rates)
  - Source attribution and PMIDs

### Literature
- Semantic Scholar + PubMed search
- LLM relevance scoring with signal type extraction
- Resistance/sensitivity signal extraction
- Structured knowledge extraction from papers

### Gene Context
- Oncogene/TSG/DDR role detection
- Pathway-actionable TSGs (PTEN, NF1, VHL, etc.)
- BRAF mutation class detection (Class I/II/III)
- DDR gene therapeutic implications

---

## Near-Term (Q1 2025)

### Pathway Context (Lightweight)
**Goal:** Add pathway context to interpret co-mutations and resistance mechanisms

**Implementation:**
- Hard-coded pathway map for ~30 key cancer genes (no API calls)
- For each gene, define: pathway name, upstream nodes, downstream nodes, parallel pathways
- Annotate co-mutations with pathway relationships (same-pathway = mutual exclusivity expected)
- Surface resistance mechanisms via parallel pathway activation

**Example data structure:**
```python
PATHWAY_CONTEXT = {
    "BRAF": {
        "pathway": "MAPK/ERK",
        "upstream": ["NRAS", "KRAS", "EGFR"],
        "downstream": ["MAP2K1", "MAP2K2"],  # MEK1/2
        "parallel": ["PIK3CA", "PTEN", "AKT1"],
    },
    "EGFR": {
        "pathway": "EGFR/RAS/MAPK",
        "downstream": ["KRAS", "BRAF", "PIK3CA"],
        "parallel": ["MET", "ERBB2"],
    },
}
```

**Output example:**
```
Pathway: MAPK/ERK signaling
├── Upstream (same pathway): NRAS (mutually exclusive, 0.9%)
├── Downstream (drug targets): MEK1/2 → trametinib, cobimetinib
└── Parallel (resistance): PIK3CA co-mutation in 8% of cases

Research implication: Co-occurring PIK3CA mutations may predict
resistance to BRAF/MEK inhibition via parallel pathway activation.
```

**Research value:** Makes co-mutation data interpretable. Explains *why* BRAF/NRAS are mutually exclusive (same pathway) while BRAF/PIK3CA co-occur (parallel pathways, resistance mechanism).

### DepMap/CCLE Integration ✅ COMPLETE
**Goal:** Add cell line dependency and drug sensitivity data

**Implementation:** ✅
- Query DepMap API for gene effect scores
- CRISPR dependency scores (CERES) per gene across 1000+ cell lines
- PRISM drug sensitivity (IC50/AUC) for cell lines with the mutation
- Identify cell line models harboring the exact variant
- Data flows to LLM for research narrative generation

**Delivered:**
- `DepMapClient` with async API querying
- `DepMapEvidence` model with gene_dependency, drug_sensitivities, cell_line_models
- Helper methods: `is_essential()`, `get_top_sensitive_drugs()`, `has_data()`
- CLI panel showing essentiality, top drugs, and model cell lines
- Gap detection integration for preclinical evidence gaps

### GDSC/CTRP Drug Sensitivity
**Goal:** Systematic cell line drug response data

**Implementation:**
- IC50/AUC values by mutation status
- Drug sensitivity heatmaps
- Identify drugs with selective activity in mutant vs wild-type

**Research value:** "Cells with this mutation show IC50 = 12nM for trametinib vs 890nM in wild-type (GDSC2, n=23 lines)."

### Enhanced Conflict Detection
**Goal:** Explicitly surface disagreements between sources

**Implementation:**
- Compare CIViC vs CGI vs VICC evidence levels
- Flag when FDA approval exists but CIViC has resistance evidence
- Detect when functional predictions disagree (AlphaMissense vs PolyPhen2)
- Add `conflicts` field to Evidence model

**Research value:** "⚠️ Conflict: CGI reports sensitivity to drug X, but 2 recent papers (PMID: 123, 456) report acquired resistance via Y mechanism."

### Hypothesis Quality Scoring
**Goal:** Rank generated hypotheses by testability

**Implementation:**
- Score hypotheses based on:
  - Number of supporting evidence elements
  - Availability of model systems (cell lines, PDX)
  - Existing tool availability (drugs, assays)
  - Gap severity being addressed
- Surface top 3 hypotheses with feasibility assessment

---

## Mid-Term (Q2-Q3 2025)

### UniProt/InterPro Integration
**Goal:** Protein structure and function context

**Implementation:**
- Protein domain boundaries and functions
- Post-translational modification sites
- 3D structure context (AlphaFold)
- Known protein-protein interactions
- Annotate which domain the variant affects

**Research value:** "V600E is in the kinase domain activation segment (residues 596-600). This region is critical for autoinhibition."

**Gap detection enhancement:** "Structural impact not characterized" when variant is in unannotated region

### bioRxiv/medRxiv Integration
**Goal:** Access preprints for cutting-edge findings

**Implementation:**
- Search bioRxiv API for gene+variant
- Flag preprints as "not peer-reviewed"
- Highlight preprints contradicting published evidence
- Track preprint → publication progression

**Research value:** "🔬 Preprint (bioRxiv, 2024-12): Reports novel resistance mechanism via XYZ. Not yet peer-reviewed."

### Full-Text Literature Access (PMC)
**Goal:** Go beyond abstracts for richer evidence extraction

**The problem:** Abstracts miss critical details:
- Methods (cell line vs clinical trial? what concentrations?)
- Resistance mechanisms (buried in results/discussion)
- Subgroup analyses ("In the BRAF V600E cohort specifically...")
- Quantitative data (IC50 values, response rates)
- Negative results (abstracts oversell; full text has caveats)

**Implementation (Phase 1 - PMC Open Access):**
- Check if paper has PMCID (indicates PMC availability)
- Fetch full text via PMC E-utilities API (free, no key required)
- Parse XML to extract Methods, Results, Discussion sections
- Run LLM extraction on richer content
- Coverage: ~30-40% of oncology papers are open access

**API:**
```
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
    ?db=pmc&id=PMC1234567&rettype=xml
```

**Research value:** "From full text (PMID:12345678): IC50 = 12nM in V600E-mutant lines vs 890nM in wild-type. Resistance emerged at 6 months via NRAS Q61K co-mutation (3/12 patients)."

**Coverage limitations:**
| Source | Coverage | Access |
|--------|----------|--------|
| PMC Open Access | ~30-40% | Free API |
| Subscription journals | ~60-70% | Requires institutional access or partnerships |

**Future phases:**
- Unpaywall API for additional OA sources
- User-provided PDFs for key papers
- Publisher API partnerships (Elsevier, Springer — business development)

### Structural Variant Support
**Goal:** Extend beyond point mutations

**Implementation:**
- Fusions (ALK, ROS1, RET, NTRK, FGFR)
  - Partner gene detection
  - Breakpoint annotation
  - Fusion-specific drug approvals
- Gene amplifications (ERBB2, MET, MYC, EGFR)
  - Copy number thresholds
  - Amplification-specific therapies
- Deletions (CDKN2A, PTEN, RB1)
  - Homozygous vs heterozygous
  - Pathway implications

**New models:** `FusionEvidence`, `CNVEvidence`

### NIH RePORTER Integration
**Goal:** Research funding landscape

**Implementation:**
- Query NIH RePORTER API for grants studying the gene/variant
- Show active grants, PIs, institutions
- Identify funding trends
- Surface collaboration opportunities

**Research value:** "3 active R01 grants studying KRAS G12D resistance mechanisms. PI: Jane Smith (Stanford), $2.1M total."

### ClinGen Integration
**Goal:** Curated gene-disease validity

**Implementation:**
- Gene-disease validity scores
- Dosage sensitivity classifications
- Expert panel assertions
- Actionability assessments

**Research value:** "ClinGen: BRCA1 has DEFINITIVE evidence for hereditary breast/ovarian cancer (3/3 validity score)."

---

## Long-Term (2025+)

### Research Tracking & Alerts
**Goal:** Longitudinal evidence monitoring

**Implementation:**
- Watch list for variants of interest
- Daily/weekly scans for:
  - New PubMed publications
  - New clinical trials
  - FDA approval changes
  - CIViC/CGI updates
- Email/Slack alerts with diff summaries
- Evidence quality trend tracking

**Research value:** "📬 Alert: 2 new papers on EGFR C797S this week. Evidence quality upgraded from 'limited' to 'moderate'."

### Batch Analysis & Cohort Features
**Goal:** Multi-variant analysis

**Implementation:**
- Upload VCF/MAF files
- Cohort-level gap analysis ("Your cohort has 40% variants with no clinical evidence")
- Cross-variant pattern detection
- Prioritization scoring for research portfolio
- Pathway enrichment across variants
- Co-mutation network visualization

**Research value:** "Your 50-variant cohort: 12 have FDA-approved therapies, 18 have clinical trial options, 20 have critical evidence gaps. Top pathway: RAS/MAPK (60%)."

### API & Developer Platform
**Goal:** Programmatic access

**Implementation:**
- RESTful API with OpenAPI spec
- Python SDK (`pip install oncomind`)
- Jupyter notebook widgets
- LIMS/EHR integration hooks (FHIR resources)
- Webhook callbacks for watch list alerts

**Endpoints:**
```
GET  /api/v1/variants/{gene}/{variant}
POST /api/v1/variants/batch
GET  /api/v1/variants/{id}/gaps
GET  /api/v1/variants/{id}/hypotheses
POST /api/v1/watchlist
```

### Agentic Multi-Source Research
**Goal:** Autonomous deep-dive on novel variants

**Implementation:**
- Agent workflow for poorly-characterized variants:
  1. Identify gaps
  2. Query additional sources (UniProt, DepMap, bioRxiv)
  3. Find similar variants with evidence
  4. Generate extrapolation hypotheses
  5. Compile research brief
- Human-in-the-loop validation
- Confidence calibration

**Research value:** "No direct evidence for ATM R2993H. Agent found: 3 nearby pathogenic variants, similar residue in BRCA1 shows LOF, DepMap shows ATM-dependent lines. Extrapolated hypothesis: Likely LOF with PARP inhibitor sensitivity (confidence: 0.6)."

### Comparative Analysis Mode
**Goal:** Compare variants head-to-head

**Implementation:**
- Side-by-side evidence comparison
- Differential gap analysis
- Shared vs unique therapeutic options
- Mutation class comparisons (BRAF V600E vs G469A)

**Research value:** "BRAF V600E vs G469A: V600E has 5 FDA-approved drugs, G469A has 1. Both share MEK inhibitor sensitivity. V600E-specific: BRAF monomer inhibitors. G469A-specific: Requires combination therapy."

### Evidence Versioning & History
**Goal:** Track how evidence evolves

**Implementation:**
- Snapshot evidence state at query time
- Track when evidence quality changed
- Show evidence trajectory over time
- Identify rapidly evolving variants

**Research value:** "KRAS G12C evidence timeline: 2019 (preclinical only) → 2021 (Phase 2 sotorasib) → 2022 (FDA approved) → 2024 (resistance mechanisms emerging)."

### Pathway-Level Gap Analysis
**Goal:** Answer "What don't we know about the JAK-STAT pathway?"

This is the full vision beyond the lightweight pathway context in Near-Term. Instead of annotating a single variant with pathway context, this enables querying an entire pathway for research gaps.

**The challenge:**
- Pathways aren't well-defined entities (Reactome: 87 genes, KEGG: 42, MSigDB: 158 for "JAK-STAT")
- Requires aggregating gaps across many genes with different evidence types
- Need to weight by biological importance (JAK2 matters more than adapter proteins)
- LLM context limits make dumping 50+ gene summaries impractical

**Implementation:**
- Reactome API integration for pathway membership
- Gene-level gap scoring (new capability)
- Pathway aggregation logic with importance weighting
- Smart summarization: which gaps matter most?
- New CLI: `mind pathway MAPK --tumor NSCLC`

**Output example:**
```
JAK-STAT Pathway Gap Analysis

Well-characterized nodes:
├── JAK2 V617F: comprehensive (FDA drugs, 500+ papers, resistance known)
├── STAT3: 200+ papers, multiple inhibitors in trials
└── JAK1: moderate evidence, emerging therapeutics

Under-studied nodes:
├── STAT5B: limited functional data, no trials
├── SOCS1: tumor suppressor role unclear in solid tumors
└── PTPN11: preclinical only, mechanism debated

Pathway-level gaps:
├── Cross-talk with PI3K: conflicting data
├── Resistance mechanisms: JAK1 << JAK2 (less studied)
└── Biomarkers for patient selection: weak

Research opportunities:
├── "SOCS1 loss as synthetic lethal target"
├── "JAK1-specific resistance vs JAK2"
└── "STAT5B in solid tumors"
```

**Research value:** Genuinely novel capability. Nobody does "pathway gap analysis" well — it's all manual literature review today. This could help PIs identify white space for grant proposals.

**Estimated effort:** 6-12 months (not blocking current release)

---

## Technical Debt & Infrastructure

### Performance Optimization
- Pre-compute cBioPortal prevalence for common genes (cache layer)
- Async batch fetching for all API clients
- Response caching with TTL (Redis/SQLite)
- Rate limit handling with exponential backoff

### Testing & Validation
- Golden dataset of well-characterized variants for regression testing
- Automated comparison against OncoKB/CIViC for accuracy
- LLM output validation (JSON schema, citation verification)
- Evidence gap detector unit tests

### Observability
- Structured logging with correlation IDs
- API latency metrics per source
- LLM token usage tracking
- Error rate monitoring per data source

---

## Design Principles

1. **Gaps over facts** — Prioritize surfacing what's unknown
2. **Sources required** — Every claim must link to evidence
3. **Conflicts visible** — Never hide disagreements between sources
4. **Research-first** — Optimize for research workflows, not clinical decisions
5. **Hypothesis-generating** — Output should suggest next experiments
6. **Quantitative when possible** — Include sample sizes, p-values, IC50s
7. **Extrapolations labeled** — Clearly mark inferences vs direct data

---

## Success Metrics

| Metric | Current | Q1 Target | Q2 Target |
|--------|---------|-----------|-----------|
| Data sources integrated | 14 | 17 | 22 |
| Gap categories detected | 8 | 10 | 12 |
| Avg gaps identified per query | ~3 | ~5 | ~6 |
| Evidence quality "comprehensive" rate | ~20% | ~25% | ~30% |
| Hypothesis generation rate | 100% | 100% | 100% |
| User research questions answered | - | Track | Improve |

*Note: Data sources increased from 12→14 with DepMap (gene essentiality + drug sensitivity) and cBioPortal integrations.*

---

## Contributing

We welcome contributions! Priority areas:
- New data source integrations (see Near-Term goals)
- Evidence gap detection improvements
- Hypothesis generation algorithms
- Structural variant support
- Conflict detection logic

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.



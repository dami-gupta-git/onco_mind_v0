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
  - Gap categories: functional, clinical, tumor_type, drug_response, resistance, preclinical, prevalence, prognostic, discordant, validation (10 categories)
  - Severity levels: critical, significant, minor
  - Suggested studies and addressable data sources per gap
  - Context-aware enrichment (dynamic suggestions based on actual evidence)
  - Hotspot and near-hotspot detection (via cancerhotspots.org integration)
- `CharacterizedAspect` model with basis explanations for what's well-characterized
- Evidence quality assessment (comprehensive/moderate/limited/minimal) using net scoring
- Research priority scoring (very_high/high/medium/low based on gene importance + gap profile)
- Cross-source conflict detection (CGI vs CIViC vs VICC)
- LLM-powered hypothesis generation with evidence basis tags
- Cited sources for all claims (enforced in prompts)

### Therapeutic Data Model
- `TherapeuticData` model with:
  - All evidence levels (FDA → preclinical → computational)
  - Response type (sensitivity/resistance)
  - Mechanism of action
  - Quantitative data (IC50, response rates)
  - Source attribution and PMIDs
  - Match specificity tracking (variant/codon/gene level)
  - Cell lines tested for preclinical data
  - Confidence scoring

### Literature
- Semantic Scholar + PubMed search
- LLM relevance scoring with signal type extraction (resistance/sensitivity/mixed/prognostic/unclear)
- `LiteratureKnowledge` model with structured resistance/sensitivity extraction
- `DrugResistance` and `DrugSensitivity` models with predictive vs prognostic distinction
- Multi-provider LLM support via litellm (Claude, GPT-4o, GPT-4o-mini, Haiku)

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

### Robust FDA Drug Fetching Algorithm
**Goal:** Reliable, comprehensive FDA approval and label data fetching

**Current problems:**
- OpenFDA API is rate-limited and sometimes returns incomplete data
- Drug name normalization is fragile (brand vs generic, spelling variants)
- Curated `ONCOLOGY_DRUG_MAPPINGS` requires manual maintenance
- **Development codes not mapped:** CIViC/CGI use development codes (AZD5363, AMG 510, MRTX849) but OpenFDA indexes by INN (capivasertib, sotorasib, adagrasib)
- Biomarker specificity extraction (variant vs gene vs phenotype level) has edge cases
- FDA label cache can become stale
- **Combination therapy not handled:** FDA approvals for combination regimens (e.g., Adagrasib + cetuximab for KRAS G12C CRC) are not parsed or displayed as combinations. Each drug appears separately without indicating the approved combination context.

**Implementation:**
1. **Multi-source drug resolution:**
   - OpenFDA (primary)
   - DailyMed API (backup for label text)
   - RxNorm for drug name normalization
   - ChEMBL for drug-target mappings

2. **Improved drug name matching:**
   - Fuzzy matching for spelling variants
   - Brand ↔ generic name bidirectional lookup
   - Handle combination drugs (e.g., "dabrafenib + trametinib")
   - Normalize suffixes (-ib, -mab, -nib patterns)
   - **Development code → INN mapping** via ChEMBL/PubChem synonyms (AZD5363 → capivasertib)
   - **Combination therapy parsing:** Parse indication text for "in combination with" patterns and display as unified regimen (e.g., "Adagrasib + Cetuximab" instead of separate entries)

3. **Better biomarker specificity extraction:**
   - Handle negation ("not approved for", "excluding patients with")
   - Parse companion diagnostic requirements
   - Extract line-of-therapy constraints
   - Detect tumor-agnostic vs tumor-specific approvals

4. **Cache management:**
   - Nightly refresh job for all oncology drugs
   - Version tracking for FDA label updates
   - Delta detection for new approvals

**Research value:** Reliable FDA data is foundation for all regulatory evidence. False negatives (missing approvals) and false positives (wrong specificity) directly impact clinical utility.

### DGIdb Integration
**Goal:** Replace manual drug-target mappings with comprehensive drug-gene interaction data

**The problem:** The manual `BIOMARKER_SELECTION_DRUGS` mapping in `constants.py` is incomplete and requires manual curation. We need to distinguish between:
- Drugs that directly target a gene (e.g., osimertinib targets EGFR)
- Drugs where a gene mutation is a patient selection biomarker but not the drug target (e.g., datopotamab deruxtecan targets TROP2, but is indicated for EGFR-mutant NSCLC patients)

**Implementation:**
- Integrate DGIdb (Drug-Gene Interaction Database) API
- ~40k drug-gene interactions with mechanism annotations
- Cache results per gene (interactions don't change frequently)
- Use interaction_types to classify: inhibitor, agonist, antibody, etc.
- API: https://dgidb.org/api

**Research value:** Automatic, comprehensive drug-target annotation without manual curation. Better distinction between direct targets and biomarker associations.

### Enhanced Conflict Detection
**Goal:** Explicitly surface disagreements between sources

**Implementation:**
- Compare CIViC vs CGI vs VICC evidence levels
- Flag when FDA approval exists but CIViC has resistance evidence
- Detect when functional predictions disagree (AlphaMissense vs PolyPhen2)
- Add `conflicts` field to Evidence model

**Research value:** "⚠️ Conflict: CGI reports sensitivity to drug X, but 2 recent papers (PMID: 123, 456) report acquired resistance via Y mechanism."

### LLM Therapeutic Contextualization
**Goal:** Explain *why* certain drugs are being studied — connect variant biology to therapeutic rationale

**Current state:** We already send gene role (oncogene/TSG) and pathway to the LLM. The Cross-Source Drug Analysis identifies drug corroboration across sources.

**Enhancement — Mechanistic Rationale (without hallucinating):**
```
"AKT1 is an oncogene in the PI3K/AKT pathway. E17K is a known activating
mutation, making direct AKT inhibition a logical therapeutic strategy."
```

**Prompt instruction:**
```
"Explain why the identified drugs are relevant to this variant based on the
gene role and pathway information provided. Do not invent mechanisms — only
use what's stated in GENE ROLE and PATHWAY sections."
```

**Implementation consideration:** We currently tell the LLM not to describe mechanisms ("NEVER describe HOW a variant works mechanistically"). Need to relax this slightly to allow pathway-level context while still prohibiting invented molecular details.

**Research value:** Researchers understand not just *what* drugs are relevant, but *why* — enabling better experimental design and hypothesis refinement.

### Cross-Source Conflict Detection (Pre-computed)
**Goal:** Flag when evidence sources disagree on the same drug, with structured explanations

**Types of conflicts to detect:**
1. **Same drug, opposite association** — CGI says Responsive, CIViC says Resistant
2. **Same drug class, mixed signals** — One AKT inhibitor works, another doesn't
3. **Different evidence levels disagree** — FDA approved but literature shows resistance emerging

**Implementation — Pre-compute conflicts before LLM:**
Detect conflicts in code and pass a structured `conflicts` section to the LLM prompt:
```
CONFLICTS DETECTED:
- Imatinib: CGI says Responsive (gene-level), CGI also says Resistant (D816V-specific)
- Vemurafenib: FDA approved (Melanoma), Literature reports resistance in CRC context
```

**Prompt instruction:**
```
"Review the evidence for contradictions. Distinguish between:
- Expected biology (e.g., resistance mutation to one drug, sensitivity to next-gen drug)
- True conflicts (same drug, same context, opposite outcomes)
Flag true conflicts explicitly."
```

**Output enhancement:** Add `conflicts_identified` field to JSON output:
```json
{
  "conflicts_identified": [
    {
      "drug": "Imatinib",
      "conflict_type": "same_drug_opposite_signal",
      "sources": ["CGI (Responsive)", "CGI (Resistant)"],
      "likely_explanation": "D816V-specific resistance vs gene-level sensitivity",
      "clinical_implication": "Variant-level testing critical for imatinib decisions"
    }
  ]
}
```

**Research value:** Surfaces contradictions that researchers would otherwise miss by manually reviewing multiple databases. Enables smarter trial design and biomarker selection.

### Hypothesis Quality Scoring
**Goal:** Rank generated hypotheses by testability

**Implementation:**
- Score hypotheses based on:
  - Number of supporting evidence elements
  - Availability of model systems (cell lines, PDX)
  - Existing tool availability (drugs, assays)
  - Gap severity being addressed
- Surface top 3 hypotheses with feasibility assessment

### Ensemble LLM Extraction
**Goal:** Run extraction across multiple LLMs for high-stakes calls

**When to use:**
| Use Case | Ensemble Worth It? |
|----------|-------------------|
| Resistance/sensitivity extraction | Yes |
| Paper relevance scoring | Maybe |
| Narrative generation | No |
| Trial eligibility parsing | Yes |

**Implementation:**
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

### Variant-to-Phenotype (V2P) Inference
**Goal:** Rule-based phenotype inference from variants

**Well-established mappings:**
- MSI-H from MMR gene mutations (MLH1, MSH2, MSH6, PMS2)
- HRD from BRCA1/2 loss
- Hypermutation from POLE/POLD1 mutations

**Output model:**
```python
class PhenotypeInference(BaseModel):
    msi_status: Literal["MSI-H", "MSS", "unknown"]
    hrd_status: Literal["HRD", "HRP", "unknown"]
    hypermutated: bool
    actionable_phenotypes: list[str]  # ["immunotherapy_responsive", "PARP_sensitive"]
    inferred_from: list[str]
    confidence: float
```

### Synthetic Lethality Mapping
**Goal:** Surface established synthetic lethal relationships

**Well-established pairs:**
| Context | Vulnerability | Drug Class | Evidence |
|---------|--------------|------------|----------|
| BRCA1/2 loss | PARP dependency | PARP inhibitors | FDA-approved |
| MSI-H | Immune evasion | Checkpoint inhibitors | FDA-approved |
| ATM loss | DNA repair | PARP/ATR inhibitors | Emerging |
| MTAP loss | Methionine salvage | MAT2A inhibitors | Emerging |

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

### Multi-Variant Reasoning
**Goal:** "I have these 3 variants together, what does that mean?"

Current tools annotate variants independently. Real tumors have multiple variants that interact.

**Concrete examples:**
- EGFR L858R + TP53 → worse prognosis, shorter TKI response
- EGFR L858R + MET amp → potential primary resistance, consider combo
- BRAF V600E + PIK3CA → may reduce BRAF inhibitor efficacy

**Implementation:**
1. Annotate each variant individually
2. Search literature for co-occurrence patterns
3. LLM synthesis of combined implications
4. Surface known interaction patterns (can curate top combinations)

### Resistance Escape Routes
**Goal:** Most tools tell you what the variant *is*. OncoMind tells you what happens *next*.

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

### Clinical Question Answering
**Goal:** Not a report, but an answer to a specific question.

Instead of: "Here's everything about BRAF V600E"
Answer: "Should I try pembrolizumab given BRAF V600E and melanoma?"

**Question types:**
- "Should I try [drug] given [variant] and [tumor]?"
- "What's the expected response duration for [variant] on [drug]?"
- "What should I monitor for resistance?"
- "Is this patient eligible for [trial NCT#]?"

**Risk:** Close to clinical decision support. Requires strong disclaimers.

### Contextual Query with Treatment History
**Goal:** "What does this variant mean for a patient who already failed X therapy?"

**Example:**
- Input: "EGFR L858R, failed erlotinib (8mo), now failed osimertinib (14mo)"
- Reasoning:
  1. T790M likely emerged after erlotinib
  2. Osimertinib failure suggests C797S or bypass (MET, HER2)
  3. Check resistance mechanism to guide next therapy
- Output: Ranked options with rationale specific to treatment history

### Dynamic Evidence Synthesis with Temporal Weighting
**Goal:** Evidence has a temporal dimension. A 2018 FDA approval doesn't account for 2024 resistance reports.

**Features:**
- Weight by recency, evidence level, tumor specificity
- Surface temporal trends: "Older sources say sensitive, newer sources report resistance"
- Conflict resolution with rationale

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

Your roadmap mentions this—do it soon. Reactome has a clean REST API for:
Mapping genes/variants to pathways/reactions.
Disease-specific annotations (e.g., "Signaling by EGFR in Cancer").
Visual diagrams and downstream effectors.
→ For a variant like MAP2K1 P124L: "Activates MEK-ERK pathway → downstream targets include MYC, CCND1" → Hypothesis: "Test ERK-dependent transcription in models."
→ In LLM prompts: "Pathway involvement: Hyperactive MAPK signaling (Reactome R-HSA-5684996)" + gaps like "No annotated resistance bypass in this branch."

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

How You'd Use Reactome Data in OncoMind
The ContentService API is straightforward—RESTful endpoints for querying by gene/protein (UniProt ID), pathway ID, or variant. Key ways to leverage it:

Provide Mechanistic Context
For a variant like BRAF V600E in melanoma:
Query pathways involving BRAF → Returns "Signaling by BRAF and RAF fusions," "RAF/MAP kinase cascade," and disease-specific ones like "Signaling by BRAF in Cancer."
Highlights how V600E leads to constitutive RAF activation → downstream hyperactive ERK signaling.
→ In your evidence backbone: Add a "Pathway Involvement" section with pathway names, IDs, and brief descriptions.

Annotate Functional Impact at Reaction Level
Reactome has dedicated "disease reactions" for variants (e.g., "Constitutive Signaling by Mutant BRAF" or similar for MAP2K1 alterations).
Shows altered inputs/outputs (e.g., BRAF V600E dimerizes independently of RAS, bypassing negative feedback).
→ For under-studied variants (your sweet spot): If no clinical data, but pathway shows gain-of-function in MAPK → "Predicted oncogenic activation of RAF-independent ERK signaling."

Generate Richer Hypotheses (LLM Mode)
Feed structured Reactome context into prompts:textPathway Context <a href="https://reactome.org" target="_blank" rel="noopener noreferrer nofollow"></a>:
  - Involved in: RAF/MAP kinase cascade (R-HSA-5673001), Signaling by BRAF in Cancer
  - Altered reaction: Constitutive BRAF activation → hyperphosphorylation of MEK/ERK
  - Downstream effectors: MYC, CCND1 (proliferation), DUSP6 (feedback inhibition bypassed)
  - Co-occurring genes in pathway: NRAS (mutually exclusive), CDKN2A (co-loss)→ LLM outputs: "Hypothesis: Test if [variant] bypasses RAS dependency, similar to BRAF V600E—use RAS-independent cell lines or ERK readout assays."
→ Flags gaps: "No annotated bypass mechanisms in this subpathway" or "Missing preclinical data on downstream node X."
Enhance Gap Detection & Evidence Quality
High pathway annotation + low clinical data → "Mechanistically well-understood (Reactome disease pathway) but lacking therapeutic validation."
For co-mutations (your cBioPortal strength): Reactome shows if they hit the same reaction/pathway → "Potential synthetic lethality in [pathway node]."

Suggest Testable Nodes/Experiments
List upstream regulators, downstream effectors, or parallel pathways.
→ "Research implication: Target downstream ERK (if inhibitor-resistant) or combine with PI3K inhibitors (parallel survival pathway)."


Practical Integration Steps (Low Effort)

Use gene symbol/UniProt to query /data/participates/gene/{symbol}/human or /data/pathways/ids/{uniprot}.
For variants: Reactome annotates many cancer hotspots directly in disease pathways (e.g., AKT1 E17K constitutive signaling).
Cache results per gene (pathways don't change often).
Start with top pathways + disease-flagged ones.

Example Outputs in OncoMind
For BRAF V600E in Melanoma:
Pathway Insight (Reactome): Constitutively active in "RAF-independent MAPK signaling" → drives ERK without RAS. Mutual exclusivity with NRAS explained by pathway redundancy.
Knowledge Gap: Well-characterized activation mechanism, but limited data on variant-specific resistance bypasses.
Hypothesis: Co-occurring CDKN2A loss may enhance proliferative output—test in isogenic models.
For rarer MAP2K1 P124L:
If Reactome shows MEK gain-of-function in MAPK cascade → "Predicted hyperactivation downstream of BRAF; potential MEK inhibitor sensitivity, but resistance via ERK rebound."
This adds the "why" and "how to test" layer that's missing from most tools—perfect for turning gaps into hypotheses. It's exactly what elevates OncoMind from annotator to true research intelligence.
When you're ready, this is the feature that will make people go "wow." 🚀

**Research value:** Genuinely novel capability. Nobody does "pathway gap analysis" well — it's all manual literature review today. This could help PIs identify white space for grant proposals.

**Estimated effort:** 6-12 months (not blocking current release)


### GDSC (Genomics of Drug Sensitivity in Cancer)
Complementary to DepMap PRISM:
~1,000 cell lines, ~500 compounds (more targeted oncology drugs).
Public API/downloads.
→ Adds sensitivity data for drugs not in PRISM; cross-validation boosts confidence in preclinical hypotheses.

---

## Technical Debt & Infrastructure

### Code Quality
- [ ] **Tumor type matching consistency**: Check if `_tumor_type_matches` should be used everywhere for tumor matching, e.g. in `evidence.py:524` (clinical trial cancer specificity matching). The curated lists in `gene_context.py` are already normalized so simple substring matching is appropriate there, but external data sources (ClinicalTrials.gov, CIViC, CGI) use varied terminology that may benefit from alias resolution.

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

#### Test Coverage Gaps (as of Dec 2024)

**Current Coverage:** 406 tests

**Well-Covered:**
| Source | Unit Tests | Integration Tests |
|--------|-----------|-------------------|
| cBioPortal | 24 | 12 |
| DepMap | 39 | 17 |
| CGI | 21 | via evidence_fields |
| VICC | 25 | 16 |
| VEP | - | 20 |
| MyVariant | - | 16 |
| PubMed | - | 14 |
| Semantic Scholar | - | 10 |
| FDA | - | 23 |

**Missing Integration Tests (Priority Order):**
1. CIViC - Core evidence source, no dedicated API tests
2. ClinicalTrials.gov - Important for clinical relevance
3. OncoKB - If enabled, needs coverage

**Missing Unit Tests:**
- CIViC client
- ClinicalTrials client
- PubMed client (only integration)
- Semantic Scholar client (only integration)
- FDA client
- VEP client

**Other Gaps:**
- End-to-end Conductor tests (full pipeline)
- Error handling tests (network failures, rate limits)
- Edge case tests (rare variants, missing data)
- Performance tests (timeout handling, concurrent API calls)

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
| Gap categories detected | 10 | 11 | 12 |
| Avg gaps identified per query | ~4 | ~5 | ~6 |
| Evidence quality "comprehensive" rate | ~20% | ~25% | ~30% |
| Hypothesis generation rate | 100% | 100% | 100% |
| User research questions answered | - | Track | Improve |

*Note: Data sources increased from 12→14 with DepMap (gene essentiality + drug sensitivity) and cBioPortal integrations. Gap categories expanded from 8→10 with addition of DISCORDANT and VALIDATION categories.*

### Hardcoded Data Drift
Location: config/constants.py - DEPMAP_GENE_DEPENDENCIES_FALLBACK, DEPMAP_DRUG_SENSITIVITIES_FALLBACK
Problem: These constants will become stale. DepMap updates quarterly.
Suggestion:

Add version/date comments to the fallback data
Consider a script to refresh from DepMap downloads
Log when fallback data is used so users know it may be dated


#### Minor
Fix docstrings - inconsistent
More typehints
---

## Contributing

We welcome contributions! Priority areas:
- New data source integrations (see Near-Term goals)
- Evidence gap detection improvements
- Hypothesis generation algorithms
- Structural variant support
- Conflict detection logic

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.



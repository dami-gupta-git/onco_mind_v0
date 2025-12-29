# OncoMind Evidence Sources

This document describes each external data source OncoMind fetches from, what data is retrieved, and how it flows through the evidence synthesis pipeline.

---

## Table of Contents

1. [Evidence Architecture](#evidence-architecture) ← **How evidence flows through the system**
2. [Functional Annotations](#functional-annotations) ← **Computational pathogenicity predictions**
3. [MyVariant.info](#myvariantinfo)
4. [CIViC (Clinical Interpretation of Variants in Cancer)](#civic)
5. [VICC MetaKB](#vicc-metakb)
6. [CGI (Cancer Genome Interpreter)](#cgi)
7. [FDA OpenFDA](#fda-openfda)
8. [ClinicalTrials.gov](#clinicaltrialsgov)
9. [Ensembl VEP](#ensembl-vep)
10. [Semantic Scholar](#semantic-scholar)
11. [PubMed](#pubmed)
12. [OncoTree](#oncotree)
13. [cBioPortal](#cbioportal)
14. [DepMap (Cancer Dependency Map)](#depmap)
15. [Evidence Gap Detection](#evidence-gap-detection)
16. [LLM Research Synthesis](#llm-research-synthesis)
17. [Therapeutic Evidence Model](#therapeutic-evidence-model)
18. [Literature Knowledge Model](#literature-knowledge-model)

---

## Evidence Architecture

OncoMind uses a layered evidence architecture that separates deterministic annotation from LLM-powered synthesis.

### Two-Layer Design

```
┌─────────────────────────────────────────────────────────────────────┐
│                     LLM RESEARCH LAYER (optional)                   │
│  • Research-focused narrative synthesis                             │
│  • Evidence quality calibration                                     │
│  • Hypothesis generation with evidence basis tags                   │
│  • Knowledge gap identification                                     │
└─────────────────────────────────────────────────────────────────────┘
                                  ▲
                                  │ Structured context
                                  │
┌─────────────────────────────────────────────────────────────────────┐
│                    ANNOTATION BACKBONE (deterministic)              │
│  Evidence Model ◄── EvidenceAggregator ◄── API Clients              │
│       │                                         │                   │
│       ├── identifiers (COSMIC, dbSNP, ClinVar) │                   │
│       ├── functional (AlphaMissense, CADD)     │                   │
│       ├── clinical (FDA, CIViC, VICC, CGI)     │                   │
│       ├── biological (cBioPortal, DepMap)      │                   │
│       ├── literature (PubMed, Semantic Scholar)│                   │
│       └── evidence_gaps (computed)             │                   │
└─────────────────────────────────────────────────────────────────────┘
```

### Match Specificity Tracking

All evidence items track how precisely they match the query variant:

| Match Level | Meaning | Example |
|-------------|---------|---------|
| `variant` | Exact variant match | BRAF V600E → evidence specifically for V600E |
| `codon` | Same position, different change | BRAF V600K → evidence for "V600 mutations" |
| `gene` | Gene-level only | BRAF V600E → evidence for "BRAF mutations" |

This is tracked via the `match_level` field on `TherapeuticEvidence`, `CIViCAssertionEvidence`, `VICCEvidence`, `CGIBiomarkerEvidence`, and `FDAApproval`.

### Evidence Model Structure

The `Evidence` model ([evidence.py](src/oncomind/models/evidence/evidence.py)) aggregates all sources:

```python
class Evidence(BaseModel):
    # Core info
    identifiers: VariantIdentifiers  # gene, variant, COSMIC ID, dbSNP, HGVS
    functional: FunctionalScores     # AlphaMissense, CADD, PolyPhen2, gnomAD
    context: VariantContext          # tumor_type, gene_role, pathway

    # Evidence lists (one per source)
    fda_approvals: list[FDAApproval]
    civic_assertions: list[CIViCAssertionEvidence]
    civic_evidence: list[CIViCEvidence]
    vicc_evidence: list[VICCEvidence]
    cgi_biomarkers: list[CGIBiomarkerEvidence]
    clinvar_entries: list[ClinVarEvidence]
    cosmic_entries: list[COSMICEvidence]
    clinical_trials: list[ClinicalTrialEvidence]
    pubmed_articles: list[PubMedEvidence]
    preclinical_biomarkers: list[CGIBiomarkerEvidence]
    early_phase_biomarkers: list[CGIBiomarkerEvidence]

    # Rich context sources
    cbioportal_evidence: CBioPortalEvidence | None
    depmap_evidence: DepMapEvidence | None
    literature_knowledge: LiteratureKnowledge | None

    # Computed analysis
    evidence_gaps: EvidenceGaps | None
```

### Key Helper Methods

The `Evidence` model provides methods for downstream consumption:

| Method | Purpose |
|--------|---------|
| `get_therapeutic_evidence()` | Returns `list[TherapeuticEvidence]` sorted by evidence tier |
| `get_resistance_summary()` | Synthesized narrative of resistance signals |
| `get_sensitivity_summary()` | Synthesized narrative of sensitivity signals |
| `get_evidence_summary_for_llm()` | Compact text for LLM prompt |
| `get_biological_context_for_llm()` | cBioPortal + DepMap formatted for LLM |
| `get_literature_summary_for_llm()` | PubMed articles with signal types |
| `compute_evidence_gaps()` | Runs gap detection and caches result |

---

## Functional Annotations

**Purpose:** Computational pathogenicity predictions and population frequency data displayed in the UI "Functional" tab.

### Data Source

All functional annotations are fetched from **MyVariant.info** API, which aggregates data from:
- **dbNSFP** (AlphaMissense, PolyPhen2, SIFT, CADD)
- **gnomAD** (population allele frequencies)
- **SnpEff** (variant effect predictions)

Fallback: When MyVariant.info returns no data, **Ensembl VEP** provides AlphaMissense, PolyPhen2, and SIFT predictions.

### Fields Retrieved

| Field | Source | MyVariant.info Path | Description |
|-------|--------|---------------------|-------------|
| AlphaMissense Score | dbNSFP | `dbnsfp.alphamissense.score` | Pathogenicity score (0-1), higher = more pathogenic |
| AlphaMissense Prediction | dbNSFP | `dbnsfp.alphamissense.pred` | P = Pathogenic, B = Benign, A = Ambiguous |
| PolyPhen2 Prediction | dbNSFP | `dbnsfp.polyphen2.hdiv.pred` | D = Damaging, P = Possibly damaging, B = Benign |
| PolyPhen2 Score | dbNSFP | `dbnsfp.polyphen2.hdiv.score` | Score (0-1), higher = more damaging |
| SIFT Prediction | dbNSFP | `dbnsfp.sift.pred` | D = Deleterious, T = Tolerated |
| SIFT Score | dbNSFP | `dbnsfp.sift.score` | Score (0-1), lower = more deleterious |
| CADD PHRED | dbNSFP | `dbnsfp.cadd.phred` | PHRED-scaled score (>20 = top 1% deleterious) |
| gnomAD Exome AF | gnomAD | `gnomad_exome.af.af` | Population allele frequency in gnomAD exomes |
| gnomAD Genome AF | gnomAD | `gnomad_genome.af.af` | Population allele frequency in gnomAD genomes |
| SnpEff Effect | SnpEff | `snpeff.ann[].effect` | Predicted variant effect (e.g., missense_variant) |
| SnpEff Impact | SnpEff | `snpeff.ann[].impact` | Impact level: HIGH, MODERATE, LOW, MODIFIER |

### Interpretation Guide

#### AlphaMissense
- **P (Pathogenic)**: Score ≥ 0.564 — likely pathogenic
- **A (Ambiguous)**: Score 0.340–0.564 — uncertain significance
- **B (Benign)**: Score < 0.340 — likely benign

#### PolyPhen2 (HumDiv)
- **D (Probably Damaging)**: Score ≥ 0.909
- **P (Possibly Damaging)**: Score 0.447–0.909
- **B (Benign)**: Score < 0.447

#### SIFT
- **D (Deleterious)**: Score ≤ 0.05
- **T (Tolerated)**: Score > 0.05

#### CADD PHRED Score
- **≥ 30**: Top 0.1% most deleterious
- **≥ 20**: Top 1% most deleterious
- **≥ 10**: Top 10% most deleterious

#### gnomAD Allele Frequency
- **< 0.01%**: Rare (potential somatic or germline pathogenic)
- **0.01-1%**: Low frequency
- **> 1%**: Common (likely benign polymorphism)

### UI Display

The "Functional" tab displays these fields in a condensed format:

| UI Label | Field | Example Display |
|----------|-------|-----------------|
| AlphaMissense | score, prediction | `0.985, P` |
| PolyPhen2 | score, prediction | `0.99, D` |
| gnomAD AF | gnomad_exome_af | `3.98e-06, Rare` |
| SnpEff | effect | `missense_variant` |

### How Functional Data is Used

#### 1. Evidence Gap Detection

Functional scores are used to detect evidence gaps in `gap_detector.py`:

```python
has_functional = (
    evidence.functional.alphamissense_score is not None or
    evidence.functional.cadd_score is not None or
    evidence.functional.polyphen2_prediction is not None
)

if has_functional:
    well_characterized.append("computational pathogenicity")
else:
    gaps.append(EvidenceGap(
        category=GapCategory.FUNCTIONAL,
        severity=GapSeverity.SIGNIFICANT,
        description=f"No computational pathogenicity predictions for {gene} {variant}",
        suggested_studies=["Run AlphaMissense, CADD, PolyPhen2"],
        addressable_with=["MyVariant.info", "VEP"]
    ))
```

#### 2. LLM Context

Functional data informs the LLM's `functional_summary` output. The LLM uses these scores to:

- **Calibrate confidence**: Strong pathogenicity predictions (AlphaMissense P, PolyPhen2 D) support functional impact claims
- **Flag discordance**: When predictions disagree (e.g., AlphaMissense P but SIFT Tolerated), the LLM notes this as a conflict
- **Support rarity**: Rare gnomAD AF (< 0.01%) supports the variant being a somatic driver

The `FunctionalScores.get_pathogenicity_summary()` method generates a compact summary:

```
AlphaMissense: Pathogenic (0.99) | CADD: 32.5 | PolyPhen2: probably_damaging | SIFT: deleterious
```

#### 3. Evidence Quality Assessment

Functional predictions contribute to the overall evidence quality score:

| Has Functional Data | Contribution |
|---------------------|--------------|
| Yes | "computational pathogenicity" added to `well_characterized` |
| No | SIGNIFICANT gap flagged, lowers evidence quality |

### Data Model

Functional scores are stored in the `FunctionalScores` model ([evidence.py:78](src/oncomind/models/evidence/evidence.py#L78)):

```python
class FunctionalScores(BaseModel):
    # AlphaMissense
    alphamissense_score: float | None
    alphamissense_prediction: str | None  # P/B/A

    # CADD
    cadd_score: float | None  # PHRED score
    cadd_raw: float | None

    # PolyPhen2
    polyphen2_prediction: str | None
    polyphen2_score: float | None

    # SIFT
    sift_prediction: str | None
    sift_score: float | None

    # SnpEff
    snpeff_effect: str | None
    snpeff_impact: str | None  # HIGH/MODERATE/LOW/MODIFIER

    # Population frequencies
    gnomad_exome_af: float | None
    gnomad_genome_af: float | None
```

### Code References

- **FunctionalScores model**: [evidence.py](src/oncomind/models/evidence/evidence.py) - `FunctionalScores` class
- **Gap detection for functional**: [gap_detector.py](src/oncomind/insight_builder/gap_detector.py) - `_check_functional_predictions()`
- **MyVariant API client**: [myvariant.py](src/oncomind/api/myvariant.py)
- **VEP API client**: [vep.py](src/oncomind/api/vep.py)

---

## MyVariant.info

**API Endpoint:** `https://myvariant.info/v1`

**Purpose:** Primary aggregator for variant annotations from multiple databases.

### Data Retrieved

| Field | Source | Description |
|-------|--------|-------------|
| `cosmic_id` | COSMIC | COSMIC mutation ID (e.g., "COSM476") |
| `dbsnp_id` | dbSNP | dbSNP rs number (e.g., "rs113488022") |
| `clinvar_id` | ClinVar | ClinVar variation ID |
| `clinvar_clinical_significance` | ClinVar | Clinical significance (e.g., "Pathogenic") |
| `clinvar_accession` | ClinVar | RCV accession number |
| `hgvs_genomic` | Various | Genomic HGVS notation (e.g., "chr7:g.140453136A>T") |
| `hgvs_protein` | Various | Protein HGVS notation |
| `hgvs_transcript` | Various | Transcript HGVS notation |
| `polyphen2_prediction` | dbNSFP | PolyPhen-2 prediction ("probably_damaging", "possibly_damaging", "benign") |
| `cadd_score` | CADD | CADD PHRED score (>20 = top 1% deleterious) |
| `gnomad_exome_af` | gnomAD | Population allele frequency |
| `alphamissense_score` | AlphaMissense | Pathogenicity score (0-1) |
| `alphamissense_prediction` | AlphaMissense | Prediction ("P" = pathogenic, "B" = benign, "A" = ambiguous) |
| `snpeff_effect` | SnpEff | Predicted variant effect |
| `transcript_id` | Various | Transcript identifier |
| `transcript_consequence` | Various | Consequence type (e.g., "missense_variant") |

### Fallback Behavior

When MyVariant returns no data:
1. Falls back to **direct CIViC GraphQL API** for clinical evidence
2. Falls back to **NCBI E-utilities** for ClinVar data
3. Uses **VEP** for functional predictions

---

## CIViC

**API Endpoint:** `https://civicdb.org/api/graphql` (GraphQL)

**Purpose:** Curated clinical interpretations with AMP/ASCO/CAP tier classifications.

### Data Retrieved

| Field | Description |
|-------|-------------|
| `assertion_id` | CIViC assertion ID |
| `name` | Assertion name |
| `amp_level` | AMP tier classification (e.g., "TIER_I_LEVEL_A") |
| `assertion_type` | Type: PREDICTIVE, PROGNOSTIC, DIAGNOSTIC, ONCOGENIC |
| `assertion_direction` | SUPPORTS or DOES_NOT_SUPPORT |
| `significance` | SENSITIVITYRESPONSE, RESISTANCE, ONCOGENIC, etc. |
| `status` | ACCEPTED, SUBMITTED, REJECTED |
| `molecular_profile` | Molecular profile name (e.g., "BRAF V600E") |
| `disease` | Associated disease |
| `therapies` | List of associated therapies |
| `fda_companion_test` | Whether FDA companion test exists |
| `nccn_guideline` | Associated NCCN guideline |
| `description` | Assertion description |

---

## VICC MetaKB

**API Endpoint:** `https://search.cancervariants.org/api/v1`

**Purpose:** Harmonized clinical interpretations aggregated from multiple knowledgebases.

### Sources Aggregated

- CIViC (Clinical Interpretations of Variants in Cancer)
- Cancer Genome Interpreter (CGI)
- JAX-CKB (Jackson Laboratory Clinical Knowledgebase)
- OncoKB
- PMKB (Precision Medicine Knowledgebase)
- MolecularMatch

### Data Retrieved

| Field | Description |
|-------|-------------|
| `description` | Clinical interpretation description |
| `gene` | Gene symbol |
| `variant` | Variant notation |
| `disease` | Associated disease |
| `drugs` | List of associated drugs |
| `evidence_level` | Evidence level (A, B, C, D) |
| `response_type` | Response type (Responsive, Resistant, 1A, 1B, etc.) |
| `source` | Original source database (civic, cgi, jax, oncokb, pmkb) |
| `publication_url` | Publication reference |
| `oncogenic` | Oncogenicity classification |

### Special Handling

- **KIT exon-level search**: For KIT variants, also searches by exon (9, 11, 13, 17) to find evidence catalogued as "KIT exon 11 deletion"
- **Compound mutation filtering**: Filters out resistance caused by secondary mutations, not the queried variant

---

## CGI

**API Endpoint:** `https://www.cancergenomeinterpreter.org/data/biomarkers/cgi_biomarkers_latest.tsv`

**Purpose:** FDA approval and guideline information from curated biomarker database.

### Data Retrieved

| Field | Description |
|-------|-------------|
| `gene` | Gene symbol |
| `alteration` | Alteration pattern (e.g., "EGFR:G719.", "KRAS:.12.") |
| `drug` | Drug name |
| `drug_status` | "Approved", "Clinical trial", etc. |
| `association` | "Responsive" or "Resistant" |
| `evidence_level` | "FDA guidelines", "NCCN guidelines", etc. |
| `source` | Evidence source |
| `tumor_type` | Primary tumor type abbreviation |
| `tumor_type_full` | Full tumor type name |

### Pattern Matching

CGI uses wildcard patterns:
- `EGFR:G719.` matches G719S, G719A, G719C, etc.
- `KRAS:.12.` matches any mutation at position 12 (G12D, G12V, etc.)
- `EGFR:.` matches any EGFR mutation

### Caching

- Downloads TSV file once and caches locally
- Cache expires after 7 days
- Falls back to stale cache if download fails

---

## FDA OpenFDA

**API Endpoint:** `https://api.fda.gov/drug/label.json`

**Purpose:** FDA-approved drug information and indications.

### Data Retrieved

| Field | Description |
|-------|-------------|
| `drug_name` | Brand or generic name |
| `brand_name` | Brand name |
| `generic_name` | Generic name |
| `indication` | Full indication text (truncated to 2500 chars) |
| `marketing_status` | "Prescription" |
| `variant_in_indications` | Whether variant explicitly mentioned in indications |
| `variant_in_clinical_studies` | Whether variant found in clinical studies section |

### Search Strategies

1. **Gene + variant search**: Full-text search for gene and variant
2. **Codon-level patterns**: For G719S, also searches "G719X" (FDA convention)
3. **Gene-only search**: Falls back to gene-only if no variant matches
4. **Disease-based search**: For KIT, also searches "GIST"; for MPL/JAK2/CALR, searches myelofibrosis
5. **MSI-H/dMMR search**: For MMR genes (MLH1, MSH2, MSH6, PMS2), searches microsatellite instability terms

### Special Handling

- **BRCA patterns**: Recognizes "BRCA-mutated", "gBRCAm", "deleterious BRCA"
- **Exclusion detection**: Filters out drugs where variant is mentioned in exclusion context (e.g., "without the D816V mutation")

---

## ClinicalTrials.gov

**API Endpoint:** `https://clinicaltrials.gov/api/v2/studies`

**Purpose:** Active clinical trials for investigational therapies.

### Data Retrieved

| Field | Description |
|-------|-------------|
| `nct_id` | NCT identifier |
| `title` | Brief trial title |
| `status` | RECRUITING, ENROLLING_BY_INVITATION, ACTIVE_NOT_RECRUITING, etc. |
| `phase` | Trial phase (PHASE1, PHASE2, PHASE3, etc.) |
| `conditions` | List of conditions |
| `interventions` | List of drug/intervention names |
| `brief_summary` | Trial summary |
| `eligibility_criteria` | Eligibility criteria text |
| `sponsor` | Lead sponsor name |
| `url` | ClinicalTrials.gov link |

### Filtering

- **Status filter**: Defaults to recruiting trials only
- **Variant mention filter**: Checks if variant appears in title, eligibility, or summary
- **Gene specificity**: Avoids false positives (e.g., KRAS G12D trial matching NRAS G12D query)

### Rate Limiting

- ~50 requests/minute per IP
- Exponential backoff with jitter on 429/403 errors

---

## Ensembl VEP

**API Endpoint:** `https://rest.ensembl.org/vep/human/hgvs`

**Purpose:** Variant normalization and functional predictions.

### Data Retrieved

| Field | Description |
|-------|-------------|
| `hgvs_genomic` | Genomic HGVS notation |
| `hgvs_transcript` | Transcript HGVS notation |
| `chromosome` | Chromosome |
| `position` | Genomic position |
| `ref_allele` | Reference allele |
| `alt_allele` | Alternate allele |
| `polyphen_prediction` | PolyPhen-2 prediction |
| `polyphen_score` | PolyPhen-2 score |
| `sift_prediction` | SIFT prediction ("tolerated", "deleterious") |
| `sift_score` | SIFT score |
| `cadd_phred` | CADD PHRED score |
| `alphamissense_prediction` | AlphaMissense prediction |
| `alphamissense_score` | AlphaMissense score |
| `consequence_terms` | List of consequence terms |
| `impact` | Impact level (HIGH, MODERATE, LOW, MODIFIER) |
| `transcript_id` | Canonical transcript ID |
| `myvariant_query` | HGVS notation for re-querying MyVariant |

### Variant Notation Handling

Converts various formats to HGVS protein notation:
- `E1978K` -> `ATM:p.Glu1978Lys`
- `p.E1978K` -> `ATM:p.Glu1978Lys`
- `E746_A750del` -> deletion notation
- `R348*` -> nonsense notation
- `W288fs` -> frameshift notation

---

## Semantic Scholar

**API Endpoint:** `https://api.semanticscholar.org/graph/v1`

**Purpose:** Literature enrichment with citation metrics and AI summaries.

### Data Retrieved

| Field | Description |
|-------|-------------|
| `paper_id` | Semantic Scholar paper ID |
| `pmid` | PubMed ID (if available) |
| `title` | Paper title |
| `abstract` | Abstract text |
| `citation_count` | Total citation count |
| `influential_citation_count` | Influential citations count |
| `reference_count` | Number of references |
| `year` | Publication year |
| `venue` | Journal/venue name |
| `is_open_access` | Open access status |
| `open_access_pdf_url` | PDF URL if available |
| `tldr` | AI-generated paper summary |
| `fields_of_study` | Classification fields |
| `publication_types` | Publication type (review, clinical trial, etc.) |

### Search Modes

1. **Resistance literature**: Searches for gene + variant + "resistance" + cancer context
2. **Variant literature**: General search for gene + variant

### Analysis Methods

- `mentions_resistance()` - Check for resistance terms
- `mentions_sensitivity()` - Check for sensitivity/response terms
- `get_signal_type()` - Classify as "resistance", "sensitivity", "mixed", or "unknown"
- `extract_drug_mentions()` - Extract mentioned drug names
- `get_impact_score()` - Calculate 0-1 impact score from citations

### Rate Limiting

- 1 request/second without API key
- Higher limits with API key

---

## PubMed

**API Endpoint:** `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`

**Purpose:** Primary literature search via NCBI E-utilities.

### Data Retrieved

| Field | Description |
|-------|-------------|
| `pmid` | PubMed ID |
| `title` | Article title |
| `abstract` | Abstract text |
| `authors` | Author list |
| `journal` | Journal name |
| `year` | Publication year |
| `doi` | DOI |
| `keywords` | Keywords + MeSH terms |
| `url` | PubMed link |

### Search Modes

1. **Resistance search**: `gene[Title/Abstract] AND variant[Title/Abstract] AND resistance[Title/Abstract]`
2. **General search**: `gene[Title/Abstract] AND variant[Title/Abstract]`

### Analysis Methods

- `mentions_resistance()` - Check for resistance terms
- `mentions_sensitivity()` - Check for sensitivity/response terms
- `get_signal_type()` - Classify paper focus
- `extract_drug_mentions()` - Extract targeted therapy drug names

### Rate Limiting

- 3 requests/second without API key
- 10 requests/second with API key

---

## OncoTree

**API Endpoint:** `https://oncotree.mskcc.org/api`

**Purpose:** Tumor type normalization and standardization.

### Data Retrieved

| Field | Description |
|-------|-------------|
| `code` | OncoTree code (e.g., "NSCLC", "LUAD", "MEL") |
| `name` | Full tumor type name |
| `tissue` | Tissue of origin |
| `mainType` | Main tumor type category |

### Usage

Primary method: `resolve_tumor_type(user_input)` - Resolves user input to standardized tumor type name:
- `"NSCLC"` -> `"Non-Small Cell Lung Cancer"`
- `"nsclc"` -> `"Non-Small Cell Lung Cancer"` (case-insensitive)
- `"NSCLC - Non-Small Cell Lung Cancer"` -> `"Non-Small Cell Lung Cancer"`

### Caching

- In-memory cache for API responses
- Avoids repeated API calls within same session

---

## cBioPortal

**API Endpoint:** `https://www.cbioportal.org/api`

**Purpose:** Variant prevalence, co-mutation patterns, and biological context from large-scale cancer genomics studies (TCGA, MSK-IMPACT, institutional cohorts).

### Data Retrieved

OncoMind returns a `CBioPortalEvidence` model with prevalence and co-mutation data:

#### Prevalence Fields

| Field | Type | Description |
|-------|------|-------------|
| `gene` | `str` | Gene symbol (e.g., "BRAF") |
| `variant` | `str \| None` | Variant notation (e.g., "V600E") |
| `tumor_type` | `str \| None` | Tumor type used for study selection |
| `study_id` | `str` | cBioPortal study ID (e.g., "skcm_tcga_pan_can_atlas_2018") |
| `total_samples` | `int` | Total samples in the study |
| `samples_with_gene_mutation` | `int` | Samples with any mutation in the gene |
| `samples_with_exact_variant` | `int` | Samples with the exact variant |
| `gene_prevalence_pct` | `float` | Gene mutation frequency (%) |
| `variant_prevalence_pct` | `float` | Exact variant frequency (%) |

**Example:** BRAF V600E in melanoma returns:
- `gene_prevalence_pct`: 51.2% (BRAF mutations overall)
- `variant_prevalence_pct`: 48.3% (V600E specifically)
- `samples_with_exact_variant`: 217
- `total_samples`: 449

#### Co-Mutation Fields

Each entry in `co_occurring` and `mutually_exclusive` lists contains:

| Field | Type | Description |
|-------|------|-------------|
| `gene` | `str` | Co-mutated gene symbol |
| `count` | `int` | Number of samples with both mutations |
| `pct` | `float` | Percentage of queried variant samples that also have this mutation |
| `odds_ratio` | `float \| None` | Odds ratio (>1 = co-occurrence, <1 = mutual exclusivity) |

**Example:** For BRAF V600E in melanoma:
```python
co_occurring = [
    {"gene": "CDKN2A", "count": 98, "pct": 45.2, "odds_ratio": 2.34},
    {"gene": "TP53", "count": 28, "pct": 12.9, "odds_ratio": 1.56},
]
mutually_exclusive = [
    {"gene": "NRAS", "count": 2, "pct": 0.9, "odds_ratio": 0.08},
]
```

### Study Selection

cBioPortal contains hundreds of studies. OncoMind maps tumor types to specific study IDs (TCGA PanCancer Atlas preferred for broad coverage):

| Tumor Type | Study IDs |
|------------|-----------|
| NSCLC/Lung | luad_tcga_pan_can_atlas_2018, lusc_tcga_pan_can_atlas_2018, nsclc_mskcc_2018 |
| SCLC | sclc_ucologne_2015 |
| Melanoma | skcm_tcga_pan_can_atlas_2018, mel_ucla_2016, skcm_mskcc_2014 |
| Colorectal | coadread_tcga_pan_can_atlas_2018, crc_msk_2017 |
| Breast | brca_tcga_pan_can_atlas_2018, breast_msk_2018 |
| Pancreatic | paad_tcga_pan_can_atlas_2018, paad_qcmg_uq_2016 |
| Glioblastoma | gbm_tcga_pan_can_atlas_2018 |
| Low-Grade Glioma | lgg_tcga_pan_can_atlas_2018 |
| Ovarian | ov_tcga_pan_can_atlas_2018 |
| Prostate | prad_tcga_pan_can_atlas_2018, prad_mskcc_2017 |
| Bladder/Urothelial | blca_tcga_pan_can_atlas_2018 |
| Kidney (Clear Cell) | kirc_tcga_pan_can_atlas_2018 |
| Kidney (Papillary) | kirp_tcga_pan_can_atlas_2018 |
| Thyroid | thca_tcga_pan_can_atlas_2018 |
| Head and Neck | hnsc_tcga_pan_can_atlas_2018 |
| Liver/HCC | lihc_tcga_pan_can_atlas_2018 |
| Gastric/Stomach | stad_tcga_pan_can_atlas_2018 |
| Esophageal | esca_tcga_pan_can_atlas_2018 |
| Endometrial/Uterine | ucec_tcga_pan_can_atlas_2018 |
| Cervical | cesc_tcga_pan_can_atlas_2018 |
| Cholangiocarcinoma | chol_tcga_pan_can_atlas_2018 |
| Sarcoma | sarc_tcga_pan_can_atlas_2018 |
| GIST | gist_mskcc |
| Mesothelioma | meso_tcga_pan_can_atlas_2018 |
| Adrenocortical | acc_tcga_pan_can_atlas_2018 |
| Pheochromocytoma | pcpg_tcga_pan_can_atlas_2018 |
| Testicular | tgct_tcga_pan_can_atlas_2018 |
| Thymoma | thym_tcga_pan_can_atlas_2018 |
| Uveal Melanoma | uvm_tcga_pan_can_atlas_2018 |
| AML | laml_tcga_pan_can_atlas_2018 |
| DLBCL | dlbc_tcga_pan_can_atlas_2018 |

If no tumor-specific study is found, falls back to MSK-IMPACT pan-cancer cohort (`msk_impact_2017`).

### Co-Mutation Calculation

OncoMind calculates co-occurrence statistics by:

1. Fetching all mutations for the query gene from the selected study
2. Building sample sets (samples with gene mutation, samples with exact variant)
3. Fetching mutations for ~35 curated cancer genes selected for biologically meaningful co-occurrence patterns:
   - **Tumor suppressors**: TP53, PTEN, CDKN2A, RB1, SMAD4, STK11, NF1, APC, ARID1A, KEAP1
   - **DDR genes**: ATM, BRCA1, BRCA2
   - **Oncogenes**: KRAS, NRAS, BRAF, PIK3CA, EGFR, ERBB2, MET, IDH1, IDH2, CTNNB1
   - **Oxidative stress**: NFE2L2 (co-mutates with KEAP1)
   - **Chromatin/epigenetic**: KMT2D, ARID2
   - **RTKs**: FGFR1/2/3, KIT, PDGFRA, ALK, ROS1, RET
4. Calculating odds ratios: `OR = (both × neither) / (only_query × only_other)`
5. Categorizing based on odds ratio and sample count:
   - **Co-occurring**: OR > 1.5 (appear together 50%+ more than expected) AND count ≥ 3 (minimum samples to avoid noise)
   - **Mutually exclusive**: OR < 0.5 (appear together less than half as expected) AND count ≤ 2 (rarely seen together)

### How cBioPortal Data Flows to LLM

The `CBioPortalEvidence.to_prompt_context()` method formats data for LLM consumption:

```
Cohort: 449 melanoma samples

PREVALENCE (copy the source citation when quoting these numbers):
  BRAF mutations: 51.2% (230/449) - cite as: ([cBioPortal: skcm_tcga](https://...))
  BRAF V600E specifically: 48.3% (217/449) - cite as: ([cBioPortal: skcm_tcga](https://...))

CO-OCCURRING MUTATIONS (cite source when quoting these percentages):
  - CDKN2A: 45.2% of BRAF-mutant samples (98 cases, OR=2.34) - cite as: ([cBioPortal: ...])
  - TP53: 12.9% of BRAF-mutant samples (28 cases, OR=1.56) - cite as: ([cBioPortal: ...])

MUTUALLY EXCLUSIVE MUTATIONS:
  - NRAS: 0.9% co-occurrence (2 cases, OR=0.08) - cite as: ([cBioPortal: ...])

INTERPRETATION:
  - BRAF mutations frequently co-occur with CDKN2A (45.2%)
  - BRAF and NRAS mutations are mutually exclusive (0.9%)
```

The LLM uses this to:
1. **Contextualize prevalence**: Common vs rare variant in this tumor type
2. **Identify co-mutation hypotheses**: Co-occurring mutations may suggest combination therapy targets or resistance mechanisms
3. **Flag mutual exclusivity**: Mutually exclusive mutations often act in the same pathway (e.g., BRAF/NRAS in MAPK)
4. **Generate research questions**: "Does CDKN2A co-deletion affect BRAF inhibitor response?"

### Evidence Gap Detection

cBioPortal data contributes to evidence quality:

| Scenario | Impact |
|----------|--------|
| Tumor-specific data available | Counts toward "well characterized" |
| Co-mutation patterns found | Supports biological context section |
| No tumor-specific data | Gap: "No {tumor_type}-specific prevalence data" |
| No co-mutation data | Gap: "Co-mutation patterns not characterized" |

### Rate Limiting

- No strict rate limits on cBioPortal public API
- OncoMind uses 30-second timeout per request
- Gene symbol → Entrez ID mappings are cached within session

---

## DepMap

**API Endpoint:** `https://depmap.org/portal/` (via DepMap Portal API)

**Purpose:** Preclinical research context from cancer cell line screens, including gene essentiality (CRISPR knockout screens) and drug sensitivity data (PRISM drug screens).

### Data Retrieved

#### Gene Dependency (CRISPR Screens)

| Field | Description |
|-------|-------------|
| `gene` | Gene symbol |
| `mean_dependency_score` | Mean CERES score across cell lines (negative = essential) |
| `n_dependent_lines` | Number of cell lines dependent on this gene |
| `n_total_lines` | Total cell lines screened |
| `dependency_pct` | Percentage of cell lines showing dependency |

**CERES Score Interpretation:**
- Score < -0.5: Gene is **ESSENTIAL** (knockout kills the cell)
- Score ≈ 0: Gene is not essential
- Example: BRAF has CERES score of -0.80, meaning it's essential in 45% of cancer cell lines

#### Drug Sensitivities (PRISM Screens)

| Field | Description |
|-------|-------------|
| `drug_name` | Drug name |
| `ic50_nm` | IC50 in nanomolar (lower = more sensitive) |
| `auc` | Area under dose-response curve |
| `n_cell_lines` | Number of cell lines tested |

**IC50 Interpretation:**
- IC50 < 100nM: Highly sensitive
- IC50 100-1000nM: Moderately sensitive
- IC50 > 1000nM: Resistant

#### Cell Line Models

| Field | Description |
|-------|-------------|
| `name` | Cell line name (e.g., "A375", "SK-MEL-28") |
| `primary_disease` | Cancer type (e.g., "Melanoma") |
| `subtype` | Disease subtype |
| `has_mutation` | Whether cell line has the queried mutation |
| `mutation_details` | Specific mutation (e.g., "V600E") |

### How DepMap Data Flows to LLM

DepMap data is sent to the LLM as part of the biological context, influencing the narrative synthesis:

```
GENE DEPENDENCY ([DepMap](https://depmap.org/)):
  BRAF is ESSENTIAL (CERES score: -0.80)
  Dependent in 450/1000 cell lines (45.0%)

DRUG SENSITIVITIES in BRAF-mutant lines ([DepMap](https://depmap.org/)):
  - trametinib: IC50=8.0nM (n=10)
  - cobimetinib: IC50=15.0nM (n=10)
  - encorafenib: IC50=20.0nM (n=10)
  - dabrafenib: IC50=30.0nM (n=10)
  - vemurafenib: IC50=50.0nM (n=10)

AVAILABLE MODEL CELL LINES with BRAF mutation ([DepMap](https://depmap.org/)):
  - A375 (Melanoma) [V600E]
  - SK-MEL-28 (Melanoma) [V600E]
  - COLO800 (Melanoma) [V600E]
  ...
```

The LLM uses this to:
1. **Validate functional significance**: Essential genes (CERES < -0.5) support the biological importance of the mutation
2. **Identify preclinical drug candidates**: IC50 values inform the "Preclinical" section of the therapeutic landscape
3. **Suggest model systems**: Cell lines with the mutation can be recommended for research

### Example LLM Output

When DepMap shows drug sensitivity data, the LLM incorporates it into the "Therapeutic Landscape" section:

> **Therapeutic Landscape:** FDA-approved: Mekinist (trametinib), BRAFTOVI (encorafenib), ZELBORAF (vemurafenib); **Preclinical: trametinib (IC50=8nM), cobimetinib (IC50=15nM), encorafenib (IC50=20nM), dabrafenib (IC50=30nM), vemurafenib (IC50=50nM)**

The drugs listed as "Preclinical" come directly from DepMap's PRISM drug screen data, ordered by IC50 (most sensitive first).

### Evidence Gap Detection

DepMap data contributes to evidence quality assessment:

| DepMap Data Present | Contribution to "Well Characterized" |
|---------------------|--------------------------------------|
| Gene dependency score | "gene essentiality (DepMap CRISPR)" |
| Gene is essential (CERES < -0.5) | "{gene} is essential in cancer cells" |
| Drug sensitivities | "preclinical drug sensitivity (DepMap)" |
| Cell line models with mutation | "model cell lines (N with mutation)" |

If DepMap data is missing for a clinically relevant variant, gaps are flagged:
- "No cell line models identified for {gene} {variant}"
- Suggested: "Identify cell lines with mutation", "Generate isogenic model"
- Addressable with: DepMap, CCLE, Cellosaurus

---

## Evidence Gap Detection

**Purpose:** Identify what's unknown or understudied about a variant to guide research priorities.

### Gap Categories

| Category | Description | Example Gap |
|----------|-------------|-------------|
| `FUNCTIONAL` | Mechanism unknown | "Functional impact of V657F on JAK1 protein is unknown" |
| `CLINICAL` | No clinical trials/outcomes | "No curated clinical evidence for JAK1 V657F" |
| `TUMOR_TYPE` | Not studied in this tumor type | "No evidence specific to T-ALL for JAK1 V657F" |
| `DRUG_RESPONSE` | No drug sensitivity data | "No drug sensitivity/resistance data" |
| `RESISTANCE` | Resistance mechanisms unknown | "Resistance mechanisms not well characterized" |
| `PRECLINICAL` | No cell line/model data | "No cell line models identified" |
| `PREVALENCE` | Frequency unknown | "Prevalence in T-ALL unknown" |
| `PROGNOSTIC` | Survival impact unknown | "Prognostic impact not characterized" |
| `DISCORDANT` | Conflicting evidence between sources | "Conflicting drug response for osimertinib" |
| `VALIDATION` | Strong oncogenic signal but lacks therapeutic validation | "Strong oncogenic signal but limited therapeutic validation" |

### Gap Severity Levels

| Severity | Weight | Meaning |
|----------|--------|---------|
| `CRITICAL` | 3.0 | No data at all in key area |
| `SIGNIFICANT` | 2.0 | Limited data, needs more research |
| `MINOR` | 1.0 | Some data exists but could be deeper |

### Well-Characterized Aspects

The gap detector also identifies what IS well-characterized, with structured basis:

```python
class CharacterizedAspect(BaseModel):
    aspect: str      # "clinical actionability"
    basis: str       # "3 FDA approvals + 5 CIViC assertions"
    category: GapCategory | None  # For grouped display
```

Example output:
```
well_characterized_detailed = [
    CharacterizedAspect(
        aspect="Clinical Actionability",
        basis="3 FDA approvals + 5 CIViC assertions",
        category=GapCategory.CLINICAL
    ),
    CharacterizedAspect(
        aspect="Computational Pathogenicity",
        basis="AlphaMissense=0.99 | CADD=32.5 | PolyPhen2=D",
        category=GapCategory.FUNCTIONAL
    ),
]
```

### Gap Detection Checks

The `detect_evidence_gaps()` function ([gap_detector.py](src/oncomind/insight_builder/gap_detector.py)) runs these checks:

1. **Hotspot context** - Is variant at or near a known cancer hotspot?
2. **Functional predictions** - AlphaMissense, CADD, PolyPhen2 available?
3. **Gene mechanism** - Gene role, pathway, DepMap essentiality known?
4. **Clinical evidence** - FDA approvals, CIViC assertions present?
5. **Tumor-type evidence** - Evidence specific to query tumor type?
6. **Drug response** - CGI, VICC, preclinical biomarkers present?
7. **Resistance mechanisms** - Resistance signals from any source?
8. **Discordant evidence** - Cross-source conflicts detected?
9. **Prevalence** - cBioPortal data available?
10. **Clinical trials** - Active trials found?
11. **Preclinical models** - Cell lines with mutation in DepMap?
12. **Literature depth** - PubMed articles found (if searched)?
13. **Validation gap** - Strong oncogenic signal but no therapeutic validation?

### Context-Aware Suggestions

Gaps are enriched with dynamic, evidence-aware study suggestions:

```python
# If resistance gap + primary drug known:
"Test bypass mechanisms for vemurafenib resistance in BRAF-mutant models"
"ctDNA monitoring for BRAF V600E emergence under vemurafenib treatment"

# If preclinical gap + DepMap drugs known:
"Validate sensitivity to trametinib, cobimetinib in isogenic BRAF V600E models"

# If validation gap + strong co-occurrence:
"Investigate synthetic lethality with CDKN2A co-mutation"
"CRISPR screen in BRAF/CDKN2A double-mutant background"
```

### Overall Quality Scoring

Evidence quality uses net scoring (gaps minus well-characterized):

```python
gap_score = sum(
    GAP_CATEGORY_WEIGHTS[gap.category] * SEVERITY_MULTIPLIERS[gap.severity]
    for gap in gaps
)
positive_credit = well_characterized_count * 1.5
net_score = gap_score - positive_credit

# Thresholds
if net_score >= 12.0: return "minimal"
elif net_score >= 6.0: return "limited"
elif net_score >= 0.0: return "moderate"
else: return "comprehensive"
```

### Research Priority

Research priority considers gene importance and gap profile:

| Priority | Criteria |
|----------|----------|
| `very_high` | Strong oncogenic signal (pathogenic + essential) + biological gaps |
| `very_high` | Hotspot-adjacent variant with pathogenic signal + biological gaps |
| `high` | Cancer gene with critical gaps |
| `high` | Hotspot-adjacent variant in cancer gene |
| `medium` | Any critical gaps OR cancer gene with significant gaps |
| `low` | Comprehensive evidence with no significant gaps |

### EvidenceGaps Model

```python
class EvidenceGaps(BaseModel):
    gaps: list[EvidenceGap]
    overall_evidence_quality: str  # "comprehensive" | "moderate" | "limited" | "minimal"
    well_characterized: list[str]  # Simple strings (legacy)
    well_characterized_detailed: list[CharacterizedAspect]  # With basis
    poorly_characterized: list[str]
    research_priority: str  # "very_high" | "high" | "medium" | "low"

    def to_dict_for_llm(self) -> dict:
        """Format for LLM prompt with knowledge_gaps key."""

    def top_gaps(self, n: int = 3) -> list[EvidenceGap]:
        """Get N most important gaps sorted by severity."""
```

---

## LLM Research Synthesis

**Purpose:** Generate research-focused narrative synthesis calibrated to evidence quality.

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

### Paper Relevance Scoring

The LLM service also scores individual papers for relevance:

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

### Variant Knowledge Extraction

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

## Therapeutic Evidence Model

**Purpose:** Unified model for therapeutic evidence at all evidence levels (FDA → preclinical).

### TherapeuticEvidence Fields

```python
class TherapeuticEvidence(BaseModel):
    # Core fields (from legacy RecommendedTherapy)
    drug_name: str
    evidence_level: str | None   # "FDA-approved" | "Phase 3" | "Phase 2" | "Preclinical" | etc.
    approval_status: str | None  # "Approved in indication" | "Investigational" | etc.
    clinical_context: str | None # "first-line", "resistance setting", etc.

    # Research-focused fields
    response_type: str | None    # "Sensitivity" | "Resistance" | "Mixed"
    mechanism: str | None        # "Constitutive kinase activation"
    tumor_types_tested: list[str]
    cell_lines_tested: list[str]
    ic50_nm: float | None
    response_rate_pct: float | None

    # Source attribution
    source: str | None           # "CIViC" | "CGI" | "VICC" | "FDA" | "DepMap"
    source_url: str | None
    pmids: list[str]
    confidence: str              # "high" | "moderate" | "low"

    # Match specificity
    match_level: str | None      # "variant" | "codon" | "gene"
```

### Evidence Tier Ranking

```python
def get_evidence_tier(self) -> int:
    """
    1: FDA-approved
    2: Phase 3
    3: Phase 2
    4: Phase 1 / Case reports
    5: Preclinical
    6: In vitro / Computational
    7: Unknown
    """
```

### Getting Therapeutic Evidence

```python
# From Evidence model:
result.evidence.get_therapeutic_evidence(
    include_preclinical=True,  # Include DepMap/CGI preclinical
    max_results=20
) -> list[TherapeuticEvidence]

# Grouped by level:
result.evidence.get_therapeutic_evidence_by_level() -> {
    "fda_approved": [...],
    "clinical": [...],
    "preclinical": [...]
}

# Filtered by response:
result.evidence.get_resistance_evidence() -> list[TherapeuticEvidence]
result.evidence.get_sensitivity_evidence() -> list[TherapeuticEvidence]
```

### Source Priority

When building therapeutic evidence, sources are processed in order:

1. **FDA approvals** (Tier 1) - `source="FDA"`
2. **CIViC assertions** (Tier 1-2) - `source="CIViC"`
3. **CGI FDA-approved biomarkers** (Tier 1) - `source="CGI"`
4. **VICC evidence** (Tier 1-3) - `source="VICC ({original_source})"`
5. **CGI preclinical biomarkers** (Tier 5) - `source="CGI (preclinical)"`

Duplicates are removed by drug name (case-insensitive).

---

## Literature Knowledge Model

**Purpose:** Structured knowledge extracted from literature via LLM.

### LiteratureKnowledge Fields

```python
class LiteratureKnowledge(BaseModel):
    mutation_type: str           # "primary" | "secondary" | "both" | "unknown"
    is_prognostic_only: bool     # True if variant only prognostic, not predictive

    resistant_to: list[DrugResistance]
    sensitive_to: list[DrugSensitivity]

    clinical_significance: str
    evidence_level: str          # "FDA-approved" | "Phase 3" | ... | "None"
    references: list[str]        # PMIDs
    key_findings: list[str]
    confidence: float            # 0-1
```

### Predictive vs Prognostic Distinction

```python
class DrugResistance(BaseModel):
    drug: str
    evidence: str                # "in vitro" | "preclinical" | "clinical" | "FDA-labeled"
    mechanism: str | None
    is_predictive: bool          # True = affects drug selection, False = just prognostic

# Usage:
literature_knowledge.get_resistance_drugs(predictive_only=True)
literature_knowledge.is_resistance_marker(predictive_only=True)
```

---

## Evidence Flow Summary

```
User Query (gene, variant, tumor_type)
           │
           ├──► MyVariant.info ──► ClinVar, COSMIC, gnomAD, CADD, AlphaMissense
           │         │
           │         └──► VEP (fallback) ──► Functional predictions
           │
           ├──► CIViC GraphQL ──► AMP/ASCO/CAP assertions
           │
           ├──► VICC MetaKB ──► Harmonized evidence from OncoKB, CIViC, CGI, JAX
           │
           ├──► CGI Biomarkers ──► FDA/NCCN approval status
           │
           ├──► FDA OpenFDA ──► Drug approval indications
           │
           ├──► ClinicalTrials.gov ──► Active clinical trials
           │
           ├──► Semantic Scholar ──► Literature with citations, TLDR
           │
           ├──► PubMed ──► Research articles
           │
           ├──► OncoTree ──► Standardized tumor type
           │
           ├──► cBioPortal ──► Prevalence, co-mutations, biological context
           │
           └──► DepMap ──► Gene essentiality, drug sensitivity, cell line models
                   │
                   ▼
           ┌───────────────────────────────────────────┐
           │           Evidence Model                  │
           │  ├── identifiers, functional, context     │
           │  ├── clinical evidence lists              │
           │  ├── cbioportal_evidence, depmap_evidence │
           │  └── literature_knowledge                 │
           └───────────────────────────────────────────┘
                   │
                   ▼
           ┌───────────────────────────────────────────┐
           │         Gap Detection                     │
           │  detect_evidence_gaps(evidence) →         │
           │  EvidenceGaps with research_priority      │
           └───────────────────────────────────────────┘
                   │
                   ▼
           ┌───────────────────────────────────────────┐
           │      LLM Research Synthesis (optional)    │
           │  Calibrated narrative + hypotheses →      │
           │  LLMInsight with evidence_tags            │
           └───────────────────────────────────────────┘
                   │
                   ▼
           ┌───────────────────────────────────────────┐
           │              Result                       │
           │  evidence: Evidence                       │
           │  llm: LLMInsight | None                   │
           └───────────────────────────────────────────┘
```

# OncoMind Evidence Sources

This document describes each external data source OncoMind fetches from, what data is retrieved, and how it's used in variant analysis.

---

## Table of Contents

1. [MyVariant.info](#myvariantinfo)
2. [CIViC (Clinical Interpretation of Variants in Cancer)](#civic)
3. [VICC MetaKB](#vicc-metakb)
4. [CGI (Cancer Genome Interpreter)](#cgi)
5. [FDA OpenFDA](#fda-openfda)
6. [ClinicalTrials.gov](#clinicaltrialsgov)
7. [Ensembl VEP](#ensembl-vep)
8. [Semantic Scholar](#semantic-scholar)
9. [PubMed](#pubmed)
10. [OncoTree](#oncotree)
11. [cBioPortal](#cbioportal)
12. [DepMap (Cancer Dependency Map)](#depmap)

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

cBioPortal contains hundreds of studies. OncoMind maps tumor types to specific study IDs:

| Tumor Type | Study IDs |
|------------|-----------|
| Melanoma | skcm_tcga_pan_can_atlas_2018, mel_ucla_2016, skcm_mskcc_2014 |
| NSCLC | luad_tcga_pan_can_atlas_2018, lusc_tcga_pan_can_atlas_2018, nsclc_tcga_broad_2016 |
| Colorectal | coadread_tcga_pan_can_atlas_2018, crc_msk_2017 |
| Breast | brca_tcga_pan_can_atlas_2018, breast_msk_2018 |
| Pancreatic | paad_tcga_pan_can_atlas_2018 |
| GIST | gist_mskcc |

If no tumor-specific study is found, falls back to MSK-IMPACT pan-cancer cohort (`msk_impact_2017`).

### Co-Mutation Calculation

OncoMind calculates co-occurrence statistics by:

1. Fetching all mutations for the query gene from the selected study
2. Building sample sets (samples with gene mutation, samples with exact variant)
3. Fetching mutations for top 20 cancer genes (TP53, KRAS, NRAS, EGFR, PIK3CA, PTEN, etc.)
4. Calculating odds ratios: `OR = (both × neither) / (only_query × only_other)`
5. Categorizing: OR > 1.5 with count ≥ 3 → co-occurring; OR < 0.5 with count ≤ 2 → mutually exclusive

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
              Evidence Model ──► LLM Synthesis ──► Result
```

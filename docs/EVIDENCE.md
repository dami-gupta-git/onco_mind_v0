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

### Helper Methods

- `is_sensitivity()` - Check if represents sensitivity/response
- `is_resistance()` - Check if represents resistance
- `get_amp_tier()` - Extract AMP tier (I, II, III, IV)
- `get_amp_level()` - Extract level letter (A, B, C, D)

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
           └──► OncoTree ──► Standardized tumor type
                   │
                   ▼
              Evidence Model ──► LLM Synthesis ──► Result
```

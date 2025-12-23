# OncoMind API Reference

Complete reference for OncoMind's public API and data structures.

## Public API

### `get_insight()`

Primary async function for variant analysis.

```python
from oncomind import get_insight, InsightConfig

panel = await get_insight(
    variant="BRAF V600E",           # Required: gene + variant string
    tumor_type="Melanoma",          # Optional: tumor context
    config=InsightConfig(...)       # Optional: configuration
)
```

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `variant` | `str` | Required | Variant string (e.g., "BRAF V600E", "EGFR L858R") |
| `tumor_type` | `str \| None` | `None` | Tumor type for context-specific evidence |
| `config` | `InsightConfig` | Default | Configuration options |

**Returns:** `EvidencePanel`

### `get_insights()`

Batch processing for multiple variants.

```python
panels = await get_insights(
    variants=["BRAF V600E", "KRAS G12D", "EGFR L858R"],
    tumor_type="NSCLC",
    config=InsightConfig(...)
)
```

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `variants` | `list[str]` | Required | List of variant strings |
| `tumor_type` | `str \| None` | `None` | Tumor type (applied to all) |
| `config` | `InsightConfig` | Default | Configuration options |

**Returns:** `list[EvidencePanel]`

### `get_insight_sync()` / `get_insights_sync()`

Synchronous wrappers for non-async contexts.

```python
from oncomind import get_insight_sync

panel = get_insight_sync("BRAF V600E", tumor_type="Melanoma")
```

---

## InsightConfig

Configuration options for insight generation.

```python
from oncomind import InsightConfig

config = InsightConfig(
    enable_llm=True,
    enable_literature=True,
    llm_model="gpt-4o-mini",
    llm_temperature=0.1,
    max_literature_results=10,
    max_clinical_trials=5,
)
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `enable_llm` | `bool` | `False` | Enable LLM synthesis of literature |
| `enable_literature` | `bool` | `True` | Enable literature search |
| `llm_model` | `str` | `"gpt-4o-mini"` | LLM model for synthesis |
| `llm_temperature` | `float` | `0.1` | LLM temperature (0.0-1.0) |
| `max_literature_results` | `int` | `10` | Max papers to retrieve |
| `max_clinical_trials` | `int` | `5` | Max trials to retrieve |
| `max_civic_assertions` | `int` | `20` | Max CIViC assertions |

---

## EvidencePanel Structure

The `EvidencePanel` is the core output model containing all variant evidence.

```python
panel = await get_insight("BRAF V600E", tumor_type="Melanoma")
```

### `panel.identifiers`

Variant identification and normalization.

```python
panel.identifiers.gene              # "BRAF"
panel.identifiers.variant           # "V600E"
panel.identifiers.variant_input     # Original input string
panel.identifiers.cosmic_id         # "COSM476"
panel.identifiers.dbsnp_id          # "rs113488022"
panel.identifiers.hgvs_protein      # "p.Val600Glu"
panel.identifiers.hgvs_coding       # "c.1799T>A"
panel.identifiers.hgvs_genomic      # "NC_000007.14:g.140753336A>T"
panel.identifiers.transcript_id     # "ENST00000288602"
panel.identifiers.transcript_consequence  # "missense_variant"
panel.identifiers.chromosome        # "7"
panel.identifiers.position          # 140753336
panel.identifiers.ref               # "A"
panel.identifiers.alt               # "T"
```

### `panel.kb`

Knowledgebase evidence from curated sources.

```python
# CIViC curated assertions
panel.kb.civic_assertions           # list[CIViCAssertion]
for assertion in panel.kb.civic_assertions:
    assertion.id                    # "12"
    assertion.variant_name          # "V600E"
    assertion.disease               # "Melanoma"
    assertion.drugs                 # ["Dabrafenib", "Trametinib"]
    assertion.evidence_type         # "Predictive"
    assertion.clinical_significance # "Sensitivity"
    assertion.evidence_level        # "A"
    assertion.source                # "PMID:25399551"

# VICC MetaKB (aggregated from OncoKB, CIViC, MOAlmanac, etc.)
panel.kb.vicc                       # list[VICCAssociation]
for assoc in panel.kb.vicc:
    assoc.source                    # "oncokb"
    assoc.gene                      # "BRAF"
    assoc.variant                   # "V600E"
    assoc.disease                   # "Melanoma"
    assoc.drugs                     # ["Vemurafenib"]
    assoc.evidence_level            # "1"
    assoc.clinical_significance     # "Sensitive"

# CGI biomarkers
panel.kb.cgi_biomarkers             # list[CGIBiomarker]

# ClinVar clinical significance
panel.kb.clinvar                    # ClinVarEntry | None
panel.kb.clinvar.clinical_significance  # "Pathogenic"
panel.kb.clinvar.review_status     # "reviewed by expert panel"
panel.kb.clinvar.conditions        # ["Melanoma", "Colorectal cancer"]
```

### `panel.functional`

Computational predictions and population frequencies.

```python
# AlphaMissense (Google DeepMind)
panel.functional.alphamissense_score       # 0.9834 (0-1, higher = pathogenic)
panel.functional.alphamissense_prediction  # "P" (P=Pathogenic, B=Benign, A=Ambiguous)

# CADD (Combined Annotation Dependent Depletion)
panel.functional.cadd_score                # 32.0 (PHRED-scaled, >20 = top 1%)
panel.functional.cadd_raw                  # 6.23

# Other predictors
panel.functional.polyphen2_score           # 0.999
panel.functional.polyphen2_prediction      # "probably_damaging"
panel.functional.sift_score                # 0.001
panel.functional.sift_prediction           # "deleterious"
panel.functional.revel_score               # 0.92

# Population frequencies (gnomAD)
panel.functional.gnomad_exome_af           # 0.00001 (allele frequency)
panel.functional.gnomad_genome_af          # 0.000008

# Helper methods
panel.functional.is_predicted_pathogenic() # True if AlphaMissense = "P"
panel.functional.is_rare(threshold=0.01)   # True if gnomAD AF < threshold
panel.functional.get_prediction_summary()  # "Pathogenic (AM=0.98, CADD=32)"
```

### `panel.clinical`

Clinical context including FDA approvals and trials.

```python
# FDA-approved therapies
panel.clinical.fda_approvals        # list[FDAApproval]
for approval in panel.clinical.fda_approvals:
    approval.drug_name              # "Dabrafenib"
    approval.indication             # "BRAF V600E mutant melanoma"
    approval.approval_date          # "2013-05-29"
    approval.label_url              # FDA label URL

# Clinical trials
panel.clinical.clinical_trials      # list[ClinicalTrial]
for trial in panel.clinical.clinical_trials:
    trial.nct_id                    # "NCT04543188"
    trial.title                     # "Study of..."
    trial.status                    # "RECRUITING"
    trial.phase                     # "Phase 3"
    trial.interventions             # ["Encorafenib", "Binimetinib"]
    trial.url                       # ClinicalTrials.gov link

# Gene-level context
panel.clinical.gene_role            # "oncogene" | "tumor_suppressor" | "unknown"
panel.clinical.gene_summary         # Gene function summary

# Helper methods
panel.clinical.get_approved_drugs() # ["Dabrafenib", "Vemurafenib", ...]
panel.clinical.has_fda_approval()   # True
panel.clinical.get_recruiting_trials()  # Trials with status="RECRUITING"
```

### `panel.literature`

Literature evidence and LLM-synthesized insights.

```python
# Retrieved articles
panel.literature.pubmed_articles         # list[PubMedArticle]
for article in panel.literature.pubmed_articles:
    article.pmid                         # "22735384"
    article.title                        # "BRAF inhibitor resistance..."
    article.abstract                     # Full abstract text
    article.authors                      # ["Chapman PB", "Hauschild A", ...]
    article.journal                      # "N Engl J Med"
    article.year                         # "2012"
    article.url                          # PubMed link

# Semantic Scholar enrichment (if available)
panel.literature.semantic_papers         # list[SemanticPaperInfo]
for paper in panel.literature.semantic_papers:
    paper.citation_count                 # 1542
    paper.influential_citation_count     # 89
    paper.tldr                           # AI-generated summary

# LLM-synthesized knowledge (requires enable_llm=True)
panel.literature.literature_knowledge    # str | None
# Example: "BRAF V600E shows initial response to vemurafenib but
#           develops resistance in 6-12 months via NRAS mutations (20%),
#           MEK1/2 mutations (5-10%), or BRAF amplification..."

# Literature source
panel.literature.literature_source       # "semantic_scholar" | "pubmed" | None

# Helper methods
panel.literature.get_resistance_articles()    # Articles mentioning resistance
panel.literature.get_sensitivity_articles()   # Articles mentioning sensitivity
```

### `panel.meta`

Processing metadata and trust signals.

```python
# Source tracking
panel.meta.sources_queried          # ["CIViC", "VICC", "ClinVar", "COSMIC", ...]
panel.meta.sources_with_data        # ["CIViC", "VICC", "COSMIC"]
panel.meta.sources_failed           # ["ClinicalTrials.gov"]  # If any errored

# Conflicts between sources
panel.meta.conflicts                # list[str]
# Example: ["CIViC reports sensitivity, CGI reports resistance to Drug X"]

# Processing info
panel.meta.processing_time_ms       # 1234
panel.meta.timestamp                # "2024-01-15T10:30:00Z"

# Helper methods
panel.meta.has_conflicts()          # True if any cross-source disagreements
panel.meta.get_evidence_summary()   # "5 sources, 12 assertions, 2 conflicts"
```

---

## Output Methods

### `to_knowledge_header()`

Generate a dense, LLM-ready context block.

```python
header = panel.to_knowledge_header()
# Returns:
# "BRAF V600E in melanoma. Oncogenic driver via constitutive MAPK activation.
#  FDA-approved: dabrafenib + trametinib, vemurafenib + cobimetinib.
#  Resistance typical at 6-12 months via NRAS (20%), MEK1/2 (5-10%), or BRAF amp.
#  Sources: CIViC:assertion:12, FDA label, PMID:22735384"
```

### `model_dump()`

Export to dictionary (Pydantic v2).

```python
data = panel.model_dump()           # Full dict
data = panel.model_dump(mode="json")  # JSON-serializable dict
```

### `model_dump_json()`

Export to JSON string.

```python
json_str = panel.model_dump_json(indent=2)
```

---

## Variant Input Formats

OncoMind accepts flexible variant notation:

```python
# Standard formats
await get_insight("BRAF V600E")           # Gene + protein change
await get_insight("BRAF p.V600E")         # With p. prefix
await get_insight("BRAF p.Val600Glu")     # Three-letter amino acids

# HGVS notation
await get_insight("BRAF c.1799T>A")       # Coding DNA
await get_insight("NM_004333.4:c.1799T>A")  # With transcript

# Deletions/insertions
await get_insight("EGFR E746_A750del")    # Deletion
await get_insight("EGFR T790M")           # Point mutation
await get_insight("ERBB2 A775_G776insYVMA")  # Insertion

# Frameshift/nonsense
await get_insight("TP53 R248*")           # Nonsense (stop)
await get_insight("TP53 R248X")           # Nonsense (X notation)
await get_insight("APC K1462fs")          # Frameshift
```

---

## Error Handling

```python
from oncomind import get_insight
from oncomind.exceptions import VariantParseError, APIError

try:
    panel = await get_insight("INVALID_INPUT")
except VariantParseError as e:
    print(f"Could not parse variant: {e}")
except APIError as e:
    print(f"API error: {e}")
```

---

## CLI Reference

```bash
# Basic insight
mind insight BRAF V600E --tumor Melanoma

# With LLM synthesis
mind insight BRAF V600E --tumor Melanoma --llm

# Output to file
mind insight KRAS G12C --tumor NSCLC --output result.json

# Full LLM narrative mode
mind insight-llm BRAF V600E --tumor Melanoma

# Batch processing
mind batch variants.csv --output results.json

# Help
mind --help
mind insight --help
```

**CLI Options:**
| Option | Short | Description |
|--------|-------|-------------|
| `--tumor` | `-t` | Tumor type context |
| `--llm` | | Enable LLM synthesis |
| `--output` | `-o` | Output file path |
| `--format` | `-f` | Output format (json, csv, markdown) |

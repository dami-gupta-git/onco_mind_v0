# OncoMind API Reference

Complete reference for OncoMind's public API and data structures.

## Public API

### `get_insight()`

Primary async function for variant analysis.

```python
from oncomind import get_insight, InsightConfig, Result

result = await get_insight(
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

**Returns:** `Result`

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

**Returns:** `list[Result]`

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

## Result Structure

The `Result` is the core output model containing structured evidence and optional LLM insight.

```python
result = await get_insight("BRAF V600E", tumor_type="Melanoma")

# Result contains:
# - evidence: Evidence (structured data from databases)
# - llm: LLMInsight | None (LLM narrative when enabled)

# Property shortcuts on Result delegate to evidence:
result.identifiers  # → result.evidence.identifiers
result.kb           # → result.evidence.kb
result.functional   # → result.evidence.functional
result.clinical     # → result.evidence.clinical
result.literature   # → result.evidence.literature
```

### `result.identifiers`

Variant identification and normalization.

```python
result.identifiers.gene              # "BRAF"
result.identifiers.variant           # "V600E"
result.identifiers.variant_input     # Original input string
result.identifiers.cosmic_id         # "COSM476"
result.identifiers.dbsnp_id          # "rs113488022"
result.identifiers.hgvs_protein      # "p.Val600Glu"
result.identifiers.hgvs_coding       # "c.1799T>A"
result.identifiers.hgvs_genomic      # "NC_000007.14:g.140753336A>T"
result.identifiers.transcript_id     # "ENST00000288602"
result.identifiers.transcript_consequence  # "missense_variant"
result.identifiers.chromosome        # "7"
result.identifiers.position          # 140753336
result.identifiers.ref               # "A"
result.identifiers.alt               # "T"
```

### `result.kb`

Knowledgebase evidence from curated sources.

```python
# CIViC curated assertions
result.kb.civic_assertions           # list[CIViCAssertion]
for assertion in result.kb.civic_assertions:
    assertion.id                    # "12"
    assertion.variant_name          # "V600E"
    assertion.disease               # "Melanoma"
    assertion.drugs                 # ["Dabrafenib", "Trametinib"]
    assertion.evidence_type         # "Predictive"
    assertion.clinical_significance # "Sensitivity"
    assertion.evidence_level        # "A"
    assertion.source                # "PMID:25399551"

# VICC MetaKB (aggregated from OncoKB, CIViC, MOAlmanac, etc.)
result.kb.vicc                       # list[VICCAssociation]
for assoc in result.kb.vicc:
    assoc.source                    # "oncokb"
    assoc.gene                      # "BRAF"
    assoc.variant                   # "V600E"
    assoc.disease                   # "Melanoma"
    assoc.drugs                     # ["Vemurafenib"]
    assoc.evidence_level            # "1"
    assoc.clinical_significance     # "Sensitive"

# CGI biomarkers
result.kb.cgi_biomarkers             # list[CGIBiomarker]

# ClinVar clinical significance
result.kb.clinvar                    # ClinVarEntry | None
result.kb.clinvar.clinical_significance  # "Pathogenic"
result.kb.clinvar.review_status     # "reviewed by expert panel"
result.kb.clinvar.conditions        # ["Melanoma", "Colorectal cancer"]
```

### `result.functional`

Computational predictions and population frequencies.

```python
# AlphaMissense (Google DeepMind)
result.functional.alphamissense_score       # 0.9834 (0-1, higher = pathogenic)
result.functional.alphamissense_prediction  # "P" (P=Pathogenic, B=Benign, A=Ambiguous)

# CADD (Combined Annotation Dependent Depletion)
result.functional.cadd_score                # 32.0 (PHRED-scaled, >20 = top 1%)
result.functional.cadd_raw                  # 6.23

# Other predictors
result.functional.polyphen2_score           # 0.999
result.functional.polyphen2_prediction      # "probably_damaging"
result.functional.sift_score                # 0.001
result.functional.sift_prediction           # "deleterious"
result.functional.revel_score               # 0.92

# Population frequencies (gnomAD)
result.functional.gnomad_exome_af           # 0.00001 (allele frequency)
result.functional.gnomad_genome_af          # 0.000008

# Helper methods
result.functional.is_predicted_pathogenic() # True if AlphaMissense = "P"
result.functional.is_rare(threshold=0.01)   # True if gnomAD AF < threshold
result.functional.get_prediction_summary()  # "Pathogenic (AM=0.98, CADD=32)"
```

### `result.clinical`

Clinical context including FDA approvals and trials.

```python
# FDA-approved therapies
result.clinical.fda_approvals        # list[FDAApproval]
for approval in result.clinical.fda_approvals:
    approval.drug_name              # "Dabrafenib"
    approval.indication             # "BRAF V600E mutant melanoma"
    approval.approval_date          # "2013-05-29"
    approval.label_url              # FDA label URL

# Clinical trials
result.clinical.clinical_trials      # list[ClinicalTrial]
for trial in result.clinical.clinical_trials:
    trial.nct_id                    # "NCT04543188"
    trial.title                     # "Study of..."
    trial.status                    # "RECRUITING"
    trial.phase                     # "Phase 3"
    trial.interventions             # ["Encorafenib", "Binimetinib"]
    trial.url                       # ClinicalTrials.gov link

# Gene-level context
result.clinical.gene_role            # "oncogene" | "tumor_suppressor" | "unknown"
result.clinical.gene_summary         # Gene function summary

# Helper methods
result.clinical.get_approved_drugs() # ["Dabrafenib", "Vemurafenib", ...]
result.clinical.has_fda_approval()   # True
result.clinical.get_recruiting_trials()  # Trials with status="RECRUITING"
```

### `result.literature`

Literature evidence and LLM-synthesized insights.

```python
# Retrieved articles
result.literature.pubmed_articles         # list[PubMedArticle]
for article in result.literature.pubmed_articles:
    article.pmid                         # "22735384"
    article.title                        # "BRAF inhibitor resistance..."
    article.abstract                     # Full abstract text
    article.authors                      # ["Chapman PB", "Hauschild A", ...]
    article.journal                      # "N Engl J Med"
    article.year                         # "2012"
    article.url                          # PubMed link

# Semantic Scholar enrichment (if available)
result.literature.semantic_papers         # list[SemanticPaperInfo]
for paper in result.literature.semantic_papers:
    paper.citation_count                 # 1542
    paper.influential_citation_count     # 89
    paper.tldr                           # AI-generated summary

# LLM-synthesized knowledge (requires enable_llm=True)
result.literature.literature_knowledge    # str | None
# Example: "BRAF V600E shows initial response to vemurafenib but
#           develops resistance in 6-12 months via NRAS mutations (20%),
#           MEK1/2 mutations (5-10%), or BRAF amplification..."

# Literature source
result.literature.literature_source       # "semantic_scholar" | "pubmed" | None

# Helper methods
result.literature.get_resistance_articles()    # Articles mentioning resistance
result.literature.get_sensitivity_articles()   # Articles mentioning sensitivity
```

### `result.llm`

LLM-generated research-focused narrative (when `enable_llm=True`).

```python
# LLM insight is None when LLM is disabled
if result.llm:
    # Core summary
    result.llm.llm_summary              # Research-focused narrative summary
    result.llm.rationale                # Research implications and reasoning

    # Therapeutic evidence
    result.llm.therapeutic_evidence     # list[TherapeuticEvidence]
    result.llm.recommended_therapies    # Alias for therapeutic_evidence
    result.llm.clinical_trials_available # True if trials exist

    # Evidence assessment (from structured gap detection)
    result.llm.evidence_quality         # "comprehensive" | "moderate" | "limited" | "minimal"
    result.llm.well_characterized       # list[str] - aspects with strong evidence
    result.llm.knowledge_gaps           # list[str] - identified gaps
    result.llm.conflicting_evidence     # list[str] - areas where sources disagree

    # Research context
    result.llm.research_implications    # Future research directions
    result.llm.evidence_tags            # list[str] - e.g., "direct clinical data", "preclinical only"
    result.llm.references               # list[str] - cited sources
```

---

## Evidence Gaps Detection

OncoMind provides structured evidence gap detection to identify what's missing or understudied about a variant.

### `EvidenceGaps`

The `EvidenceGaps` model aggregates detected gaps with severity ratings and research recommendations.

```python
from oncomind.models.evidence import EvidenceGaps, GapCategory, GapSeverity

# Access via Evidence model
evidence = result.evidence
gaps = evidence.compute_evidence_gaps()

# Overall assessment
gaps.overall_evidence_quality    # "comprehensive" | "moderate" | "limited" | "minimal"
gaps.research_priority           # "high" | "medium" | "low"
gaps.well_characterized          # ["clinical actionability", "published literature"]
gaps.poorly_characterized        # ["resistance mechanisms", "tumor-specific data"]

# Individual gaps
for gap in gaps.gaps:
    gap.category                 # GapCategory enum
    gap.severity                 # GapSeverity enum
    gap.description              # Human-readable description
    gap.suggested_studies        # ["Case series", "Retrospective cohort", ...]
    gap.addressable_with         # ["CIViC", "Literature search", ...]

# Helper methods
gaps.has_critical_gaps()         # True if any critical gaps
gaps.get_gaps_by_category(GapCategory.CLINICAL)
gaps.get_gaps_by_severity(GapSeverity.CRITICAL)
gaps.to_summary()                # Human-readable text summary
gaps.to_dict_for_llm()           # Dict optimized for LLM prompts
```

### `GapCategory`

Categories of evidence gaps:

| Category | Description |
|----------|-------------|
| `FUNCTIONAL` | Mechanism of variant unknown |
| `CLINICAL` | No clinical trials/outcomes data |
| `TUMOR_TYPE` | Not studied in this specific tumor type |
| `DRUG_RESPONSE` | No drug sensitivity/resistance data |
| `RESISTANCE` | Resistance mechanisms unknown |
| `PRECLINICAL` | No cell line/model data |
| `PREVALENCE` | Frequency in population unknown |
| `PROGNOSTIC` | Survival impact unknown |

### `GapSeverity`

Severity levels for gaps:

| Severity | Description |
|----------|-------------|
| `CRITICAL` | No data at all in key area |
| `SIGNIFICANT` | Limited data, needs more research |
| `MINOR` | Some data exists but could be deeper |

### Gap Detection Logic

The `detect_evidence_gaps()` function analyzes evidence across multiple dimensions:

```python
from oncomind.insight_builder.gap_detector import detect_evidence_gaps

# Checks performed:
# 1. Functional characterization (AlphaMissense, CADD, PolyPhen2)
# 2. Mechanism/functional studies (gene role, pathway)
# 3. Clinical evidence (CIViC assertions, FDA approvals)
# 4. Tumor-type-specific evidence
# 5. Drug response data (CGI, VICC, preclinical)
# 6. Resistance mechanisms
# 7. Prevalence/epidemiology (cBioPortal)
# 8. Clinical trials
# 9. Literature depth (publication count)

gaps = detect_evidence_gaps(evidence)
```

### Evidence Quality Computation

Overall evidence quality is computed from gap severity:

| Quality | Criteria |
|---------|----------|
| `minimal` | 2+ critical gaps |
| `limited` | 1 critical gap OR 2+ significant gaps |
| `moderate` | 1 significant gap OR only minor gaps |
| `comprehensive` | No gaps detected |

### Research Priority

Priority is computed based on gene importance and gaps:

```python
# High priority: clinically important gene with critical gaps
# Medium priority: critical gaps OR important gene with significant gaps
# Low priority: only minor gaps
```

### LLM Integration

Evidence gaps are passed to the LLM as a structured assessment dict:

```python
assessment = gaps.to_dict_for_llm()
# Returns:
# {
#     "overall_quality": "moderate",
#     "research_priority": "high",
#     "well_characterized": ["clinical actionability", ...],
#     "knowledge_gaps": ["resistance mechanisms", ...],
#     "conflicting_evidence": [],
#     "critical_gaps": [
#         {"description": "...", "suggested_studies": [...]}
#     ],
#     "significant_gaps": [
#         {"description": "...", "suggested_studies": [...]}
#     ]
# }
```

The LLM uses this to:
- Calibrate response verbosity (brief for minimal evidence)
- Echo well-characterized vs gap areas accurately
- Generate research implications based on suggested studies

---

## Output Methods

### `get_summary()`

Get a one-line summary of the variant evidence.

```python
summary = result.get_summary()
# Returns: "BRAF V600E: 5 sources, FDA-approved therapies available"
```

### `has_evidence()`

Check if any evidence was found.

```python
if result.has_evidence():
    print("Evidence found from databases")
```

### `model_dump()`

Export to dictionary (Pydantic v2).

```python
data = result.model_dump()           # Full dict
data = result.model_dump(mode="json")  # JSON-serializable dict

# JSON structure:
# {
#   "evidence": {
#     "identifiers": {...},
#     "kb": {...},
#     "functional": {...},
#     "clinical": {...},
#     "literature": {...}
#   },
#   "llm": {...} | null
# }
```

### `model_dump_json()`

Export to JSON string.

```python
json_str = result.model_dump_json(indent=2)
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
    result = await get_insight("INVALID_INPUT")
except VariantParseError as e:
    print(f"Could not parse variant: {e}")
except APIError as e:
    print(f"API error: {e}")
```

---

## CLI Reference

### Main Command

```bash
mind insight GENE VARIANT [OPTIONS]
```

### CLI Modes

| Mode | Flag | Speed | Output |
|------|------|-------|--------|
| **Default** | (none) | ~12s | Structured evidence + LLM clinical summary |
| **Lite** | `--lite` | ~7s | Structured evidence only (no LLM) |
| **Full** | `--full` | ~25s | + Literature search + Literature panel + enhanced narrative |

### Output Panels by Mode

| Panel | Lite | Default | Full |
|-------|------|---------|------|
| Evidence Summary | ✓ | ✓ | ✓ |
| Recommended Therapies | ✓ (FDA) | ✓ (LLM) | ✓ (LLM) |
| Clinical Evidence | ✓ | ✓ | ✓ |
| Literature | - | - | ✓ |
| Variant Insight | - | ✓ | ✓ |

### Examples

```bash
# Default: structured evidence + LLM narrative (~12s)
mind insight BRAF V600E --tumor Melanoma
mind insight PIK3CA H1047R -t "Breast Cancer"

# Lite mode: structured evidence only, no LLM (~7s)
mind insight EGFR L858R -t NSCLC --lite

# Full mode: + literature search + enhanced narrative (~25s)
mind insight KRAS G12C -t NSCLC --full

# Save to JSON
mind insight BRAF V600E -t Melanoma --output result.json

# Batch processing
mind batch variants.json --output results.json
mind batch variants.json --lite              # Fastest: no LLM
mind batch variants.json --full              # Slowest: with literature

# Help
mind --help
mind insight --help
```

### CLI Options

| Option | Short | Description |
|--------|-------|-------------|
| `--tumor` | `-t` | Tumor type context |
| `--lite` | | Lite mode: structured evidence only, no LLM |
| `--full` | | Full mode: include literature search + enhanced narrative |
| `--model` | `-m` | LLM model (default: gpt-4o-mini) |
| `--output` | `-o` | Output file path (JSON) |

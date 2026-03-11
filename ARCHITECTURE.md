# Architecture

OncoMind is a two-layer system: a **deterministic evidence backbone** that aggregates data from 12+ cancer databases, and an **optional LLM synthesis layer** that generates research dossiers grounded in that evidence.

## System Overview

```
User Input (gene + variant + tumor type)
         |
    [Normalization]  -- input_parser.py, variant_normalization.py
         |
    [Conductor]  -- orchestrator (insight_builder/conductor.py)
         |
    +----+----+
    |         |
    v         v
[EvidenceAggregator]    [LLMService] (optional)
    |                       |
    v                       v
 Evidence              LLMInsight
 (structured)          (narrative)
    |                       |
    +-----------+-----------+
                |
             Result
```

## Entry Points

| Entry Point | Location | Description |
|-------------|----------|-------------|
| CLI (`mind`) | `src/oncomind/cli.py` | Typer-based CLI with `insight`, `batch`, `annotate`, `version` commands |
| Streamlit UI | `streamlit/app.py` | Interactive web interface via `streamlit/backend.py` |
| Python API | `oncomind.insight_builder.Conductor` | Async context manager for programmatic use |

## Core Modules

### Conductor (`insight_builder/conductor.py`)

Orchestrates the full pipeline. Accepts `ConductorConfig` to toggle evidence sources, LLM, literature, and processing options. Runs evidence aggregation first, then gap computation and LLM calls in parallel.

### EvidenceAggregator (`insight_builder/evidence_aggregator.py`)

Fetches evidence from all API clients in parallel via `asyncio.gather`. Handles per-source failures gracefully (logs warnings, continues). Assembles results into strongly-typed `Evidence` models. No LLM dependency.

### LLM Service (`llm/service.py`)

All LLM calls go through this service. Uses `litellm` for multi-provider support. Three LLM operations:
- **Synthesis** -- Research dossier generation (5-section structured output)
- **Paper relevance scoring** -- Scores individual papers against variant context (uses fast model)
- **Variant knowledge extraction** -- Extracts structured knowledge from paper abstracts (uses fast model)

Prompts live in `llm/prompts.py`. LLM output model is `LLMInsight` in `llm/llm_insight.py`.

### Normalization (`normalization/`)

Parses free-text variant inputs (e.g., "BRAF V600E") into `ParsedVariant`. Handles protein notation, HGVS, variant type classification, and amino acid code conversion.

### Annotator (`annotator/`)

VCF file annotation pipeline. Parses VCF records and fetches annotations from MyVariant.info. Separate from the insight pipeline.

## API Clients (`api/`)

All external data fetching is isolated in individual client modules:

| Client | Source | Data |
|--------|--------|------|
| `myvariant.py` | MyVariant.info | ClinVar, COSMIC, gnomAD, AlphaMissense, CADD |
| `civic.py` | CIViC | Curated variant-drug assertions |
| `vicc.py` | VICC MetaKB | Aggregated knowledge bases |
| `cgi.py` | CGI | Biomarker annotations (local TSV) |
| `fda_drugs.py` / `fda_label_parser.py` | FDA | Drug label biomarker evidence |
| `cbioportal.py` | cBioPortal | Prevalence, co-mutations |
| `depmap.py` | DepMap | Gene essentiality (CRISPR), drug sensitivity (PRISM), cell line models |
| `hotspots.py` | Cancer Hotspots | Known hotspot positions |
| `clinicaltrials.py` | ClinicalTrials.gov | Active and historical trials |
| `pubmed.py` | PubMed | Literature search |
| `semantic_scholar.py` | Semantic Scholar | AI-powered literature search |
| `oncotree.py` | OncoTree | Tumor type ontology |

## Data Models (`models/`)

- `Evidence` (`models/evidence/evidence.py`) -- Container for all structured evidence: identifiers, functional scores, variant context, and per-source evidence lists
- `Result` (`models/result.py`) -- Wraps `Evidence` + optional `LLMInsight` + optional cross-source analysis
- `LLMInsight` (`llm/llm_insight.py`) -- LLM-generated research dossier with functional summary, biological context, therapeutic landscape, research program
- Evidence sub-models: one per source (CIViC, VICC, CGI, FDA, DepMap, cBioPortal, Hotspots, ClinVar, COSMIC, PubMed, ClinicalTrials)
- `EvidenceGaps` (`models/extracted/evidence_gaps.py`) -- Gap detection with severity levels (CRITICAL, HIGH, SIGNIFICANT, MODERATE, INFORMATIONAL)

## Configuration

- `config/constants.py` -- All hardcoded mappings, thresholds, API endpoints, LLM parameters
- `config/debug.py` -- Logging configuration (ONCOMIND_LOG_LEVEL env var or --log-level flag)
- `.env` -- API keys (GOOGLE_API_KEY, ANTHROPIC_API_KEY, OPENAI_API_KEY, SEMANTIC_SCHOLAR_API_KEY)

## Key Design Patterns

- **Async throughout** -- All API calls and orchestration use async/await
- **Parallel fetching** -- Evidence sources are queried concurrently via `asyncio.gather`
- **Graceful degradation** -- Individual source failures do not block the pipeline
- **Separation of concerns** -- Evidence aggregation is LLM-independent; LLM synthesis is a separate optional layer
- **Match specificity tracking** -- Evidence is tagged as variant-level, codon-level, or gene-level
- **Pydantic models** -- All data structures use Pydantic for validation and serialization

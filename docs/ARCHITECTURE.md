# OncoMind Architecture

This document describes the architecture of OncoMind, an AI-powered cancer variant annotation and evidence synthesis tool.

## Overview

OncoMind follows a layered architecture that separates concerns into distinct modules:

```
┌───────────────────────────────────────────────────────────────────┐
│                       User Interfaces                             │
│  ┌─────────────┐  ┌─────────────┐  ┌───────────────────────────┐  │
│  │    CLI      │  │  Streamlit  │  │      Python API           │  │
│  │   (mind)    │  │    App      │  │   (get_insight)           │  │
│  └──────┬──────┘  └──────┬──────┘  └─────────────┬─────────────┘  │
└─────────┼────────────────┼───────────────────────┼────────────────┘
          │                │                       │
          ▼                ▼                       ▼
┌───────────────────────────────────────────────────────────────────┐
│                      Public API Layer                             │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │  api_public/insight.py                                      │  │
│  │  - get_insight(variant_str, tumor_type, config) → Insight   │  │
│  │  - get_insights(variants, tumor_type, config) → [Insight]   │  │
│  │  - InsightConfig                                            │  │
│  └─────────────────────────────────────────────────────────────┘  │
└───────────────────────────────────────────────────────────────────┘
          │
          ▼
┌───────────────────────────────────────────────────────────────────┐
│                     Normalization Layer                           │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │  normalization/                                             │  │
│  │  - input_parser.py: parse_variant_input("BRAF V600E")       │  │
│  │  - hgvs_utils.py: normalize_variant, classify_variant_type  │  │
│  │  → Output: ParsedVariant                                    │  │
│  └─────────────────────────────────────────────────────────────┘  │
└───────────────────────────────────────────────────────────────────┘
          │
          ▼
┌───────────────────────────────────────────────────────────────────┐
│                   Evidence Aggregation Layer                      │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │  evidence/builder.py                                        │  │
│  │  - EvidenceBuilder (async context manager)                  │  │
│  │  - build_insight(parsed_variant, tumor_type) → Insight      │  │
│  │  - Parallel API fetching with asyncio.gather()              │  │
│  └─────────────────────────────────────────────────────────────┘  │
└───────────────────────────────────────────────────────────────────┘
          │
          ├────────────────────────────────────┐
          ▼                                    ▼
┌─────────────────────────────┐  ┌─────────────────────────────────┐
│     API Clients Layer       │  │      LLM Layer (Optional)       │
│  ┌───────────────────────┐  │  │  ┌───────────────────────────┐  │
│  │ api/myvariant.py      │  │  │  │ llm/service.py            │  │
│  │ api/civic.py          │  │  │  │ - get_llm_insight()       │  │
│  │ api/vicc.py           │  │  │  │ - score_paper_relevance   │  │
│  │ api/fda.py            │  │  │  │ - extract_variant_knowledge│ │
│  │ api/cgi.py            │  │  │  └───────────────────────────┘  │
│  │ api/pubmed.py         │  │  └─────────────────────────────────┘
│  │ api/semantic_scholar.py│ │
│  │ api/clinicaltrials.py │  │
│  │ api/oncotree.py       │  │
│  └───────────────────────┘  │
└─────────────────────────────┘
          │
          ▼
┌───────────────────────────────────────────────────────────────────┐
│                        Models Layer                               │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │  models/evidence/insight.py                                 │  │
│  │  - Insight (top-level output)                               │  │
│  │    ├── identifiers: VariantIdentifiers                      │  │
│  │    ├── kb: KnowledgebaseEvidence                            │  │
│  │    ├── functional: FunctionalScores                         │  │
│  │    ├── clinical: ClinicalContext                            │  │
│  │    ├── literature: LiteratureEvidence                       │  │
│  │    └── llm: LLMInsight | None  ← Optional LLM narrative     │  │
│  └─────────────────────────────────────────────────────────────┘  │
└───────────────────────────────────────────────────────────────────┘
```

## Core Design Principles

### 1. Separation of Concerns

- **Normalization** is separate from evidence fetching
- **Evidence aggregation** is separate from LLM synthesis
- **API clients** are independent and composable
- **Models** are strongly-typed with Pydantic

### 2. Unified Output Model

All code paths return a single `Insight` object that contains:
- Structured evidence from databases (always populated)
- Optional LLM-generated narrative (when LLM mode is enabled)

```python
# The Insight model embeds LLM narrative when available
class Insight(BaseModel):
    identifiers: VariantIdentifiers    # Gene, variant, IDs
    kb: KnowledgebaseEvidence          # CIViC, ClinVar, VICC, etc.
    functional: FunctionalScores       # AlphaMissense, CADD, etc.
    clinical: ClinicalContext          # FDA, trials, gene role
    literature: LiteratureEvidence     # PubMed, extracted knowledge
    llm: LLMInsight | None             # LLM narrative (when enabled)
```

### 3. LLM as Optional Enhancement

The core annotation pipeline is deterministic. LLM can be enabled for clinical narrative synthesis:

```python
# Default: structured evidence + LLM narrative (~12s)
insight = await get_insight("BRAF V600E", tumor_type="Melanoma")
print(insight.llm.llm_summary)  # LLM narrative

# Lite: fast, no LLM (~7s)
config = InsightConfig(enable_llm=False)
insight = await get_insight("BRAF V600E", config=config)
print(insight.llm)  # None

# Full: + literature search + enhanced narrative (~25s)
config = InsightConfig(enable_llm=True, enable_literature=True)
insight = await get_insight("BRAF V600E", config=config)
```

**CLI equivalent:**
```bash
mind insight BRAF V600E -t Melanoma           # Default
mind insight BRAF V600E -t Melanoma --lite    # Lite
mind insight BRAF V600E -t Melanoma --full    # Full
```

### 4. Async-First Design

All I/O operations are async for efficient parallel fetching:

```python
# Evidence builder fetches from 8+ sources in parallel
results = await asyncio.gather(
    self.myvariant_client.fetch_evidence(...),
    self.fda_client.fetch_drug_approvals(...),
    self.vicc_client.fetch_associations(...),
    self.civic_client.fetch_assertions(...),
    # ...
    return_exceptions=True,
)
```

### 5. Graceful Degradation

Individual API failures don't break the pipeline:

```python
if isinstance(vicc_result, Exception):
    print(f"Warning: VICC API failed: {vicc_result}")
    sources_failed.append("VICC")
else:
    # Process successful result
    sources_with_data.append("VICC")
```

## Module Details

### Public API (`api_public/insight.py`)

The single entry point for users:

```python
async def get_insight(
    variant_str: str,              # "BRAF V600E" or "EGFR L858R in NSCLC"
    tumor_type: str | None = None, # Optional tumor context
    config: InsightConfig | None = None,
) -> Insight:
```

**Responsibilities:**
- Parse variant input
- Validate variant type (SNP/indel only)
- Orchestrate evidence building
- Apply optional LLM enhancement
- Return strongly-typed Insight

### Normalization (`normalization/`)

Converts various input formats to canonical representation:

```python
# Supported input formats
parse_variant_input("BRAF V600E")           # Gene + variant
parse_variant_input("EGFR L858R in NSCLC")  # With tumor type
parse_variant_input("BRAF:V600E")           # Colon separator
parse_variant_input("TP53 p.R248W")         # HGVS notation

# Output: ParsedVariant
@dataclass
class ParsedVariant:
    gene: str                    # "BRAF"
    variant: str                 # "V600E"
    variant_normalized: str      # "V600E"
    variant_type: str            # "missense"
    tumor_type: str | None       # "NSCLC"
    parse_confidence: float      # 1.0
    parse_warnings: list[str]    # []
```

**Key Features:**
- Regex-based pattern matching
- Gene alias resolution (HER2 → ERBB2)
- Tumor type extraction from free text
- Variant type classification

### Evidence Builder (`evidence/builder.py`)

Orchestrates parallel API calls and assembles results:

```python
class EvidenceBuilder:
    """Async context manager for evidence aggregation."""

    async def __aenter__(self):
        # Initialize all API client sessions
        await self.myvariant_client.__aenter__()
        await self.vicc_client.__aenter__()
        # ...
        return self

    async def build_insight(
        self,
        variant: ParsedVariant,
        tumor_type: str | None,
    ) -> Insight:
        # Parallel fetch from all sources
        results = await asyncio.gather(
            self.myvariant_client.fetch_evidence(...),
            self.fda_client.fetch_drug_approvals(...),
            asyncio.to_thread(self.cgi_client.fetch_biomarkers, ...),
            fetch_vicc(),            # local async function
            fetch_civic_assertions(), # local async function
            fetch_clinical_trials(),  # local async function
            fetch_literature(),       # local async function
            return_exceptions=True,
        )

        # Assemble into Insight
        return Insight(
            identifiers=...,
            kb=...,
            functional=...,
            clinical=...,
            literature=...,
            meta=...,
            llm=None,  # Populated later if LLM enabled
        )
```

### API Clients (`api/`)

Each client follows a consistent pattern:

```python
class VICCClient:
    """Client for VICC MetaKB API."""

    def __init__(self):
        self.base_url = "https://search.cancervariants.org/api/v1"
        self.session: httpx.AsyncClient | None = None

    async def __aenter__(self):
        self.session = httpx.AsyncClient(timeout=30.0)
        return self

    async def __aexit__(self, *args):
        if self.session:
            await self.session.aclose()

    async def fetch_associations(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None = None,
    ) -> list[VICCAssociation]:
        # API call and response parsing
        ...
```

**Available Clients:**

| Client | API | Data |
|--------|-----|------|
| `MyVariantClient` | myvariant.info | ClinVar, COSMIC, gnomAD, CADD |
| `CIViCClient` | CIViC GraphQL | Curated variant-drug evidence |
| `VICCClient` | VICC MetaKB | Aggregated knowledgebases |
| `FDAClient` | OpenFDA | Drug approvals |
| `CGIClient` | Local TSV | Cancer biomarkers |
| `PubMedClient` | NCBI E-utils | Literature search |
| `SemanticScholarClient` | S2 API | Literature with citations |
| `ClinicalTrialsClient` | ClinicalTrials.gov | Active trials |
| `OncoTreeClient` | OncoTree API | Tumor type resolution |

### Models (`models/`)

Pydantic models for type safety and validation:

```python
class Insight(BaseModel):
    """Top-level evidence aggregation with optional LLM narrative."""

    identifiers: VariantIdentifiers    # Gene, variant, IDs
    kb: KnowledgebaseEvidence          # CIViC, ClinVar, VICC, etc.
    functional: FunctionalScores       # AlphaMissense, CADD, etc.
    clinical: ClinicalContext          # FDA, trials, gene role
    literature: LiteratureEvidence     # PubMed, extracted knowledge
    llm: LLMInsight | None             # Optional LLM narrative

    def get_summary(self) -> str:
        """Generate one-line summary of evidence."""
        ...

    def get_evidence_summary_for_llm(self) -> str:
        """Generate compact evidence summary for LLM prompt."""
        ...
```

**Model Hierarchy:**

```
Insight
├── VariantIdentifiers
│   ├── gene, variant, variant_type
│   ├── cosmic_id, clinvar_id, dbsnp_id
│   └── hgvs_protein, hgvs_genomic
├── KnowledgebaseEvidence
│   ├── civic: list[CIViCEvidence]
│   ├── civic_assertions: list[CIViCAssertionEvidence]
│   ├── clinvar: list[ClinVarEvidence]
│   ├── vicc: list[VICCEvidence]
│   └── cgi_biomarkers: list[CGIBiomarkerEvidence]
├── FunctionalScores
│   ├── alphamissense_score, alphamissense_prediction
│   ├── cadd_score, polyphen2_prediction, sift_prediction
│   └── gnomad_exome_af, gnomad_genome_af
├── ClinicalContext
│   ├── tumor_type, tumor_type_resolved
│   ├── fda_approvals: list[FDAApproval]
│   ├── clinical_trials: list[ClinicalTrialEvidence]
│   └── gene_role, gene_class, pathway
├── LiteratureEvidence
│   ├── pubmed_articles: list[PubMedEvidence]
│   └── literature_knowledge: LiteratureKnowledge | None
└── LLMInsight | None  ← Optional, embedded when LLM mode enabled
    ├── llm_summary: str
    ├── rationale: str
    ├── recommended_therapies: list[RecommendedTherapy]
    ├── clinical_trials_available: bool
    └── references: list[str]
```

### LLM Service (`llm/service.py`)

Optional LLM-powered analysis:

```python
class LLMService:
    """LLM service for generating variant narratives."""

    async def get_llm_insight(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
        evidence_summary: str,           # From Insight.get_evidence_summary_for_llm()
        has_clinical_trials: bool,
    ) -> LLMInsight:
        """Generate LLM narrative from evidence summary."""

    async def score_paper_relevance(
        self,
        title: str,
        abstract: str | None,
        tldr: str | None,
        gene: str,
        variant: str,
        tumor_type: str | None,
    ) -> dict:
        """Score paper relevance and extract signals."""

    async def extract_variant_knowledge(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
        paper_contents: list[dict],
    ) -> dict:
        """Extract structured knowledge from papers."""
```

**LLM Integration Flow:**

```
                EvidenceBuilder
                      │
                      ▼
                   Insight (with structured evidence)
                      │
                      ▼
        insight.get_evidence_summary_for_llm()
                      │
                      ▼
        ┌─────────────────────────────────┐
        │ LLMService.get_llm_insight()    │
        │   - evidence_summary: str       │
        │   - has_clinical_trials: bool   │
        └─────────────────────────────────┘
                      │
                      ▼
                  LLMInsight
                      │
                      ▼
        insight.llm = LLMInsight  ← Embedded in parent
```

## Data Flow

### Single Variant Annotation

```
User Input: "BRAF V600E in Melanoma"
              │
              ▼
┌─────────────────────────────────┐
│ 1. Parse Input                  │
│    parse_variant_input()        │
│    → ParsedVariant              │
│      gene="BRAF"                │
│      variant="V600E"            │
│      tumor_type="Melanoma"      │
└───────────────┬─────────────────┘
                │
                ▼
┌─────────────────────────────────┐
│ 2. Validate Variant Type        │
│    - Check: missense ✓          │
│    - Reject: fusion, amp ✗      │
└───────────────┬─────────────────┘
                │
                ▼
┌─────────────────────────────────┐
│ 3. Build Insight                │
│    EvidenceBuilder              │
│    ┌────────────────────────┐   │
│    │ Parallel API Calls:    │   │
│    │ • MyVariant.info       │   │
│    │ • FDA OpenAPI          │   │
│    │ • CGI Biomarkers       │   │
│    │ • VICC MetaKB          │   │
│    │ • CIViC Assertions     │   │
│    │ • ClinicalTrials.gov   │   │
│    │ • Semantic Scholar     │   │
│    └────────────────────────┘   │
└───────────────┬─────────────────┘
                │
                ▼
┌─────────────────────────────────┐
│ 4. Assemble Insight             │
│    - Merge results              │
│    - Track failures             │
│    - llm = None                 │
└───────────────┬─────────────────┘
                │
                ▼ (if enable_llm=True)
┌─────────────────────────────────┐
│ 5. LLM Enhancement              │
│    - Get evidence summary       │
│    - Call LLMService            │
│    - Embed LLMInsight           │
│    → insight.llm = LLMInsight   │
└───────────────┬─────────────────┘
                │
                ▼
           Insight
```

### Batch Processing

```
User Input: ["BRAF V600E", "EGFR L858R", "KRAS G12C"]
              │
              ▼
┌─────────────────────────────────┐
│ For each variant:               │
│   1. Parse                      │
│   2. Validate                   │
│   3. Queue for processing       │
└───────────────┬─────────────────┘
                │
                ▼
┌─────────────────────────────────┐
│ Sequential Processing           │
│ (to respect API rate limits)    │
│                                 │
│ Variant 1 ──► Insight 1         │
│ Variant 2 ──► Insight 2         │
│ Variant 3 ──► Insight 3         │
│                                 │
│ Progress callback:              │
│   progress_callback(i, total)   │
└───────────────┬─────────────────┘
                │
                ▼
         list[Insight]
```

## Configuration

### InsightConfig

```python
@dataclass
class InsightConfig:
    # Evidence sources
    enable_vicc: bool = True
    enable_civic_assertions: bool = True
    enable_clinical_trials: bool = True
    enable_literature: bool = True

    # LLM
    enable_llm: bool = False
    llm_model: str = "gpt-4o-mini"
    llm_temperature: float = 0.1

    # Limits
    max_vicc_results: int = 15
    max_clinical_trials: int = 10
    max_literature_results: int = 6

    # Validation
    validate_variant_type: bool = True
```

### EvidenceBuilderConfig

```python
@dataclass
class EvidenceBuilderConfig:
    # Source toggles
    enable_vicc: bool = True
    enable_civic_assertions: bool = True
    enable_clinical_trials: bool = True
    enable_literature: bool = True

    # Result limits
    max_vicc_results: int = 15
    max_civic_assertions: int = 20
    max_clinical_trials: int = 10
    max_literature_results: int = 6

    # API keys
    semantic_scholar_api_key: str | None = None
```

## Error Handling

### API Failures

Each API call is wrapped with exception handling:

```python
results = await asyncio.gather(
    self._fetch_vicc(),
    self._fetch_civic(),
    return_exceptions=True,  # Don't fail on individual errors
)

for i, result in enumerate(results):
    if isinstance(result, Exception):
        sources_failed.append(source_names[i])
        print(f"Warning: {source_names[i]} failed: {result}")
    else:
        sources_with_data.append(source_names[i])
```

### Rate Limiting & Retry Strategy

All API clients use **tenacity** for standardized retry with exponential backoff and jitter:

```python
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential_jitter,
)

class ClinicalTrialsClient:
    @retry(
        retry=retry_if_exception_type(ClinicalTrialsRateLimitError),
        stop=stop_after_attempt(3),
        wait=wait_exponential_jitter(initial=2, max=15, jitter=1),
        reraise=True,
    )
    async def _make_request(self, params):
        # Request logic...
        if response.status_code == 429:
            raise ClinicalTrialsRateLimitError("Rate limited")
```

**Retry Configuration by Client:**

| Client | Trigger | Attempts | Backoff (initial → max) | Jitter |
|--------|---------|----------|-------------------------|--------|
| MyVariant | `HTTPError`, `TimeoutException` | 3 | 2s → 10s | None |
| FDA | `HTTPError`, `TimeoutException` | 3 | 2s → 10s | None |
| OncoTree | `HTTPError`, `TimeoutException` | 3 | 2s → 10s | None |
| Semantic Scholar | `SemanticScholarRateLimitError` | 3 | 1s → 10s | ±1s |
| PubMed | `PubMedRateLimitError` | 3 | 0.5s → 5s | ±0.5s |
| ClinicalTrials | `ClinicalTrialsRateLimitError` | 3 | 2s → 15s | ±1s |

### Variant Validation

Invalid variant types are rejected early:

```python
if config.validate_variant_type:
    if parsed.variant_type not in ALLOWED_TYPES:
        raise ValueError(
            f"Variant type '{parsed.variant_type}' not supported. "
            f"Only SNPs and small indels are allowed."
        )
```

## Experimental Features

### Tiering (`experimental/tiering.py`)

Non-authoritative AMP/ASCO/CAP tier computation:

```python
from oncomind.experimental import compute_experimental_tier

result = compute_experimental_tier(insight)
# TierResult(
#     tier="Tier I",
#     level="A",
#     confidence=0.8,
#     rationale="FDA-approved targeted therapy exists...",
#     caveats=["This tier is NOT authoritative..."]
# )
```

### Embeddings (`embeddings/features.py`)

Feature extraction for ML applications:

```python
from oncomind.embeddings import extract_features

features = extract_features(insight)
# {
#     "alphamissense_score": 0.98,
#     "cadd_score": 32.0,
#     "civic_evidence_count": 12,
#     "has_fda_approval": True,
#     ...
# }
```

## Performance Considerations

### Parallel Fetching

Evidence sources are queried in parallel:
- **Lite mode** (`--lite`): ~7 seconds — structured evidence only, no LLM
- **Default mode**: ~12 seconds — structured evidence + LLM clinical narrative
- **Full mode** (`--full`): ~25 seconds — + literature search + enhanced narrative

### Caching (Future)

Planned caching layers:
- HTTP response caching (httpx)
- Literature pre-fetching for common genes
- Redis/SQLite for persistent cache

### Rate Limits

| API | Rate Limit | Strategy |
|-----|------------|----------|
| MyVariant.info | 1000/day free | Tenacity retry (3 attempts, 2-10s backoff) |
| VICC MetaKB | Unlimited | None needed |
| CIViC | Unlimited | None needed |
| Semantic Scholar | 1 RPS (free), 10 RPS (key) | Tenacity retry + PubMed fallback |
| PubMed | 3/sec, 10/sec with key | Tenacity retry (3 attempts, 0.5-5s backoff) |
| ClinicalTrials.gov | ~50/min | Tenacity retry (3 attempts, 2-15s backoff) |
| FDA OpenFDA | Reasonable use | Tenacity retry (3 attempts, 2-10s backoff) |
| OncoTree | Unlimited | Tenacity retry for transient failures |

## Testing Strategy

### Unit Tests (`tests/unit/`)

- Model validation
- Parser edge cases
- API client mocking
- Normalization accuracy

```bash
pytest tests/unit/ -v
```

### Integration Tests (`tests/integration/`)

- Real API calls (marked with `@pytest.mark.integration`)
- End-to-end annotation
- Rate limit handling

```bash
pytest tests/integration/ -v -m integration
```

## Future Architecture

### Planned Enhancements

1. **VCF Pipeline**: Batch VCF annotation with streaming
2. **Structural Variants**: Fusion/amplification support
3. **Agent Workflow**: Multi-agent analysis with LangGraph
4. **Caching Layer**: Redis/SQLite for performance
5. **Webhooks**: Async result delivery for long-running jobs

### Extension Points

- New API clients: Implement base client pattern
- New evidence types: Add to models/evidence/
- New normalizers: Add patterns to input_parser.py
- New LLM tasks: Add methods to llm/service.py

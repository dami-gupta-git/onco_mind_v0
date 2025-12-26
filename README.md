# OncoMind

**Research intelligence for cancer variants. Find the gaps, not just the facts.**

For BRAF V600E, databases have the answers. For the next 10,000 variants, the interesting question isn't "what do we know?" — it's "what don't we know yet?"

OncoMind aggregates evidence across databases, identifies knowledge gaps, and generates research hypotheses — turning variant annotation into research intelligence.

> **Note**: This tool is for research purposes only. It is not intended for clinical decision-making.

## The Problem

Existing variant annotation tools excel at telling you what's known. But researchers need more:

- **What's well-characterized vs. under-studied?**
- **Where do databases disagree?**
- **What co-mutation patterns suggest testable hypotheses?**
- **Which functional studies are missing?**

OncoMind synthesizes evidence from 10+ sources and explicitly surfaces what we *don't* know — the gaps that represent research opportunities.

## What Makes OncoMind Different

| Feature | Annotation Tools | OncoMind |
|---------|------------------|----------|
| **Focus** | What the variant *is* | What we *don't yet know* |
| **Knowledge gaps** | Not addressed | Explicitly identified |
| **Co-mutation patterns** | Data only | Hypothesis generation |
| **Source conflicts** | Pick one, hide the rest | Surface and explain |
| **Output** | Static annotations | Research-ready insights |
| **Evidence quality** | Implied | Explicitly assessed |

## Features

### Evidence Aggregation
- **Clinical databases**: CIViC, ClinVar, COSMIC, VICC MetaKB, CGI, FDA labels
- **Functional predictions**: AlphaMissense, CADD, PolyPhen2, gnomAD frequencies
- **Co-mutation patterns**: cBioPortal prevalence and co-occurrence data
- **Preclinical research**: DepMap gene essentiality (CRISPR screens), drug sensitivity (PRISM), cell line models
- **Literature**: Semantic Scholar and PubMed with LLM relevance scoring
- **Clinical trials**: ClinicalTrials.gov integration

### Research Intelligence
- **Evidence gap detection**: Explicitly identifies what's missing (no functional data, limited clinical evidence, etc.)
- **Conflict surfacing**: Flags when databases disagree and explains why
- **Hypothesis generation**: Suggests testable research directions based on evidence patterns
- **Evidence quality assessment**: Rates overall evidence strength (comprehensive/moderate/limited/minimal)
- **Source attribution**: Every claim links to a PMID, database entry, or cBioPortal study

## Installation

```bash
git clone https://github.com/yourusername/oncomind.git
cd oncomind/onco_mind_v0

python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

pip install -e ".[dev]"
```

> **Note:** All dependencies are managed in `pyproject.toml`. If you encounter `ModuleNotFoundError` after pulling updates, reinstall with:
> ```bash
> pip install -e . --force-reinstall
> ```

## Quick Start

### Python API

```python
import asyncio
from oncomind import get_insight, InsightConfig, Result

async def main():
    # Default: Fast annotation mode (no LLM, ~7s)
    result = await get_insight("BRAF V600E", tumor_type="Melanoma")

    print(f"Gene: {result.identifiers.gene}")
    print(f"Therapies: {result.evidence.get_recommended_therapies()}")
    print(f"Clinical Trials: {len(result.clinical.clinical_trials)}")
    print(f"Summary: {result.get_summary()}")

    # LLM mode: literature search + AI synthesis (~25s)
    config = InsightConfig(
        enable_llm=True,           # Enables LLM synthesis
        enable_literature=True,    # Enables literature search (recommended with LLM)
        llm_model="gpt-4o-mini"
    )
    result = await get_insight("EGFR S768I", tumor_type="NSCLC", config=config)

    # Research-focused insights (only available with enable_llm=True)
    if result.llm:
        print(result.llm.llm_summary)
        print(f"Evidence quality: {result.llm.evidence_quality}")
        print(f"Knowledge gaps: {result.llm.knowledge_gaps}")
        print(f"Research implications: {result.llm.research_implications}")

asyncio.run(main())
```

### CLI

```bash
# Default: fast annotation mode (~7s)
mind insight BRAF V600E --tumor Melanoma
mind insight PIK3CA H1047R -t "Breast Cancer"

# LLM mode: literature search + AI synthesis (~25s)
mind insight KRAS G12C -t NSCLC --llm

# Save to JSON
mind insight BRAF V600E -t Melanoma --output result.json
mind insight EGFR L858R -t NSCLC --llm --output result.json
```

#### CLI Modes

| Mode | Flag | Speed | Output |
|------|------|-------|--------|
| **Annotation** | (none) | ~7s | Fast structured evidence from all databases |
| **LLM** | `--llm` | ~25s | + Literature search + AI synthesis + hypothesis generation |

**Output by mode:**

| Feature | Annotation | LLM |
|---------|------------|-----|
| Evidence Summary | ✓ | ✓ |
| Database Annotations | ✓ | ✓ |
| DepMap Research Context | ✓ | ✓ |
| Literature Search | - | ✓ |
| LLM Research Insight | - | ✓ |
| Knowledge Gaps | - | ✓ |
| Hypothesis Generation | - | ✓ |

#### Example Output

```
$ mind insight BRAF V600E -t Melanoma

╔══════════════════════════════════════════════════════════════════════════════╗
║                            BRAF V600E in Melanoma                            ║
╚══════════════════════════════════════════════════════════════════════════════╝
╭────────────────────────────── Evidence Summary ──────────────────────────────╮
│  COSMIC:COSM476                                                              │
│  ClinVar:            Pathogenic                                              │
│  Pathogenicity:      AlphaMissense: Pathogenic (0.99) | PolyPhen2: D         │
│  Gene Role:          oncogene                                                │
│  Evidence Sources:   CIViC, CIViC Assertions, ClinVar, COSMIC, CGI, VICC     │
╰──────────────────────────────────────────────────────────────────────────────╯
╭───────────────────────────── FDA Approved Drugs ─────────────────────────────╮
│  Mekinist, KEYTRUDA QLEX, BRAFTOVI, ZELBORAF, MEKTOVI                        │
╰──────────────────────────────────────────────────────────────────────────────╯
╭───────────────────────────── Clinical Evidence ──────────────────────────────╮
│  CIViC Assertions:                                                           │
│    • Trametinib, Dabrafenib → SENSITIVITYRESPONSE                            │
│    • Vemurafenib, Cobimetinib → SENSITIVITYRESPONSE                          │
│                                                                              │
│  CGI Biomarkers:                                                             │
│    • PLX4720: Responsive                                                     │
│    • Vemurafenib: Responsive                                                 │
╰──────────────────────────────────────────────────────────────────────────────╯
```

### LLM-Ready Context

OncoMind outputs context blocks designed for injection into downstream LLM systems:

```python
# Get a dense, grounded context block
context = panel.to_knowledge_header()

# Returns something like:
# "BRAF V600E in melanoma. Oncogenic driver via constitutive MAPK activation.
#  FDA-approved: dabrafenib + trametinib, vemurafenib + cobimetinib.
#  Resistance typical at 6-12 months via NRAS (20%), MEK1/2 (5-10%), or BRAF amp.
#  Sources: CIViC:assertion:12, FDA label, PMID:22735384"
```

## Configuration

Create a `.env` file:

```bash
# Required for LLM features
OPENAI_API_KEY=your-openai-key

# Optional: enhanced literature search (recommended)
SEMANTIC_SCHOLAR_API_KEY=your-s2-key

# Optional: other LLM providers
ANTHROPIC_API_KEY=your-anthropic-key
```

## Result Model

The `Result` is the core output model containing:

| Field | Contents |
|-------|----------|
| `evidence` | `Evidence` model with all structured data |
| `llm` | `LLMInsight` with LLM narrative (optional, only when LLM enabled) |

The `Evidence` model contains:

| Section | Contents |
|---------|----------|
| `identifiers` | Gene, variant, HGVS notation, COSMIC/dbSNP IDs |
| `kb` | CIViC assertions, VICC MetaKB, CGI biomarkers |
| `functional` | AlphaMissense, CADD, PolyPhen2, gnomAD frequencies |
| `clinical` | FDA approvals, clinical trials, gene role, ClinVar significance |
| `literature` | PubMed articles, literature knowledge |

```python
result = await get_insight("BRAF V600E", tumor_type="Melanoma")

# Access via shortcut properties on Result
result.identifiers.gene                # "BRAF"
result.kb.civic_assertions             # Curated drug-variant associations
result.functional.alphamissense_score  # 0.98 (pathogenicity)
result.clinical.fda_approvals          # FDA-approved therapies

# Helper methods
result.get_summary()                   # One-line summary
result.has_evidence()                  # True if any evidence found
result.evidence.get_recommended_therapies()  # FDA + CGI therapies
```

See [docs/API_REFERENCE.md](docs/API_REFERENCE.md) for complete field documentation.

## Supported Variant Types

**Input formats:**
Accepts protein-level variants (V600E, Val600Glu, p.V600E). For cDNA or genomic coordinates, use HGVS notation or let MyVariant.info handle the lookup.

**Currently supported:**
- Missense mutations (V600E, L858R)
- Nonsense mutations (R248*)
- Small insertions/deletions (E746_A750del)
- Frameshift mutations (K132fs)

**Coming soon:** Fusions, amplifications, copy number variants (see [ROADMAP.md](docs/ROADMAP.md))

## Data Sources

### Clinical & Therapeutic
| Source | Data Type | Access |
|--------|-----------|--------|
| [CIViC](https://civicdb.org/) | Curated variant-drug associations | Free API |
| [VICC MetaKB](https://search.cancervariants.org/) | Aggregated knowledgebases | Free API |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinical significance | Via MyVariant.info |
| [COSMIC](https://cancer.sanger.ac.uk/cosmic) | Somatic mutation catalog | Via MyVariant.info |
| [CGI](https://www.cancergenomeinterpreter.org/) | Biomarker annotations | Local database |
| [FDA](https://www.fda.gov/) | Drug approvals | OpenFDA API |
| [ClinicalTrials.gov](https://clinicaltrials.gov/) | Active trials | Free API |

### Functional & Biological
| Source | Data Type | Access |
|--------|-----------|--------|
| [cBioPortal](https://www.cbioportal.org/) | Co-mutation patterns, prevalence | Free API |
| [AlphaMissense](https://alphamissense.hegelab.org/) | Pathogenicity predictions | Via MyVariant.info |
| [gnomAD](https://gnomad.broadinstitute.org/) | Population frequencies | Via MyVariant.info |

### Literature
| Source | Data Type | Access |
|--------|-----------|--------|
| [Semantic Scholar](https://www.semanticscholar.org/) | AI-powered literature search | Free API |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov/) | Biomedical literature | Free E-utilities |

### Preclinical Research
| Source | Data Type | Access |
|--------|-----------|--------|
| [DepMap](https://depmap.org/) | Gene essentiality (CRISPR), drug sensitivity (PRISM), cell line models | Free API |

### Coming Soon
| Source | Data Type | Status |
|--------|-----------|--------|
| [Reactome](https://reactome.org/) | Pathway context | Planned |

See [ROADMAP.md](docs/ROADMAP.md) for the full development roadmap.

## Development

```bash
# Run tests
pytest tests/unit/ -v

# Run with coverage
pytest tests/unit/ --cov=src/oncomind --cov-report=html

# Type checking
mypy src/oncomind

# Linting
ruff check src/oncomind
ruff format src/oncomind
```

## Streamlit App

```bash
cd streamlit
streamlit run app.py
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [CIViC](https://civicdb.org/) - Clinical Interpretation of Variants in Cancer
- [VICC](https://cancervariants.org/) - Variant Interpretation for Cancer Consortium
- [MyVariant.info](https://myvariant.info/) - Variant annotation aggregator
- [Semantic Scholar](https://www.semanticscholar.org/) - AI-powered research tool
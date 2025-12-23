# OncoMind

**AI-powered cancer variant annotation and evidence synthesis**

OncoMind aggregates evidence from multiple cancer databases and uses LLMs to synthesize actionable insights for somatic variants. It provides a unified interface for variant annotation, combining structured database queries with intelligent literature analysis.

> **Disclaimer**: This tool is for research purposes only. Clinical decisions should always be made by qualified healthcare professionals.

## Features

- **Multi-source evidence aggregation**: CIViC, ClinVar, COSMIC, VICC MetaKB, CGI, FDA labels, ClinicalTrials.gov
- **Functional predictions**: AlphaMissense, CADD, PolyPhen2, SIFT, gnomAD frequencies
- **Literature search**: Semantic Scholar and PubMed with relevance scoring
- **LLM-powered synthesis**: Extract resistance/sensitivity signals from literature
- **Flexible input**: Free-text variants, CSV files, VCF files, or programmatic API
- **Strongly-typed output**: `EvidencePanel` model with organized evidence sections

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/oncomind.git
cd oncomind/onco_mind_v0

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"
```

## Quick Start

### Python API

```python
import asyncio
from oncomind import process_variant, AnnotationConfig

async def main():
    # Simple annotation (no LLM, fast)
    panel = await process_variant("BRAF V600E", tumor_type="Melanoma")

    # Access structured evidence
    print(f"Gene: {panel.identifiers.gene}")
    print(f"FDA Approved Drugs: {panel.clinical.get_approved_drugs()}")
    print(f"Clinical Trials: {len(panel.clinical.clinical_trials)}")
    print(f"Evidence Sources: {panel.kb.get_evidence_sources()}")

    # With LLM enhancement for literature analysis
    config = AnnotationConfig(enable_llm=True, llm_model="gpt-4o-mini")
    panel = await process_variant("EGFR L858R", tumor_type="NSCLC", config=config)

asyncio.run(main())
```

### Synchronous API

```python
from oncomind import process_variant_sync

# For scripts without async support
panel = process_variant_sync("BRAF V600E", tumor_type="Melanoma")
print(panel.clinical.get_approved_drugs())
```

### Batch Processing

```python
from oncomind import process_variants

# From list of strings
panels = await process_variants(["BRAF V600E", "EGFR L858R", "KRAS G12C"])

# From CSV file
panels = await process_variants("variants.csv", tumor_type="NSCLC")

# From pandas DataFrame
import pandas as pd
df = pd.DataFrame({"gene": ["BRAF", "EGFR"], "variant": ["V600E", "L858R"]})
panels = await process_variants(df)
```

### CLI

```bash
# Annotate a single variant
mind annotate "BRAF V600E" --tumor Melanoma

# Save to JSON
mind annotate "EGFR L858R in NSCLC" --output result.json

# With LLM enhancement
mind annotate "KRAS G12C" --tumor NSCLC --llm

# Legacy command with full LLM narrative
mind process BRAF V600E --tumor Melanoma --model gpt-4o-mini
```

## Configuration

Create a `.env` file with your API keys:

```bash
# Required for LLM features
OPENAI_API_KEY=your-openai-key

# Optional: for enhanced literature search
SEMANTIC_SCHOLAR_API_KEY=your-s2-key

# Optional: other LLM providers
ANTHROPIC_API_KEY=your-anthropic-key
GROQ_API_KEY=your-groq-key
```

## EvidencePanel Structure

The `EvidencePanel` is the primary output, organizing evidence into logical sections:

```python
panel = await process_variant("BRAF V600E")

# Identifiers - variant IDs and HGVS notation
panel.identifiers.gene          # "BRAF"
panel.identifiers.variant       # "V600E"
panel.identifiers.cosmic_id     # "COSM476"
panel.identifiers.hgvs_protein  # "p.V600E"

# Knowledgebase evidence
panel.kb.civic              # CIViC evidence entries
panel.kb.civic_assertions   # CIViC curated assertions
panel.kb.clinvar            # ClinVar entries
panel.kb.vicc               # VICC MetaKB associations
panel.kb.cgi_biomarkers     # CGI biomarker evidence

# Functional predictions
panel.functional.alphamissense_score       # 0.98
panel.functional.alphamissense_prediction  # "P" (Pathogenic)
panel.functional.cadd_score                # 32.0
panel.functional.gnomad_exome_af           # 0.00001

# Clinical context
panel.clinical.tumor_type           # "Melanoma"
panel.clinical.fda_approvals        # FDA-approved drugs
panel.clinical.clinical_trials      # Matching trials
panel.clinical.gene_role            # "oncogene"
panel.clinical.get_approved_drugs() # ["Dabrafenib", "Vemurafenib"]

# Literature
panel.literature.pubmed_articles      # PubMed articles
panel.literature.literature_knowledge # LLM-extracted insights

# Metadata
panel.meta.sources_queried    # ["MyVariant", "FDA", "CIViC", ...]
panel.meta.sources_with_data  # Sources that returned evidence
panel.meta.evidence_strength  # "Strong" / "Moderate" / "Weak"
```

## Supported Variant Types

OncoMind currently supports **SNPs and small indels**:

- Missense mutations (e.g., V600E, L858R)
- Nonsense mutations (e.g., R248*)
- Small insertions/deletions (e.g., E746_A750del)
- Frameshift mutations (e.g., K132fs)

**Not yet supported**: Fusions, amplifications, copy number variants, large structural variants.

## Data Sources

| Source | Data Type | Access |
|--------|-----------|--------|
| [CIViC](https://civicdb.org/) | Curated variant-drug associations | Free API |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinical significance | Via MyVariant.info |
| [COSMIC](https://cancer.sanger.ac.uk/cosmic) | Somatic mutation catalog | Via MyVariant.info |
| [VICC MetaKB](https://search.cancervariants.org/) | Aggregated knowledgebases | Free API |
| [CGI](https://www.cancergenomeinterpreter.org/) | Biomarker annotations | Local database |
| [FDA](https://www.fda.gov/) | Drug approvals | OpenFDA API |
| [ClinicalTrials.gov](https://clinicaltrials.gov/) | Active trials | Free API |
| [Semantic Scholar](https://www.semanticscholar.org/) | Literature | Free API (key recommended) |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov/) | Literature | Free E-utilities |

## Project Structure

```
src/oncomind/
├── __init__.py           # Public API exports
├── cli.py                # Command-line interface
├── engine.py             # Legacy InsightEngine
├── api/                  # API clients
│   ├── civic.py          # CIViC GraphQL client
│   ├── myvariant.py      # MyVariant.info client
│   ├── vicc.py           # VICC MetaKB client
│   ├── fda.py            # OpenFDA client
│   ├── pubmed.py         # PubMed E-utilities
│   └── ...
├── api_public/           # Public API
│   └── annotate.py       # process_variant, process_variants
├── models/               # Pydantic models
│   └── evidence/
│       ├── evidence_panel.py  # EvidencePanel model
│       ├── civic.py           # CIViC evidence models
│       └── ...
├── normalization/        # Input parsing
│   ├── input_parser.py   # parse_variant_input
│   └── hgvs_utils.py     # Variant normalization
├── evidence/             # Evidence aggregation
│   └── builder.py        # EvidenceBuilder
├── embeddings/           # Feature extraction (experimental)
│   └── features.py       # extract_features
├── experimental/         # Experimental features
│   └── tiering.py        # Non-authoritative tier computation
└── llm/                  # LLM services
    └── service.py        # LLM-based analysis
```

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

A web interface is available in the `streamlit/` directory:

```bash
cd streamlit
streamlit run app.py
```

## Roadmap

- [ ] Fusion/amplification support
- [ ] VCF annotation pipeline
- [ ] ESMFold protein structure visualization
- [ ] SpliceAI integration
- [ ] Multi-agent workflow (LangGraph)
- [ ] Pre-fetched literature cache

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [CIViC](https://civicdb.org/) - Clinical Interpretation of Variants in Cancer
- [VICC](https://cancervariants.org/) - Variant Interpretation for Cancer Consortium
- [MyVariant.info](https://myvariant.info/) - Variant annotation aggregator
- [Semantic Scholar](https://www.semanticscholar.org/) - AI-powered research tool

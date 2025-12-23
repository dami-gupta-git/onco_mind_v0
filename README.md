# OncoMind

**Grounded context for cancer variant reasoning. What happens next, with receipts.**

For common variants like BRAF V600E, databases have the answers. For everything else, there's a 2-hour PubMed rabbit hole. OncoMind does that search for you — and shows its work.

> **Disclaimer**: This tool is for research purposes only. Clinical decisions should always be made by qualified healthcare professionals.

## The Problem

Variant interpretation tools like CIViC, OncoKB, and CancerVar are excellent for well-characterized mutations. But for less common variants, researchers still spend hours manually searching PubMed, ClinicalTrials.gov, and Google Scholar for case reports, functional studies, and resistance mechanisms.

OncoMind automates that workflow: aggregate databases, search literature, synthesize findings, and cite sources — in seconds instead of hours.

## What Makes OncoMind Different

| Feature | Legacy Tools | OncoMind |
|---------|--------------|----------|
| **Output format** | 50-page PDFs, 10K-row CSVs | LLM-ready context blocks |
| **Source disagreement** | Pick one, hide the rest | Surface conflicts explicitly |
| **Uncertainty** | Implied | Quantified (ensemble agreement) |
| **Claims without sources** | Common | Forbidden by design |
| **Resistance mechanisms** | What the variant *is* | What happens *next* |

## Features

- **Multi-source aggregation**: CIViC, ClinVar, COSMIC, VICC MetaKB, CGI, FDA labels, ClinicalTrials.gov
- **Functional predictions**: AlphaMissense, CADD, PolyPhen2, SIFT, gnomAD frequencies
- **Literature search**: Semantic Scholar and PubMed with LLM relevance scoring
- **Evidence synthesis**: Extract resistance/sensitivity signals from abstracts
- **Conflict detection**: Flag when databases disagree
- **Source attribution**: Every claim links to a PMID, FDA label, or database entry
- **Strongly-typed output**: Pydantic `EvidencePanel` model for programmatic use

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
from oncomind import get_insight, AnnotationConfig

async def main():
    # Fast annotation (no LLM)
    panel = await get_insight("BRAF V600E", tumor_type="Melanoma")

    print(f"Gene: {panel.identifiers.gene}")
    print(f"FDA Approved: {panel.clinical.get_approved_drugs()}")
    print(f"Clinical Trials: {len(panel.clinical.clinical_trials)}")
    print(f"Sources: {panel.meta.sources_with_data}")

    # With LLM synthesis for literature analysis
    config = AnnotationConfig(enable_llm=True, llm_model="gpt-4o-mini")
    panel = await get_insight("EGFR S768I", tumor_type="NSCLC", config=config)

    # LLM-extracted insights from literature
    print(panel.literature.literature_knowledge)

asyncio.run(main())
```

### CLI

```bash
# Basic insight
mind insight PIK3CA H1047R --tumor "Breast Cancer"
mind insight PIK3CA H1047R -t "Breast Cancer"


# With LLM synthesis
mind insight PIK3CA H1047R --tumor "Breast Cancer" --llm

# Save to JSON
mind insight KRAS G12C --tumor NSCLC --output result.json

# Full LLM narrative (separate gene and variant args)
mind insight-llm BRAF V600E --tumor Melanoma
mind insight-llm IDH1 R132H -t "Acute Myeloid Leukemia"
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

## EvidencePanel Structure

```python
panel = await get_insight("BRAF V600E", tumor_type="Melanoma")

# Identifiers
panel.identifiers.gene              # "BRAF"
panel.identifiers.variant           # "V600E"
panel.identifiers.cosmic_id         # "COSM476"
panel.identifiers.hgvs_protein      # "p.V600E"

# Knowledgebase evidence
panel.kb.civic_assertions           # CIViC curated assertions
panel.kb.vicc                       # VICC MetaKB (OncoKB, CIViC, MOAlmanac)
panel.kb.cgi_biomarkers             # CGI biomarker evidence
panel.kb.clinvar                    # ClinVar entries

# Functional scores
panel.functional.alphamissense_score       # 0.98
panel.functional.alphamissense_prediction  # "P" (Pathogenic)
panel.functional.cadd_score                # 32.0
panel.functional.gnomad_exome_af           # 0.00001

# Clinical context
panel.clinical.fda_approvals        # FDA-approved therapies
panel.clinical.clinical_trials      # Matching active trials
panel.clinical.gene_role            # "oncogene"
panel.clinical.get_approved_drugs() # ["Dabrafenib", "Vemurafenib"]

# Literature
panel.literature.pubmed_articles         # Retrieved articles
panel.literature.literature_knowledge    # LLM-synthesized insights

# Metadata & trust
panel.meta.sources_queried          # All sources attempted
panel.meta.sources_with_data        # Sources that returned evidence
panel.meta.sources_failed           # Sources that errored
panel.meta.conflicts                # Cross-source disagreements
```

## Supported Variant Types

**Currently supported:**
- Missense mutations (V600E, L858R)
- Nonsense mutations (R248*)
- Small insertions/deletions (E746_A750del)
- Frameshift mutations (K132fs)

**Coming soon:** Fusions, amplifications, copy number variants (see [ROADMAP.md](docs/ROADMAP.md))

## Data Sources

| Source | Data Type | Access |
|--------|-----------|--------|
| [CIViC](https://civicdb.org/) | Curated variant-drug associations | Free API |
| [VICC MetaKB](https://search.cancervariants.org/) | Aggregated knowledgebases | Free API |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinical significance | Via MyVariant.info |
| [COSMIC](https://cancer.sanger.ac.uk/cosmic) | Somatic mutation catalog | Via MyVariant.info |
| [CGI](https://www.cancergenomeinterpreter.org/) | Biomarker annotations | Local database |
| [FDA](https://www.fda.gov/) | Drug approvals | OpenFDA API |
| [ClinicalTrials.gov](https://clinicaltrials.gov/) | Active trials | Free API |
| [Semantic Scholar](https://www.semanticscholar.org/) | Literature | Free API |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov/) | Literature | Free E-utilities |

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
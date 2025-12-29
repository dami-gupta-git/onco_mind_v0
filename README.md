# OncoMind

**Research intelligence for cancer variants. Find the gaps, not just the facts.**

For BRAF V600E, databases have the answers. For the next 10,000 variants, the real question is **"what don't we know yet?"**

OncoMind aggregates evidence across 10+ sources, identifies knowledge gaps, and generates research hypotheses — turning variant annotation into **research intelligence**, not clinical advice.

> **Research only.** Not for clinical decision-making or patient care.

---

## Why OncoMind?

Most tools tell you **what's known** about a variant. Researchers also need:

- What's **well-characterized vs under-studied**?
- Where do **databases disagree**?
- Which **co-mutation patterns** suggest testable hypotheses?
- What **functional / preclinical work** is missing?

OncoMind has two layers:

1. A deterministic **annotation backbone** (KBs, trials, cBioPortal, DepMap, literature).
2. An optional **LLM research layer** that surfaces gaps and proposes hypotheses, grounded in that evidence.

---

## What Makes It Different

| Feature                  | Typical Tools               | OncoMind                               |
|--------------------------|-----------------------------|----------------------------------------|
| Focus                    | "What is this variant?"     | "What don't we know (yet)?"           |
| Knowledge gaps           | Not explicit                | First-class citizens                   |
| Co-mutation patterns     | Raw data                    | Used to generate hypotheses            |
| Source conflicts         | Often hidden                | Surfaced and explained                 |
| Output                   | Static clinical annotations | Research-ready, gap-focused insights   |
| Evidence quality         | Implicit                    | Explicit rating per variant            |

---

## Quick Start

### Install

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

### Python API

```python
import asyncio
from oncomind import get_insight, InsightConfig

async def main():
    # Annotation mode: fast, deterministic evidence (~7s)
    result = await get_insight("BRAF V600E", tumor_type="Melanoma")
    print(result.get_summary())
    print(result.evidence.get_recommended_therapies())

    # LLM mode: annotation backbone + literature + research insight (~25s)
    config = InsightConfig(enable_llm=True, llm_model="claude-sonnet-4-20250514")
    result = await get_insight("MAP2K1 P124L", tumor_type="Melanoma", config=config)

    if result.llm:
        print(result.llm.llm_summary)
        print("Evidence quality:", result.llm.evidence_quality)
        print("Knowledge gaps:", result.llm.knowledge_gaps)
        print("Research implications:", result.llm.research_implications)

asyncio.run(main())
```

### CLI

```bash
# Annotation mode (~7s)
mind insight BRAF V600E --tumor Melanoma

# LLM research mode (~25s)
mind insight MAP2K1 P124L -t Melanoma --llm

# Save to JSON
mind insight EGFR L858R -t NSCLC --llm --output result.json
```

**Modes:**

| Mode           | Flag   | Output                                           |
|----------------|--------|--------------------------------------------------|
| Annotation     | (none) | Structured evidence from all databases          |
| LLM Research   | `--llm`| + Literature + research narrative + hypotheses  |

---

## What You Get

### 1) Evidence Backbone (Annotation)

OncoMind builds a rich, structured evidence model:

- **Clinical / KB**: CIViC, VICC MetaKB, ClinVar, COSMIC, CGI, FDA labels
- **Functional**: AlphaMissense, CADD, PolyPhen2, gnomAD
- **Biological**: cBioPortal prevalence and co-mutation patterns
- **Preclinical**: DepMap CRISPR essentiality, PRISM drug sensitivity
- **Trials**: ClinicalTrials.gov
- **Literature**: PubMed / Semantic Scholar entries in a structured `LiteratureEvidence` model

Example (BRAF V600E in melanoma):

```
Evidence Summary
  COSMIC: COSM476
  ClinVar: Pathogenic
  Pathogenicity: AlphaMissense 0.99 (Pathogenic) | PolyPhen2: D
  Gene Role: oncogene
FDA Approved Drugs
  Mekinist, BRAFTOVI, ZELBORAF, MEKTOVI, KEYTRUDA QLEX
Clinical Evidence
  CIViC: Trametinib, Dabrafenib → SENSITIVITY/RESPONSE
  CGI: PLX4720, Vemurafenib → Responsive
```

This backbone is available even when LLM is disabled and can be consumed directly via Python, CLI, or JSON.

#### Match Specificity

OncoMind tracks how precisely each piece of evidence matches your query:

| Match Level | Meaning | Example |
|-------------|---------|---------|
| `variant` | Exact variant match | BRAF V600E → evidence specifically for V600E |
| `codon` | Same position, different change | BRAF V600K → evidence for "V600 mutations" |
| `gene` | Gene-level only | BRAF V600E → evidence for "BRAF mutations" |

This helps distinguish targeted therapies (variant-specific) from broader biomarkers (gene-level).

#### Evidence Gap Analysis

OncoMind automatically identifies what's well-characterized vs under-studied:

```python
gaps = result.evidence.compute_evidence_gaps()

gaps.overall_evidence_quality  # "comprehensive" | "moderate" | "limited" | "minimal"
gaps.well_characterized        # ["FDA-approved therapies", "Functional impact"]
gaps.poorly_characterized      # ["Resistance mechanisms", "Preclinical models"]
gaps.gaps                      # Structured list with severity and suggested studies
```

---

### 2) Research Insight (LLM Mode)

On top of the evidence backbone, OncoMind can generate a **research-focused card** for each variant, constrained to the evidence it sees.

Example (MAP2K1 P124L in melanoma):

> **Insight Ready: MAP2K1 P124L in Melanoma**
> Functional impact: gain-of-function MAP2K1 variant activating MAPK signaling.
> Biological context: rare (~1%) but recurrent in melanoma, with frequent TP53 co-mutations in TCGA SKCM.
> Therapeutic landscape: limited clinical data; case reports describe trametinib response and resistance.
> Evidence quality: **moderate**.
> Knowledge gaps: resistance mechanisms, MAP2K1-specific preclinical models.
> Hypothesis: *Co-occurring MAP2K1 P124L and TP53 mutations may alter MEK-inhibitor sensitivity compared to MAP2K1 alone, motivating co-mutation–stratified studies in MAPK-driven melanoma.*

Structured fields on `result.llm`:

- `llm_summary` – human-readable synthesis
- `evidence_quality` – `"comprehensive" | "moderate" | "limited" | "minimal"`
- `knowledge_gaps` / `well_characterized` – what's missing vs solid
- `research_implications` – 2–3 sentence, testable hypothesis
- `key_references` – PMIDs, trials, and database IDs it relied on

---

## LLM-Ready Context

OncoMind can emit dense, grounded context strings for use in your own LLM workflows:

```python
context = panel.to_knowledge_header()
# "BRAF V600E in melanoma. Oncogenic driver via constitutive MAPK activation.
#  FDA-approved: dabrafenib + trametinib, vemurafenib + cobimetinib.
#  Resistance via NRAS, MEK1/2, BRAF amplification. Sources: CIViC, FDA, PMID:22735384."
```

---

## Configuration

Create a `.env` file:

```bash
# Required for LLM research mode (Claude Sonnet 4 is the default)
ANTHROPIC_API_KEY=your-anthropic-key

# Optional: use OpenAI models instead
OPENAI_API_KEY=your-openai-key

# Optional: better literature context (used in LLM mode)
SEMANTIC_SCHOLAR_API_KEY=your-s2-key
```

**Supported LLM models:**
- `claude-sonnet-4-20250514` (default, recommended)
- `claude-3-5-haiku-20241022` (faster, lower cost)
- `gpt-4o-mini`, `gpt-4o`, `gpt-4-turbo`

---

## Result Model

The top-level `Result` object contains:

- `evidence`: all structured data (identifiers, KBs, functional, clinical, literature)
- `llm`: optional `LLMInsight` with research narrative and gaps (only when `enable_llm=True`)

```python
result = await get_insight("BRAF V600E", tumor_type="Melanoma")

result.identifiers.gene                 # "BRAF"
result.kb.civic_assertions              # CIViC drug–variant assertions
result.clinical.fda_approvals           # FDA therapies
result.evidence.get_recommended_therapies()

if result.llm:
    result.llm.llm_summary
    result.llm.evidence_quality
    result.llm.knowledge_gaps
    result.llm.research_implications
```

See [docs/API_REFERENCE.md](docs/API_REFERENCE.md) for complete field documentation.

---

## Supported Variant Types

Currently supports:

- Missense (e.g., `V600E`, `L858R`)
- Nonsense (e.g., `R248*`)
- Small indels (e.g., `E746_A750del`)
- Frameshift (e.g., `K132fs`)

Variants can be provided as simple protein changes (`V600E`, `p.V600E`) or in HGVS; MyVariant.info is used to normalize where needed.

Planned: fusions, amplifications, copy-number variants (see [docs/ROADMAP.md](docs/ROADMAP.md)).

---

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

See [docs/ROADMAP.md](docs/ROADMAP.md) for the full development roadmap.

---

## Development

```bash
pytest tests/unit/ -v
pytest tests/unit/ --cov=src/oncomind --cov-report=html

mypy src/oncomind

ruff check src/oncomind
ruff format src/oncomind
```

Run the Streamlit app:

```bash
cd streamlit
streamlit run app.py
```

---

## License

MIT License – see [LICENSE](LICENSE).

## Acknowledgments

CIViC, VICC, MyVariant.info, DepMap, Semantic Scholar, and all the open-data projects that make this kind of research tooling possible.

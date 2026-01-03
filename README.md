**OncoMind**

**Research intelligence for cancer variants. Find the gaps, not just the facts.**

For BRAF V600E, databases already agree. For the next 10,000 variants, the key question is **“what don’t we know yet?”**

OncoMind aggregates evidence from 10+ public resources, scores knowledge gaps, and proposes testable hypotheses — turning variant annotation into **research intelligence** for discovery and translational work, not medical advice.

> **Research use only.** Not for diagnosis, treatment selection, or any clinical decision-making.

***

## Why OncoMind?

Most platforms summarize **what’s known** about a variant. Researchers also need to see:

- Which areas are **well-characterized vs under-studied**  
- Where **databases conflict** or leave holes  
- Which **co-mutation patterns** hint at mechanism or resistance  
- What **functional and preclinical experiments** are still missing  

It is a **rapid hypothesis‑generation and gap‑spotting tool** around variants, not as a clinical decision support system.  

OncoMind is built in two layers:

1. A deterministic **annotation backbone** (knowledge bases, trials, cBioPortal, DepMap, literature).  
2. An optional **LLM research layer** that highlights gaps and drafts hypotheses, constrained by that backbone.  

OncoMind is designed for **translational teams and small biotechs** that need to decide *which variants are worth a project*, not how to treat an individual patient. It integrates public knowledge bases, clinical trials, prevalence data, and mechanistic literature into an interpretable, research‑grade summary plus gap analysis that a discovery team can act on.  

Typical users include:

- Small biotechs prioritizing targets, models, or combination strategies  
- Academic labs planning functional studies or resistance screens  
- Platform teams triaging large variant lists from NGS or CRISPR screens  


Using 2–3 contrasted variants (for example, AKT1 E17K in breast cancer, IDH1 R132H in glioma, and ERBB2 S310F in   
bladder cancer) makes it clear that OncoMind handles both rich‑evidence and sparse‑evidence settings and that its   
differentiator is surfacing **where to push next in the biology**, rather than merely restating known facts.  

***

## What Makes It Different

| Feature                  | Typical tools                 | OncoMind                                       |
|--------------------------|-------------------------------|------------------------------------------------|
| Primary question         | “What is this variant?”       | “What don’t we know yet?”                      |
| Knowledge gaps           | Rarely explicit               | First-class outputs                            |
| Co-mutation patterns     | Raw prevalence tables         | Used to drive mechanistic hypotheses           |
| Source conflicts         | Buried in details             | Detected, surfaced, and explained              |
| Output                   | Static clinical-style notes   | Research-ready, gap-focused variant briefs     |
| Evidence quality         | Implicit                      | Explicit per-variant grading and rationale     |

***

## Quick Start

### Install

```bash
git clone https://github.com/yourusername/oncomind.git
cd oncomind/onco_mind_v0

python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

pip install -e ".[dev]"
```

If you hit a `ModuleNotFoundError` after pulling updates, reinstall with:

```bash
pip install -e . --force-reinstall
```

### Python API

```python
import asyncio
from oncomind import get_insight, InsightConfig

async def main():
    # Fast, deterministic evidence (~7s)
    result = await get_insight("BRAF V600E", tumor_type="Melanoma")
    print(result.get_summary())
    print(result.evidence.get_recommended_therapies())

    # Evidence backbone + literature + LLM research layer (~25s)
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
# Annotation backbone (~7s)
mind insight BRAF V600E --tumor Melanoma

# Literature search only (~15s)
mind insight EGFR L858R -t NSCLC --lit

# Research mode: backbone + LLM (~20s)
mind insight MAP2K1 P124L -t Melanoma --llm

# Full: backbone + literature + LLM (~25s)
mind insight KRAS G12D -t CRC --full

# Save to JSON
mind insight EGFR L858R -t NSCLC --llm --output result.json

# Debug logging
mind insight BRAF V600E --log-level DEBUG
```

**Modes**

| Mode       | Flag     | Output                                             |
|------------|----------|----------------------------------------------------|
| Annotation | (none)   | Structured evidence from all data sources          |
| Literature | `--lit`  | + PubMed / Semantic Scholar hits                   |
| LLM        | `--llm`  | + Research narrative and gap analysis              |
| Full       | `--full` | Annotation + literature + LLM layer                |

Logging can be set via CLI or environment:

```bash
mind insight BRAF V600E --log-level DEBUG
mind insight BRAF V600E -l DEBUG  # shorthand

ONCOMIND_LOG_LEVEL=DEBUG mind insight BRAF V600E
mind batch variants.json --log-level DEBUG
```

***

## What You Get

### 1) Evidence Backbone

OncoMind constructs a structured, variant‑centric evidence model:

- **Clinical / KB:** CIViC, VICC MetaKB, ClinVar, COSMIC, CGI, FDA labels  
- **Functional:** AlphaMissense, CADD, PolyPhen2, gnomAD  
- **Biological:** cBioPortal prevalence and co‑mutation structure  
- **Preclinical:** DepMap CRISPR essentiality and PRISM drug response  
- **Trials:** ClinicalTrials.gov  
- **Literature:** PubMed / Semantic Scholar summarized in a `LiteratureEvidence` model  

Match specificity is tracked so you can separate variant‑level from gene‑level signals:

| Match level | Meaning                            | Example                          |
|------------|-------------------------------------|----------------------------------|
| `variant`  | Exact amino-acid change            | BRAF V600E specific data         |
| `codon`    | Same residue, different change     | BRAF V600K in “V600 variants”    |
| `gene`     | Gene-level-only evidence           | “BRAF mutation” basket trials    |

Gap analysis is available programmatically:

```python
gaps = result.evidence.compute_evidence_gaps()

gaps.overall_evidence_quality  # "comprehensive" | "moderate" | "limited" | "minimal"
gaps.well_characterized        # e.g. ["FDA-approved therapies", "Functional impact"]
gaps.poorly_characterized      # e.g. ["Resistance mechanisms", "Preclinical models"]
gaps.gaps                      # structured list with severity + suggested work
```

### 2) Research Insight (LLM Layer)

When enabled, OncoMind adds a research card on top of the evidence backbone:

- `llm_summary` – concise synthesis of function, biology, and therapeutic landscape  
- `evidence_quality` – comprehensive / moderate / limited / minimal  
- `knowledge_gaps` / `well_characterized` – structured view of what’s missing vs solid  
- `research_implications` – short, testable hypotheses  
- `key_references` – PMIDs, trials, and KB IDs supporting the card  

Example (MAP2K1 P124L in melanoma) might include: gain‑of‑function MAPK activation, rarity and co‑mutation context, limited MEK‑inhibitor data, and hypotheses around resistance mechanisms or TP53‑stratified sensitivity.

***

## LLM-Ready Context

OncoMind can emit compact, grounded context strings for use in your own LLM workflows:

```python
context = panel.to_knowledge_header()
# "BRAF V600E in melanoma. Oncogenic driver via constitutive MAPK activation.
#  FDA-approved: dabrafenib + trametinib, vemurafenib + cobimetinib.
#  Resistance via NRAS, MEK1/2, BRAF amplification. Sources: CIViC, FDA, PMID:22735384."
```

***

## Configuration

Create a `.env` file:

```bash
# Required for LLM research mode
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

***

## Result Model

The top-level `Result` object contains:

- `evidence` – all structured data (identifiers, KBs, functional, clinical, literature)  
- `llm` – optional `LLMInsight` with research narrative and gaps (only when `enable_llm=True`)  

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

See `docs/API_REFERENCE.md` for full field documentation.

***

## Supported Variant Types

Currently supports:

- Missense (e.g., `V600E`, `L858R`)  
- Nonsense (e.g., `R248*`)  
- Small indels (e.g., `E746_A750del`)  
- Frameshift (e.g., `K132fs`)  

Variants can be provided as simple protein changes (`V600E`, `p.V600E`) or in HGVS notation; normalization is handled under the hood.

Planned: fusions, amplifications, and copy‑number variants (see `docs/ROADMAP.md`).

***

## Data Sources

### Clinical & Therapeutic

| Source | Data Type | Access |
|--------|-----------|--------|
| CIViC | Curated variant–drug associations | API / dump |
| VICC MetaKB | Aggregated knowledge bases | API |
| ClinVar | Clinical significance | Via aggregation layer |
| COSMIC | Somatic mutation catalog | Via aggregation layer |
| CGI | Biomarker annotations | Local DB |
| FDA | Drug approvals | Public APIs / labels |
| ClinicalTrials.gov | Active and historical trials | Public API |

### Functional & Biological

| Source | Data Type | Access |
|--------|-----------|--------|
| cBioPortal | Co-mutation patterns, prevalence | API |
| AlphaMissense | Pathogenicity predictions | Precomputed scores |
| gnomAD | Population frequencies | Via aggregation layer |

### Literature

| Source | Data Type | Access |
|--------|-----------|--------|
| Semantic Scholar | AI-powered literature search | API |
| PubMed | Biomedical literature | E-utilities |

### Preclinical Research

| Source | Data Type | Access |
|--------|-----------|--------|
| DepMap | Gene essentiality (CRISPR), drug sensitivity (PRISM), cell line models | API / downloads |

### Coming Soon

| Source | Data Type | Status |
|--------|-----------|--------|
| Reactome | Pathway context | Planned |

***

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

***

## License

MIT License – see `LICENSE`.

## Acknowledgments

Built on the work of CIViC, VICC, MyVariant.info, DepMap, Semantic Scholar, cBioPortal, and the broader open‑data oncology community.

[1](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/images/150567280/813a06e7-3022-4ad0-b979-01ad4311d148/image.jpg)
[2](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/images/150567280/a7656279-d3ae-457f-b821-43bc10216bd9/image.jpg)
[3](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/images/150567280/f16c56cc-f47f-4b88-baba-f3aa8c964980/image.jpg)
[4](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/images/150567280/7ede33b8-62a5-4df4-8e6e-07f8593fabbd/image.jpg)
[5](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/images/150567280/f68319f9-791b-4200-a419-4cdaf76e7ab0/image.jpg)
[6](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/images/150567280/f689d40e-3225-40b6-a5e5-954dcfab3b03/image.jpg)
[7](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/images/150567280/5c3f51a7-f556-4043-932e-8e135adf640d/image.jpg)
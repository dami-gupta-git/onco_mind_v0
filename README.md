---
title: OncoMind
emoji: 🚀
colorFrom: blue
colorTo: red
sdk: docker
app_file: app.py
app_port: 8501
pinned: false
tags:
- oncology
- cancer
- genomics
- precision-medicine
- bioinformatics
- llm
- biomedical
- genetics
- variant-analysis
- variant-interpretation
- structured-llm
- evidence-synthesis
- research
- streamlit
- demo
thumbnail: https://huggingface.co/spaces/damigupta/onco_mind/resolve/main/images/starter.png
short_description: Cancer Research Copilot - Gap Analysis
---

# OncoMind
**Cancer Research Copilot for Gap Analysis**

### IN PROGRESS 

---
**Research intelligence for cancer variants. Find the gaps, not just the facts.**

For BRAF V600E, databases already agree. For the next 10,000 variants, the key question is **"what don't we know yet?"**

OncoMind is a research intelligence platform that identifies evidence gaps in cancer variant knowledge—surfacing where research is thin, conflicting, or missing entirely. It's built for translational teams and small biotechs deciding *which variants are worth a project*, not for treating individual patients.

> **Status: Proof-of-Concept / Architectural Demo**
> 
> This demonstrates an approach to systematic evidence gap detection. It is **research use only** — not for diagnosis, treatment selection, or any clinical decision-making.  
>   
> ⚠️ **SNPs and small indels only.** Fusions, amplifications, and copy-number variants are not yet supported.
---

## What Makes It Different

| Feature                        | Typical tools | OncoMind |
|--------------------------------|---------------|----------|
| Primary question               | "What is this variant?" | "What don't we know yet?" |
| Knowledge gaps                 | Rarely explicit | First-class outputs with severity scoring |
| Source conflicts               | Buried in details | Detected, surfaced, and explained |
| Match specificity              | Not tracked | Variant vs codon vs gene-level evidence labeled |
| Source attribution             | Often missing | Every claim linked to PMID, FDA label, or DB entry |
| Cancer hotspots                | Binary yes/no | + Adjacent hotspot detection (±5 codons) |
| Research hypotheses            | None | Generated with evidence basis tags |
| LLM synthesis                  | Generic summaries | Grounded in structured evidence backbone |
| LLM Cross-source drug analysis | Manual comparison | Corroboration, conflicts, and emerging targets surfaced |
| Output                         | Static clinical-style notes | LLM-ready context blocks with receipts |

---

## Screenshots

![Title](https://huggingface.co/spaces/damigupta/onco_mind/resolve/main/images/title.png)

### Gap Analysis
![Gap Analysis](https://huggingface.co/spaces/damigupta/onco_mind/resolve/main/images/gap_large2.png)

<details>

<summary>📸 Click to see more screenshots</summary>

### LLM Research Synthesis
![LLM Synthesis](https://huggingface.co/spaces/damigupta/onco_mind/resolve/main/images/synthesis.png)


</details>

---

## What This Demonstrates

### Architecture & Integration
- Multi-source evidence aggregation (CIViC, ClinVar, CGI, VICC, DepMap, cBioPortal) with conflict detection
- Evidence hierarchy (variant > codon > gene level) with match specificity tracking
- Resistance mechanism annotation using cross-database validation
- Gap detection with severity scoring (CRITICAL → INFORMATIONAL)

### Technical Decisions
- Structured data extraction before LLM synthesis (reduces hallucination risk)
- Evidence provenance tracking across 6+ databases
- Context-aware research hypothesis generation
- Deterministic annotation backbone + optional LLM research layer

---

## Known Limitations (By Design)

This is a proof-of-concept, not production-ready software. Known issues include:

- **Validation**: Needs systematic validation, at the technical level, as well as will as by a SME expert.
- **SNPs and small indels only**: Fusions, amplifications, and copy-number variants are not yet supported.
- **Negation detection:** FDA label parsing may miss negative indicators ("not demonstrated", "not approved")
- **Edge cases:** Rare variant-disease pairings may have inconsistent evidence grading
- **Display formatting:** Some compound identifiers (CAS numbers) may appear in clinical evidence sections
- **LLM variability:** Research hypothesis quality varies; some may be speculative

**Production deployment would require:**
- Systematic validation pipeline with domain expert review
- Robust regulatory text parsing (negation detection, contraindications)
- Automated testing across representative variant sets
- Human-in-the-loop review for high-stakes clinical use cases

---

## Use Cases

### Good for:
- Exploring architectural approaches to clinical evidence synthesis
- Understanding systematic challenges in multi-database genomics integration
- Generating research directions for understudied variants
- Portfolio demonstration of domain expertise + technical execution
- Prioritizing targets, models, or combination strategies for small biotechs
- Planning functional studies or resistance screens for academic labs
- Triaging large variant lists from NGS or CRISPR screens

### Not suitable for:
- Clinical decision-making (use validated tools like OncoKB, CIViC)
- Regulatory submissions
- Production therapeutic recommendations without expert review

---

## Why This Approach?

Traditional variant knowledgebases focus on *summarizing what's known*. OncoMind inverts this: it systematically identifies *what's unknown* to guide research prioritization. The gap detection architecture could inform:
- Research funding decisions
- Clinical trial design
- Functional validation studies

The platform works in two layers:
1. **Deterministic annotation backbone** – structured evidence from knowledge bases, trials, cBioPortal, DepMap, Hotspots, literature
2. **Optional LLM research layer** – highlights gaps and drafts hypotheses, constrained by that backbone

---

## Quick Start

### Install

```bash
git clone https://github.com/dami-gupta-git/onco_mind_v0.git
cd oncomind/onco_mind_v0

python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

pip install -e ".[dev]"
```

If you hit a `ModuleNotFoundError` after pulling updates, reinstall with:

```bash
pip install -e . --force-reinstall
```

### Streamlit UI (Recommended)

The interactive web interface is the easiest way to explore variants:

```bash
cd streamlit
streamlit run app.py
```

Features:
- Enter any gene + variant + tumor type
- Browse evidence by source (CIViC, VICC, CGI, FDA, DepMap, etc.)
- View gap analysis with severity scoring
- See match specificity (variant vs codon vs gene level)
- Cross-source drug analysis identifying corroboration and conflicts
- Export results to JSON

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

| Mode | Flag | Output |
|------|------|--------|
| Annotation | (none) | Structured evidence from all data sources |
| Literature | `--lit` | + PubMed / Semantic Scholar hits |
| LLM | `--llm` | + Research narrative and gap analysis |
| Full | `--full` | Annotation + literature + LLM layer |

---

## What You Get

### Evidence Backbone

OncoMind constructs a structured, variant-centric evidence model:

- **Clinical / KB:** CIViC, VICC MetaKB, ClinVar, COSMIC, CGI, FDA labels
- **Functional:** AlphaMissense, CADD, PolyPhen2, gnomAD
- **Biological:** cBioPortal prevalence and co-mutation structure
- **Preclinical:** DepMap CRISPR essentiality and PRISM drug response
- **Trials:** ClinicalTrials.gov
- **Literature:** PubMed / Semantic Scholar summarized in a `LiteratureEvidence` model

Match specificity is tracked so you can separate variant-level from gene-level signals:

| Match level | Meaning | Example |
|-------------|---------|---------|
| `variant` | Exact amino-acid change | BRAF V600E specific data |
| `codon` | Same residue, different change | BRAF V600K in "V600 variants" |
| `gene` | Gene-level-only evidence | "BRAF mutation" basket trials |

### Research Insight (LLM Layer)

When enabled, OncoMind adds a research card on top of the evidence backbone:

- `llm_summary` – concise synthesis of function, biology, and therapeutic landscape
- `evidence_quality` – comprehensive / moderate / limited / minimal
- `knowledge_gaps` / `well_characterized` – structured view of what's missing vs solid
- `research_implications` – short, testable hypotheses
- `key_references` – PMIDs, trials, and KB IDs supporting the card

### Cross-Source Drug Analysis (LLM Layer)

A separate LLM analysis that synthesizes therapeutic evidence across CGI, CIViC, VICC, and Literature sources:

- **Strongest Evidence** – Drugs with corroboration across multiple independent sources, with biological rationale
- **Conflicting Signals** – Drugs where sources disagree, with likely explanations (tumor type differences, sequential therapy, acquired mutations)
- **Emerging Targets** – Single-source preclinical or early-phase evidence worth investigating
- **Key Gaps** – Expected drugs not found, tumor type extrapolation concerns

This runs in parallel with the main LLM synthesis for faster response times.

---

## Configuration

Create a `.env` file (or set as Hugging Face Secrets):

```bash
# Required for LLM research mode (Gemini 2.0 Flash is the default)
GOOGLE_API_KEY=your-google-api-key
# Optional: use Anthropic or OpenAI models instead
ANTHROPIC_API_KEY=your-anthropic-key
# Optional: use OpenAI models instead
OPENAI_API_KEY=your-openai-key
# Optional: better literature context 
SEMANTIC_SCHOLAR_API_KEY=your-s2-key
```

**Supported LLM models:**


- `gemini/gemini-2.0-flash` (default, fast)
- `gemini/gemini-1.5-pro`
- `claude-sonnet-4-20250514`
- `claude-3-5-haiku-20241022`
- `claude-sonnet-4-20250514` 
- `claude-3-5-haiku-20241022` 
- `gpt-4o-mini`, `gpt-4o`, `gpt-4-turbo`

---

## Supported Variant Types

Currently supports:

- Missense (e.g., `V600E`, `L858R`)
- Nonsense (e.g., `R248*`)
- Small indels (e.g., `E746_A750del`)
- Frameshift (e.g., `K132fs`)

Variants can be provided as simple protein changes (`V600E`, `p.V600E`) or in HGVS notation; normalization is handled under the hood.

**Planned:** fusions, amplifications, and copy-number variants.

---

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

---

## Development

```bash
pytest tests/unit/ -v
pytest tests/unit/ --cov=src/oncomind --cov-report=html

mypy src/oncomind
ruff check src/oncomind
ruff format src/oncomind
```

---

## Validation

| Variant       | Tumor Type | Source | Expected Data | OncoMind Match | Error Description |
|---------------|------------|--------|---------------|----------------|-------------------|
| BRAF V600E    | Melanoma   | FDA Labels | Trametinib Mekinist, Trametinib + Dabrafenib Mekinist + Tafinlar, Encorafenib + Binimetinib Braftovi + Mektovi, Vemurafenib Zelboraf, Cobimetinib + Vemurafenib Cotellic + Zelboraf, Dabrafenib Tafinlar, Atezolizumab + Cobimetinib + Vemurafenib Tecentriq + Cotellic + Zelboraf, Atezolizumab and Hyaluronidase-tqjs + Cobimetinib + Vemurafenib Tecentriq Hybreza + Cotellic + Zelboraf | ✓              |                   |
| BRAF V600E    | Melanoma   | CIViC | Level A, 5 evidence items | ✓              |                   |
| BRAF V600E    | Melanoma   | ClinVar | Pathogenic | ✓              |                   |
| BRAF V600E    | Melanoma   | cBioPortal | Melanoma (MSK, Clin Cancer Res 2021) | ✓              |                   |
| EGFR L858R    | NSCLC      | FDA Labels | Gefitinib Iressa, Erlotinib Tarceva, Afatinib Gilotrif, Dacomitinib Vizimpro, Osimertinib Tagrisso, Osimertinib + pemetrexed and platinum-based chemotherapy Tagrisso, Amivantamab Rybrevant + lazertinib Lazcluze, Amivantamab Rybrevant + carboplatin + pemetrexed, Amivantamab and hyaluronidase-lpuj Rybrevant Faspro | ✓              |                   |
| EGFR L858R    | NSCLC      | CIViC | Level A, 14 evidence items | ✓              |                   |
| EGFR L858R    | NSCLC      | DepMap | Gene essentiality - EGFR is not essential | ✓              |                   |
| KRAS G12C     | NSCLC      | FDA Labels | Sotorasib, adagrasib | ✓              |                   |
| KRAS G12C     | NSCLC      | ClinVar | Pathogenic | ✓              |                   |
| PIK3CA H1047R | Breast     | FDA Labels | Inavolisib + palbociclib + fulvestrant, Alpelisib + fulvestrant, Capivasertib + fulvestrant | ✓              |                   |
| PIK3CA H1047R | Breast     | CIViC | 3 Level A variant specific items | ✓              |                   |
| IDH1 R132H    | Glioma     | FDA Labels | Vorasidenib (Voranigo) | ✓              |                   |
| IDH1 R132H    | Glioma     | ClinVar | Pathogenic | ✓              |                   |
| IDH1 R132H    | Glioma     | Hotspots | This variant is at known cancer hotspot | ✓              |                   |
| ERBB2 S310F   | Breast     | CIViC | Level B, 3 evidence items | ✓              |                   |


**Summary:**

- Total checks: 15
- Matches: 15 (100%)

---

## License

MIT License – see `LICENSE`.

---

## Acknowledgments

Built on the work of CIViC, VICC, MyVariant.info, DepMap, Semantic Scholar, cBioPortal, and the broader open-data oncology community.
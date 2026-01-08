# OncoMind Streamlit Application

Clean, single-container Streamlit implementation of OncoMind variant insight generation.

> ⚠️ **SNPs and small indels only.** Fusions, amplifications, and copy-number variants are not yet supported.

## Quick Start

1. **Set up environment variables:**
```bash
cp .env.example .env
# Edit .env with your API keys
```

2. **Start the application:**
```bash
cd streamlit
docker compose up --build
```

3. **Open in browser:**
```
http://localhost:8502
```

## Features

### 🔬 Single Variant Insight

Two modes available:

| Mode | Speed | Features |
|------|-------|----------|
| **Annotation** (default) | ~7s | Fast structured evidence from all databases |
| **LLM** (toggle on) | ~25s | + Literature search + AI synthesis + hypothesis generation |

**Annotation Mode** provides:
- Evidence summary with pathogenicity scores
- Recommended therapies (FDA-approved, clinical, preclinical)
- Database identifiers (COSMIC, ClinVar, dbSNP)
- HGVS notations
- Functional annotations (AlphaMissense, CADD, PolyPhen2)
- Cancer hotspot detection (MSK Cancer Hotspots) with adjacent hotspot flagging
- cBioPortal co-mutation patterns
- DepMap preclinical data (gene essentiality, drug sensitivity, cell line models)
- Evidence from CIViC, VICC MetaKB, CGI with match specificity tracking

**LLM Mode** adds:
- Literature search via Semantic Scholar
- AI-powered research synthesis
- Evidence gap analysis with severity scoring
- Well-characterized aspects with basis explanations
- Research implications and hypothesis generation
- Cross-source drug analysis (corroboration, conflicts, emerging targets)

### 📊 Batch Upload
- Upload CSV with variant data
- Process multiple variants concurrently
- Download results as CSV or JSON
- Real-time progress tracking

### ✅ Validation
- Run against gold standard dataset
- Get accuracy metrics
- Per-tier precision/recall/F1 scores
- Failure analysis

## CSV Format for Batch Upload

```csv
gene,variant,tumor_type
BRAF,V600E,Melanoma
EGFR,L858R,Lung Adenocarcinoma
KRAS,G12D,Colorectal Cancer
```

- Required columns: `gene`, `variant`
- Optional column: `tumor_type`

## Development

To run locally without Docker:

```bash
# Install the oncomind package
cd ..
pip install -e .

# Install streamlit dependencies
cd streamlit
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

## Environment Variables

- `OPENAI_API_KEY`: OpenAI API key
- `ANTHROPIC_API_KEY`: Anthropic API key
- `GOOGLE_API_KEY`: Google API key
- `GROQ_API_KEY`: Groq API key

## Volume Mounts

- `../benchmarks:/app/benchmarks`: Gold standard datasets for validation
- `../data:/app/data`: Additional data files

## Notes

- App runs on port 8502
- All scientific functionality from the original oncomind package is preserved
- Backend uses async/await for efficient concurrent processing
- LiteLLM provides unified interface to multiple LLM providers
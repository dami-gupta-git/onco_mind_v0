# OncoMind Streamlit Application

Clean, single-container Streamlit implementation of OncoMind variant insight generation.

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
- Input gene, variant, and tumor type (The tumor type should exactly match values from the OncoTree ontology or CIViC database)
- Select LLM model (OpenAI, Anthropic, Google, Groq)
- Get comprehensive insight with:
  - Evidence summary with pathogenicity scores
  - Recommended therapies (FDA-approved, clinical, preclinical)
  - Database identifiers (COSMIC, ClinVar, dbSNP)
  - HGVS notations
  - Functional annotations (AlphaMissense, CADD, PolyPhen2)
  - cBioPortal co-mutation patterns
  - 🧬 DepMap preclinical data (gene essentiality, drug sensitivity, cell line models)
  - Evidence gap analysis and research implications

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
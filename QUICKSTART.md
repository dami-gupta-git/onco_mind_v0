# OncoMind Quick Start Guide

This guide will get you up and running with OncoMind in 5 minutes.

## Choose Your Interface

OncoMind is available as:
- **🐳 Docker**: One-command deployment (easiest!)
- **🌐 Web Application**: Modern Angular UI with Flask REST API
- **💻 Command-Line Interface**: Python CLI tool for batch processing

## 🐳 Docker Quick Start (Recommended for POC)

### Prerequisites
- Docker
- OpenAI API key

### 3-Command Setup

```bash
# 1. Set API key
echo "OPENAI_API_KEY=your-key-here" > .env

# 2. Start container
docker compose up

# 3. Open in browser (wait ~30s for build)
open http://localhost:4200
```

Done! See [DOCKER.md](DOCKER.md) for advanced setup.

---

## Web Application Quick Start (Manual)

### Prerequisites
- Python 3.11+
- Node.js 18+ and npm
- OpenAI API key

### 5-Minute Setup

1. **Setup Python environment**:
```bash
cd tumor_board_v0
python -m venv venv
source venv/bin/activate
pip install -e .
pip install -r backend/requirements.txt
```

2. **Configure API key**:
```bash
cp .env.example .env
# Edit .env and add: OPENAI_API_KEY=your-key-here
```

3. **Start both servers** (easiest method):
```bash
chmod +x start-dev.sh
./start-dev.sh
```

4. **Open in browser**: http://localhost:4200

That's it! See the [streamlit README](streamlit/README.md) for detailed web app documentation.

---

## CLI Quick Start

### Prerequisites

- Python 3.11 or higher
- An OpenAI API key (or Anthropic/other LLM provider)

## Installation

1. Clone and navigate to the directory:
```bash
cd tumor_board
```

2. Create a virtual environment with Python 3.11:
```bash
python3.11 -m venv venv
source venv/bin/activate
```

3. Install the package:
```bash
pip install -e ".[dev]"
```

4. Create a `.env` file in the project root with your API key:
```bash
echo "OPENAI_API_KEY=sk-..." > .env
```

The `.env` file is automatically loaded when you run `oncomind` commands!

## Your First Insight

Get insight on the BRAF V600E mutation in melanoma:

```bash
mind insight BRAF V600E --tumor Melanoma
```

### CLI Modes

| Mode | Flag | Speed | Output |
|------|------|-------|--------|
| **Default** | (none) | ~12s | Structured evidence + LLM clinical summary |
| **Lite** | `--lite` | ~7s | Structured evidence only (no LLM) |
| **Full** | `--full` | ~25s | + Literature search + enhanced narrative |

```bash
# Fast mode - just the structured evidence
mind insight EGFR L858R -t NSCLC --lite

# Full mode - includes literature search
mind insight KRAS G12C -t NSCLC --full

# Save output to JSON
mind insight BRAF V600E -t Melanoma --output result.json
```

## Batch Processing

1. Create a file `my_variants.json`:
```json
[
  {
    "gene": "BRAF",
    "variant": "V600E",
    "tumor_type": "Melanoma"
  },
  {
    "gene": "EGFR",
    "variant": "L858R",
    "tumor_type": "Lung Adenocarcinoma"
  },
  {
    "gene": "KRAS",
    "variant": "G12C",
    "tumor_type": "Non-Small Cell Lung Cancer"
  }
]
```

2. Run batch processing:
```bash
mind batch my_variants.json --output my_results.json

# Fast mode (no LLM)
mind batch my_variants.json --output my_results.json --lite

# Full mode (with literature)
mind batch my_variants.json --output my_results.json --full
```

3. View results:
```bash
cat my_results.json
```

## Try Different Models

Use a different LLM model for the narrative:

```bash
# Use GPT-4o for better accuracy
mind insight BRAF V600E --tumor Melanoma --model gpt-4o

# Use Claude (requires ANTHROPIC_API_KEY)
export ANTHROPIC_API_KEY="sk-ant-..."
mind insight BRAF V600E --tumor Melanoma --model claude-3-5-sonnet-20241022
```

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Explore the [gold_standard.json](benchmarks/bak/gold_standard.json) to understand the benchmark dataset
- Try the [sample_batch.json](benchmarks/sample_batch.json) for more examples
- Run the test suite: `pytest`
- Experiment with different prompts by modifying [src/oncomind/llm/prompts.py](src/oncomind/llm/prompts.py)

## Common Issues

**Issue**: `ModuleNotFoundError: No module named 'oncomind'`
**Solution**: Make sure you installed with `pip install -e .` from the project root

**Issue**: `litellm.exceptions.AuthenticationError`
**Solution**: Check that your API key is set correctly: `echo $OPENAI_API_KEY`

**Issue**: `MyVariantAPIError`
**Solution**: This is usually a network issue or the variant isn't in the database. Check your internet connection.

## Getting Help

- Check the [README.md](README.md) for full documentation
- Look at test files in `tests/` for usage examples
- Open an issue on GitHub if you find bugs

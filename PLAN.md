# Plan

## Progress Overview

- Phase 1: Core Infrastructure -- Complete
- Phase 2: Evidence Sources -- Complete
- Phase 3: LLM Synthesis Layer -- Complete
- Phase 4: Gap Detection & Analysis -- Complete
- Phase 5: User Interfaces -- Complete
- Phase 6: Validation & Hardening -- In Progress
- Phase 7: Expansion -- Not Started

## Phase 1: Core Infrastructure

- [x] Project scaffolding (pyproject.toml, hatchling build)
- [x] Variant normalization and input parsing
- [x] Amino acid code conversion (3-letter/1-letter)
- [x] Variant type classification (missense, nonsense, indel, frameshift)
- [x] Async HTTP client infrastructure (httpx)
- [x] Centralized constants and configuration
- [x] Logging module (configurable via env var or CLI flag)
- [x] Pydantic data models for all evidence types

## Phase 2: Evidence Sources

- [x] MyVariant.info integration (ClinVar, COSMIC, gnomAD, AlphaMissense, CADD)
- [x] CIViC API client (variant-drug assertions by evidence level)
- [x] VICC MetaKB client (aggregated knowledge bases)
- [x] CGI biomarker annotations (local TSV lookup)
- [x] FDA drug label parsing (biomarker-drug associations)
- [x] cBioPortal integration (prevalence, co-mutations, tumor-specific studies)
- [x] DepMap integration (CRISPR essentiality, PRISM drug sensitivity, cell line models)
- [x] Cancer Hotspots API
- [x] ClinicalTrials.gov client
- [x] PubMed literature search
- [x] Semantic Scholar literature search
- [x] OncoTree tumor type ontology
- [x] EvidenceAggregator with parallel fetching and graceful degradation

## Phase 3: LLM Synthesis Layer

- [x] LLM service with multi-provider support (litellm)
- [x] Research dossier synthesis (5-section structured output)
- [x] Paper relevance scoring (fast model)
- [x] Variant knowledge extraction from literature
- [x] Cross-source drug analysis (separate LLM call)
- [x] Prompt engineering for grounded synthesis (no hallucination of drugs/data)
- [x] Match specificity awareness in prompts (variant vs codon vs gene level)
- [x] JSON response parsing with repair logic

## Phase 4: Gap Detection & Analysis

- [x] Evidence gap detection with severity scoring (CRITICAL through INFORMATIONAL)
- [x] Gap categories (therapeutic, functional, biological, literature)
- [x] Overall evidence quality assessment
- [x] Research priority computation
- [x] Well-characterized vs knowledge gap identification
- [x] Conflicting evidence detection
- [x] Acquired resistance mutation handling (context-aware)
- [x] Targetable sensitizing variant tracking (avoid false conflict flags)

## Phase 5: User Interfaces

- [x] CLI with `insight` command (single variant)
- [x] CLI with `batch` command (JSON input, multiple variants)
- [x] CLI with `annotate` command (VCF file annotation)
- [x] CLI with `version` command
- [x] Rich terminal output (panels, color coding, gap severity indicators)
- [x] Streamlit web UI with tabbed evidence display
- [x] LLM synthesis rendering in Streamlit
- [x] JSON export from both CLI and UI
- [x] Docker deployment (Dockerfile for HuggingFace Spaces)

## Phase 6: Validation & Hardening

- [x] Manual validation across 15 variant/source checks (100% match)
- [x] Unit tests for LLM service (mocked)
- [ ] Systematic validation pipeline with domain expert review
- [ ] Automated testing across representative variant sets
- [ ] Negation detection in FDA label parsing
- [ ] Edge case handling for rare variant-disease pairings

## Phase 7: Expansion

- [ ] Structural variant support (fusions, amplifications, copy-number variants)
- [ ] Additional tumor types beyond current coverage
- [ ] Pre-fetching and caching for top cancer genes
- [ ] CSV and Markdown output formats (CLI currently outputs JSON only)
- [ ] Per-variant output splitting in batch mode

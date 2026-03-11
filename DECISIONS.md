# Decisions

## LLM as Synthesis Layer, Not Classifier

- **Date**: 2024 (pre-project, from TumorBoard Lite experience)
- **Status**: Accepted
- **Context**: TumorBoard Lite showed that LLMs are poor at rules-based variant classification but good at evidence synthesis and narrative generation.
- **Decision**: Use LLM only for synthesis on top of structured evidence, never for classification or data generation.
- **Rationale**: Deterministic annotation backbone ensures correctness; LLM adds interpretive value without inventing data. All LLM outputs are grounded in provided evidence with explicit provenance.

## Multi-Provider LLM via litellm

- **Date**: 2024
- **Status**: Accepted
- **Context**: Need to support multiple LLM providers (Google, Anthropic, OpenAI) without vendor lock-in and with different models for different tasks.
- **Decision**: Use litellm as the LLM abstraction layer. Default model is `gemini/gemini-2.0-flash` for synthesis; `gpt-4o-mini` for high-volume paper scoring.
- **Rationale**: litellm provides a unified API across providers. Using a cheaper/faster model for paper scoring keeps costs manageable while the primary synthesis model can be swapped freely.

## Single-Stage LLM Pipeline

- **Date**: 2025
- **Status**: Accepted
- **Context**: Earlier design considered multi-stage LLM pipelines (extract then synthesize). Single-stage research dossier proved sufficient.
- **Decision**: One LLM call produces a 5-section research dossier (functional impact, tumor biology, therapeutic landscape, evidence quality, research program). Cross-source drug analysis is a separate parallel call.
- **Rationale**: Reduces latency, cost, and complexity. The structured evidence input is rich enough that a single well-prompted call produces high-quality output.

## Parallel Evidence Fetching with Graceful Degradation

- **Date**: 2024
- **Status**: Accepted
- **Context**: 12+ external APIs must be queried per variant. Sequential fetching would be unacceptably slow; any single source failure should not block results.
- **Decision**: Use `asyncio.gather` with per-source exception handling. Failed sources log warnings but do not prevent result assembly.
- **Rationale**: Keeps response time at ~7s for annotation mode. Individual API outages are common and should not block the user.

## Match Specificity Tracking (Variant vs Codon vs Gene Level)

- **Date**: 2025
- **Status**: Accepted
- **Context**: Evidence from databases often matches at different specificity levels. A CIViC entry for "BRAF V600" applies differently than one for "BRAF V600E" specifically.
- **Decision**: Tag all evidence with match level (variant, codon, gene). Surface this in gap analysis and LLM prompts.
- **Rationale**: Prevents overconfident claims based on gene-level evidence when variant-specific data is absent. Critical for clinical genomics accuracy.

## Gap Detection as First-Class Feature

- **Date**: 2025
- **Status**: Accepted
- **Context**: Existing tools focus on summarizing what is known. For research prioritization, knowing what is unknown is equally valuable.
- **Decision**: Build a rule-based gap detection system with severity levels (CRITICAL through INFORMATIONAL). Compute gaps deterministically before LLM sees anything.
- **Rationale**: Gap detection is the core differentiator. Making it deterministic ensures reproducibility and avoids LLM-generated gap assessments that could miss or hallucinate issues.

## Acquired Resistance Mutation Handling

- **Date**: 2025
- **Status**: Accepted
- **Context**: VICC may report "resistance" for variants like EGFR T790M + erlotinib. But T790M is an acquired resistance mutation -- erlotinib works initially, and resistance develops later. Flagging this as "conflicting evidence" (FDA approval vs VICC resistance) would be misleading.
- **Decision**: Maintain a curated `ACQUIRED_RESISTANCE_MUTATIONS` mapping in constants.py, with optional tumor-type restriction. Similarly, `TARGETABLE_SENSITIZING_VARIANTS` maps drugs to their known sensitizing targets.
- **Rationale**: Prevents false conflict detection. Clinical nuance (acquired vs intrinsic resistance) cannot be derived from database labels alone.

## SNPs and Small Indels Only (Initial Scope)

- **Date**: 2024
- **Status**: Accepted (pending expansion)
- **Context**: Full variant type support (fusions, amplifications, CNVs) requires significant additional normalization and API integration work.
- **Decision**: Ship with SNPs and small indels first. Structural variants are tracked as a known limitation and planned expansion.
- **Rationale**: Covers the majority of actionable variants in solid tumors. Allows shipping a working tool while building structural variant support incrementally.

## Pydantic for All Data Models

- **Date**: 2024
- **Status**: Accepted
- **Context**: Evidence data flows through many layers (API -> aggregator -> models -> LLM prompt -> CLI/UI). Type safety and serialization are critical.
- **Decision**: Use Pydantic BaseModel for all evidence types, Result, Evidence, LLMInsight.
- **Rationale**: Pydantic provides validation, JSON serialization (`model_dump`), and clear field documentation. Essential for a clinical genomics tool where data integrity matters.

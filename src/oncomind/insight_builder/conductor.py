"""Conductor for orchestrating evidence aggregation and LLM insight generation.

This module provides the Conductor class that coordinates:
1. Evidence aggregation via EvidenceAggregator
2. LLM narrative synthesis via LLMService

ARCHITECTURE:
    Conductor
        ├── EvidenceAggregator → Evidence (structured data)
        └── LLMService → LLMInsight (narrative)
                ↓
            Result(evidence=Evidence, llm=LLMInsight)
"""

from dataclasses import dataclass
from typing import Callable

from oncomind.insight_builder.evidence_aggregator import (
    EvidenceAggregator,
    EvidenceAggregatorConfig,
)
from oncomind.models.evidence import Evidence
from oncomind.models.result import Result
from oncomind.normalization import ParsedVariant, parse_variant_input


@dataclass
class ConductorConfig:
    """Configuration for the Conductor.

    Controls evidence sources, LLM integration, and processing options.
    """

    # Evidence source toggles
    enable_vicc: bool = True
    enable_civic_assertions: bool = True
    enable_clinical_trials: bool = True
    enable_literature: bool = True

    # LLM configuration
    enable_llm: bool = False  # Disabled by default for fast annotation
    llm_model: str = "gpt-4o-mini"
    llm_temperature: float = 0.1

    # Result limits
    max_vicc_results: int = 15
    max_civic_assertions: int = 20
    max_clinical_trials: int = 10
    max_literature_results: int = 6

    def to_aggregator_config(self) -> EvidenceAggregatorConfig:
        """Convert to EvidenceAggregatorConfig."""
        return EvidenceAggregatorConfig(
            enable_vicc=self.enable_vicc,
            enable_civic_assertions=self.enable_civic_assertions,
            enable_clinical_trials=self.enable_clinical_trials,
            enable_literature=self.enable_literature,
            max_vicc_results=self.max_vicc_results,
            max_civic_assertions=self.max_civic_assertions,
            max_clinical_trials=self.max_clinical_trials,
            max_literature_results=self.max_literature_results,
        )


class Conductor:
    """Orchestrates evidence aggregation and LLM insight generation.

    Use as an async context manager:

        async with Conductor(config) as conductor:
            insight = await conductor.run("BRAF V600E", tumor_type="Melanoma")

    Or for batch processing:

        async with Conductor(config) as conductor:
            insights = await conductor.run_batch(
                ["BRAF V600E", "EGFR L858R"],
                tumor_type="NSCLC"
            )
    """

    def __init__(self, config: ConductorConfig | None = None):
        self.config = config or ConductorConfig()
        self._aggregator: EvidenceAggregator | None = None

    async def __aenter__(self) -> "Conductor":
        aggregator_config = self.config.to_aggregator_config()
        self._aggregator = EvidenceAggregator(aggregator_config)
        await self._aggregator.__aenter__()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb) -> None:
        if self._aggregator:
            await self._aggregator.__aexit__(exc_type, exc_val, exc_tb)
            self._aggregator = None

    async def run(
        self,
        variant: ParsedVariant | str,
        tumor_type: str | None = None,
    ) -> Result:
        """Run the full pipeline for a single variant.

        Args:
            variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
            tumor_type: Optional tumor type for clinical context

        Returns:
            Result with evidence and optional LLM narrative
        """
        if self._aggregator is None:
            raise RuntimeError("Conductor must be used as async context manager")

        # Parse string input if needed
        if isinstance(variant, str):
            parsed = parse_variant_input(variant, tumor_type)
        else:
            parsed = variant
            if tumor_type:
                parsed.tumor_type = tumor_type

        # Step 1: Aggregate evidence
        evidence = await self._aggregator.build_evidence(parsed, parsed.tumor_type)

        # Step 2: Generate LLM narrative if enabled
        llm_insight = None
        if self.config.enable_llm:
            llm_insight = await self._generate_llm_insight(evidence)

        return Result(evidence=evidence, llm=llm_insight)

    async def run_batch(
        self,
        variants: list[ParsedVariant | str],
        tumor_type: str | None = None,
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> list[Result]:
        """Run the full pipeline for multiple variants.

        Args:
            variants: List of ParsedVariant objects or variant strings
            tumor_type: Optional tumor type (applied to all variants)
            progress_callback: Optional callback(current, total) for progress

        Returns:
            List of Result objects
        """
        if self._aggregator is None:
            raise RuntimeError("Conductor must be used as async context manager")

        results = []
        total = len(variants)

        for i, variant in enumerate(variants):
            if progress_callback:
                progress_callback(i, total)

            try:
                result = await self.run(variant, tumor_type)
                results.append(result)
            except Exception as e:
                gene = variant if isinstance(variant, str) else f"{variant.gene} {variant.variant}"
                print(f"  Warning: Failed to process {gene}: {str(e)}")

        if progress_callback:
            progress_callback(total, total)

        return results

    async def _generate_llm_insight(self, evidence: Evidence):
        """Generate LLM narrative insight from evidence.

        Args:
            evidence: Evidence with aggregated data

        Returns:
            LLMInsight with narrative summary
        """
        from oncomind.llm.service import LLMService

        llm_service = LLMService(
            model=self.config.llm_model,
            temperature=self.config.llm_temperature,
        )

        # Get evidence summary for LLM
        evidence_summary = evidence.get_evidence_summary_for_llm()

        # Generate LLM insight
        return await llm_service.get_llm_insight(
            gene=evidence.identifiers.gene,
            variant=evidence.identifiers.variant,
            tumor_type=evidence.context.tumor_type,
            evidence_summary=evidence_summary,
            has_clinical_trials=bool(evidence.clinical_trials),
        )


async def conduct(
    variant: ParsedVariant | str,
    tumor_type: str | None = None,
    config: ConductorConfig | None = None,
) -> Result:
    """Convenience function to run the full pipeline for a single variant.

    This is the recommended entry point for single-variant analysis.

    Args:
        variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
        tumor_type: Optional tumor type for clinical context
        config: Optional configuration for the pipeline

    Returns:
        Result with evidence and optional LLM narrative

    Example:
        >>> result = await conduct("BRAF V600E", tumor_type="Melanoma")
        >>> print(result.identifiers.gene)
        BRAF
    """
    async with Conductor(config) as conductor:
        return await conductor.run(variant, tumor_type)


async def conduct_batch(
    variants: list[ParsedVariant | str],
    tumor_type: str | None = None,
    config: ConductorConfig | None = None,
    progress_callback: Callable[[int, int], None] | None = None,
) -> list[Result]:
    """Convenience function to run the full pipeline for multiple variants.

    Args:
        variants: List of ParsedVariant objects or variant strings
        tumor_type: Optional tumor type (applied to all variants)
        config: Optional configuration for the pipeline
        progress_callback: Optional callback(current, total) for progress

    Returns:
        List of Result objects
    """
    async with Conductor(config) as conductor:
        return await conductor.run_batch(variants, tumor_type, progress_callback)


__all__ = [
    "Conductor",
    "ConductorConfig",
    "conduct",
    "conduct_batch",
]

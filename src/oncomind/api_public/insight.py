"""Public API for variant insight.

This module provides the main entry points for OncoMind:
- get_insight(): Get insight for a single variant
- get_insights(): Batch get insights for multiple variants

ARCHITECTURE:
    Input (str/VCF/DataFrame) → parse → build_insight → Insight
                                              ↓ (optional)
                                       LLM summarization

The API separates:
1. Core evidence aggregation (deterministic, always runs)
2. LLM narrative synthesis (optional, controlled by config)
"""

import asyncio
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from oncomind.insight_builder import (
    InsightBuilder,
    InsightBuilderConfig,
    build_insight as _build_insight,
    build_insights as _build_insights,
)
from oncomind.models.insight import Insight
from oncomind.normalization import (
    parse_variant_input,
    parse_variant_row,
    ParsedVariant,
)


@dataclass
class InsightConfig:
    """Configuration for the insight pipeline.

    Controls data sources, LLM integration, and output options.

    Example:
        >>> config = InsightConfig(
        ...     enable_llm=True,
        ...     llm_model="gpt-4o-mini",
        ...     enable_literature=True,
        ... )
        >>> panel = await get_insight("BRAF V600E", config=config)
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

    # Processing options
    validate_variant_type: bool = True  # Reject fusions/amplifications

    def to_builder_config(self) -> InsightBuilderConfig:
        """Convert to InsightBuilderConfig."""
        return InsightBuilderConfig(
            enable_vicc=self.enable_vicc,
            enable_civic_assertions=self.enable_civic_assertions,
            enable_clinical_trials=self.enable_clinical_trials,
            enable_literature=self.enable_literature,
            max_vicc_results=self.max_vicc_results,
            max_civic_assertions=self.max_civic_assertions,
            max_clinical_trials=self.max_clinical_trials,
            max_literature_results=self.max_literature_results,
        )


async def get_insight(
    variant_str: str,
    tumor_type: str | None = None,
    config: InsightConfig | None = None,
) -> Insight:
    """Get insight for a single variant and return aggregated evidence.

    This is the main entry point for single-variant analysis.

    Args:
        variant_str: Variant string in any supported format:
            - "BRAF V600E" (gene + variant)
            - "EGFR L858R in NSCLC" (with tumor type)
            - "BRAF:V600E" (colon-separated)
            - "TP53 p.R248W" (with p. notation)
        tumor_type: Optional tumor type for clinical context.
            Overrides any tumor type extracted from the string.
        config: Optional configuration for the insight pipeline.

    Returns:
        Insight with all aggregated evidence.

    Raises:
        ValueError: If variant string cannot be parsed or variant type
            is not supported (when validate_variant_type=True).

    Example:
        >>> insight = await get_insight("BRAF V600E", tumor_type="Melanoma")
        >>> print(insight.identifiers.gene, insight.identifiers.variant)
        BRAF V600E
        >>> print(insight.clinical.get_approved_drugs())
        ['Dabrafenib', 'Vemurafenib', 'Encorafenib']

        >>> insight = await get_insight("EGFR L858R in lung cancer")
        >>> print(insight.clinical.tumor_type)
        lung cancer
    """
    config = config or InsightConfig()

    # Parse the variant string
    parsed = parse_variant_input(variant_str, tumor_type)

    # Validate variant type if enabled
    if config.validate_variant_type:
        from oncomind.utils.variant_normalization import VariantNormalizer
        if parsed.variant_type and parsed.variant_type not in VariantNormalizer.ALLOWED_VARIANT_TYPES:
            raise ValueError(
                f"Variant type '{parsed.variant_type}' is not supported. "
                f"Only SNPs and small indels are allowed (missense, nonsense, insertion, deletion, frameshift). "
                f"Got variant: {variant_str}"
            )

    # Build the insight
    builder_config = config.to_builder_config()
    insight = await _build_insight(parsed, parsed.tumor_type, builder_config)

    # Apply LLM enhancement if enabled
    if config.enable_llm:
        insight = await _apply_llm_enhancement(insight, config)

    return insight


async def get_insights(
    variants: list[str] | list[dict] | pd.DataFrame | Path | str,
    tumor_type: str | None = None,
    config: InsightConfig | None = None,
    progress_callback: Callable[[int, int], None] | None = None,
) -> list[Insight]:
    """Get insights for multiple variants and return aggregated evidence.

    Supports various input formats:
    - List of variant strings: ["BRAF V600E", "EGFR L858R"]
    - List of dictionaries: [{"gene": "BRAF", "variant": "V600E"}, ...]
    - pandas DataFrame with gene/variant columns
    - Path to CSV file
    - Path to VCF file (basic support)

    Args:
        variants: Variants in any supported format.
        tumor_type: Optional tumor type (applied to all variants unless
            the variant has its own tumor_type).
        config: Optional configuration for the insight pipeline.
        progress_callback: Optional callback(current, total) for progress updates.

    Returns:
        List of Insight objects.

    Example:
        >>> insights = await get_insights(["BRAF V600E", "EGFR L858R"])
        >>> for insight in insights:
        ...     print(f"{insight.identifiers.gene} {insight.identifiers.variant}")

        >>> insights = await get_insights(
        ...     "variants.csv",
        ...     progress_callback=lambda i, t: print(f"{i}/{t}")
        ... )
    """
    config = config or InsightConfig()
    parsed_variants: list[ParsedVariant] = []

    # Parse input based on type
    if isinstance(variants, (str, Path)):
        path = Path(variants)
        if path.suffix.lower() == ".csv":
            parsed_variants = _parse_csv(path, tumor_type)
        elif path.suffix.lower() == ".vcf":
            parsed_variants = _parse_vcf(path)
        else:
            raise ValueError(f"Unsupported file type: {path.suffix}")

    elif isinstance(variants, pd.DataFrame):
        parsed_variants = _parse_dataframe(variants, tumor_type)

    elif isinstance(variants, list):
        if not variants:
            return []

        if isinstance(variants[0], str):
            for v in variants:
                parsed = parse_variant_input(v, tumor_type)
                parsed_variants.append(parsed)
        elif isinstance(variants[0], dict):
            for v in variants:
                parsed = parse_variant_row(v)
                if tumor_type and not parsed.tumor_type:
                    parsed.tumor_type = tumor_type
                parsed_variants.append(parsed)
        else:
            raise ValueError(f"Unsupported variant type: {type(variants[0])}")

    else:
        raise ValueError(f"Unsupported input type: {type(variants)}")

    # Validate variant types if enabled
    if config.validate_variant_type:
        from oncomind.utils.variant_normalization import VariantNormalizer
        valid_variants = []
        for parsed in parsed_variants:
            if parsed.variant_type and parsed.variant_type not in VariantNormalizer.ALLOWED_VARIANT_TYPES:
                print(f"  Skipping unsupported variant type: {parsed.gene} {parsed.variant} ({parsed.variant_type})")
                continue
            valid_variants.append(parsed)
        parsed_variants = valid_variants

    # Build insights
    builder_config = config.to_builder_config()

    async with InsightBuilder(builder_config) as builder:
        insights = []
        total = len(parsed_variants)

        for i, parsed in enumerate(parsed_variants):
            if progress_callback:
                progress_callback(i, total)

            try:
                insight = await builder.build_insight(
                    parsed, parsed.tumor_type or tumor_type
                )
                insights.append(insight)
            except Exception as e:
                print(f"  Warning: Failed to process {parsed.gene} {parsed.variant}: {str(e)}")

        if progress_callback:
            progress_callback(total, total)

    # Apply LLM enhancement if enabled
    if config.enable_llm:
        insights = await asyncio.gather(*[
            _apply_llm_enhancement(insight, config) for insight in insights
        ])

    return insights


def _parse_csv(path: Path, tumor_type: str | None = None) -> list[ParsedVariant]:
    """Parse variants from a CSV file."""
    df = pd.read_csv(path)
    return _parse_dataframe(df, tumor_type)


def _parse_dataframe(
    df: pd.DataFrame,
    tumor_type: str | None = None,
) -> list[ParsedVariant]:
    """Parse variants from a DataFrame."""
    # Find column names (case-insensitive)
    col_map = {}
    for col in df.columns:
        col_lower = col.lower()
        if col_lower in ("gene", "gene_symbol", "hugo_symbol"):
            col_map["gene"] = col
        elif col_lower in ("variant", "mutation", "alteration", "protein_change"):
            col_map["variant"] = col
        elif col_lower in ("tumor_type", "cancer_type", "disease", "tumor"):
            col_map["tumor_type"] = col

    if "gene" not in col_map or "variant" not in col_map:
        raise ValueError(
            "DataFrame must have 'gene' and 'variant' columns. "
            f"Found columns: {list(df.columns)}"
        )

    parsed_variants = []
    for _, row in df.iterrows():
        row_dict = {
            "gene": row[col_map["gene"]],
            "variant": row[col_map["variant"]],
        }
        if "tumor_type" in col_map:
            row_dict["tumor_type"] = row[col_map["tumor_type"]]

        parsed = parse_variant_row(row_dict)
        if tumor_type and not parsed.tumor_type:
            parsed.tumor_type = tumor_type
        parsed_variants.append(parsed)

    return parsed_variants


def _parse_vcf(path: Path) -> list[ParsedVariant]:
    """Parse variants from a VCF file (basic support)."""
    from oncomind.normalization.input_parser import parse_vcf_variant

    parsed_variants = []

    with open(path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            chrom, pos, _, ref, alt = parts[:5]

            # Parse INFO field if available
            info_dict = {}
            if len(parts) > 7:
                info_str = parts[7]
                for item in info_str.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info_dict[key] = value

            try:
                parsed = parse_vcf_variant(
                    chrom=chrom,
                    pos=int(pos),
                    ref=ref,
                    alt=alt,
                    info=info_dict if info_dict else None,
                )
                parsed_variants.append(parsed)
            except Exception as e:
                print(f"  Warning: Failed to parse VCF line: {str(e)}")

    return parsed_variants


async def _apply_llm_enhancement(
    insight: Insight,
    config: InsightConfig,
) -> Insight:
    """Apply LLM enhancement to an Insight.

    This adds:
    - Paper relevance scoring
    - Literature knowledge extraction
    - Evidence strength assessment
    """
    from oncomind.llm.service import LLMService

    llm_service = LLMService(
        model=config.llm_model,
        temperature=config.llm_temperature,
        enable_logging=True,
    )

    # Score paper relevance and extract knowledge
    if insight.literature.pubmed_articles:
        # Score each paper
        scored_articles = []
        for article in insight.literature.pubmed_articles:
            try:
                relevance = await llm_service.score_paper_relevance(
                    title=article.title,
                    abstract=article.abstract,
                    tldr=article.tldr,
                    gene=insight.identifiers.gene,
                    variant=insight.identifiers.variant,
                    tumor_type=insight.clinical.tumor_type,
                )

                if relevance["is_relevant"]:
                    # Update article with LLM-extracted info
                    article.signal_type = relevance.get("signal_type")
                    article.drugs_mentioned = relevance.get("drugs_mentioned", [])
                    if relevance.get("key_finding"):
                        article.tldr = relevance["key_finding"]
                    scored_articles.append(article)
            except Exception as e:
                print(f"  Warning: Failed to score paper {article.pmid}: {str(e)}")
                scored_articles.append(article)

        insight.literature.pubmed_articles = scored_articles

        # Extract structured knowledge from relevant papers
        if scored_articles:
            paper_contents = [
                {
                    "title": p.title,
                    "abstract": p.abstract,
                    "tldr": p.tldr,
                    "pmid": p.pmid,
                    "url": p.url,
                }
                for p in scored_articles
            ]

            try:
                knowledge_data = await llm_service.extract_variant_knowledge(
                    gene=insight.identifiers.gene,
                    variant=insight.identifiers.variant,
                    tumor_type=insight.clinical.tumor_type,
                    paper_contents=paper_contents,
                )

                from oncomind.models.insight.literature_knowledge import (
                    LiteratureKnowledge, DrugResistance, DrugSensitivity
                )

                resistant_to = []
                for r in knowledge_data.get("resistant_to", []):
                    if isinstance(r, dict):
                        if "is_predictive" not in r:
                            r["is_predictive"] = True
                        resistant_to.append(DrugResistance(**r))
                    else:
                        resistant_to.append(DrugResistance(drug=str(r), is_predictive=True))

                insight.literature.literature_knowledge = LiteratureKnowledge(
                    mutation_type=knowledge_data.get("mutation_type", "unknown"),
                    is_prognostic_only=knowledge_data.get("is_prognostic_only", False),
                    resistant_to=resistant_to,
                    sensitive_to=[
                        DrugSensitivity(**s) if isinstance(s, dict) else DrugSensitivity(drug=str(s))
                        for s in knowledge_data.get("sensitive_to", [])
                    ],
                    clinical_significance=knowledge_data.get("clinical_significance", ""),
                    evidence_level=knowledge_data.get("evidence_level", "None"),
                    references=knowledge_data.get("references", []),
                    key_findings=knowledge_data.get("key_findings", []),
                    confidence=knowledge_data.get("confidence", 0.0),
                )
            except Exception as e:
                print(f"  Warning: Failed to extract literature knowledge: {str(e)}")

    return insight


# Synchronous wrappers for convenience

def get_insight_sync(
    variant_str: str,
    tumor_type: str | None = None,
    config: InsightConfig | None = None,
) -> Insight:
    """Synchronous wrapper for get_insight.

    Use this when you're not in an async context:

        >>> insight = get_insight_sync("BRAF V600E")
    """
    return asyncio.run(get_insight(variant_str, tumor_type, config))


def get_insights_sync(
    variants: list[str] | list[dict] | pd.DataFrame | Path | str,
    tumor_type: str | None = None,
    config: InsightConfig | None = None,
    progress_callback: Callable[[int, int], None] | None = None,
) -> list[Insight]:
    """Synchronous wrapper for get_insights.

    Use this when you're not in an async context:

        >>> insights = get_insights_sync(["BRAF V600E", "EGFR L858R"])
    """
    return asyncio.run(get_insights(variants, tumor_type, config, progress_callback))


__all__ = [
    "get_insight",
    "get_insights",
    "get_insight_sync",
    "get_insights_sync",
    "InsightConfig",
]

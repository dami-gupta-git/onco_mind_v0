"""Public API for variant annotation.

This module provides the main entry points for OncoMind:
- process_variant(): Annotate a single variant
- process_variants(): Batch annotate multiple variants

ARCHITECTURE:
    Input (str/VCF/DataFrame) → parse → build_evidence_panel → EvidencePanel
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

from oncomind.evidence import (
    EvidenceBuilder,
    EvidenceBuilderConfig,
    build_evidence_panel as _build_panel,
    build_evidence_panels as _build_panels,
)
from oncomind.models.evidence.evidence_panel import EvidencePanel
from oncomind.normalization import (
    parse_variant_input,
    parse_variant_row,
    ParsedVariant,
)


@dataclass
class AnnotationConfig:
    """Configuration for the annotation pipeline.

    Controls data sources, LLM integration, and output options.

    Example:
        >>> config = AnnotationConfig(
        ...     enable_llm=True,
        ...     llm_model="gpt-4o-mini",
        ...     enable_literature=True,
        ... )
        >>> panel = await process_variant("BRAF V600E", config=config)
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

    def to_builder_config(self) -> EvidenceBuilderConfig:
        """Convert to EvidenceBuilderConfig."""
        return EvidenceBuilderConfig(
            enable_vicc=self.enable_vicc,
            enable_civic_assertions=self.enable_civic_assertions,
            enable_clinical_trials=self.enable_clinical_trials,
            enable_literature=self.enable_literature,
            max_vicc_results=self.max_vicc_results,
            max_civic_assertions=self.max_civic_assertions,
            max_clinical_trials=self.max_clinical_trials,
            max_literature_results=self.max_literature_results,
        )


async def process_variant(
    variant_str: str,
    tumor_type: str | None = None,
    config: AnnotationConfig | None = None,
) -> EvidencePanel:
    """Process a single variant and return aggregated evidence.

    This is the main entry point for single-variant annotation.

    Args:
        variant_str: Variant string in any supported format:
            - "BRAF V600E" (gene + variant)
            - "EGFR L858R in NSCLC" (with tumor type)
            - "BRAF:V600E" (colon-separated)
            - "TP53 p.R248W" (with p. notation)
        tumor_type: Optional tumor type for clinical context.
            Overrides any tumor type extracted from the string.
        config: Optional configuration for the annotation pipeline.

    Returns:
        EvidencePanel with all aggregated evidence.

    Raises:
        ValueError: If variant string cannot be parsed or variant type
            is not supported (when validate_variant_type=True).

    Example:
        >>> panel = await process_variant("BRAF V600E", tumor_type="Melanoma")
        >>> print(panel.identifiers.gene, panel.identifiers.variant)
        BRAF V600E
        >>> print(panel.clinical.get_approved_drugs())
        ['Dabrafenib', 'Vemurafenib', 'Encorafenib']

        >>> panel = await process_variant("EGFR L858R in lung cancer")
        >>> print(panel.clinical.tumor_type)
        lung cancer
    """
    config = config or AnnotationConfig()

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

    # Build the evidence panel
    builder_config = config.to_builder_config()
    panel = await _build_panel(parsed, parsed.tumor_type, builder_config)

    # Apply LLM enhancement if enabled
    if config.enable_llm:
        panel = await _apply_llm_enhancement(panel, config)

    return panel


async def process_variants(
    variants: list[str] | list[dict] | pd.DataFrame | Path | str,
    tumor_type: str | None = None,
    config: AnnotationConfig | None = None,
    progress_callback: Callable[[int, int], None] | None = None,
) -> list[EvidencePanel]:
    """Process multiple variants and return aggregated evidence.

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
        config: Optional configuration for the annotation pipeline.
        progress_callback: Optional callback(current, total) for progress updates.

    Returns:
        List of EvidencePanel objects.

    Example:
        >>> panels = await process_variants(["BRAF V600E", "EGFR L858R"])
        >>> for panel in panels:
        ...     print(f"{panel.identifiers.gene} {panel.identifiers.variant}")

        >>> panels = await process_variants(
        ...     "variants.csv",
        ...     progress_callback=lambda i, t: print(f"{i}/{t}")
        ... )
    """
    config = config or AnnotationConfig()
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

    # Build evidence panels
    builder_config = config.to_builder_config()

    async with EvidenceBuilder(builder_config) as builder:
        panels = []
        total = len(parsed_variants)

        for i, parsed in enumerate(parsed_variants):
            if progress_callback:
                progress_callback(i, total)

            try:
                panel = await builder.build_evidence_panel(
                    parsed, parsed.tumor_type or tumor_type
                )
                panels.append(panel)
            except Exception as e:
                print(f"  Warning: Failed to process {parsed.gene} {parsed.variant}: {str(e)}")

        if progress_callback:
            progress_callback(total, total)

    # Apply LLM enhancement if enabled
    if config.enable_llm:
        panels = await asyncio.gather(*[
            _apply_llm_enhancement(panel, config) for panel in panels
        ])

    return panels


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
    panel: EvidencePanel,
    config: AnnotationConfig,
) -> EvidencePanel:
    """Apply LLM enhancement to an EvidencePanel.

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
    if panel.literature.pubmed_articles:
        # Score each paper
        scored_articles = []
        for article in panel.literature.pubmed_articles:
            try:
                relevance = await llm_service.score_paper_relevance(
                    title=article.title,
                    abstract=article.abstract,
                    tldr=article.tldr,
                    gene=panel.identifiers.gene,
                    variant=panel.identifiers.variant,
                    tumor_type=panel.clinical.tumor_type,
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

        panel.literature.pubmed_articles = scored_articles

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
                    gene=panel.identifiers.gene,
                    variant=panel.identifiers.variant,
                    tumor_type=panel.clinical.tumor_type,
                    paper_contents=paper_contents,
                )

                from oncomind.models.evidence.literature_knowledge import (
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

                panel.literature.literature_knowledge = LiteratureKnowledge(
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

    return panel


# Synchronous wrappers for convenience

def process_variant_sync(
    variant_str: str,
    tumor_type: str | None = None,
    config: AnnotationConfig | None = None,
) -> EvidencePanel:
    """Synchronous wrapper for process_variant.

    Use this when you're not in an async context:

        >>> panel = process_variant_sync("BRAF V600E")
    """
    return asyncio.run(process_variant(variant_str, tumor_type, config))


def process_variants_sync(
    variants: list[str] | list[dict] | pd.DataFrame | Path | str,
    tumor_type: str | None = None,
    config: AnnotationConfig | None = None,
    progress_callback: Callable[[int, int], None] | None = None,
) -> list[EvidencePanel]:
    """Synchronous wrapper for process_variants.

    Use this when you're not in an async context:

        >>> panels = process_variants_sync(["BRAF V600E", "EGFR L858R"])
    """
    return asyncio.run(process_variants(variants, tumor_type, config, progress_callback))


__all__ = [
    "process_variant",
    "process_variants",
    "process_variant_sync",
    "process_variants_sync",
    "AnnotationConfig",
]

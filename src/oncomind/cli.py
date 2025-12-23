"""Command-line interface for OncoMind.

ARCHITECTURE:
    CLI Commands → get_insight/get_insights → EvidencePanel/VariantInsight

Two main workflows:
- insight: Uses new public API (EvidencePanel output)
- insight-llm: Uses InsightEngine (VariantInsight output with LLM narrative)

Key Design:
- Typer framework for auto-help and type validation
- asyncio.run() bridges sync CLI → async API
- Flexible I/O: stdout or JSON file output
"""

import asyncio
import json
import warnings
from pathlib import Path
from typing import Optional
import typer
from dotenv import load_dotenv


from oncomind import get_insight, AnnotationConfig
from oncomind.engine import InsightEngine
from oncomind.models import VariantInput

# Suppress litellm's async cleanup warnings (harmless internal warnings)
warnings.filterwarnings("ignore", message=".*async_success_handler.*")
warnings.filterwarnings("ignore", message=".*coroutine.*was never awaited.*")

load_dotenv()

app = typer.Typer(
    name="mind",
    help="AI-powered cancer variant insight and evidence synthesis",
    add_completion=False,
)


@app.command()
def insight(
    gene: str = typer.Argument(..., help="Gene symbol (e.g., BRAF)"),
    variant: str = typer.Argument(..., help="Variant notation (e.g., V600E)"),
    tumor: Optional[str] = typer.Option(None, "--tumor", "-t", help="Tumor type"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output JSON file"),
    llm: bool = typer.Option(False, "--llm/--no-llm", help="Enable LLM enhancement"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model (if --llm enabled)"),
    fast: bool = typer.Option(False, "--fast", "-f", help="Fast mode: skip literature search"),
) -> None:
    """Get insight for a single variant and return evidence panel.

    This is the recommended command for variant insight.
    Returns structured EvidencePanel output.

    Examples:
        mind insight BRAF V600E --tumor Melanoma
        mind insight EGFR L858R -t NSCLC
        mind insight TP53 R248W --output result.json
        mind insight KRAS G12D -t NSCLC --fast  # Skip literature for speed
    """
    variant_str = f"{gene} {variant}"

    async def run_insight() -> None:
        print(f"\nGetting insight for {gene} {variant}...")
        if tumor:
            print(f"  Tumor type: {tumor}")
        if fast:
            print("  Mode: fast (skipping literature)")

        config = AnnotationConfig(
            enable_llm=llm,
            llm_model=model,
            enable_literature=not fast,
        )

        panel = await get_insight(variant_str, tumor_type=tumor, config=config)

        # Print summary
        print(f"\n{'='*60}")
        print(f"EVIDENCE PANEL: {panel.identifiers.gene} {panel.identifiers.variant}")
        print(f"{'='*60}")

        stats = panel.get_summary_stats()
        print(f"Tumor Type: {stats['tumor_type'] or 'Not specified'}")
        print(f"Evidence Sources: {', '.join(stats['evidence_sources']) or 'None'}")
        print(f"FDA Approved Drugs: {', '.join(stats['fda_approved_drugs']) or 'None'}")
        print(f"Clinical Trials: {stats['clinical_trials_count']}")
        print(f"PubMed Articles: {stats['pubmed_articles_count']}")

        # Functional scores
        func_summary = panel.functional.get_pathogenicity_summary()
        if func_summary != "No functional predictions available":
            print(f"\nFunctional Predictions:")
            print(f"  {func_summary}")

        # Gene context
        if panel.clinical.gene_role:
            print(f"\nGene Context:")
            print(f"  Role: {panel.clinical.gene_role}")
            if panel.clinical.pathway:
                print(f"  Pathway: {panel.clinical.pathway}")

        print(f"{'='*60}")

        # Save JSON output if requested
        if output:
            output_data = panel.model_dump(mode="json")
            with open(output, "w") as f:
                json.dump(output_data, f, indent=2)
            print(f"\nSaved to {output}")

    asyncio.run(run_insight())


@app.command(name="insight-llm")
def insight_llm(
    gene: str = typer.Argument(..., help="Gene symbol (e.g., BRAF)"),
    variant: str = typer.Argument(..., help="Variant notation (e.g., V600E)"),
    tumor: Optional[str] = typer.Option(None, "--tumor", "-t", help="Tumor type"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model"),
    temperature: float = typer.Option(0.1, "--temperature", help="LLM temperature (0.0-1.0)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output JSON file"),
    log: bool = typer.Option(True, "--log/--no-log", help="Enable LLM decision logging"),
    vicc: bool = typer.Option(True, "--vicc/--no-vicc", help="Enable VICC MetaKB integration"),
) -> None:
    """Get insight for a single variant with LLM-generated narrative.

    This uses the InsightEngine for full LLM narrative generation.
    For faster insight without LLM, use 'mind insight' instead.
    """
    # Import engine only when needed
    from oncomind.engine import InsightEngine
    from oncomind.models import VariantInput

    async def run_llm_insight() -> None:
        variant_input = VariantInput(gene=gene, variant=variant, tumor_type=tumor)

        if tumor:
            print(f"\nGetting insight for {gene} {variant} in {tumor}...")
        else:
            print(f"\nGetting insight for {gene} {variant}...")

        async with InsightEngine(llm_model=model, llm_temperature=temperature, enable_logging=log, enable_vicc=vicc) as engine:
            result = await engine.get_insight(variant_input)

            print(result.get_insight())

            if output:
                output_data = result.model_dump(mode="json")
                with open(output, "w") as f:
                    json.dump(output_data, f, indent=2)
                print(f"Saved to {output}")

    asyncio.run(run_llm_insight())


@app.command()
def batch(
    input_file: Path = typer.Argument(..., help="Input JSON file with variants"),
    output: Path = typer.Option("results.json", "--output", "-o", help="Output file"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model"),
    temperature: float = typer.Option(0.1, "--temperature", help="LLM temperature (0.0-1.0)"),
    fast: bool = typer.Option(False, "--fast", "-f", help="Fast mode: skip literature search"),
    llm: bool = typer.Option(True, "--llm/--no-llm", help="Enable LLM synthesis (use --no-llm for faster results)"),
) -> None:
    """Batch process multiple variants.

    Examples:
        mind batch variants.json --output results.json
        mind batch variants.json --fast --no-llm  # Fastest mode
        mind batch variants.json --fast            # Skip literature, keep LLM
    """

    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        raise typer.Exit(1)

    async def run_batch() -> None:
        from oncomind import get_insights, AnnotationConfig

        with open(input_file, "r") as f:
            data = json.load(f)

        # Build variant strings and tumor types
        variant_strs = []
        tumor_types = []
        for item in data:
            gene = item.get('gene', '')
            variant = item.get('variant', '')
            variant_strs.append(f"{gene} {variant}")
            tumor_types.append(item.get('tumor_type'))

        print(f"\nLoaded {len(variant_strs)} variants from {input_file}")
        if fast:
            print("  Mode: fast (skipping literature)")
        if not llm:
            print("  LLM: disabled")

        config = AnnotationConfig(
            enable_llm=llm,
            llm_model=model,
            llm_temperature=temperature,
            enable_literature=not fast,
        )

        def progress_callback(current: int, total: int) -> None:
            print(f"  Processing {current}/{total}...", end='\r')

        panels = await get_insights(variant_strs, config=config, progress_callback=progress_callback)
        print()  # Clear progress line

        # Apply tumor types and build output
        output_data = []
        for i, panel in enumerate(panels):
            if tumor_types[i]:
                panel.clinical.tumor_type = tumor_types[i]
            output_data.append(panel.model_dump(mode="json"))

        with open(output, "w") as f:
            json.dump(output_data, f, indent=2)

        print(f"\nSuccessfully processed {len(panels)}/{len(variant_strs)} variants")
        print(f"Results saved to {output}")

        # Summary of evidence strength
        strength_counts: dict[str, int] = {}
        for panel in panels:
            strength = panel.meta.evidence_strength or "Unknown"
            strength_counts[strength] = strength_counts.get(strength, 0) + 1

        print("\nEvidence Strength Distribution:")
        for strength, count in sorted(strength_counts.items()):
            print(f"  {strength}: {count}")

    asyncio.run(run_batch())


@app.command()
def version() -> None:
    """Show version information."""
    from oncomind import __version__
    print(f"OncoMind version {__version__}")


if __name__ == "__main__":
    app()

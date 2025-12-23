"""Command-line interface for OncoMind.

ARCHITECTURE:
    CLI Commands → process_variant/process_variants → EvidencePanel/VariantInsight

Two main workflows:
- annotate: Uses new public API (EvidencePanel output)
- process: Uses legacy InsightEngine (VariantInsight output with LLM narrative)

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


from oncomind import process_variant, AnnotationConfig

# Suppress litellm's async cleanup warnings (harmless internal warnings)
warnings.filterwarnings("ignore", message=".*async_success_handler.*")
warnings.filterwarnings("ignore", message=".*coroutine.*was never awaited.*")

load_dotenv()

app = typer.Typer(
    name="mind",
    help="LLM-powered cancer variant annotation and evidence synthesis",
    add_completion=False,
)


@app.command()
def annotate(
    variant_str: str = typer.Argument(..., help="Variant (e.g., 'BRAF V600E' or 'EGFR L858R in NSCLC')"),
    tumor: Optional[str] = typer.Option(None, "--tumor", "-t", help="Tumor type"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output JSON file"),
    llm: bool = typer.Option(False, "--llm/--no-llm", help="Enable LLM enhancement"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model (if --llm enabled)"),
) -> None:
    """Annotate a single variant and return evidence panel.

    This is the recommended command for variant annotation.
    Returns structured EvidencePanel output.

    Examples:
        mind annotate "BRAF V600E" --tumor Melanoma
        mind annotate "EGFR L858R in NSCLC"
        mind annotate "TP53 R248W" --output result.json
    """
    async def run_annotation() -> None:
        print(f"\nAnnotating {variant_str}...")
        if tumor:
            print(f"  Tumor type: {tumor}")

        config = AnnotationConfig(
            enable_llm=llm,
            llm_model=model,
        )

        panel = await process_variant(variant_str, tumor_type=tumor, config=config)

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

    asyncio.run(run_annotation())


@app.command()
def process(
    gene: str = typer.Argument(..., help="Gene symbol (e.g., BRAF)"),
    variant: str = typer.Argument(..., help="Variant notation (e.g., V600E)"),
    tumor: Optional[str] = typer.Option(None, "--tumor", "-t", help="Tumor type"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model"),
    temperature: float = typer.Option(0.1, "--temperature", help="LLM temperature (0.0-1.0)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output JSON file"),
    log: bool = typer.Option(True, "--log/--no-log", help="Enable LLM decision logging"),
    vicc: bool = typer.Option(True, "--vicc/--no-vicc", help="Enable VICC MetaKB integration"),
) -> None:
    """Process a single variant with LLM-generated insight narrative.

    This uses the legacy InsightEngine for full LLM narrative generation.
    For faster annotation without LLM, use 'mind annotate' instead.
    """
    # Import legacy engine only when needed
    from oncomind.engine import InsightEngine
    from oncomind.models import VariantInput

    async def get_variant_insight() -> None:
        variant_input = VariantInput(gene=gene, variant=variant, tumor_type=tumor)

        if tumor:
            print(f"\nProcessing {gene} {variant} in {tumor}...")
        else:
            print(f"\nProcessing {gene} {variant}...")

        async with InsightEngine(llm_model=model, llm_temperature=temperature, enable_logging=log, enable_vicc=vicc) as engine:
            result = await engine.get_insight(variant_input)

            print(result.get_insight())

            if output:
                output_data = result.model_dump(mode="json")
                with open(output, "w") as f:
                    json.dump(output_data, f, indent=2)
                print(f"Saved to {output}")

    asyncio.run(get_variant_insight())


@app.command()
def batch(
    input_file: Path = typer.Argument(..., help="Input JSON file with variants"),
    output: Path = typer.Option("results.json", "--output", "-o", help="Output file"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model"),
    temperature: float = typer.Option(0.1, "--temperature", help="LLM temperature (0.0-1.0)"),
    log: bool = typer.Option(True, "--log/--no-log", help="Enable LLM decision logging"),
) -> None:
    """Batch process multiple variants."""

    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        raise typer.Exit(1)

    async def run_batch() -> None:
        with open(input_file, "r") as f:
            data = json.load(f)

        variants = [VariantInput(**item) for item in data]
        print(f"\nLoaded {len(variants)} variants from {input_file}")

        async with InsightEngine(llm_model=model, llm_temperature=temperature, enable_logging=log) as engine:
            print(f"Processing {len(variants)} variants...")
            insights = await engine.batch_report(variants)

            output_data = [insight.model_dump(mode="json") for insight in insights]
            with open(output, "w") as f:
                json.dump(output_data, f, indent=2)

            print(f"\nSuccessfully processed {len(insights)}/{len(variants)} variants")
            print(f"Results saved to {output}")

            # Summary of evidence strength
            strength_counts: dict[str, int] = {}
            for insight in insights:
                strength = insight.evidence_strength or "Unknown"
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

"""Command-line interface for OncoMind.

ARCHITECTURE:
    CLI Commands → Conductor → Result (evidence + optional LLM narrative)

Main command:
    mind insight GENE VARIANT [--tumor] [--llm]

Modes:
    (default)  Structured evidence only, fast annotation (~7s)
    --llm      + Literature search + LLM synthesis (~25s)
"""

import asyncio
import json
import warnings
from pathlib import Path
from typing import Optional
import typer
from dotenv import load_dotenv


from oncomind import get_insight, InsightConfig

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
    llm: bool = typer.Option(False, "--llm", help="Enable LLM mode: literature search + AI synthesis"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model (only used with --llm)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output JSON file"),
) -> None:
    """Get variant insight with structured evidence and optional LLM narrative.

    By default, fetches evidence from all databases (fast annotation mode).
    Use --llm to enable literature search and AI-powered synthesis.

    Modes:
        (default)  Structured evidence only, fast annotation (~7s)
        --llm      + Literature search + LLM synthesis (~25s)

    Examples:
        mind insight BRAF V600E --tumor Melanoma
        mind insight EGFR L858R -t NSCLC
        mind insight KRAS G12D -t CRC --llm
        mind insight TP53 R248W --llm --output result.json
    """

    from rich.console import Console
    from rich.panel import Panel
    from oncomind.insight_builder import Conductor, ConductorConfig
    import textwrap

    console = Console(width=80)

    async def run_insight() -> None:
        # Generate Result with evidence + optional LLM narrative

        console.print(f"\n[dim]Generating insight ...[/dim]\n", highlight=False)

        # Configure and run the Conductor
        config = ConductorConfig(
            enable_literature=llm,  # Literature search only with --llm
            enable_llm=llm,         # LLM synthesis only with --llm
            llm_model=model,
        )
        async with Conductor(config) as conductor:
            result = await conductor.run(f"{gene} {variant}", tumor_type=tumor)

        # === RENDER OUTPUT ===

        # Variant header with metrics
        variant_title = f"{gene} {variant}"
        if tumor:
            variant_title += f" [dim]in[/dim] {tumor}"

        # Build metrics line
        therapies_count = len(result.llm.recommended_therapies if result.llm else result.evidence.get_recommended_therapies())
        clinvar_sig = result.evidence.clinvar_significance or "N/A"
        am_score = result.functional.alphamissense_score
        am_display = f"{am_score:.2f}" if am_score else "N/A"

        metrics_line = f"[dim]Therapies:[/dim] {therapies_count} [dim]|[/dim] [dim]ClinVar:[/dim] {clinvar_sig} [dim]|[/dim] [dim]AlphaMissense:[/dim] {am_display}"

        from rich.align import Align
        from rich.box import DOUBLE
        from rich.console import Group

        # Center each line individually
        title_line = Align.center(f"[bold bright_white]{variant_title}[/bold bright_white]")
        metrics_centered = Align.center(metrics_line)
        header_content = Group(title_line, metrics_centered)

        console.print(Panel(
            header_content,
            border_style="bold bright_white",
            box=DOUBLE,
            padding=(0, 2),
        ))

        # Summary panel
        summary_text = result.get_summary()
        wrapped_summary = textwrap.fill(summary_text, width=74)
        console.print(Panel(
            f"[cyan]{wrapped_summary}[/cyan]",
            title="[bold]Summary[/bold]",
            border_style="cyan",
            padding=(0, 2),
        ))



        # Build IDs and Scores content
        ids_lines = []

        # Database Identifiers
        if result.identifiers.cosmic_id:
            ids_lines.append(f"[dim]COSMIC:[/dim]       {result.identifiers.cosmic_id}")
        if result.identifiers.dbsnp_id:
            ids_lines.append(f"[dim]dbSNP:[/dim]        {result.identifiers.dbsnp_id}")
        if result.identifiers.clinvar_id:
            ids_lines.append(f"[dim]ClinVar ID:[/dim]   {result.identifiers.clinvar_id}")
        if result.identifiers.ncbi_gene_id:
            ids_lines.append(f"[dim]NCBI Gene:[/dim]    {result.identifiers.ncbi_gene_id}")

        # HGVS notations
        if result.identifiers.hgvs_protein:
            ids_lines.append(f"[dim]HGVS.p:[/dim]       {result.identifiers.hgvs_protein}")
        if result.identifiers.hgvs_transcript:
            ids_lines.append(f"[dim]HGVS.c:[/dim]       {result.identifiers.hgvs_transcript}")
        if result.identifiers.hgvs_genomic:
            ids_lines.append(f"[dim]HGVS.g:[/dim]       {result.identifiers.hgvs_genomic}")

        # Add separator before scores if we have IDs
        if ids_lines:
            ids_lines.append("")

        # Gene role
        if result.context.gene_role:
            ids_lines.append(f"[dim]Gene Role:[/dim]    {result.context.gene_role}")

        # Pathogenicity scores
        func_summary = result.functional.get_pathogenicity_summary()
        if func_summary != "No functional predictions available":
            ids_lines.append(f"[dim]Pathogenicity:[/dim] {func_summary}")

        # gnomAD frequency
        if result.functional.gnomad_exome_af is not None:
            af = result.functional.gnomad_exome_af
            ids_lines.append(f"[dim]gnomAD AF:[/dim]    {af:.2e}" if af > 0 else "[dim]gnomAD AF:[/dim]    0")

        # Show panel only if we have content
        if ids_lines:
            console.print(Panel(
                "\n".join(ids_lines),
                title="[bold]IDs and Scores[/bold]",
                border_style="blue",
                padding=(0, 2),
            ))


        # LLM Insight (when LLM mode is enabled)
        if result.llm:
            # Build formatted LLM insight from raw components
            llm_parts = []
            if result.llm.functional_summary:
                wrapped = textwrap.fill(result.llm.functional_summary, width=70)
                llm_parts.append(f"[bold]Functional Impact:[/bold] {wrapped}")
            if result.llm.biological_context:
                wrapped = textwrap.fill(result.llm.biological_context, width=70)
                llm_parts.append(f"[bold]Biological Context:[/bold] {wrapped}")
            if result.llm.therapeutic_landscape:
                tl = result.llm.therapeutic_landscape
                therapy_parts = []
                if tl.get("fda_approved"):
                    therapy_parts.append(f"FDA-approved: {', '.join(tl['fda_approved'])}")
                if tl.get("clinical_evidence"):
                    therapy_parts.append(f"Clinical evidence: {', '.join(tl['clinical_evidence'])}")
                if tl.get("preclinical"):
                    therapy_parts.append(f"Preclinical: {', '.join(tl['preclinical'])}")
                if tl.get("resistance_mechanisms"):
                    therapy_parts.append(f"Resistance: {', '.join(tl['resistance_mechanisms'])}")
                if therapy_parts:
                    wrapped = textwrap.fill("; ".join(therapy_parts), width=70)
                    llm_parts.append(f"[bold]Therapeutic Landscape:[/bold] {wrapped}")
            if result.llm.research_implications:
                wrapped = textwrap.fill(result.llm.research_implications, width=70)
                llm_parts.append(f"[bold]Research Implications:[/bold] {wrapped}")

            # Fall back to plain summary if no structured components
            llm_content = "\n\n".join(llm_parts) if llm_parts else textwrap.fill(result.llm.llm_summary, width=74)

            console.print(Panel(
                llm_content,
                title="[bold]LLM Insight[/bold]",
                border_style="magenta",
                padding=(1, 2),
            ))

        # Save JSON if requested
        if output:
            # The result object now contains everything, including llm if present
            with open(output, "w") as f:
                json.dump(result.model_dump(mode="json"), f, indent=2)
            console.print(f"\n[dim]Saved to {output}[/dim]")

    asyncio.run(run_insight())


@app.command()
def batch(
    input_file: Path = typer.Argument(..., help="Input JSON file with variants"),
    output: Path = typer.Option("results.json", "--output", "-o", help="Output file"),
    llm: bool = typer.Option(False, "--llm", help="Enable LLM mode: literature search + AI synthesis"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model (only used with --llm)"),
    temperature: float = typer.Option(0.1, "--temperature", help="LLM temperature (0.0-1.0)"),
) -> None:
    """Batch process multiple variants.

    By default, fetches evidence from all databases (fast annotation mode).
    Use --llm to enable literature search and AI-powered synthesis.

    Examples:
        mind batch variants.json --output results.json
        mind batch variants.json --llm               # With literature + LLM
        mind batch variants.json --llm --model gpt-4o
    """
    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        raise typer.Exit(1)

    async def run_batch() -> None:
        from oncomind import get_insights, InsightConfig

        with open(input_file, "r") as f:
            data = json.load(f)

        # Build variant strings and tumor types
        variant_strs = []
        tumor_types = []
        for item in data:
            g = item.get('gene', '')
            v = item.get('variant', '')
            variant_strs.append(f"{g} {v}")
            tumor_types.append(item.get('tumor_type'))

        mode_str = "llm" if llm else "annotation"
        print(f"\nLoaded {len(variant_strs)} variants from {input_file}")
        print(f"  Mode: {mode_str}")

        config = InsightConfig(
            enable_llm=llm,
            llm_model=model,
            llm_temperature=temperature,
            enable_literature=llm,
        )

        def progress_callback(current: int, total: int) -> None:
            print(f"  Processing {current}/{total}...", end='\r')

        results = await get_insights(variant_strs, config=config, progress_callback=progress_callback)
        print()  # Clear progress line

        # Apply tumor types and build output
        output_data = []
        for i, result in enumerate(results):
            if tumor_types[i]:
                result.context.tumor_type = tumor_types[i]
            output_data.append(result.model_dump(mode="json"))

        with open(output, "w") as f:
            json.dump(output_data, f, indent=2)

        print(f"\nSuccessfully processed {len(results)}/{len(variant_strs)} variants")
        print(f"Results saved to {output}")

    asyncio.run(run_batch())


@app.command()
def version() -> None:
    """Show version information."""
    from oncomind import __version__
    print(f"OncoMind version {__version__}")


if __name__ == "__main__":
    app()

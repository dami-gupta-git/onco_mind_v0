"""Command-line interface for OncoMind.

ARCHITECTURE:
    CLI Commands → Conductor → Result (evidence + optional LLM narrative)

Main command:
    mind insight GENE VARIANT [--tumor] [--lit] [--llm] [--full]

Modes:
    (default)  Structured evidence only, fast annotation (~7s)
    --lit      + Literature search (PubMed/Semantic Scholar) (~15s)
    --llm      + LLM synthesis (~20s)
    --full     Both --lit and --llm (~25s)

Logging:
    --log-level  Set log level (DEBUG, INFO, WARN, ERROR). Default: INFO
    Environment: ONCOMIND_LOG_LEVEL=DEBUG|INFO|WARN|ERROR
"""

import asyncio
import json
import warnings
from pathlib import Path
from typing import Optional
import typer
from dotenv import load_dotenv


from oncomind import get_insight, InsightConfig
from oncomind.config.debug import set_log_level, get_logger

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
    lit: bool = typer.Option(False, "--lit", help="Enable literature search (PubMed/Semantic Scholar)"),
    llm: bool = typer.Option(False, "--llm", help="Enable LLM synthesis"),
    full: bool = typer.Option(False, "--full", help="Enable both --lit and --llm"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model (only used with --llm)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output JSON file"),
    log_level: str = typer.Option("INFO", "--log-level", "-l", help="Log level: DEBUG, INFO, WARN, ERROR"),
) -> None:
    """Get variant insight with structured evidence and optional LLM narrative.

    By default, fetches evidence from all databases (fast annotation mode).
    Use --lit for literature search, --llm for AI synthesis, --full for both.

    Modes:
        (default)  Structured evidence only, fast annotation (~7s)
        --lit      + Literature search (PubMed/Semantic Scholar) (~15s)
        --llm      + LLM synthesis (~20s)
        --full     Both --lit and --llm (~25s)

    Logging:
        --log-level  DEBUG, INFO, WARN, ERROR (default: INFO)
        Or set ONCOMIND_LOG_LEVEL environment variable

    Examples:
        mind insight BRAF V600E --tumor Melanoma
        mind insight EGFR L858R -t NSCLC --lit
        mind insight KRAS G12D -t CRC --llm --log-level DEBUG
        mind insight TP53 R248W --full --output result.json
    """
    # Initialize logging
    set_log_level(log_level)
    logger = get_logger(__name__)

    from rich.console import Console
    from rich.panel import Panel
    from oncomind.insight_builder import Conductor, ConductorConfig
    import textwrap

    console = Console(width=80)

    async def run_insight() -> None:
        # Generate Result with evidence + optional LLM narrative

        console.print(f"\n[dim]Generating insight ...[/dim]\n", highlight=False)

        # Handle flag combinations:
        # --full implies both --lit and --llm
        # --lit and --llm are independent
        enable_lit = lit or full
        enable_llm_flag = llm or full

        logger.debug(f"Processing {gene} {variant} (tumor={tumor})")
        logger.debug(f"Flags: lit={enable_lit}, llm={enable_llm_flag}, model={model}")

        # Configure and run the Conductor
        config = ConductorConfig(
            enable_literature=enable_lit,
            enable_llm=enable_llm_flag,
            llm_model=model,
        )
        logger.debug(f"ConductorConfig: {config}")

        async with Conductor(config) as conductor:
            result = await conductor.run(f"{gene} {variant}", tumor_type=tumor)

        logger.debug(f"Result received: evidence_sources={len(result.evidence.fda_approvals)} FDA, "
                    f"{len(result.evidence.civic_assertions)} CIViC, llm={'yes' if result.llm else 'no'}")

        # LLM Insight (debug logging only)
        logger.debug("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        if result.llm:
            logger.debug(f"LLM Insight functional_summary: {result.llm.functional_summary}")
            logger.debug(f"LLM Insight biological_context: {result.llm.biological_context}")
            logger.debug(f"LLM Insight therapeutic_landscape: {result.llm.therapeutic_landscape}")
            logger.debug(f"LLM Insight research_implications: {result.llm.research_implications}")
            logger.debug(f"LLM Insight llm_summary: {result.llm.llm_summary}")
        logger.debug("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

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

        # Check for LLM error and display to user
        if result.llm and result.llm.rationale and result.llm.rationale.startswith("LLM narrative generation failed:"):
            error_msg = result.llm.rationale.replace("LLM narrative generation failed: ", "")
            console.print(Panel(
                f"[yellow]LLM synthesis failed: {error_msg}[/yellow]\n\n[dim]Showing evidence-only results below.[/dim]",
                title="[bold yellow]⚠ LLM Error[/bold yellow]",
                border_style="yellow",
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

        # Gap Analysis panel (shown in both LLM and annotation modes)
        evidence_gaps = result.evidence.evidence_gaps
        if evidence_gaps:
            gap_lines = []

            # Overall quality and research priority
            quality = evidence_gaps.overall_evidence_quality or "unknown"
            priority = evidence_gaps.research_priority or "unknown"
            quality_colors = {"comprehensive": "green", "moderate": "yellow", "limited": "red", "minimal": "red"}
            priority_colors = {"very_high": "red", "high": "red", "medium": "yellow", "low": "green"}
            q_color = quality_colors.get(quality.lower(), "white")
            p_color = priority_colors.get(priority.lower(), "white")

            gap_lines.append(f"[dim]Evidence Quality:[/dim]   [{q_color}]{quality.capitalize()}[/{q_color}]")
            gap_lines.append(f"[dim]Research Priority:[/dim]  [{p_color}]{priority.replace('_', ' ').title()}[/{p_color}]")

            # Well characterized aspects (compact)
            if evidence_gaps.well_characterized:
                gap_lines.append("")
                gap_lines.append("[dim]✅ Well Characterized:[/dim]")
                for aspect in evidence_gaps.well_characterized[:5]:
                    gap_lines.append(f"   • {aspect}")

            # Evidence gaps by severity
            from oncomind.models.evidence.evidence_gaps import GapSeverity
            critical_gaps = evidence_gaps.get_gaps_by_severity(GapSeverity.CRITICAL)
            significant_gaps = evidence_gaps.get_gaps_by_severity(GapSeverity.SIGNIFICANT)

            if critical_gaps:
                gap_lines.append("")
                gap_lines.append(f"[red]🔴 Critical Gaps ({len(critical_gaps)}):[/red]")
                for g in critical_gaps[:3]:
                    gap_lines.append(f"   • {g.description}")

            if significant_gaps:
                gap_lines.append("")
                gap_lines.append(f"[yellow]🟠 Significant Gaps ({len(significant_gaps)}):[/yellow]")
                for g in significant_gaps[:3]:
                    gap_lines.append(f"   • {g.description}")

            if gap_lines:
                console.print(Panel(
                    "\n".join(gap_lines),
                    title="[bold]Gap Analysis[/bold]",
                    border_style="magenta",
                    padding=(0, 2),
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
    lit: bool = typer.Option(False, "--lit", help="Enable literature search (PubMed/Semantic Scholar)"),
    llm: bool = typer.Option(False, "--llm", help="Enable LLM synthesis"),
    full: bool = typer.Option(False, "--full", help="Enable both --lit and --llm"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model (only used with --llm)"),
    temperature: float = typer.Option(0.1, "--temperature", help="LLM temperature (0.0-1.0)"),
    log_level: str = typer.Option("INFO", "--log-level", "-l", help="Log level: DEBUG, INFO, WARN, ERROR"),
) -> None:
    """Batch process multiple variants.

    By default, fetches evidence from all databases (fast annotation mode).
    Use --lit for literature, --llm for AI synthesis, --full for both.

    Examples:
        mind batch variants.json --output results.json
        mind batch variants.json --lit               # With literature search
        mind batch variants.json --llm               # With literature + LLM
        mind batch variants.json --full --model gpt-4o
        mind batch variants.json --log-level DEBUG   # Verbose logging
    """
    # Initialize logging
    set_log_level(log_level)
    logger = get_logger(__name__)

    logger.debug(f"Batch processing: input={input_file}, output={output}")

    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
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

        # Handle flag combinations:
        # --full implies both --lit and --llm
        # --lit and --llm are independent
        enable_lit = lit or full
        enable_llm_mode = llm or full

        mode_str = "llm" if enable_llm_mode else ("literature" if enable_lit else "annotation")
        logger.info(f"Loaded {len(variant_strs)} variants from {input_file}")
        logger.info(f"Mode: {mode_str}")
        print(f"\nLoaded {len(variant_strs)} variants from {input_file}")
        print(f"  Mode: {mode_str}")

        config = InsightConfig(
            enable_llm=enable_llm_mode,
            llm_model=model,
            llm_temperature=temperature,
            enable_literature=enable_lit,
        )

        def progress_callback(current: int, total: int) -> None:
            logger.info(f"Processing {current}/{total}...")
            print(f"  Processing {current}/{total}...", end='\r')

        results = await get_insights(variant_strs, config=config, progress_callback=progress_callback)
        print()  # Clear progress line

        # Apply tumor types and build output, tracking LLM errors
        output_data = []
        llm_errors = []
        for i, result in enumerate(results):
            if tumor_types[i]:
                result.context.tumor_type = tumor_types[i]
            output_data.append(result.model_dump(mode="json"))
            # Track LLM errors
            if result.llm and result.llm.rationale and result.llm.rationale.startswith("LLM narrative generation failed:"):
                error_msg = result.llm.rationale.replace("LLM narrative generation failed: ", "")
                llm_errors.append(f"{variant_strs[i]}: {error_msg}")

        with open(output, "w") as f:
            json.dump(output_data, f, indent=2)

        logger.info(f"Successfully processed {len(results)}/{len(variant_strs)} variants")
        print(f"\nSuccessfully processed {len(results)}/{len(variant_strs)} variants")

        # Show LLM errors if any occurred
        if llm_errors:
            logger.warning(f"LLM errors ({len(llm_errors)}):")
            print(f"\n⚠ LLM errors ({len(llm_errors)}):")
            for err in llm_errors:
                logger.warning(f"  - {err}")
                print(f"  - {err}")

        logger.info(f"Results saved to {output}")
        print(f"Results saved to {output}")

    asyncio.run(run_batch())


@app.command()
def version() -> None:
    """Show version information."""
    from oncomind import __version__
    print(f"OncoMind version {__version__}")


if __name__ == "__main__":
    app()

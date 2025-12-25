"""Command-line interface for OncoMind.

ARCHITECTURE:
    CLI Commands → InsightBuilder → Insight
                 → LLMService → LLMInsight (optional, embedded in Insight.llm)

Main command:
    mind insight GENE VARIANT [--tumor] [--lite|--full]

Modes:
    (default)  Structured evidence + LLM narrative (~12s)
    --lite     Structured evidence only, no LLM (~7s)
    --full     + Literature search + enhanced narrative (~25s)
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
    lite: bool = typer.Option(False, "--lite", help="Lite mode: structured evidence only, no LLM"),
    full: bool = typer.Option(False, "--full", help="Full mode: include literature search + enhanced narrative"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output JSON file"),
) -> None:
    """Get variant insight with structured evidence and optional LLM narrative.

    By default, fetches evidence from all databases and generates an LLM summary.

    Modes:
        (default)  Structured evidence + LLM narrative (~12s)
        --lite     Structured evidence only, no LLM (~7s)
        --full     + Literature search + enhanced narrative (~25s)

    Examples:
        mind insight BRAF V600E --tumor Melanoma
        mind insight EGFR L858R -t NSCLC --lite
        mind insight KRAS G12D -t CRC --full
        mind insight TP53 R248W --output result.json
    """
    if lite and full:
        print("Error: Cannot use both --lite and --full")
        raise typer.Exit(1)

    from rich.console import Console
    from rich.panel import Panel
    from oncomind.insight_builder import InsightBuilder, InsightBuilderConfig
    from oncomind.llm.service import LLMService
    import textwrap

    console = Console(width=80)

    async def run_insight() -> None:
        # Generate Insight with optional LLM narrative embedded

        # Determine mode
        if lite:
            mode_str = "lite"
        elif full:
            mode_str = "full"
        else:
            mode_str = "default"

        console.print(f"\n[dim]Generating insight ...[/dim]\n", highlight=False)

        # Step 1: Fetch evidence
        enable_literature = full  # Only fetch literature in full mode
        config = InsightBuilderConfig(enable_literature=enable_literature)
        async with InsightBuilder(config) as builder:
            insight = await builder.build_insight(f"{gene} {variant}", tumor_type=tumor)

        # Step 2: Generate LLM narrative (unless --lite)
        if not lite:
            # Get evidence summary for LLM
            evidence_summary = insight.get_evidence_summary_for_llm()

            # Get LLM insight
            llm_service = LLMService(model=model, temperature=0.1)
            llm_insight = await llm_service.get_llm_insight(
                gene=gene,
                variant=variant,
                tumor_type=tumor,
                evidence_summary=evidence_summary,
                has_clinical_trials=bool(insight.clinical.clinical_trials),
            )

            # Embed LLM insight in the main insight object
            insight.llm = llm_insight

        # === RENDER OUTPUT ===

        # Variant header with metrics
        variant_title = f"{gene} {variant}"
        if tumor:
            variant_title += f" [dim]in[/dim] {tumor}"

        # Build metrics line
        therapies_count = len(insight.llm.recommended_therapies) if insight.llm else len(insight.clinical.fda_approvals)
        clinvar_sig = insight.clinical.clinvar_clinical_significance or "N/A"
        am_score = insight.functional.alphamissense_score
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
        summary_text = insight.get_summary()
        wrapped_summary = textwrap.fill(summary_text, width=74)
        console.print(Panel(
            f"[cyan]{wrapped_summary}[/cyan]",
            title="[bold]Summary[/bold]",
            border_style="cyan",
            padding=(0, 2),
        ))

        # LLM Insight (when LLM mode is enabled)
        if insight.llm:
            wrapped_llm = textwrap.fill(insight.llm.llm_summary, width=74)
            console.print(Panel(
                wrapped_llm,
                title="[bold]LLM Insight[/bold]",
                border_style="magenta",
                padding=(1, 2),
            ))

        # Build header content with all variant info
        header_lines = []

        # IDs row
        ids = []
        if insight.identifiers.cosmic_id:
            ids.append(f"COSMIC:{insight.identifiers.cosmic_id}")
        if insight.identifiers.dbsnp_id:
            ids.append(f"dbSNP:{insight.identifiers.dbsnp_id}")
        if insight.identifiers.clinvar_id:
            ids.append(f"ClinVar:{insight.identifiers.clinvar_id}")
        if ids:
            header_lines.append(f"[dim]{' | '.join(ids)}[/dim]")

        # ClinVar significance
        if insight.clinical.clinvar_clinical_significance:
            header_lines.append(f"[dim]ClinVar:[/dim]            {insight.clinical.clinvar_clinical_significance}")

        # Functional predictions
        func_summary = insight.functional.get_pathogenicity_summary()
        if func_summary != "No functional predictions available":
            header_lines.append(f"[dim]Pathogenicity:[/dim]      {func_summary}")

        # Gene role
        if insight.clinical.gene_role:
            header_lines.append(f"[dim]Gene Role:[/dim]          {insight.clinical.gene_role}")

        # Evidence sources
        sources = insight.kb.get_evidence_sources()
        if sources:
            header_lines.append(f"[dim]Evidence Sources:[/dim]   {', '.join(sources)}")

        console.print(Panel(
            "\n".join(header_lines),
            title="[bold]Evidence Overview[/bold]",
            border_style="blue",
            padding=(0, 2),
        ))

        # FDA Approved Drugs (always show if available)
        if insight.clinical.fda_approvals:
            drugs = insight.clinical.get_approved_drugs()
            if drugs:
                console.print(Panel(
                    "[bold green]" + ", ".join(drugs) + "[/bold green]",
                    title="[bold]FDA Approved Drugs[/bold]",
                    border_style="green",
                    padding=(0, 2),
                ))

        # Recommended Therapies (from LLM insight or extracted from insight)
        therapy_lines = []
        if insight.llm and insight.llm.recommended_therapies:
            # Use LLM-recommended therapies
            for t in insight.llm.recommended_therapies:
                level = f"Level {t.evidence_level}" if t.evidence_level else ""
                status = t.approval_status or ""
                parts = [p for p in [level, status] if p]
                suffix = f" ({', '.join(parts)})" if parts else ""
                therapy_lines.append(f"  • {t.drug_name}{suffix}")
        else:
            # Extract from insight (lite mode)
            seen_drugs = set()
            # From FDA approvals
            for approval in insight.clinical.fda_approvals:
                drug_name = approval.brand_name or approval.generic_name or approval.drug_name
                if drug_name and drug_name not in seen_drugs:
                    seen_drugs.add(drug_name)
                    therapy_lines.append(f"  • {drug_name} (Level A, FDA Approved)")
            # From CGI biomarkers
            for biomarker in insight.kb.cgi_biomarkers:
                if biomarker.fda_approved and biomarker.drug and biomarker.drug not in seen_drugs:
                    seen_drugs.add(biomarker.drug)
                    level = biomarker.evidence_level or "A"
                    therapy_lines.append(f"  • {biomarker.drug} (Level {level}, FDA Approved)")

        if therapy_lines:
            console.print(Panel(
                "\n".join(therapy_lines),
                title="[bold]Recommended Therapies[/bold]",
                border_style="cyan",
                padding=(0, 2),
            ))

        # Clinical evidence details (CIViC/CGI/Trials in one box)
        has_clinical_evidence = insight.kb.civic_assertions or insight.kb.cgi_biomarkers or insight.clinical.clinical_trials
        if has_clinical_evidence:
            evidence_lines = []

            if insight.kb.civic_assertions:
                evidence_lines.append("[bold]CIViC Assertions:[/bold]")
                for a in insight.kb.civic_assertions[:3]:
                    therapies = ", ".join(a.therapies) if a.therapies else "N/A"
                    evidence_lines.append(f"  • {therapies} → {a.significance}")

            if insight.kb.cgi_biomarkers:
                if evidence_lines:
                    evidence_lines.append("")
                evidence_lines.append("[bold]CGI Biomarkers:[/bold]")
                for b in insight.kb.cgi_biomarkers[:3]:
                    fda_tag = " [green](FDA)[/green]" if b.fda_approved else ""
                    evidence_lines.append(f"  • {b.drug}: {b.association}{fda_tag}")

            if insight.clinical.clinical_trials:
                if evidence_lines:
                    evidence_lines.append("")
                trials_count = len(insight.clinical.clinical_trials)
                variant_specific = sum(1 for t in insight.clinical.clinical_trials if t.variant_specific)
                evidence_lines.append("[bold]Clinical Trials:[/bold]")
                trials_text = f"  {trials_count} recruiting"
                if variant_specific:
                    trials_text += f" ({variant_specific} variant-specific)"
                evidence_lines.append(trials_text)

            console.print(Panel(
                "\n".join(evidence_lines),
                title="[bold]Knowledge Base Evidence[/bold]",
                border_style="cyan",
                padding=(0, 2),
            ))

        # Literature section (only in full mode)
        if full and insight.literature.pubmed_articles:
            lit_lines = []
            lit_lines.append(f"[bold]PubMed Articles ({len(insight.literature.pubmed_articles)}):[/bold]")
            for article in insight.literature.pubmed_articles[:3]:
                title = article.title[:55] + "..." if len(article.title) > 55 else article.title
                lit_lines.append(f"  • PMID:{article.pmid} {title}")
            if len(insight.literature.pubmed_articles) > 3:
                lit_lines.append("  ...")

            if insight.literature.literature_knowledge:
                lk = insight.literature.literature_knowledge
                if lk.resistance_mechanisms:
                    lit_lines.append("")
                    lit_lines.append("[bold]Resistance Mechanisms:[/bold]")
                    for r in lk.resistance_mechanisms[:2]:
                        lit_lines.append(f"  • {r.drug}: {r.mechanism}")
                if lk.sensitivity_markers:
                    lit_lines.append("")
                    lit_lines.append("[bold]Sensitivity Markers:[/bold]")
                    for s in lk.sensitivity_markers[:2]:
                        lit_lines.append(f"  • {s.drug}: {s.marker}")

            console.print(Panel(
                "\n".join(lit_lines),
                title="[bold]Literature[/bold]",
                border_style="yellow",
                padding=(0, 2),
            ))

        # Save JSON if requested
        if output:
            # The insight object now contains everything, including llm if present
            with open(output, "w") as f:
                json.dump(insight.model_dump(mode="json"), f, indent=2)
            console.print(f"\n[dim]Saved to {output}[/dim]")

    asyncio.run(run_insight())


@app.command()
def batch(
    input_file: Path = typer.Argument(..., help="Input JSON file with variants"),
    output: Path = typer.Option("results.json", "--output", "-o", help="Output file"),
    model: str = typer.Option("gpt-4o-mini", "--model", "-m", help="LLM model"),
    temperature: float = typer.Option(0.1, "--temperature", help="LLM temperature (0.0-1.0)"),
    lite: bool = typer.Option(False, "--lite", help="Lite mode: skip LLM synthesis"),
    full: bool = typer.Option(False, "--full", help="Full mode: include literature search"),
) -> None:
    """Batch process multiple variants.

    Examples:
        mind batch variants.json --output results.json
        mind batch variants.json --lite              # Fastest: no LLM
        mind batch variants.json --full              # Slowest: with literature
    """
    if lite and full:
        print("Error: Cannot use both --lite and --full")
        raise typer.Exit(1)

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

        mode_str = "lite" if lite else ("full" if full else "default")
        print(f"\nLoaded {len(variant_strs)} variants from {input_file}")
        print(f"  Mode: {mode_str}")

        config = InsightConfig(
            enable_llm=not lite,
            llm_model=model,
            llm_temperature=temperature,
            enable_literature=full,
        )

        def progress_callback(current: int, total: int) -> None:
            print(f"  Processing {current}/{total}...", end='\r')

        insights = await get_insights(variant_strs, config=config, progress_callback=progress_callback)
        print()  # Clear progress line

        # Apply tumor types and build output
        output_data = []
        for i, insight in enumerate(insights):
            if tumor_types[i]:
                insight.clinical.tumor_type = tumor_types[i]
            output_data.append(insight.model_dump(mode="json"))

        with open(output, "w") as f:
            json.dump(output_data, f, indent=2)

        print(f"\nSuccessfully processed {len(insights)}/{len(variant_strs)} variants")
        print(f"Results saved to {output}")

    asyncio.run(run_batch())


@app.command()
def version() -> None:
    """Show version information."""
    from oncomind import __version__
    print(f"OncoMind version {__version__}")


if __name__ == "__main__":
    app()

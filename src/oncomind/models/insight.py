"""Variant insight and annotation models."""

import textwrap

from pydantic import BaseModel, Field
from rich.console import Console
from rich.panel import Panel

from oncomind.models.annotations import VariantAnnotations


class RecommendedTherapy(BaseModel):
    """Recommended therapy based on variant."""

    drug_name: str = Field(..., description="Name of the therapeutic agent")
    evidence_level: str | None = Field(None, description="Level of supporting evidence")
    approval_status: str | None = Field(None, description="FDA approval status for this indication")
    clinical_context: str | None = Field(
        None, description="Clinical context (e.g., first-line, resistant)"
    )


class LLMInsight(VariantAnnotations):
    """LLM-generated narrative insight for a variant.

    This contains the LLM's synthesized interpretation of the evidence,
    including clinical summary, recommended therapies, and rationale.
    """

    gene: str
    variant: str
    tumor_type: str | None
    llm_summary: str = Field(..., description="LLM-generated narrative summary of the variant")
    recommended_therapies: list[RecommendedTherapy] = Field(default_factory=list)
    rationale: str = Field(..., description="Detailed rationale for clinical interpretation")
    evidence_strength: str | None = Field(
        None, description="Overall strength of evidence (Strong/Moderate/Weak)"
    )
    clinical_trials_available: bool = Field(
        default=False, description="Whether relevant clinical trials exist"
    )
    references: list[str] = Field(
        default_factory=list, description="Key references supporting the insight"
    )

    def get_insight(self) -> str:
        """Pretty report output with Rich formatting and soft-wrapping."""
        console = Console(width=80, force_terminal=True)

        # Build content sections
        tumor_display = self.tumor_type if self.tumor_type else "Not specified"

        # Header line
        header = f"[bold cyan]{self.gene} {self.variant}[/bold cyan]  |  Tumor: [italic]{tumor_display}[/italic]"

        content_lines = [header, ""]

        # Add identifiers if available
        identifiers = []
        if self.cosmic_id:
            identifiers.append(f"COSMIC: {self.cosmic_id}")
        if self.ncbi_gene_id:
            identifiers.append(f"NCBI: {self.ncbi_gene_id}")
        if self.dbsnp_id:
            identifiers.append(f"dbSNP: {self.dbsnp_id}")
        if self.clinvar_id:
            identifiers.append(f"ClinVar: {self.clinvar_id}")
        if identifiers:
            content_lines.append(f"[dim]IDs:[/dim] {' | '.join(identifiers)}")

        # Add HGVS notations if available (compact)
        if self.hgvs_protein:
            content_lines.append(f"[dim]HGVS:[/dim] {self.hgvs_protein}")

        # Add ClinVar significance if available
        if self.clinvar_clinical_significance:
            content_lines.append(f"[dim]ClinVar:[/dim] {self.clinvar_clinical_significance}")

        # Add key functional annotations if available
        annotations = []
        if self.alphamissense_prediction:
            am_display = {"P": "Pathogenic", "B": "Benign", "A": "Ambiguous"}.get(
                self.alphamissense_prediction, self.alphamissense_prediction
            )
            score_str = f" ({self.alphamissense_score:.2f})" if self.alphamissense_score else ""
            annotations.append(f"AlphaMissense: {am_display}{score_str}")
        if self.cadd_score is not None:
            annotations.append(f"CADD: {self.cadd_score:.2f}")
        if annotations:
            content_lines.append(f"[dim]Scores:[/dim] {' | '.join(annotations)}")

        # Evidence strength if available
        if self.evidence_strength:
            content_lines.append(f"[dim]Evidence:[/dim] {self.evidence_strength}")

        # Clinical narrative - soft-wrapped
        content_lines.append("")
        wrapped_summary = textwrap.fill(self.llm_summary, width=74)
        content_lines.append(wrapped_summary)

        # Therapies section
        if self.recommended_therapies:
            content_lines.append("")
            therapy_names = ", ".join([t.drug_name for t in self.recommended_therapies])
            wrapped_therapies = textwrap.fill(f"Therapies: {therapy_names}", width=74)
            content_lines.append(f"[bold green]{wrapped_therapies}[/bold green]")

        # Join all content
        content = "\n".join(content_lines)

        # Create panel with box styling
        panel = Panel(
            content,
            title="[bold white]Variant Insight[/bold white]",
            border_style="blue",
            padding=(1, 2),
        )

        # Render to string
        with console.capture() as capture:
            console.print(panel)

        return capture.get()

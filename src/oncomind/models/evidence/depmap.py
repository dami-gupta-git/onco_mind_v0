"""DepMap/CCLE evidence model for cell line dependency and drug sensitivity data."""

from pydantic import BaseModel, Field

from oncomind.models.evidence.base import EvidenceItemBase


class DrugSensitivity(BaseModel):
    """Drug sensitivity data from PRISM screens.

    Based on mean log2 fold change across cell lines with the mutation:
    - log2fc <= -1.7: sensitive (drug kills cells)
    - log2fc > -1.7: resistant (drug doesn't kill cells)
    """

    drug_name: str = Field(..., description="Drug name")
    mean_log2fc: float | None = Field(None, description="Mean log2 fold change across cell lines (more negative = more sensitive)")
    n_cell_lines: int = Field(0, description="Number of cell lines tested")
    sensitive_lines: list[str] = Field(default_factory=list, description="Names of sensitive cell lines")


class GeneDependency(BaseModel):
    """CRISPR gene dependency data."""

    gene: str = Field(..., description="Gene symbol")
    mean_dependency_score: float | None = Field(
        None, description="Mean CERES score (negative = essential, <-0.5 typically essential)"
    )
    n_dependent_lines: int = Field(0, description="Number of cell lines dependent on this gene")
    n_total_lines: int = Field(0, description="Total cell lines tested")
    dependency_pct: float = Field(0.0, description="Percent of lines showing dependency")
    top_dependent_lines: list[str] = Field(
        default_factory=list, description="Cell lines most dependent on this gene"
    )


class CellLineModel(BaseModel):
    """A cell line that can be used as a model system."""

    name: str = Field(..., description="Cell line name (e.g., A375)")
    ccle_name: str | None = Field(None, description="CCLE formatted name")
    primary_disease: str | None = Field(None, description="Primary disease/cancer type")
    subtype: str | None = Field(None, description="Cancer subtype")
    has_mutation: bool = Field(False, description="Whether this line has the query mutation")
    mutation_details: str | None = Field(None, description="Specific mutation if known")


class DepMapEvidence(EvidenceItemBase):
    """Evidence from DepMap/CCLE showing gene dependencies and drug sensitivities.

    DepMap provides:
    - CRISPR gene dependency scores (is this gene essential?)
    - PRISM drug sensitivity (which drugs kill cells with this mutation?)
    - Cell line models (which cell lines have this mutation?)

    This data is critical for:
    - Identifying synthetic lethality opportunities
    - Finding drugs with selective activity in mutant cells
    - Selecting model systems for functional studies
    """

    gene: str = Field(..., description="Query gene")
    variant: str | None = Field(None, description="Query variant if applicable")

    # Gene dependency (is this gene essential?)
    gene_dependency: GeneDependency | None = Field(
        None, description="CRISPR dependency data for the query gene"
    )

    # Co-dependencies (what other genes are essential when this gene is mutated?)
    co_dependencies: list[GeneDependency] = Field(
        default_factory=list,
        description="Genes that are co-essential with the query gene (synthetic lethality candidates)"
    )

    # Drug sensitivities in mutant cell lines
    drug_sensitivities: list[DrugSensitivity] = Field(
        default_factory=list,
        description="Drugs showing differential sensitivity in cell lines with this mutation"
    )

    # Available cell line models
    cell_line_models: list[CellLineModel] = Field(
        default_factory=list,
        description="Cell lines harboring this gene mutation"
    )

    # Metadata
    data_version: str | None = Field(None, description="DepMap data release version")
    n_cell_lines_screened: int = Field(0, description="Total cell lines in dataset")

    def has_data(self) -> bool:
        """Check if there is meaningful data."""
        return bool(
            self.gene_dependency or
            self.co_dependencies or
            self.drug_sensitivities or
            self.cell_line_models
        )

    def get_essential_score(self) -> float | None:
        """Get the gene essentiality score (negative = more essential)."""
        if self.gene_dependency:
            return self.gene_dependency.mean_dependency_score
        return None

    def is_essential(self, threshold: float = -0.5) -> bool:
        """Check if gene is essential (CERES score < threshold)."""
        score = self.get_essential_score()
        return score is not None and score < threshold

    def get_top_sensitive_drugs(self, n: int = 5) -> list[DrugSensitivity]:
        """Get top drugs by sensitivity (most negative log2fc = most effective)."""
        def sort_key(d: DrugSensitivity) -> float:
            if d.mean_log2fc is not None:
                return d.mean_log2fc  # More negative = more sensitive
            return float('inf')

        return sorted(self.drug_sensitivities, key=sort_key)[:n]

    def get_model_cell_lines(self, with_mutation_only: bool = True) -> list[str]:
        """Get cell line names for experimental use."""
        if with_mutation_only:
            return [cl.name for cl in self.cell_line_models if cl.has_mutation]
        return [cl.name for cl in self.cell_line_models]

    def to_prompt_context(self) -> str:
        """Format DepMap evidence for LLM prompt.

        Returns:
            Formatted string with dependency and drug sensitivity context
        """
        if not self.has_data():
            return "No DepMap/CCLE data available for this gene."

        lines = []
        source_cite = "[DepMap](https://depmap.org/)"

        # Gene dependency
        if self.gene_dependency:
            gd = self.gene_dependency
            score = gd.mean_dependency_score
            if score is not None:
                essential_str = "ESSENTIAL" if score < -0.5 else "not essential"
                lines.append(f"GENE DEPENDENCY ({source_cite}):")
                lines.append(f"  {self.gene} is {essential_str} (CERES score: {score:.2f})")
                lines.append(f"  Dependent in {gd.n_dependent_lines}/{gd.n_total_lines} cell lines ({gd.dependency_pct:.1f}%)")
                if gd.top_dependent_lines:
                    lines.append(f"  Most dependent lines: {', '.join(gd.top_dependent_lines[:5])}")
                lines.append("")

        # Co-dependencies (synthetic lethality candidates)
        if self.co_dependencies:
            lines.append(f"CO-DEPENDENCIES (synthetic lethality candidates, {source_cite}):")
            for cd in self.co_dependencies[:5]:
                if cd.mean_dependency_score is not None:
                    lines.append(f"  - {cd.gene}: CERES={cd.mean_dependency_score:.2f}, "
                               f"{cd.n_dependent_lines} dependent lines")
            lines.append("")

        # Drug sensitivities
        if self.drug_sensitivities:
            lines.append(f"DRUG SENSITIVITIES in {self.gene}-mutant lines ({source_cite}):")
            for ds in self.get_top_sensitive_drugs(5):
                parts = [f"  - {ds.drug_name}:"]
                if ds.mean_log2fc is not None:
                    parts.append(f"log2FC={ds.mean_log2fc:.2f}")
                parts.append(f"(n={ds.n_cell_lines})")
                lines.append(" ".join(parts))
            lines.append("")

        # Available cell line models
        if self.cell_line_models:
            mutant_lines = [cl for cl in self.cell_line_models if cl.has_mutation]
            if mutant_lines:
                lines.append(f"AVAILABLE MODEL CELL LINES with {self.gene} mutation ({source_cite}):")
                for cl in mutant_lines[:8]:
                    disease_str = f" ({cl.primary_disease})" if cl.primary_disease else ""
                    mutation_str = f" [{cl.mutation_details}]" if cl.mutation_details else ""
                    lines.append(f"  - {cl.name}{disease_str}{mutation_str}")
                if len(mutant_lines) > 8:
                    lines.append(f"  ... and {len(mutant_lines) - 8} more")
                lines.append("")

        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "gene": self.gene,
            "variant": self.variant,
            "gene_dependency": self.gene_dependency.model_dump() if self.gene_dependency else None,
            "co_dependencies": [cd.model_dump() for cd in self.co_dependencies],
            "drug_sensitivities": [ds.model_dump() for ds in self.drug_sensitivities],
            "cell_line_models": [cl.model_dump() for cl in self.cell_line_models],
            "data_version": self.data_version,
            "n_cell_lines_screened": self.n_cell_lines_screened,
        }

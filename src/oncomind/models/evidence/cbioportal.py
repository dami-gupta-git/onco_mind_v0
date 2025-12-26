"""cBioPortal co-mutation evidence model."""

from pydantic import BaseModel, Field
from typing import Any


class CoMutationEntry(BaseModel):
    """A single co-mutation or mutual exclusivity entry."""

    gene: str
    count: int = 0
    pct: float = 0.0
    odds_ratio: float | None = None


class CBioPortalEvidence(BaseModel):
    """Evidence from cBioPortal showing mutation prevalence and co-occurrence patterns."""

    gene: str
    variant: str | None = None
    tumor_type: str | None = None
    study_id: str | None = None

    # Prevalence in tumor type
    total_samples: int = 0
    samples_with_gene_mutation: int = 0
    samples_with_exact_variant: int = 0
    gene_prevalence_pct: float = 0.0
    variant_prevalence_pct: float = 0.0

    # Co-occurring mutations (genes frequently mutated together)
    co_occurring: list[CoMutationEntry] = Field(default_factory=list)

    # Mutually exclusive mutations (genes rarely mutated together)
    mutually_exclusive: list[CoMutationEntry] = Field(default_factory=list)

    def has_data(self) -> bool:
        """Check if there is meaningful data."""
        return self.total_samples > 0

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "gene": self.gene,
            "variant": self.variant,
            "tumor_type": self.tumor_type,
            "study_id": self.study_id,
            "total_samples": self.total_samples,
            "samples_with_gene_mutation": self.samples_with_gene_mutation,
            "samples_with_exact_variant": self.samples_with_exact_variant,
            "gene_prevalence_pct": self.gene_prevalence_pct,
            "variant_prevalence_pct": self.variant_prevalence_pct,
            "co_occurring": [c.model_dump() for c in self.co_occurring],
            "mutually_exclusive": [m.model_dump() for m in self.mutually_exclusive],
        }

    def get_study_url(self) -> str | None:
        """Get the cBioPortal study URL."""
        if self.study_id:
            return f"https://www.cbioportal.org/study/summary?id={self.study_id}"
        return None

    def to_prompt_context(self) -> str:
        """Format cBioPortal evidence for LLM prompt.

        Returns:
            Formatted string with prevalence and co-mutation context
        """
        if not self.has_data():
            return "No cBioPortal data available for this variant."

        lines = []

        # Build source citation that must be used inline with statistics
        tumor_str = self.tumor_type or "pan-cancer"
        study_url = self.get_study_url()
        if study_url:
            source_cite = f"[cBioPortal: {self.study_id}]({study_url})"
        else:
            source_cite = f"cBioPortal ({self.study_id})"

        lines.append(f"Cohort: {self.total_samples} {tumor_str} samples")
        lines.append("")

        # Prevalence - embed source citation directly with each statistic
        lines.append("PREVALENCE (copy the source citation when quoting these numbers):")
        lines.append(f"  {self.gene} mutations: {self.gene_prevalence_pct:.1f}% ({self.samples_with_gene_mutation}/{self.total_samples}) - cite as: ({source_cite})")
        if self.variant and self.samples_with_exact_variant > 0:
            lines.append(f"  {self.gene} {self.variant} specifically: {self.variant_prevalence_pct:.1f}% ({self.samples_with_exact_variant}/{self.total_samples}) - cite as: ({source_cite})")
        elif self.variant:
            lines.append(f"  {self.gene} {self.variant} specifically: Not observed in this cohort")
        lines.append("")

        # Co-occurring mutations
        if self.co_occurring:
            lines.append("CO-OCCURRING MUTATIONS (cite source when quoting these percentages):")
            for co in self.co_occurring[:5]:
                or_str = f", OR={co.odds_ratio:.2f}" if co.odds_ratio else ""
                lines.append(f"  - {co.gene}: {co.pct:.1f}% of {self.gene}-mutant samples ({co.count} cases{or_str}) - cite as: ({source_cite})")
            if len(self.co_occurring) > 5:
                lines.append(f"  ... and {len(self.co_occurring) - 5} more")
            lines.append("")

        # Mutually exclusive mutations
        if self.mutually_exclusive:
            lines.append("MUTUALLY EXCLUSIVE MUTATIONS (cite source when quoting these percentages):")
            for me in self.mutually_exclusive[:5]:
                or_str = f", OR={me.odds_ratio:.2f}" if me.odds_ratio else ""
                lines.append(f"  - {me.gene}: {me.pct:.1f}% co-occurrence ({me.count} cases{or_str}) - cite as: ({source_cite})")
            if len(self.mutually_exclusive) > 5:
                lines.append(f"  ... and {len(self.mutually_exclusive) - 5} more")
            lines.append("")

        # Interpretation hints for LLM
        if self.co_occurring or self.mutually_exclusive:
            lines.append("INTERPRETATION (cite source for all statistics):")
            if self.co_occurring:
                top_co = self.co_occurring[0]
                lines.append(f"  - {self.gene} mutations frequently co-occur with {top_co.gene} ({top_co.pct:.1f}%) - cite as: ({source_cite})")
            if self.mutually_exclusive:
                top_me = self.mutually_exclusive[0]
                lines.append(f"  - {self.gene} and {top_me.gene} mutations are mutually exclusive ({top_me.pct:.1f}%) - cite as: ({source_cite})")

        return "\n".join(lines)

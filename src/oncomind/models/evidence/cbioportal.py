"""cBioPortal co-mutation evidence model."""

from pydantic import BaseModel, Field
from typing import Any

from oncomind.models.evidence.base import EvidenceItemBase


class CoMutationEntry(BaseModel):
    """A single co-mutation or mutual exclusivity entry."""

    gene: str
    count: int = 0
    pct: float = 0.0
    odds_ratio: float | None = None


class CBioPortalEvidence(EvidenceItemBase):
    """Evidence from cBioPortal showing mutation prevalence and co-occurrence patterns."""

    gene: str
    variant: str | None = None
    tumor_type: str | None = None
    study_id: str | None = None
    study_name: str | None = None

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
            "study_name": self.study_name,
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
        """Format cBioPortal evidence for LLM prompt (compact version).

        Returns:
            Formatted string with prevalence and co-mutation context
        """
        if not self.has_data():
            return "No cBioPortal data available for this variant."

        lines = []

        # Build source citation once at the end
        tumor_str = self.tumor_type or "pan-cancer"
        study_url = self.get_study_url()
        study_display = self.study_name or self.study_id
        if study_url:
            source_cite = f"[cBioPortal: {study_display}]({study_url})"
        else:
            source_cite = f"cBioPortal ({study_display})"

        # Compact header
        lines.append(f"Study: {self.study_id} ({self.total_samples} {tumor_str} samples)")

        # Prevalence - compact format
        lines.append(f"Prevalence: {self.gene} {self.gene_prevalence_pct:.1f}% ({self.samples_with_gene_mutation}/{self.total_samples})")
        if self.variant and self.samples_with_exact_variant > 0:
            lines.append(f"  {self.gene} {self.variant}: {self.variant_prevalence_pct:.1f}% ({self.samples_with_exact_variant}/{self.total_samples})")

        # Co-occurring mutations - top 3 only, compact format
        if self.co_occurring:
            co_strs = []
            for co in self.co_occurring[:3]:
                co_strs.append(f"{co.gene} {co.pct:.1f}%")
            lines.append(f"Co-occurring: {', '.join(co_strs)}")

        # Mutually exclusive - top 2 only
        if self.mutually_exclusive:
            me_strs = []
            for me in self.mutually_exclusive[:2]:
                me_strs.append(f"{me.gene}")
            lines.append(f"Mutually exclusive: {', '.join(me_strs)}")

        # Single citation at end
        lines.append(f"Source: {source_cite}")

        return "\n".join(lines)

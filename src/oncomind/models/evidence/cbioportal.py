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

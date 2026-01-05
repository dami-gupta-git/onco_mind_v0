"""Tumor-specific evidence match aggregation.

Tracks which evidence sources have tumor-specific data and at what match level.
Used by gap_detector to populate the "Evidence In [Tumor]" well-characterized aspect.
"""

from pydantic import BaseModel, Field
from typing import Literal


MatchLevel = Literal["variant", "codon", "gene"]


class SourceMatch(BaseModel):
    """Match info from a single evidence source."""

    source: str = Field(..., description="Source name (e.g., 'CIViC', 'FDA', 'VICC', 'CGI')")
    count: int = Field(default=0, description="Number of matching items from this source")
    variant_matches: int = Field(default=0, description="Items matching at variant level")
    codon_matches: int = Field(default=0, description="Items matching at codon level")
    gene_matches: int = Field(default=0, description="Items matching at gene level only")


class TumorEvidenceMatch(BaseModel):
    """Aggregated tumor-specific evidence matches across sources.

    Tracks which sources have evidence matching the queried tumor type,
    and at what locus level (variant vs gene).

    Example:
        For BRAF V600E in Melanoma, this might show:
        - 3 CIViC items (2 variant-level, 1 gene-level)
        - 1 FDA approval (variant-level)
        - 2 VICC items (variant-level)
    """

    tumor_type: str = Field(..., description="The queried tumor type")
    has_tumor_evidence: bool = Field(default=False, description="Whether any tumor-specific evidence exists")

    # Per-source match info
    civic_assertions: SourceMatch | None = Field(default=None)
    civic_evidence: SourceMatch | None = Field(default=None)
    fda_approvals: SourceMatch | None = Field(default=None)
    vicc_evidence: SourceMatch | None = Field(default=None)
    cgi_biomarkers: SourceMatch | None = Field(default=None)

    # Aggregate counts
    total_matches: int = Field(default=0, description="Total items matching tumor type")
    total_variant_level: int = Field(default=0, description="Items at variant match level")
    total_codon_level: int = Field(default=0, description="Items at codon match level")
    total_gene_level: int = Field(default=0, description="Items at gene match level only")

    @property
    def matches_on_str(self) -> str | None:
        """Build matches_on string (e.g., '2 variant, 1 codon, 1 gene')."""
        if self.total_matches == 0:
            return None

        parts = []
        if self.total_variant_level > 0:
            parts.append(f"{self.total_variant_level} variant")
        if self.total_codon_level > 0:
            parts.append(f"{self.total_codon_level} codon")
        if self.total_gene_level > 0:
            parts.append(f"{self.total_gene_level} gene")

        return ", ".join(parts) if parts else None

    @property
    def tumor_count_str(self) -> str | None:
        """Build tumor count string (e.g., '6 tumor')."""
        if self.total_matches == 0:
            return None
        return f"{self.total_matches} tumor"

    @property
    def sources_str(self) -> str:
        """Build source summary string (e.g., '3 CIViC, 1 FDA')."""
        parts = []
        for source_match in [
            self.civic_assertions,
            self.civic_evidence,
            self.fda_approvals,
            self.vicc_evidence,
            self.cgi_biomarkers,
        ]:
            if source_match and source_match.count > 0:
                parts.append(f"{source_match.count} {source_match.source}")

        return ", ".join(parts) if parts else "No tumor-specific evidence"

    def add_source_match(
        self,
        source: str,
        count: int,
        variant_matches: int,
        codon_matches: int,
        gene_matches: int
    ) -> None:
        """Add match info for a source and update totals."""
        source_match = SourceMatch(
            source=source,
            count=count,
            variant_matches=variant_matches,
            codon_matches=codon_matches,
            gene_matches=gene_matches,
        )

        # Store in appropriate field
        source_lower = source.lower()
        if source_lower == "civic_assertions" or source == "CIViC Assertions":
            self.civic_assertions = source_match
        elif source_lower == "civic_evidence" or source == "CIViC":
            self.civic_evidence = source_match
        elif source_lower == "fda" or source == "FDA":
            self.fda_approvals = source_match
        elif source_lower == "vicc" or source == "VICC":
            self.vicc_evidence = source_match
        elif source_lower == "cgi" or source == "CGI":
            self.cgi_biomarkers = source_match

        # Update totals
        self.total_matches += count
        self.total_variant_level += variant_matches
        self.total_codon_level += codon_matches
        self.total_gene_level += gene_matches

        if count > 0:
            self.has_tumor_evidence = True

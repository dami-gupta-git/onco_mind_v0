"""MSK Cancer Hotspots evidence model.

Data source: cancerhotspots.org (Memorial Sloan Kettering)
File format: Tab-separated values with columns:
    Gene, Residue, Type, Variants, Q-value, Samples, Tumor Type Composition
"""

from pydantic import BaseModel, Field
from typing import Any

from oncomind.models.evidence.base import EvidenceItemBase


class HotspotVariant(BaseModel):
    """A specific amino acid change at a hotspot residue."""

    amino_acid: str  # e.g., "E", "K", "R"
    count: int  # Number of samples with this specific change


class HotspotTumorType(BaseModel):
    """Tumor type composition for a hotspot."""

    tumor_type: str  # e.g., "skin", "lung", "bowel"
    count: int  # Number of samples from this tumor type


class HotspotEntry(BaseModel):
    """A single hotspot entry from MSK Cancer Hotspots.

    Represents a recurrent mutation site (residue) in a gene with statistical
    significance for being a cancer hotspot.
    """

    gene: str  # Gene symbol (e.g., "BRAF", "KRAS")
    residue: str  # Hotspot residue (e.g., "V600", "G12")
    hotspot_type: str  # "single residue" or "in-frame indel"
    q_value: float  # Statistical significance (smaller = more significant)
    total_samples: int  # Total samples with mutations at this residue

    # Breakdown of specific amino acid changes at this residue
    # e.g., V600 -> E:833, M:29, K:24
    variants: list[HotspotVariant] = Field(default_factory=list)

    # Tumor type distribution
    tumor_type_composition: list[HotspotTumorType] = Field(default_factory=list)

    def get_variant_count(self, amino_acid: str) -> int:
        """Get the count for a specific amino acid change.

        Args:
            amino_acid: Single letter amino acid code (e.g., "E", "K")

        Returns:
            Count of samples with this change, or 0 if not found
        """
        for v in self.variants:
            if v.amino_acid.upper() == amino_acid.upper():
                return v.count
        return 0

    def get_tumor_count(self, tumor_type: str) -> int:
        """Get the count for a specific tumor type.

        Args:
            tumor_type: Tumor type string (e.g., "skin", "lung")

        Returns:
            Count of samples from this tumor type, or 0 if not found
        """
        tumor_lower = tumor_type.lower()
        for t in self.tumor_type_composition:
            if t.tumor_type.lower() == tumor_lower:
                return t.count
        return 0

    def get_top_variants(self, n: int = 3) -> list[HotspotVariant]:
        """Get the top N most frequent variants."""
        return sorted(self.variants, key=lambda v: v.count, reverse=True)[:n]

    def get_top_tumor_types(self, n: int = 3) -> list[HotspotTumorType]:
        """Get the top N most frequent tumor types."""
        return sorted(
            self.tumor_type_composition, key=lambda t: t.count, reverse=True
        )[:n]


class HotspotsEvidence(EvidenceItemBase):
    """Evidence from MSK Cancer Hotspots database.

    Cancer hotspots are recurrent mutation sites that are statistically enriched
    across cancer samples, indicating likely driver mutations.

    Source: https://www.cancerhotspots.org/
    Publication: Chang et al., Cancer Discovery 2017
    """

    gene: str
    variant: str | None = None  # The queried variant (e.g., "V600E")
    queried_residue: str | None = None  # Extracted residue (e.g., "V600")

    # The matching hotspot entry (if found)
    hotspot: HotspotEntry | None = None

    # Match metadata
    is_hotspot: bool = False  # True if variant is at a known hotspot
    is_exact_variant_match: bool = False  # True if exact AA change is in hotspot data

    def has_data(self) -> bool:
        """Check if there is meaningful hotspot data."""
        return self.hotspot is not None

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        result = {
            "gene": self.gene,
            "variant": self.variant,
            "queried_residue": self.queried_residue,
            "is_hotspot": self.is_hotspot,
            "is_exact_variant_match": self.is_exact_variant_match,
        }
        if self.hotspot:
            result["hotspot"] = {
                "residue": self.hotspot.residue,
                "hotspot_type": self.hotspot.hotspot_type,
                "q_value": self.hotspot.q_value,
                "total_samples": self.hotspot.total_samples,
                "top_variants": [
                    {"amino_acid": v.amino_acid, "count": v.count}
                    for v in self.hotspot.get_top_variants(5)
                ],
                "top_tumor_types": [
                    {"tumor_type": t.tumor_type, "count": t.count}
                    for t in self.hotspot.get_top_tumor_types(5)
                ],
            }
        return result

    def get_source_url(self) -> str:
        """Get the Cancer Hotspots URL."""
        return "https://www.cancerhotspots.org/"

    def to_prompt_context(self) -> str:
        """Format hotspot evidence for LLM prompt.

        Returns:
            Formatted string with hotspot context
        """
        if not self.has_data() or not self.hotspot:
            return ""

        lines = []
        h = self.hotspot

        # Header with match status
        if self.is_exact_variant_match:
            lines.append(
                f"CANCER HOTSPOT: {self.gene} {h.residue} is a known cancer hotspot "
                f"(exact variant match, q={h.q_value:.2e})"
            )
        else:
            lines.append(
                f"CANCER HOTSPOT: {self.gene} {h.residue} is a known cancer hotspot "
                f"(codon-level match, q={h.q_value:.2e})"
            )

        # Sample count
        lines.append(f"  Observed in {h.total_samples} cancer samples")

        # Top variants at this residue
        top_vars = h.get_top_variants(4)
        if top_vars:
            var_strs = [f"{h.residue[0]}{h.residue[1:]}{v.amino_acid}:{v.count}" for v in top_vars]
            lines.append(f"  Common changes: {', '.join(var_strs)}")

        # Top tumor types
        top_tumors = h.get_top_tumor_types(4)
        if top_tumors:
            tumor_strs = [f"{t.tumor_type}:{t.count}" for t in top_tumors]
            lines.append(f"  Tumor distribution: {', '.join(tumor_strs)}")

        lines.append(f"  Source: [cancerhotspots.org]({self.get_source_url()})")

        return "\n".join(lines)

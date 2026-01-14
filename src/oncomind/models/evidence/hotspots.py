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

    # Adjacent hotspot data (for variants within ±5 codons of a hotspot)
    adjacent_hotspot: HotspotEntry | None = None  # Nearest hotspot if variant is not at one
    adjacent_distance: int | None = None  # Distance in codons to the adjacent hotspot

    def has_data(self) -> bool:
        """Check if there is meaningful hotspot data."""
        return self.hotspot is not None

    def is_adjacent_to_hotspot(self) -> bool:
        """Check if variant is adjacent to (but not at) a hotspot."""
        return self.adjacent_hotspot is not None and not self.is_hotspot

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        result = {
            "gene": self.gene,
            "variant": self.variant,
            "queried_residue": self.queried_residue,
            "is_hotspot": self.is_hotspot,
            "is_exact_variant_match": self.is_exact_variant_match,
            "is_adjacent_to_hotspot": self.is_adjacent_to_hotspot(),
            "adjacent_distance": self.adjacent_distance,
        }
        if self.hotspot:
            result["hotspot"] = {
                "residue": self.hotspot.residue,
                "hotspot_type": self.hotspot.hotspot_type,
                "q_value": self.hotspot.q_value,
                "total_samples": self.hotspot.total_samples,
                "variants": [
                    {"amino_acid": v.amino_acid, "count": v.count}
                    for v in self.hotspot.variants
                ],
                "tumor_type_composition": [
                    {"tumor_type": t.tumor_type, "count": t.count}
                    for t in self.hotspot.tumor_type_composition
                ],
            }
        if self.adjacent_hotspot:
            result["adjacent_hotspot"] = {
                "residue": self.adjacent_hotspot.residue,
                "hotspot_type": self.adjacent_hotspot.hotspot_type,
                "q_value": self.adjacent_hotspot.q_value,
                "total_samples": self.adjacent_hotspot.total_samples,
                "variants": [
                    {"amino_acid": v.amino_acid, "count": v.count}
                    for v in self.adjacent_hotspot.variants
                ],
                "tumor_type_composition": [
                    {"tumor_type": t.tumor_type, "count": t.count}
                    for t in self.adjacent_hotspot.tumor_type_composition
                ],
            }
        return result

    def get_source_url(self) -> str:
        """Get the Cancer Hotspots URL."""
        return "https://www.cancerhotspots.org/"

    def _tumor_matches_hotspot(self, tumor_type: str | None, hotspot: HotspotEntry) -> bool:
        """Check if the queried tumor type matches the hotspot's tumor distribution.

        Args:
            tumor_type: Clinical tumor type (e.g., "Glioma", "NSCLC", "Melanoma")
            hotspot: The hotspot entry to check

        Returns:
            True if tumor type is None (pan-cancer query) or matches hotspot tumor distribution
        """
        if not tumor_type:
            # No tumor type specified - include all hotspot data
            return True

        from oncomind.config.constants import HOTSPOT_TUMOR_MAPPING

        # Get the hotspot organ terms for this clinical tumor type
        tumor_lower = tumor_type.lower()
        expected_organs = HOTSPOT_TUMOR_MAPPING.get(tumor_lower, [])

        if not expected_organs:
            # Unknown tumor type - be conservative and exclude from LLM
            # (hotspot data will still be shown in UI, just not sent to LLM)
            return False

        # Check if any of the expected organs appear in hotspot tumor composition
        hotspot_organs = {t.tumor_type.lower() for t in hotspot.tumor_type_composition}

        for organ in expected_organs:
            if organ.lower() in hotspot_organs:
                return True

        return False

    def to_prompt_context(self, tumor_type: str | None = None) -> str:
        """Format hotspot evidence for LLM prompt.

        Args:
            tumor_type: Optional tumor type to filter relevance. If the queried tumor
                       type doesn't match the hotspot's tumor distribution, returns empty
                       string to avoid sending irrelevant data to LLM synthesis.

        Returns:
            Formatted string with hotspot context, or empty string if tumor doesn't match
        """
        lines = []

        # Direct hotspot match
        if self.hotspot:
            h = self.hotspot

            # Check if this hotspot is relevant to the queried tumor type
            if not self._tumor_matches_hotspot(tumor_type, h):
                # Hotspot exists but not relevant to this tumor - skip for LLM
                return ""

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

        # Adjacent hotspot (not at hotspot, but near one)
        elif self.is_adjacent_to_hotspot() and self.adjacent_hotspot:
            h = self.adjacent_hotspot

            # Check if the adjacent hotspot is relevant to the queried tumor type
            if not self._tumor_matches_hotspot(tumor_type, h):
                # Adjacent hotspot not relevant to this tumor - skip for LLM
                return ""

            lines.append(
                f"NEAR HOTSPOT: {self.gene} {self.variant} is {self.adjacent_distance} codon(s) "
                f"from known hotspot {h.residue}"
            )
            lines.append(
                f"  The nearby {h.residue} hotspot is observed in {h.total_samples} cancer samples (q={h.q_value:.2e})"
            )

            # Top tumor types for the nearby hotspot
            top_tumors = h.get_top_tumor_types(4)
            if top_tumors:
                tumor_strs = [f"{t.tumor_type}:{t.count}" for t in top_tumors]
                lines.append(f"  Hotspot tumor distribution: {', '.join(tumor_strs)}")

            lines.append(
                f"  Note: Proximity to hotspot suggests possible functional relevance, but requires validation"
            )
            lines.append(f"  Source: [cancerhotspots.org]({self.get_source_url()})")

        return "\n".join(lines)

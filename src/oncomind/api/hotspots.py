"""MSK Cancer Hotspots data loader.

ARCHITECTURE:
    Gene + Variant → Local TSV file → Hotspot evidence

Cancer Hotspots is a resource for identifying statistically significant
recurrent mutations in cancer from large-scale cancer genomics data.

Source: https://www.cancerhotspots.org/
Publication: Chang et al., Cancer Discovery 2017, 2018

Data file: data/hotspots_msk/hotspots.txt (TSV format)
"""

import re
from functools import lru_cache
from pathlib import Path

from oncomind.models.evidence.hotspots import (
    HotspotsEvidence,
    HotspotEntry,
    HotspotVariant,
    HotspotTumorType,
)
from oncomind.config.debug import get_logger

logger = get_logger(__name__)

# Path to the hotspots data file
HOTSPOTS_FILE = Path(__file__).parent.parent.parent.parent / "data" / "hotspots_msk" / "hotspots.txt"

# Regex to extract residue position and reference AA from variant notation
# Matches patterns like: V600E, G12D, p.V600E, L858R
_VARIANT_PATTERN = re.compile(r"^p?\.?([A-Z])(\d+)([A-Z*]?)$", re.IGNORECASE)


class HotspotsError(Exception):
    """Exception raised for Hotspots data errors."""
    pass


def _parse_variants_field(variants_str: str) -> list[HotspotVariant]:
    """Parse the Variants field from the TSV.

    Format: "E:833|M:29|K:24|R:4|V:4|G:3"

    Args:
        variants_str: Pipe-separated variant:count pairs

    Returns:
        List of HotspotVariant objects
    """
    if not variants_str:
        return []

    result = []
    for part in variants_str.split("|"):
        if ":" not in part:
            continue
        try:
            aa, count_str = part.split(":", 1)
            result.append(HotspotVariant(
                amino_acid=aa.strip(),
                count=int(count_str.strip()),
            ))
        except (ValueError, IndexError):
            continue

    return result


def _parse_tumor_composition(composition_str: str) -> list[HotspotTumorType]:
    """Parse the Tumor Type Composition field from the TSV.

    Format: "skin:357|thyroid:316|bowel:113|lung:33"

    Args:
        composition_str: Pipe-separated tumor:count pairs

    Returns:
        List of HotspotTumorType objects
    """
    if not composition_str:
        return []

    result = []
    for part in composition_str.split("|"):
        if ":" not in part:
            continue
        try:
            tumor, count_str = part.split(":", 1)
            result.append(HotspotTumorType(
                tumor_type=tumor.strip(),
                count=int(count_str.strip()),
            ))
        except (ValueError, IndexError):
            continue

    return result


def _parse_hotspot_row(row: dict) -> HotspotEntry | None:
    """Parse a single row from the hotspots TSV.

    Args:
        row: Dictionary with column names as keys

    Returns:
        HotspotEntry or None if parsing fails
    """
    try:
        gene = row.get("Gene", "").strip()
        residue = row.get("Residue", "").strip()

        if not gene or not residue:
            return None

        # Parse Q-value (scientific notation)
        q_value_str = row.get("Q-value", "0")
        try:
            q_value = float(q_value_str)
        except ValueError:
            q_value = 0.0

        # Parse total samples
        samples_str = row.get("Samples", "0")
        try:
            total_samples = int(samples_str)
        except ValueError:
            total_samples = 0

        return HotspotEntry(
            gene=gene.upper(),
            residue=residue,
            hotspot_type=row.get("Type", "single residue"),
            q_value=q_value,
            total_samples=total_samples,
            variants=_parse_variants_field(row.get("Variants", "")),
            tumor_type_composition=_parse_tumor_composition(row.get("Tumor Type Composition", "")),
        )
    except Exception as e:
        logger.debug(f"Failed to parse hotspot row: {e}")
        return None


@lru_cache(maxsize=1)
def _load_hotspots_data() -> dict[str, list[HotspotEntry]]:
    """Load and cache the hotspots data from the TSV file.

    Returns:
        Dictionary mapping gene symbols to list of HotspotEntry objects
    """
    if not HOTSPOTS_FILE.exists():
        logger.warning(f"Hotspots file not found: {HOTSPOTS_FILE}")
        return {}

    hotspots_by_gene: dict[str, list[HotspotEntry]] = {}

    try:
        with open(HOTSPOTS_FILE, "r", encoding="utf-8-sig") as f:
            # Read header (utf-8-sig handles BOM automatically)
            header_line = f.readline().strip()
            headers = header_line.split("\t")

            # Parse each row
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < len(headers):
                    continue

                row = dict(zip(headers, parts))
                entry = _parse_hotspot_row(row)

                if entry:
                    gene = entry.gene.upper()
                    if gene not in hotspots_by_gene:
                        hotspots_by_gene[gene] = []
                    hotspots_by_gene[gene].append(entry)

        logger.info(f"Loaded {sum(len(v) for v in hotspots_by_gene.values())} hotspots for {len(hotspots_by_gene)} genes")
        return hotspots_by_gene

    except Exception as e:
        logger.error(f"Failed to load hotspots file: {e}")
        return {}


def _extract_residue(variant: str) -> tuple[str | None, str | None, str | None]:
    """Extract reference AA, position, and alternate AA from variant notation.

    Args:
        variant: Variant notation (e.g., "V600E", "G12D", "p.L858R")

    Returns:
        Tuple of (ref_aa, position, alt_aa) or (None, None, None) if not parseable
    """
    if not variant:
        return None, None, None

    # Clean up variant string
    clean = variant.strip()
    if clean.lower().startswith("p."):
        clean = clean[2:]

    match = _VARIANT_PATTERN.match(clean)
    if match:
        ref_aa = match.group(1).upper()
        position = match.group(2)
        alt_aa = match.group(3).upper() if match.group(3) else None
        return ref_aa, position, alt_aa

    return None, None, None


class HotspotsClient:
    """Client for Cancer Hotspots data.

    This is a file-based client that loads data from the local TSV file.
    No network requests are made.
    """

    def __init__(self):
        """Initialize the Hotspots client."""
        self._data: dict[str, list[HotspotEntry]] | None = None

    def _ensure_loaded(self) -> dict[str, list[HotspotEntry]]:
        """Ensure data is loaded."""
        if self._data is None:
            self._data = _load_hotspots_data()
        return self._data

    def get_hotspots_for_gene(self, gene: str) -> list[HotspotEntry]:
        """Get all hotspot entries for a gene.

        Args:
            gene: Gene symbol (case-insensitive)

        Returns:
            List of HotspotEntry objects for this gene
        """
        data = self._ensure_loaded()
        return data.get(gene.upper(), [])

    def get_hotspot_codons(self, gene: str) -> list[int]:
        """Get list of hotspot codon positions for a gene.

        This replaces the hardcoded CANCER_HOTSPOTS dict.

        Args:
            gene: Gene symbol

        Returns:
            List of codon positions that are hotspots
        """
        entries = self.get_hotspots_for_gene(gene)
        codons = []

        for entry in entries:
            # Extract position from residue (e.g., "V600" -> 600)
            _, position, _ = _extract_residue(entry.residue)
            if position:
                try:
                    codons.append(int(position))
                except ValueError:
                    pass

        return sorted(set(codons))

    def is_hotspot(self, gene: str, variant: str) -> bool:
        """Check if a variant is at a known hotspot.

        Args:
            gene: Gene symbol
            variant: Variant notation (e.g., "V600E")

        Returns:
            True if variant position is a hotspot
        """
        _, query_pos, _ = _extract_residue(variant)
        if not query_pos:
            return False

        hotspot_codons = self.get_hotspot_codons(gene)
        try:
            return int(query_pos) in hotspot_codons
        except ValueError:
            return False

    def fetch_hotspot_evidence(
        self,
        gene: str,
        variant: str | None = None,
    ) -> HotspotsEvidence | None:
        """Fetch hotspot evidence for a gene/variant.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional variant notation (e.g., "V600E")

        Returns:
            HotspotsEvidence if hotspot data exists, None otherwise
        """
        gene_upper = gene.upper()
        entries = self.get_hotspots_for_gene(gene_upper)

        if not entries:
            return None

        # If no variant, just return gene-level evidence
        if not variant:
            return HotspotsEvidence(
                gene=gene_upper,
                variant=None,
                queried_residue=None,
                hotspot=None,
                is_hotspot=False,
                is_exact_variant_match=False,
            )

        # Extract residue info from variant
        ref_aa, position, alt_aa = _extract_residue(variant)
        if not position:
            return HotspotsEvidence(
                gene=gene_upper,
                variant=variant,
                queried_residue=None,
                hotspot=None,
                is_hotspot=False,
                is_exact_variant_match=False,
            )

        queried_residue = f"{ref_aa}{position}" if ref_aa else position

        # Find matching hotspot entry
        matching_hotspot: HotspotEntry | None = None
        is_exact_match = False

        for entry in entries:
            # Extract position from hotspot residue
            hotspot_ref, hotspot_pos, _ = _extract_residue(entry.residue)
            if not hotspot_pos:
                continue

            # Check if position matches
            if hotspot_pos == position:
                matching_hotspot = entry

                # Check for exact variant match (same AA change)
                if alt_aa:
                    for v in entry.variants:
                        if v.amino_acid.upper() == alt_aa:
                            is_exact_match = True
                            break

                break  # Found the hotspot at this position

        return HotspotsEvidence(
            gene=gene_upper,
            variant=variant,
            queried_residue=queried_residue,
            hotspot=matching_hotspot,
            is_hotspot=matching_hotspot is not None,
            is_exact_variant_match=is_exact_match,
        )


# Singleton instance for convenience
_client: HotspotsClient | None = None


def get_hotspots_client() -> HotspotsClient:
    """Get or create the singleton HotspotsClient."""
    global _client
    if _client is None:
        _client = HotspotsClient()
    return _client


# Convenience functions that mirror gene_context.py interface


def is_hotspot_variant(gene: str, variant: str) -> bool:
    """Check if a variant is at a known cancer hotspot.

    This replaces the function in gene_context.py.

    Args:
        gene: Gene symbol (case-insensitive)
        variant: Variant notation (e.g., V600E, G12D)

    Returns:
        True if variant is at a known hotspot codon
    """
    client = get_hotspots_client()
    return client.is_hotspot(gene, variant)


def get_cancer_hotspots(gene: str) -> list[int]:
    """Get list of hotspot codon positions for a gene.

    This replaces CANCER_HOTSPOTS dict in gene_context.py.

    Args:
        gene: Gene symbol

    Returns:
        List of codon positions that are hotspots
    """
    client = get_hotspots_client()
    return client.get_hotspot_codons(gene)

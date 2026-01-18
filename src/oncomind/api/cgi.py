"""Cancer Genome Interpreter (CGI) biomarkers client.

ARCHITECTURE:
    Gene + Variant → CGI Biomarkers TSV → FDA/NCCN approval status

Fetches FDA approval and guideline information from the Cancer Genome Interpreter
biomarkers database, which provides curated variant-drug associations with
explicit regulatory approval status.

Key Design:
- Downloads and caches the biomarkers TSV file
- Parses variant patterns to match specific mutations (e.g., G719S matches G719.)
- Returns structured CGIBiomarker objects with approval status
- Complements FDA label search which uses generic text (e.g., "non-resistant mutations")
"""

import csv
import re
from datetime import datetime, timedelta
from typing import Any

import httpx

from oncomind.config.constants import CGI_DATA_DIR, CGI_BIOMARKERS_FILE
from oncomind.models.evidence.base import (
    is_pan_cancer_term,
    tumor_types_match,
    determine_locus_match,
    extract_variant_position,
)


class CGIError(Exception):
    """Exception raised for CGI-related errors."""

    pass


class CGIBiomarker:
    """A biomarker entry from CGI database."""

    def __init__(
        self,
        gene: str,
        alteration: str,
        drug: str,
        drug_status: str,
        association: str,
        evidence_level: str,
        source: str,
        tumor_type: str,
        tumor_type_full: str,
        comments: str | None = None,
        drug_full_name: str | None = None,
    ):
        self.gene = gene
        self.alteration = alteration
        self.drug = drug
        self.drug_status = drug_status  # "Approved", "Clinical trial", etc.
        self.association = association  # "Responsive", "Resistant"
        self.evidence_level = evidence_level  # "FDA guidelines", "NCCN guidelines", etc.
        self.source = source
        self.tumor_type = tumor_type
        self.tumor_type_full = tumor_type_full
        self.comments = comments  # FDA label text or description
        self.drug_full_name = drug_full_name  # Drug name with mechanism (e.g., "Tazemetostat (EZH2 inhibitor)")

    def is_fda_approved(self) -> bool:
        """Check if this biomarker represents FDA-approved therapy."""
        return self.drug_status == "Approved" and (
            "FDA" in self.evidence_level.upper()
            or self.evidence_level.upper() in ["NCCN GUIDELINES", "NCCN/CGC GUIDELINES"]
        )

    def get_fda_url(self) -> str | None:
        """Extract FDA URL from source field if present.

        Source field format: "FDA:https://www.fda.gov/...;PMID:12345"
        """
        if not self.source:
            return None
        # Source can contain multiple references separated by semicolons
        # Look for FDA URL pattern
        for part in self.source.split(";"):
            part = part.strip()
            if part.startswith("FDA:"):
                url = part[4:]  # Remove "FDA:" prefix
                if url.startswith("http"):
                    return url
            elif "fda.gov" in part.lower() and part.startswith("http"):
                return part
        return None

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "gene": self.gene,
            "alteration": self.alteration,
            "drug": self.drug,
            "drug_status": self.drug_status,
            "association": self.association,
            "evidence_level": self.evidence_level,
            "source": self.source,
            "tumor_type": self.tumor_type,
            "tumor_type_full": self.tumor_type_full,
            "fda_approved": self.is_fda_approved(),
            "fda_url": self.get_fda_url(),
            "comments": self.comments,
            "drug_full_name": self.drug_full_name,
        }


class CGIClient:
    """Client for Cancer Genome Interpreter biomarkers database.

    CGI provides curated biomarker-drug associations with explicit
    FDA approval status, complementing the FDA label API which uses
    generic text descriptions.

    Data is downloaded once and cached locally in data/cgi/.
    """

    BIOMARKERS_URL = "https://www.cancergenomeinterpreter.org/data/biomarkers/cgi_biomarkers_latest.tsv"
    CACHE_DIR = CGI_DATA_DIR
    CACHE_FILE = CGI_BIOMARKERS_FILE
    CACHE_MAX_AGE = timedelta(days=7)  # Re-download after 7 days

    def __init__(self, timeout: float = 30.0):
        """Initialize the CGI client.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout
        self._biomarkers: list[dict[str, str]] | None = None

    def _cache_is_valid(self) -> bool:
        """Check if the cached file exists and is recent enough."""
        if not self.CACHE_FILE.exists():
            return False
        mtime = datetime.fromtimestamp(self.CACHE_FILE.stat().st_mtime)
        return datetime.now() - mtime < self.CACHE_MAX_AGE

    def _download_biomarkers(self) -> None:
        """Download the biomarkers TSV file."""
        self.CACHE_DIR.mkdir(parents=True, exist_ok=True)

        with httpx.Client(timeout=self.timeout) as client:
            response = client.get(self.BIOMARKERS_URL)
            response.raise_for_status()
            self.CACHE_FILE.write_text(response.text)

    def _load_biomarkers(self) -> list[dict[str, str]]:
        """Load and parse the biomarkers TSV file."""
        if not self._cache_is_valid():
            try:
                self._download_biomarkers()
            except Exception as e:
                if not self.CACHE_FILE.exists():
                    raise CGIError(f"Failed to download CGI biomarkers: {e}")
                # Use stale cache if download fails

        with open(self.CACHE_FILE, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            return list(reader)

    def _get_biomarkers(self) -> list[dict[str, str]]:
        """Get biomarkers, loading from cache if needed."""
        if self._biomarkers is None:
            self._biomarkers = self._load_biomarkers()
        return self._biomarkers

    def _variant_matches(self, cgi_alteration: str, gene: str, variant: str) -> bool:
        """Check if a CGI alteration pattern matches a specific variant.

        CGI uses patterns like:
        - "EGFR:G719." matches G719S, G719A, G719C, etc.
        - "EGFR:G719A,G719S,G719C" matches any of these specific variants
        - "EGFR:L858R" matches exactly L858R
        - "KRAS:.12.,.13." matches any mutation at position 12 or 13 (G12D, G13D, etc.)
        - "KRAS:." matches any KRAS mutation
        - "KIT:449-514,788-828" matches any mutation in exon ranges (e.g., D816V is in 788-828)

        Consequence-qualified alterations (e.g., "::consequence::inframe_insertion:762-823")
        only match variants of that type, not point mutations.

        Args:
            cgi_alteration: CGI alteration string (e.g., "EGFR:G719.,L858R")
            gene: Gene symbol (e.g., "EGFR")
            variant: Variant notation (e.g., "G719S")

        Returns:
            True if the variant matches the CGI pattern
        """
        if not cgi_alteration:
            return False

        # Check if this is a consequence-qualified alteration (e.g., inframe_insertion, inframe_deletion)
        # These should NOT match point mutations - only variants of the specified type
        cgi_alt_lower = cgi_alteration.lower()
        has_consequence_qualifier = "::consequence::" in cgi_alt_lower
        is_insertion_specific = "inframe_insertion" in cgi_alt_lower
        is_deletion_specific = "inframe_deletion" in cgi_alt_lower

        # Parse the CGI alteration string: "GENE:variant1,variant2,..."
        # Handle both formats: "EGFR:V600E" and just "V600E" after gene match
        gene_upper = gene.upper()
        variant_upper = variant.upper().replace("P.", "")  # Remove p. prefix if present

        # Point mutations (like C797S) should not match insertion/deletion-specific entries
        # Point mutations have format: single letter + digits + single letter (e.g., C797S, V600E)
        is_point_mutation = bool(re.match(r'^[A-Z]\d+[A-Z]$', variant_upper))

        if has_consequence_qualifier and is_point_mutation:
            # Don't match point mutations to insertion/deletion-specific biomarkers
            if is_insertion_specific or is_deletion_specific:
                return False

        # Extract variant position for range matching (e.g., 816 from D816V)
        variant_position = None
        variant_match = re.match(r'^([A-Z]?)(\d+)([A-Z]?)$', variant_upper)
        if variant_match:
            variant_position = int(variant_match.group(2))

        # Split by comma to get individual variants
        parts = cgi_alteration.replace(f"{gene_upper}:", "").split(",")

        for part in parts:
            part = part.strip().upper()
            if not part:
                continue

            # Remove gene prefix if still present
            if ":" in part:
                part = part.split(":")[-1]

            # Exact match
            if part == variant_upper:
                return True

            # Wildcard for any mutation in gene: "." alone matches any variant
            if part == ".":
                return True

            # Range pattern: "449-514" or "788-828" matches positions in that range
            # Used for exon-based patterns like "KIT:449-514,550-592,627-664,664-714,788-828"
            range_match = re.match(r'^(\d+)-(\d+)$', part)
            if range_match and variant_position is not None:
                range_start = int(range_match.group(1))
                range_end = int(range_match.group(2))
                if range_start <= variant_position <= range_end:
                    return True

            # Pattern match: "G719." matches G719S, G719A, etc.
            # The dot represents a wildcard for any amino acid
            if part.endswith(".") and not part.startswith("."):
                # Extract base pattern (e.g., "G719" from "G719.")
                base_pattern = part[:-1]
                # Check if variant starts with base and has exactly one more character
                if variant_upper.startswith(base_pattern) and len(variant_upper) == len(
                    base_pattern
                ) + 1:
                    return True

            # Position-based wildcard: ".13." matches any mutation at position 13
            # Format: .{position}. where position is a number
            # This matches variants like G13D, G13C, G13V, etc.
            if part.startswith(".") and part.endswith(".") and len(part) > 2:
                position_str = part[1:-1]  # Extract position number (e.g., "13" from ".13.")
                if position_str.isdigit():
                    position = position_str
                    # Extract position from variant (e.g., "13" from "G13D")
                    # Variant format: {ref_aa}{position}{alt_aa} like G13D
                    pos_match = re.match(r'^([A-Z])(\d+)([A-Z])$', variant_upper)
                    if pos_match:
                        var_pos = pos_match.group(2)
                        if var_pos == position:
                            return True

        return False

    def _determine_locus_match(self, cgi_alteration: str, gene: str, variant: str) -> str:
        """Determine the match specificity level for a CGI alteration.

        CGI-specific parsing handles:
        - Comma-separated patterns (e.g., "EGFR:G719.,L858R")
        - Wildcard patterns: "G719." for codon-level, "." for gene-level
        - Range patterns: "449-514,788-828" for exon-based matching
        - Consequence-qualified ranges (e.g., "::consequence::inframe_insertion:762-823")
          only match variants of that consequence type, not point mutations

        After parsing, delegates to the core determine_locus_match() function.

        Args:
            cgi_alteration: CGI alteration string (e.g., "EGFR:G719.,L858R")
            gene: Gene symbol (e.g., "EGFR")
            variant: Variant notation (e.g., "G719S")

        Returns:
            Match level: 'variant' (exact), 'codon' (same position/pattern), or 'gene' (gene-only)
        """
        if not cgi_alteration:
            return "gene"

        # Check if this is a consequence-qualified alteration (e.g., inframe_insertion, inframe_deletion)
        # These should NOT match point mutations - only variants of the specified type
        cgi_alt_lower = cgi_alteration.lower()
        has_consequence_qualifier = "::consequence::" in cgi_alt_lower
        is_insertion_specific = "inframe_insertion" in cgi_alt_lower
        is_deletion_specific = "inframe_deletion" in cgi_alt_lower

        # Point mutations (like C797S) should not match insertion/deletion-specific entries
        # Point mutations have format: single letter + digits + single letter (e.g., C797S, V600E)
        variant_upper = variant.upper().replace("P.", "")
        is_point_mutation = bool(re.match(r'^[A-Z]\d+[A-Z]$', variant_upper))

        if has_consequence_qualifier and is_point_mutation:
            # Don't match point mutations to insertion/deletion-specific biomarkers
            if is_insertion_specific or is_deletion_specific:
                return "gene"

        gene_upper = gene.upper()
        queried_position = extract_variant_position(variant)

        # Extract numeric position for range matching
        variant_position_int = None
        if queried_position and queried_position.isdigit():
            variant_position_int = int(queried_position)

        # Split by comma to get individual variants
        parts = cgi_alteration.replace(f"{gene_upper}:", "").split(",")

        best_match = "gene"

        for part in parts:
            part = part.strip().upper()
            if not part:
                continue

            if ":" in part:
                part = part.split(":")[-1]

            # Wildcard for any mutation: "." = gene level (skip, already default)
            if part == ".":
                continue

            # Range pattern: "449-514" or "788-828" = codon level (exon-based)
            # Matches if variant position falls within the range
            range_match = re.match(r'^(\d+)-(\d+)$', part)
            if range_match and variant_position_int is not None:
                range_start = int(range_match.group(1))
                range_end = int(range_match.group(2))
                if range_start <= variant_position_int <= range_end:
                    if best_match != "variant":
                        best_match = "codon"
                    continue

            # CGI-specific: Pattern match "G719." = codon level (same position, any AA)
            if part.endswith(".") and not part.startswith("."):
                base_pattern = part[:-1]
                if variant_upper.startswith(base_pattern):
                    if best_match != "variant":
                        best_match = "codon"
                    continue

            # CGI-specific: Position wildcard ".13." = codon level
            if part.startswith(".") and part.endswith(".") and len(part) > 2:
                position_str = part[1:-1]
                if position_str.isdigit() and queried_position == position_str:
                    if best_match != "variant":
                        best_match = "codon"
                    continue

            # For non-wildcard patterns, use the core matching function
            match_result = determine_locus_match(part, variant_upper)
            if match_result == "variant":
                return "variant"  # Exact match, return immediately
            if match_result == "codon" and best_match != "variant":
                best_match = "codon"

        return best_match

    def fetch_biomarkers(
        self, gene: str, variant: str, tumor_type: str | None = None
    ) -> list[CGIBiomarker]:
        """Fetch CGI biomarkers for a gene/variant combination.

        Args:
            gene: Gene symbol (e.g., "EGFR")
            variant: Variant notation (e.g., "G719S")
            tumor_type: Optional tumor type to filter results

        Returns:
            List of matching CGIBiomarker objects
        """
        biomarkers = self._get_biomarkers()
        matches = []
        gene_upper = gene.upper()

        for row in biomarkers:
            # Check gene match
            if row.get("Gene", "").upper() != gene_upper:
                continue

            # Check variant match
            alteration = row.get("Alteration", "")
            if not self._variant_matches(alteration, gene, variant):
                continue

            # Check tumor type match if specified
            cgi_tumor = row.get("Primary Tumor type", "")
            if tumor_type and not is_pan_cancer_term(cgi_tumor) and not tumor_types_match(cgi_tumor, tumor_type):
                continue

            # Get drug name, falling back to drug family if Drug is empty/placeholder
            drug = row.get("Drug", "")
            if not drug or drug == "[]":
                # Use drug family as fallback, removing brackets if present
                drug_family = row.get("Drug family", "")
                if drug_family:
                    # Clean up brackets: "[EGFR TK inhibitor]" -> "EGFR TK inhibitor"
                    drug = drug_family.strip("[]")
                else:
                    # Skip entries with no drug information
                    continue

            # Create biomarker object
            # Get comments (FDA label text) and drug_full_name if available
            comments = row.get("Comments", "")
            drug_full_name = row.get("Drug full name", "")

            matches.append(
                CGIBiomarker(
                    gene=row.get("Gene", ""),
                    alteration=alteration,
                    drug=drug,
                    drug_status=row.get("Drug status", ""),
                    association=row.get("Association", ""),
                    evidence_level=row.get("Evidence level", ""),
                    source=row.get("Source", ""),
                    tumor_type=row.get("Primary Tumor type", ""),
                    tumor_type_full=row.get("Primary Tumor type full name", ""),
                    comments=comments if comments else None,
                    drug_full_name=drug_full_name if drug_full_name else None,
                )
            )

        return matches

    def fetch_fda_approved(
        self, gene: str, variant: str, tumor_type: str | None = None
    ) -> list[CGIBiomarker]:
        """Fetch only FDA-approved biomarkers for a gene/variant.

        This is a convenience method that filters for:
        - Drug status = "Approved"
        - Evidence level contains "FDA" or is "NCCN guidelines"
        - Association = "Responsive" (excludes resistance markers)

        Args:
            gene: Gene symbol (e.g., "EGFR")
            variant: Variant notation (e.g., "G719S")
            tumor_type: Optional tumor type to filter results

        Returns:
            List of FDA-approved CGIBiomarker objects
        """
        all_biomarkers = self.fetch_biomarkers(gene, variant, tumor_type)
        return [
            b
            for b in all_biomarkers
            if b.is_fda_approved() and b.association.upper() == "RESPONSIVE"
        ]

    def fetch_biomarker_evidence(
        self, gene: str, variant: str, tumor_type: str | None = None
    ) -> list["CGIBiomarkerEvidence"]:
        """Fetch CGI biomarkers and convert to evidence model.

        This method fetches biomarkers and converts them to the
        CGIBiomarkerEvidence model for use in Insight.

        Args:
            gene: Gene symbol (e.g., "EGFR")
            variant: Variant notation (e.g., "G719S")
            tumor_type: Optional tumor type to filter results

        Returns:
            List of CGIBiomarkerEvidence objects
        """
        from oncomind.models.evidence.cgi import CGIBiomarkerEvidence

        biomarkers = self.fetch_biomarkers(gene, variant, tumor_type)
        evidence_list = []

        for biomarker in biomarkers:
            # Determine match specificity
            locus_match = self._determine_locus_match(biomarker.alteration, gene, variant)

            # Determine tumor match using centralized function
            # Args: source_disease (from biomarker), queried_tumor (user's query)
            tumor_matches = tumor_types_match(biomarker.tumor_type, tumor_type) if tumor_type else False

            # Build EvidenceLevel objects for consistency with other models
            from oncomind.models.evidence.base import EvidenceLevel
            locus_variant_match = EvidenceLevel(
                level=locus_match,
                scope="specific" if locus_match == "variant" else "unspecified",
                origin="kb",
            )
            cancer_type_match = None
            if tumor_type:
                cancer_type_match = EvidenceLevel(
                    level="cancer_specific" if tumor_matches else biomarker.tumor_type,
                    scope="specific" if tumor_matches else "unspecified",
                    origin="kb",
                )

            evidence_list.append(CGIBiomarkerEvidence(
                gene=biomarker.gene,
                alteration=biomarker.alteration,
                drug=biomarker.drug,
                drug_status=biomarker.drug_status,
                association=biomarker.association,
                evidence_level=biomarker.evidence_level,
                source=biomarker.source,
                tumor_type=biomarker.tumor_type,
                fda_approved=biomarker.is_fda_approved(),
                fda_url=biomarker.get_fda_url(),
                matched_alteration=biomarker.alteration,
                locus_variant_match=locus_variant_match,
                cancer_type_match=cancer_type_match,
            ))

        return evidence_list

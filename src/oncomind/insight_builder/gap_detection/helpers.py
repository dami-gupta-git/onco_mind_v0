"""Helper classes and utility functions for gap detection.

Contains MatchCounts dataclass and counting/scoring functions used
across gap detection modules.
"""

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from oncomind.models.extracted.evidence_gaps import (
    EvidenceGap,
    GapCategory,
    GapSeverity,
    CharacterizedAspect,
)
from oncomind.models.evidence.base import tumor_types_match, is_pan_cancer_term
from oncomind.config.constants import CADD_DELETERIOUS_THRESHOLD

if TYPE_CHECKING:
    from oncomind.models.evidence import Evidence


# =============================================================================
# MATCH COUNTING HELPERS
# =============================================================================


@dataclass
class MatchCounts:
    """Counts of evidence items by match level and tumor match.

    Reduces code duplication across gap detection functions that need to
    track variant/codon/gene match levels and tumor-specific counts.
    """

    total: int = 0
    variant: int = 0
    codon: int = 0
    gene: int = 0
    tumor: int = 0
    other_cancers: set = field(default_factory=set)  # Cancer types that didn't match

    @property
    def matches_on_str(self) -> str:
        """Build matches_on string (e.g., '2 variant, 1 codon, 3 gene').

        Only includes non-zero counts.
        """
        parts = []
        if self.variant > 0:
            parts.append(f"{self.variant} variant")
        if self.codon > 0:
            parts.append(f"{self.codon} codon")
        if self.gene > 0:
            parts.append(f"{self.gene} gene")
        return ", ".join(parts) if parts else "0 matches"

    @property
    def tumor_breakdown_str(self) -> str:
        """Build tumor breakdown string (e.g., '2 tumor, 4 other')."""
        other = self.total - self.tumor
        parts = [f"{self.tumor} tumor"]
        if other > 0:
            parts.append(f"{other} other")
        return ", ".join(parts)

    def add(self, other: "MatchCounts") -> "MatchCounts":
        """Add counts from another MatchCounts and return self (for chaining)."""
        self.total += other.total
        self.variant += other.variant
        self.codon += other.codon
        self.gene += other.gene
        self.tumor += other.tumor
        self.other_cancers.update(other.other_cancers)
        return self


def _get_locus_match(item) -> str:
    """Extract locus match level from an evidence item.

    Returns 'variant', 'codon', or 'gene'.
    """
    # Use locus_match computed property (returns string: 'variant', 'codon', or 'gene')
    if hasattr(item, "locus_match") and item.locus_match:
        return item.locus_match
    return "gene"


def count_with_levels(
    items: list, tumor_type: str | None = None, tumor_check_fn: callable = None
) -> MatchCounts:
    """Count items with match level breakdown and tumor matching.

    Args:
        items: List of evidence items (must have locus_match attribute or be passable to _get_locus_match)
        tumor_type: The tumor type to match against (optional)
        tumor_check_fn: Function(item) -> bool to check tumor match. If None, uses
            tumor_types_match with item.tumor_type or item.disease attribute.

    Returns:
        MatchCounts with all counts populated
    """
    counts = MatchCounts()

    for item in items:
        counts.total += 1
        level = _get_locus_match(item)
        if level == "variant":
            counts.variant += 1
        elif level == "codon":
            counts.codon += 1
        else:
            counts.gene += 1

        # Tumor matching
        if tumor_type:
            # Get the disease/tumor from the item
            tumor_attr = getattr(item, "tumor_type", None)
            disease_attr = getattr(item, "disease", None)
            tissue = (
                tumor_attr
                if isinstance(tumor_attr, str) and tumor_attr
                else (
                    disease_attr
                    if isinstance(disease_attr, str) and disease_attr
                    else None
                )
            )

            if tumor_check_fn:
                if tumor_check_fn(item):
                    counts.tumor += 1
                elif tissue and not is_pan_cancer_term(tissue):
                    counts.other_cancers.add(tissue)
            else:
                # Default: use the tumor_match property if available (single source of truth)
                # This ensures consistency with evidence tables in the UI
                item_tumor_match = getattr(item, "tumor_match", None)
                if item_tumor_match is True:
                    counts.tumor += 1
                elif (
                    item_tumor_match is False
                    and tissue
                    and not is_pan_cancer_term(tissue)
                ):
                    counts.other_cancers.add(tissue)
                elif item_tumor_match is None:
                    # Fallback for items without tumor_match property
                    # Args: source_disease (tissue from evidence), queried_tumor (user's query)
                    if tumor_types_match(tissue, tumor_type):
                        counts.tumor += 1
                    elif tissue and not is_pan_cancer_term(tissue):
                        counts.other_cancers.add(tissue)

    return counts


def count_fda_match_levels(
    filtered_fda: list,
    tumor_type: str | None = None,
) -> MatchCounts:
    """Count FDA evidence by variant match level and tumor match.

    Uses FDABiomarkerEvidence.variant_match_result for locus level and
    checks tumor_types for tumor matching.

    Args:
        filtered_fda: List of FDABiomarkerEvidence items (already filtered)
        tumor_type: The tumor type to match against (optional)

    Returns:
        MatchCounts with locus and tumor breakdown
    """
    counts = MatchCounts()

    for fda_ev in filtered_fda:
        counts.total += 1

        # Get locus match from variant_match_result (set by get_filtered_fda_evidence)
        # Values: "exact" -> variant, "codon" -> codon, "gene" -> gene
        match_result = fda_ev.variant_match_result
        if match_result == "exact":
            counts.variant += 1
        elif match_result == "codon":
            counts.codon += 1
        elif match_result == "gene":
            counts.gene += 1
        else:
            # Fallback based on specificity level
            from oncomind.models.evidence.fda_biomarker import SpecificityLevel

            if fda_ev.specificity == SpecificityLevel.VARIANT:
                counts.variant += 1
            elif fda_ev.specificity == SpecificityLevel.CODON:
                counts.codon += 1
            else:
                counts.gene += 1

        # Tumor matching - check if this FDA indication covers the queried tumor
        if tumor_type and fda_ev.tumor_types:
            # Check if any tumor type in the indication matches
            tumor_matched = any(
                tumor_types_match(t, tumor_type) or is_pan_cancer_term(t)
                for t in fda_ev.tumor_types
            )
            if tumor_matched:
                counts.tumor += 1
            else:
                # Track which cancers this drug IS approved for
                for t in fda_ev.tumor_types:
                    if not is_pan_cancer_term(t):
                        counts.other_cancers.add(t)
        elif not tumor_type:
            # No tumor specified - count all as tumor-matched
            counts.tumor += 1

    return counts


# =============================================================================
# PATHOGENICITY SIGNAL DETECTION
# =============================================================================


def has_pathogenic_signal(evidence: "Evidence") -> bool:
    """Check if variant has any signal suggesting pathogenicity.

    Used to avoid overcalling gaps on clearly benign variants (common polymorphisms).

    Returns True if any of:
    - AlphaMissense predicts pathogenic (P or likely_pathogenic)
    - CADD score >= 20 (predicted deleterious)
    - PolyPhen2 predicts damaging (D or probably_damaging)
    - SIFT predicts deleterious (D)
    - Has any clinical assertions or FDA approvals
    - Has any ClinVar pathogenic/likely pathogenic entries
    - Is a truncating variant (nonsense, frameshift)
    """
    func = evidence.functional

    # AlphaMissense pathogenic prediction
    if func.alphamissense_prediction and func.alphamissense_prediction.lower() in (
        "p",
        "pathogenic",
        "likely_pathogenic",
        "lp",
    ):
        return True

    # CADD score >= 20 suggests deleteriousness
    if func.cadd_score is not None and func.cadd_score >= 20:
        return True

    # PolyPhen2 damaging prediction
    if func.polyphen2_prediction and func.polyphen2_prediction.lower() in (
        "d",
        "damaging",
        "probably_damaging",
        "possibly_damaging",
    ):
        return True

    # SIFT deleterious prediction (D = deleterious, T = tolerated)
    if func.sift_prediction and func.sift_prediction.lower() in ("d", "deleterious"):
        return True

    # Has clinical evidence (strongest signal)
    if evidence.civic_assertions or evidence.fda_biomarker_evidence:
        return True

    # ClinVar pathogenic entries
    for entry in evidence.clinvar_entries:
        if (
            entry.clinical_significance
            and "pathogenic" in entry.clinical_significance.lower()
        ):
            return True

    # Check if overall ClinVar significance is pathogenic
    if (
        evidence.clinvar_significance
        and "pathogenic" in evidence.clinvar_significance.lower()
    ):
        return True

    # Truncating variants (nonsense, frameshift) are generally pathogenic
    if func.snpeff_effect:
        effect_lower = func.snpeff_effect.lower()
        if any(
            term in effect_lower
            for term in [
                "stop_gained",
                "frameshift",
                "splice_donor",
                "splice_acceptor",
                "start_lost",
                "nonsense",
            ]
        ):
            return True

    return False


# =============================================================================
# SOURCE NORMALIZATION
# =============================================================================


def normalize_source(source: str) -> str:
    """Normalize source names to detect duplicates.

    VICC/civic is essentially the same as CIViC, VICC/oncokb is OncoKB, etc.
    """
    source_lower = source.lower()
    if source_lower in ("civic", "vicc/civic"):
        return "CIViC"
    if source_lower in ("oncokb", "vicc/oncokb"):
        return "OncoKB"
    if source_lower in ("cgi", "vicc/cgi"):
        return "CGI"
    if source_lower in ("molecularmatch", "vicc/molecularmatch"):
        return "MolecularMatch"
    if source_lower.startswith("vicc/"):
        return source[5:].title()  # Strip "vicc/" prefix
    return source


# =============================================================================
# CATEGORY ORDERING AND SORTING
# =============================================================================

# Define category ordering for sorting well-characterized aspects
CATEGORY_ORDER: dict[GapCategory | None, int] = {
    GapCategory.CLINICAL: 0,
    GapCategory.TUMOR_TYPE: 1,
    GapCategory.DRUG_RESPONSE: 2,
    GapCategory.RESISTANCE: 3,
    GapCategory.FUNCTIONAL: 4,
    GapCategory.PRECLINICAL: 5,
    GapCategory.PREVALENCE: 6,
    GapCategory.VALIDATION: 7,
    GapCategory.PROGNOSTIC: 8,
    GapCategory.DISCORDANT: 9,
    None: 99,  # Uncategorized items go last
}


def sort_characterized_by_category(
    items: list[CharacterizedAspect],
) -> list[CharacterizedAspect]:
    """Sort well-characterized aspects by category for grouped display.

    Args:
        items: List of CharacterizedAspect items

    Returns:
        Sorted list with items grouped by category
    """
    return sorted(items, key=lambda x: CATEGORY_ORDER.get(x.category, 99))


# =============================================================================
# SCORING FUNCTIONS
# =============================================================================

# Research-oriented gap weights: biological unknowns weighted higher than clinical gaps
GAP_CATEGORY_WEIGHTS: dict[GapCategory, float] = {
    GapCategory.VALIDATION: 3.5,  # Strong signal + no validation = high research value
    GapCategory.FUNCTIONAL: 3.0,  # Mechanism unknown
    GapCategory.PRECLINICAL: 2.5,  # No models to test hypotheses
    GapCategory.RESISTANCE: 2.0,  # Resistance mechanisms unknown
    GapCategory.DISCORDANT: 2.0,  # Conflicting evidence needs resolution
    GapCategory.DRUG_RESPONSE: 1.5,  # Drug sensitivity unknown
    GapCategory.TUMOR_TYPE: 1.5,  # Not studied in this tumor
    GapCategory.PREVALENCE: 1.0,  # Epidemiology unknown
    GapCategory.CLINICAL: 1.0,  # Lower weight for research context
    GapCategory.PROGNOSTIC: 1.0,  # Prognostic impact unknown
}

SEVERITY_MULTIPLIERS: dict[GapSeverity, float] = {
    GapSeverity.CRITICAL: 4.0,
    GapSeverity.HIGH: 3.0,
    GapSeverity.SIGNIFICANT: 2.0,
    GapSeverity.MODERATE: 1.5,
    GapSeverity.MINOR: 1.0,
    GapSeverity.INFORMATIONAL: 0.5,
}


def compute_overall_quality(
    gaps: list[EvidenceGap],
    well_characterized_count: int,
) -> str:
    """Compute overall evidence quality using net scoring (gaps vs well-characterized).

    A variant with many well-characterized aspects and few gaps scores better than
    one with few well-characterized aspects and the same gaps.

    FLOOR RULES (severity-based caps):
    - Any CRITICAL gap -> cannot be better than "limited"
    - Any HIGH gap -> cannot be better than "moderate"
    - Any SIGNIFICANT gap -> cannot be better than "moderate"

    Args:
        gaps: List of evidence gaps found
        well_characterized_count: Number of well-characterized aspects

    Returns:
        Quality rating: "comprehensive" | "moderate" | "limited" | "minimal"
    """
    # Check for severity-based floor caps FIRST
    has_critical = any(g.severity == GapSeverity.CRITICAL for g in gaps)
    has_high = any(g.severity == GapSeverity.HIGH for g in gaps)
    has_significant = any(g.severity == GapSeverity.SIGNIFICANT for g in gaps)

    # Calculate gap penalty score
    gap_score = 0.0
    for gap in gaps:
        category_weight = GAP_CATEGORY_WEIGHTS.get(gap.category, 1.0)
        severity_mult = SEVERITY_MULTIPLIERS.get(gap.severity, 1.0)
        gap_score += category_weight * severity_mult

    # Give credit for well-characterized aspects (each worth 1.5 points of offset)
    positive_credit = well_characterized_count * 1.5

    # Net score: higher gap_score is worse, positive_credit offsets it
    net_score = gap_score - positive_credit

    # Compute base quality from net score
    if net_score >= 12.0:
        base_quality = "minimal"
    elif net_score >= 6.0:
        base_quality = "limited"
    elif net_score >= 0.0:
        base_quality = "moderate"
    else:
        base_quality = "comprehensive"

    # Apply floor caps: high-severity gaps prevent "comprehensive" classification
    # Order matters: quality levels are minimal < limited < moderate < comprehensive
    quality_order = ["minimal", "limited", "moderate", "comprehensive"]
    base_idx = quality_order.index(base_quality)

    if has_critical:
        # CRITICAL gaps -> cannot be better than "limited" (index 1)
        max_idx = 1
    elif has_high or has_significant:
        # HIGH or SIGNIFICANT gaps -> cannot be better than "moderate" (index 2)
        max_idx = 2
    else:
        max_idx = 3  # No cap

    final_idx = min(base_idx, max_idx)
    return quality_order[final_idx]


def compute_research_priority(
    evidence: "Evidence",
    gaps: list[EvidenceGap],
    overall_quality: str,
    is_cancer_gene: bool,
    has_pathogenic: bool,
) -> str:
    """Compute research priority based on gene importance and gap profile.

    Args:
        evidence: The aggregated evidence
        gaps: List of identified gaps
        overall_quality: The computed overall quality ("comprehensive", "moderate", etc.)
        is_cancer_gene: Whether the gene is a known cancer gene (from ctx.is_cancer_gene)
        has_pathogenic: Whether evidence shows pathogenic signals (from ctx.has_pathogenic_signal)

    Returns: "very_high" | "high" | "medium" | "low"
    """
    from oncomind.models.gene_context import is_hotspot_adjacent
    from oncomind.config.constants import HOTSPOT_ADJACENCY_WINDOW

    # Count gaps first
    critical_count = sum(1 for g in gaps if g.severity == GapSeverity.CRITICAL)
    high_count = sum(1 for g in gaps if g.severity == GapSeverity.HIGH)
    significant_count = sum(1 for g in gaps if g.severity == GapSeverity.SIGNIFICANT)

    # Only return low if comprehensive AND no significant/critical/high gaps
    if (
        overall_quality == "comprehensive"
        and critical_count == 0
        and high_count == 0
        and significant_count == 0
    ):
        return "low"

    gene = evidence.identifiers.gene
    variant = evidence.identifiers.variant

    has_strong_oncogenic_signal = (
        has_pathogenic
        and evidence.depmap_evidence is not None
        and evidence.depmap_evidence.is_essential()
    )

    has_biological_gaps = any(
        g.category
        in (GapCategory.VALIDATION, GapCategory.FUNCTIONAL, GapCategory.PRECLINICAL)
        for g in gaps
    )

    is_adjacent, _ = is_hotspot_adjacent(gene, variant, window=HOTSPOT_ADJACENCY_WINDOW)

    # Very high: strong oncogenic signal + biological gaps = prime research target
    if has_strong_oncogenic_signal and has_biological_gaps:
        return "very_high"

    # Very high: hotspot-adjacent variant with pathogenic signal
    if is_adjacent and has_pathogenic and has_biological_gaps:
        return "very_high"

    # High: cancer gene with critical gaps
    if is_cancer_gene and critical_count > 0:
        return "high"

    # High: hotspot-adjacent variant in cancer gene
    if is_adjacent and is_cancer_gene:
        return "high"

    # High: any HIGH severity gaps (conflicting evidence)
    if high_count > 0:
        return "high"

    # Medium: any critical gaps OR cancer gene with significant gaps
    if critical_count > 0 or (is_cancer_gene and significant_count > 0):
        return "medium"

    return "low"


# =============================================================================
# CONTEXT-AWARE ENRICHMENT HELPERS
# =============================================================================


def get_primary_drug(evidence: "Evidence") -> str | None:
    """Get the primary approved drug for this variant."""
    # Use FDA biomarker evidence (only REQUIRED_POSITIVE)
    from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement

    for ev in evidence.fda_biomarker_evidence:
        if ev.requirement != BiomarkerRequirement.REQUIRED_NEGATIVE and ev.drug_name:
            return ev.brand_name or ev.drug_name
    # Fall back to CIViC assertion therapies
    for assertion in evidence.civic_assertions:
        if assertion.therapies:
            return (
                assertion.therapies[0]
                if isinstance(assertion.therapies, list)
                else assertion.therapies
            )
    return None


def get_top_sensitive_drugs(
    evidence: "Evidence", tumor_type: str | None = None
) -> list[str]:
    """Get top sensitive drugs from DepMap data.

    Only returns drugs if there are tumor-matched cell lines with the mutation,
    since sensitivity data from unrelated tumor types is not applicable.
    """
    if not evidence.depmap_evidence or not evidence.depmap_evidence.drug_sensitivities:
        return []

    # Check for tumor-matched cell lines
    if tumor_type and evidence.depmap_evidence.cell_line_models:
        mutant_models = [
            cl for cl in evidence.depmap_evidence.cell_line_models if cl.has_mutation
        ]
        # Args: source_disease (from model), queried_tumor (user's query)
        tumor_matched = [
            m
            for m in mutant_models
            if m.primary_disease and tumor_types_match(m.primary_disease, tumor_type)
        ]
        if not tumor_matched:
            # No tumor-matched cell lines - don't suggest drugs from unrelated tumor types
            return []

    top_drugs = evidence.depmap_evidence.get_top_sensitive_drugs(5)
    return [ds.drug_name for ds in top_drugs]


def get_top_cooccurring_gene(evidence: "Evidence") -> str | None:
    """Get the top co-occurring gene from cBioPortal data."""
    if evidence.cbioportal_evidence and evidence.cbioportal_evidence.co_occurring:
        return evidence.cbioportal_evidence.co_occurring[0].gene
    return None


def has_strong_cooccurrence(evidence: "Evidence", threshold_pct: float = None) -> bool:
    """Check if there's a strong co-occurrence signal (>threshold% co-mutation rate)."""
    if threshold_pct is None:
        from oncomind.config.constants import COOCCURRENCE_STRONG_THRESHOLD_PCT

        threshold_pct = COOCCURRENCE_STRONG_THRESHOLD_PCT

    if evidence.cbioportal_evidence and evidence.cbioportal_evidence.co_occurring:
        top_cooc = evidence.cbioportal_evidence.co_occurring[0]
        return top_cooc.pct is not None and top_cooc.pct >= threshold_pct
    return False


def get_top_codependencies(evidence: "Evidence") -> list[str]:
    """Get top co-dependent genes from DepMap data."""
    if evidence.depmap_evidence and evidence.depmap_evidence.co_dependencies:
        return [cd.gene for cd in evidence.depmap_evidence.co_dependencies[:3]]
    return []

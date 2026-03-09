"""Individual gap detection check functions.

Each function checks for a specific type of evidence gap and updates
the GapDetectionContext with findings.
"""

from typing import TYPE_CHECKING

from oncomind.models.extracted.evidence_gaps import GapCategory, GapSeverity
from oncomind.models.evidence.tumor_evidence import TumorEvidenceMatch
from oncomind.models.gene_context import (
    is_hotspot_variant,
    is_hotspot_adjacent,
    _extract_codon_position,
)
from oncomind.config.constants import (
    GNOMAD_COMMON_AF_THRESHOLD,
    LITERATURE_WELL_CHARACTERIZED_THRESHOLD,
    HOTSPOT_ADJACENCY_WINDOW,
    MAX_LITERATURE_RESULTS,
    CADD_DELETERIOUS_THRESHOLD,
)
from oncomind.models.evidence.base import tumor_types_match, is_pan_cancer_term

from .helpers import (
    MatchCounts,
    count_with_levels,
    count_fda_match_levels,
    _get_locus_match,
)

if TYPE_CHECKING:
    from oncomind.models.evidence import Evidence
    from .context import GapDetectionContext


# =============================================================================
# HOTSPOT CONTEXT
# =============================================================================


def check_hotspot_context(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check if variant is at or near a known cancer hotspot."""
    is_hotspot = is_hotspot_variant(ctx.gene, ctx.variant)
    is_adjacent, nearest_hotspot = is_hotspot_adjacent(
        ctx.gene, ctx.variant, window=HOTSPOT_ADJACENCY_WINDOW
    )

    if is_hotspot:
        ctx.add_well_characterized(
            "known cancer hotspot",
            f"Codon {_extract_codon_position(ctx.variant)} is in cancerhotspots.org",
            category=GapCategory.FUNCTIONAL,
            matches_on="codon",
        )
    elif is_adjacent and nearest_hotspot:
        # Rare variant near a hotspot - high research value
        ctx.add_well_characterized(
            f"near hotspot codon {nearest_hotspot}",
            f"Within {HOTSPOT_ADJACENCY_WINDOW} codons of known hotspot — structural similarity likely",
            category=GapCategory.FUNCTIONAL,
        )
        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Variant near known hotspot (codon {nearest_hotspot}) — functional characterization needed",
            suggested_studies=[
                f"Compare to nearby hotspot {ctx.gene} codon {nearest_hotspot}",
                "Structural modeling to assess activation mechanism",
                "Functional assay (transformation, signaling)",
            ],
            addressable_with=[
                "AlphaFold",
                "Literature on nearby hotspot",
                "Isogenic models",
            ],
        )
        ctx.add_poorly_characterized("rare-near-hotspot variant function")


# =============================================================================
# FUNCTIONAL PREDICTIONS
# =============================================================================


def check_functional_predictions(
    evidence: "Evidence", ctx: "GapDetectionContext"
) -> None:
    """Check for computational pathogenicity predictions."""
    has_functional = (
        evidence.functional.alphamissense_score is not None
        or evidence.functional.cadd_score is not None
        or evidence.functional.polyphen2_prediction is not None
        or evidence.functional.sift_prediction is not None
    )

    if has_functional:
        func_sources = []
        if evidence.functional.alphamissense_score is not None:
            func_sources.append(
                f"AlphaMissense={evidence.functional.alphamissense_score:.2f}"
            )
        if evidence.functional.cadd_score is not None:
            func_sources.append(f"CADD={evidence.functional.cadd_score:.1f}")
        if evidence.functional.polyphen2_prediction:
            func_sources.append(f"PolyPhen2={evidence.functional.polyphen2_prediction}")
        if evidence.functional.sift_prediction:
            func_sources.append(f"SIFT={evidence.functional.sift_prediction}")
        ctx.add_well_characterized(
            "computational pathogenicity",
            " | ".join(func_sources) if func_sources else "Predictions available",
            category=GapCategory.FUNCTIONAL,
            matches_on="variant",
        )
    else:
        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"No computational pathogenicity predictions for {ctx.gene} {ctx.variant}",
            suggested_studies=["Run AlphaMissense, CADD, PolyPhen2, SIFT"],
            addressable_with=["MyVariant.info", "VEP"],
        )
        ctx.add_poorly_characterized("pathogenicity predictions")

    # Check for conflicting functional predictions
    # PolyPhen2 has known limitations for somatic cancer driver mutations, especially RAS/RAF hotspots
    # Flag when PolyPhen2 says Benign but other predictors indicate pathogenicity
    check_conflicting_predictions(evidence, ctx)

    # Check for gnomAD population frequency data (informational, not a penalty)
    # If gnomAD AF > 0.01% (0.0001), note that the variant is observed in the general population
    gnomad_af = (
        evidence.functional.gnomad_exome_af or evidence.functional.gnomad_genome_af
    )
    if gnomad_af is not None and gnomad_af > GNOMAD_COMMON_AF_THRESHOLD:
        af_pct = gnomad_af * 100
        ctx.add_well_characterized(
            "population frequency",
            f"Observed in general population (gnomAD AF: {af_pct:.3f}%)",
            category=GapCategory.PREVALENCE,
        )


def check_conflicting_predictions(
    evidence: "Evidence", ctx: "GapDetectionContext"
) -> None:
    """Check for conflicting functional predictions, especially PolyPhen2 vs others.

    PolyPhen2 was trained primarily on germline disease variants and has known
    limitations for somatic cancer driver mutations, particularly RAS/RAF hotspots.
    When PolyPhen2 predicts Benign but other tools (AlphaMissense, CADD, SIFT)
    predict deleterious/pathogenic, flag this discordance.

    This is particularly important for known cancer hotspots where PolyPhen2's
    benign call could mislead clinical interpretation.
    """
    func = evidence.functional
    polyphen2_pred = func.polyphen2_prediction

    # Only check if PolyPhen2 predicts Benign (B)
    if not polyphen2_pred or polyphen2_pred.upper() != "B":
        return

    # Check if other predictors indicate pathogenicity
    alphamissense_pathogenic = (
        func.alphamissense_prediction is not None
        and func.alphamissense_prediction.upper() in ("P", "PATHOGENIC")
    )
    cadd_deleterious = (
        func.cadd_score is not None and func.cadd_score > CADD_DELETERIOUS_THRESHOLD
    )
    sift_deleterious = (
        func.sift_prediction is not None
        and func.sift_prediction.upper() in ("D", "DELETERIOUS")
    )

    # Count how many other predictors disagree with PolyPhen2
    disagreeing_predictors = []
    if alphamissense_pathogenic:
        score = func.alphamissense_score
        disagreeing_predictors.append(
            f"AlphaMissense={score:.2f}" if score else "AlphaMissense=P"
        )
    if cadd_deleterious:
        disagreeing_predictors.append(f"CADD={func.cadd_score:.1f}")
    if sift_deleterious:
        disagreeing_predictors.append("SIFT=D")

    # Only flag if at least 2 other predictors disagree, or if it's a known hotspot
    is_hotspot = evidence.hotspots_evidence and evidence.hotspots_evidence.is_hotspot
    min_disagreements = 1 if is_hotspot else 2

    if len(disagreeing_predictors) >= min_disagreements:
        pp2_score = func.polyphen2_score
        pp2_str = f"PolyPhen2={pp2_score:.2f}" if pp2_score else "PolyPhen2=B"

        # Build description based on whether it's a hotspot
        if is_hotspot:
            description = (
                f"PolyPhen2 predicts Benign ({pp2_str}) but {', '.join(disagreeing_predictors)} "
                f"indicate pathogenicity — PolyPhen2 has known limitations for cancer hotspot mutations"
            )
        else:
            description = (
                f"Conflicting predictions: {pp2_str} vs {', '.join(disagreeing_predictors)} — "
                f"consider functional validation"
            )

        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.MODERATE,  # Conflicting data in specific context
            description=description,
            suggested_studies=[
                "Prioritize AlphaMissense/CADD over PolyPhen2 for cancer variants",
                "Check ClinVar/COSMIC for clinical evidence",
                "Consider functional assay if interpretation is critical",
            ],
            addressable_with=["Literature review", "ClinVar", "Functional assay"],
        )


# =============================================================================
# GENE MECHANISM
# =============================================================================


def check_gene_mechanism(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for gene function and essentiality data."""
    has_mechanism = bool(evidence.context.gene_role) or bool(evidence.context.pathway)
    has_depmap_essentiality = (
        evidence.depmap_evidence is not None
        and evidence.depmap_evidence.gene_dependency is not None
    )

    if has_mechanism:
        role = evidence.context.gene_role or "unknown"
        pathway = evidence.context.pathway or ""
        basis = f"Role: {role}" + (f", Pathway: {pathway}" if pathway else "")
        ctx.add_well_characterized(
            "gene function", basis, category=GapCategory.FUNCTIONAL, matches_on="gene"
        )

    # Only show gene essentiality if the gene IS essential (CERES < -0.5)
    # Gene essentiality is pan-cancer data, so don't show non-essential scores
    if has_depmap_essentiality and evidence.depmap_evidence.is_essential():
        dep = evidence.depmap_evidence.gene_dependency
        score = dep.mean_dependency_score if dep else None
        pct = dep.dependency_pct if dep else 0
        score_str = f"CERES={score:.2f}, {pct:.0f}% of cell lines depend on it"
        ctx.add_well_characterized(
            "gene essentiality",
            score_str,
            category=GapCategory.FUNCTIONAL,
            matches_on="gene",
        )

    if not has_mechanism and not has_depmap_essentiality:
        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Functional impact of {ctx.variant} on {ctx.gene} protein is unknown",
            suggested_studies=[
                "Functional assay",
                "Structural modeling",
                "Cell-based reporter",
            ],
            addressable_with=["UniProt", "Literature search", "DepMap"],
        )
        ctx.add_poorly_characterized("functional mechanism")


# =============================================================================
# CLINICAL EVIDENCE
# =============================================================================


def check_clinical_evidence(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for FDA-approved therapies using fda_biomarker_evidence.

    Calls get_filtered_fda_evidence() to filter FDA drugs by gene/variant/tumor match.
    Also checks for CIViC Level A assertions (FDA-approved/guideline-included) and
    CGI biomarkers with fda_approved=True.

    Note: VICC MetaKB is intentionally NOT used for clinical actionability.
    VICC aggregates data from CIViC, CGI, OncoKB, and other sources. Using VICC
    would double-count approvals already counted via civic_assertions and cgi_biomarkers.
    The three sources used here (FDA labels, CIViC Level A, CGI FDA-approved) are
    independent and non-overlapping.
    """
    # Filter FDA evidence by gene/variant/tumor to get matching drugs
    filtered_fda = evidence.get_filtered_fda_evidence(
        queried_gene=ctx.gene,
        queried_variant=ctx.variant,
        queried_tumor=ctx.tumor_type,
    )
    fda_count = len(filtered_fda)

    # Count CIViC Level A assertions (FDA-approved companion diagnostic or guideline-included)
    # Only count Assertions (AIDs) with amp_level_letter == "A", not Evidence Items (EIDs)
    civic_assertions_level_a = [
        a
        for a in evidence.civic_assertions
        if a.amp_level_letter and a.amp_level_letter.upper() == "A"
    ]
    civic_level_a_count = len(civic_assertions_level_a)

    # Count CGI biomarkers with fda_approved=True
    cgi_fda_approved = [b for b in evidence.cgi_biomarkers if b.fda_approved]
    cgi_approved_count = len(cgi_fda_approved)

    ctx.has_clinical = (
        bool(evidence.civic_assertions)
        or bool(evidence.civic_evidence)
        or fda_count > 0
    )

    # Build clinical actionability basis from FDA, CIViC Level A, and CGI FDA-approved
    # No locus/tumor match icons for clinical actionability - it's about approvals, not match specificity
    actionability_parts = []
    if fda_count > 0:
        actionability_parts.append(
            f"{fda_count} FDA-approved indication{'s' if fda_count > 1 else ''}"
        )
    if civic_level_a_count > 0:
        actionability_parts.append(
            f"{civic_level_a_count} approval{'s' if civic_level_a_count > 1 else ''} from CIViC assertions"
        )
    if cgi_approved_count > 0:
        actionability_parts.append(
            f"{cgi_approved_count} approval{'s' if cgi_approved_count > 1 else ''} from CGI"
        )

    if actionability_parts:
        basis = ", ".join(actionability_parts)
        ctx.add_well_characterized(
            "clinical actionability",
            basis,
            category=GapCategory.CLINICAL,
            # No matches_on or tumor_match - clinical actionability doesn't show these icons
        )

    if fda_count == 0 and civic_level_a_count == 0:
        if evidence.fda_biomarker_evidence:
            # FDA drugs exist but don't match this variant/tumor (filtered count is 0 but raw list has items)
            ctx.add_gap(
                category=GapCategory.CLINICAL,
                severity=GapSeverity.SIGNIFICANT,
                description=f"FDA-approved therapies for {ctx.gene} exist but don't match {ctx.variant}",
                suggested_studies=["Basket trial", "Off-label use case series"],
                addressable_with=["ClinicalTrials.gov", "FDA label expansion studies"],
            )
        elif not ctx.has_clinical:
            # No clinical evidence at all (no FDA, CIViC assertions, or CIViC evidence)
            ctx.add_gap(
                category=GapCategory.CLINICAL,
                severity=GapSeverity.CRITICAL,
                description=f"No curated clinical evidence for {ctx.gene} {ctx.variant}",
                suggested_studies=[
                    "Case series",
                    "Retrospective cohort",
                    "Basket trial inclusion",
                ],
                addressable_with=["CIViC submission", "Literature curation"],
            )
            ctx.add_poorly_characterized("clinical evidence")


# =============================================================================
# TUMOR TYPE EVIDENCE
# =============================================================================


def check_tumor_type_evidence(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for tumor-type-specific evidence."""
    if not ctx.tumor_type:
        return

    tumor_match = check_tumor_specific_evidence(evidence, ctx.tumor_type)

    if tumor_match.has_tumor_evidence:
        # Only mark as well-characterized if we have VARIANT-level tumor-matched evidence
        if tumor_match.total_variant_level > 0:
            ctx.add_well_characterized(
                f"evidence items for {ctx.tumor_type}",
                tumor_match.sources_str,
                category=GapCategory.TUMOR_TYPE,
                matches_on=tumor_match.matches_on_str,
                tumor_match=tumor_match.tumor_count_str,
            )
        else:
            # Tumor evidence exists but only at gene/codon level — MODERATE gap
            if tumor_match.total_gene_level > 0 and tumor_match.total_codon_level == 0:
                ctx.add_gap(
                    category=GapCategory.TUMOR_TYPE,
                    severity=GapSeverity.MODERATE,
                    description=f"Evidence in {ctx.tumor_type} exists for {ctx.gene} but not specifically for {ctx.variant} (gene-level only)",
                    suggested_studies=[
                        f"Variant-specific case series in {ctx.tumor_type}",
                        "Retrospective analysis of variant outcomes",
                    ],
                    addressable_with=[
                        "ClinicalTrials.gov",
                        "Real-world evidence databases",
                    ],
                )
                ctx.add_poorly_characterized(f"variant-specific {ctx.tumor_type} data")
            elif tumor_match.total_codon_level > 0:
                from oncomind.models.evidence.base import extract_variant_position

                codon_pos = extract_variant_position(ctx.variant) or ""
                ctx.add_gap(
                    category=GapCategory.TUMOR_TYPE,
                    severity=GapSeverity.MODERATE,
                    description=f"Evidence in {ctx.tumor_type} exists for {ctx.gene} codon {codon_pos} but not specifically for {ctx.variant}",
                    suggested_studies=[
                        "Variant-specific response comparison",
                        f"Case series of {ctx.variant} in {ctx.tumor_type}",
                    ],
                    addressable_with=[
                        "Published case reports",
                        "Basket trial subgroup analyses",
                    ],
                )
                ctx.add_poorly_characterized(f"variant-specific {ctx.tumor_type} data")
    else:
        # No tumor-specific evidence at all — severity depends on gene importance
        if ctx.is_cancer_gene and not ctx.has_clinical and ctx.has_pathogenic_signal:
            severity = GapSeverity.CRITICAL
        elif ctx.is_cancer_gene or ctx.has_pathogenic_signal:
            severity = GapSeverity.SIGNIFICANT
        else:
            severity = GapSeverity.MINOR

        ctx.add_gap(
            category=GapCategory.TUMOR_TYPE,
            severity=severity,
            description=f"No evidence specific to {ctx.tumor_type} for {ctx.gene} {ctx.variant}",
            suggested_studies=[
                f"Case series in {ctx.tumor_type}",
                f"Retrospective analysis in {ctx.tumor_type} cohort",
                "Basket trial with histology-specific cohort",
            ],
            addressable_with=["ClinicalTrials.gov", "Literature search"],
        )
        ctx.add_poorly_characterized(f"{ctx.tumor_type}-specific data")


def check_tumor_specific_evidence(
    evidence: "Evidence", tumor_type: str
) -> TumorEvidenceMatch:
    """Check which evidence sources have tumor-specific data.

    Returns a TumorEvidenceMatch with per-source breakdown of:
    - How many items match the tumor type
    - Match level (variant/codon/gene) for each item

    Args:
        evidence: The aggregated evidence
        tumor_type: The queried tumor type (e.g., "Lung Cancer")

    Returns:
        TumorEvidenceMatch with aggregated match info
    """
    result = TumorEvidenceMatch(tumor_type=tumor_type)

    def count_matches(items, tumor_check_fn) -> tuple[int, int, int, int]:
        """Count items matching tumor and their match levels.

        Returns (count, variant_matches, codon_matches, gene_matches)
        """
        count = variant = codon = gene = 0
        for item in items:
            if tumor_check_fn(item):
                count += 1
                level = _get_locus_match(item)
                if level == "variant":
                    variant += 1
                elif level == "codon":
                    codon += 1
                else:
                    gene += 1
        return count, variant, codon, gene

    # CIViC assertions - use the tumor_match property already computed by API client
    counts = count_matches(
        evidence.civic_assertions, lambda a: a.disease and a.tumor_match
    )
    if counts[0] > 0:
        result.add_source_match("CIViC Assertions", *counts)

    # CIViC evidence
    counts = count_matches(
        evidence.civic_evidence, lambda c: c.disease and c.tumor_match
    )
    if counts[0] > 0:
        result.add_source_match("CIViC", *counts)

    # FDA biomarker evidence - filter for REQUIRED_POSITIVE only AND tumor match
    from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement

    fda_filtered = [
        ev
        for ev in evidence.fda_biomarker_evidence
        if ev.requirement != BiomarkerRequirement.REQUIRED_NEGATIVE
    ]
    # Only count FDA evidence that actually matches the queried tumor
    counts = count_matches(
        fda_filtered,
        lambda f: f.tumor_types
        and any(
            tumor_types_match(t, tumor_type) or is_pan_cancer_term(t)
            for t in f.tumor_types
        ),
    )
    if counts[0] > 0:
        result.add_source_match("FDA", *counts)

    # VICC evidence - use get_vicc_unique() to exclude CIViC/CGI sources (avoid double-counting)
    vicc_unique = evidence.get_vicc_unique()
    counts = count_matches(vicc_unique, lambda v: v.disease and v.tumor_match)
    if counts[0] > 0:
        result.add_source_match("VICC (meta-KB)", *counts)

    # CGI biomarkers (all tiers combined) - use tumor_match property
    all_cgi = (
        list(evidence.cgi_biomarkers)
        + list(evidence.preclinical_biomarkers)
        + list(evidence.early_phase_biomarkers)
    )
    counts = count_matches(all_cgi, lambda c: c.tumor_type and c.tumor_match)
    if counts[0] > 0:
        result.add_source_match("CGI", *counts)

    return result


# =============================================================================
# DRUG RESPONSE
# =============================================================================


def check_drug_response(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for FDA-approved therapies and drug response data.

    This function checks two separate things:
    1. FDA-approved therapies - filtered by gene/variant/tumor match
    2. Drug response data - from VICC, CGI (for general "well characterized" tracking)

    A gap is added if there's no FDA-approved therapy for this variant.
    """
    # =========================================================================
    # 1. Check for FDA-approved therapies
    # =========================================================================
    # Filter FDA evidence by gene/variant/tumor to get matching drugs
    filtered_fda = evidence.get_filtered_fda_evidence(
        queried_gene=ctx.gene,
        queried_variant=ctx.variant,
        queried_tumor=ctx.tumor_type,
    )
    fda_count = len(filtered_fda)

    if fda_count > 0:
        # Compute match breakdown for FDA evidence
        fda_match_counts = count_fda_match_levels(filtered_fda, ctx.tumor_type)

        ctx.add_well_characterized(
            "FDA-approved therapy",
            f"{fda_count} FDA approval{'s' if fda_count > 1 else ''}",
            category=GapCategory.DRUG_RESPONSE,
            matches_on=fda_match_counts.matches_on_str,
            tumor_match=(
                fda_match_counts.tumor_breakdown_str if ctx.tumor_type else None
            ),
        )
    elif not evidence.fda_biomarker_evidence:
        # No FDA-approved therapy at all - this is a significant gap
        ctx.add_gap(
            category=GapCategory.DRUG_RESPONSE,
            severity=GapSeverity.SIGNIFICANT,
            description=f"No FDA-approved therapy for {ctx.gene} {ctx.variant}",
            suggested_studies=["Clinical trial enrollment", "Off-label use evaluation"],
            addressable_with=["ClinicalTrials.gov", "NCCN guidelines", "Basket trials"],
        )
    else:
        # FDA drugs exist for the gene but none match this variant/tumor
        ctx.add_gap(
            category=GapCategory.DRUG_RESPONSE,
            severity=GapSeverity.SIGNIFICANT,
            description=f"No FDA-approved therapy for {ctx.gene} {ctx.variant}",
            suggested_studies=["Basket trial", "Off-label use case series"],
            addressable_with=["ClinicalTrials.gov", "FDA label expansion studies"],
        )

    # =========================================================================
    # 2. Set context flag for drug data (VICC, CGI used for context, not displayed)
    # =========================================================================
    # Note: VICC/CGI drug response data is shown in the Therapies tab, not Gap Analysis
    # Gap Analysis focuses only on FDA-approved therapies
    cgi_counts = count_with_levels(evidence.cgi_biomarkers, ctx.tumor_type)
    vicc_counts = count_with_levels(evidence.vicc_evidence, ctx.tumor_type)
    drug_response_counts = MatchCounts().add(cgi_counts).add(vicc_counts)

    # Set context flag for any drug data (used by other gap checks)
    ctx.has_drug_data = drug_response_counts.total > 0 or fda_count > 0


def check_preclinical_biomarkers(
    evidence: "Evidence", ctx: "GapDetectionContext"
) -> None:
    """Check for preclinical and early phase biomarker data.

    These are lower evidence tier than FDA-approved biomarkers:
    - Preclinical: cell line studies, xenografts
    - Early phase: Phase I/II trials, case reports
    """
    preclin_counts = count_with_levels(evidence.preclinical_biomarkers, ctx.tumor_type)
    early_counts = count_with_levels(evidence.early_phase_biomarkers, ctx.tumor_type)

    if not preclin_counts.total and not early_counts.total:
        return

    # Aggregate counts
    counts = MatchCounts().add(preclin_counts).add(early_counts)

    # Build source string
    sources = []
    if preclin_counts.total:
        sources.append(f"{preclin_counts.total} preclinical")
    if early_counts.total:
        sources.append(f"{early_counts.total} early phase")

    ctx.add_well_characterized(
        "preclinical/early phase biomarkers",
        " + ".join(sources),
        category=GapCategory.PRECLINICAL,
        matches_on=counts.matches_on_str,
        tumor_match=counts.tumor_breakdown_str,
    )

    # Add gap if evidence exists only in other tumors (not tumor-matched)
    if ctx.tumor_type and counts.tumor == 0 and counts.other_cancers:
        other_cancers_str = ", ".join(sorted(counts.other_cancers)[:3])
        ctx.add_gap(
            category=GapCategory.PRECLINICAL,
            severity=GapSeverity.MODERATE,  # Data exists, just in different tumor context
            description=f"Preclinical biomarker data for {ctx.gene} {ctx.variant} exists only in other cancers ({other_cancers_str}), not {ctx.tumor_type}",
            suggested_studies=["Tumor-specific cell line studies", "PDX models"],
            addressable_with=["DepMap", "GDSC", "Literature search"],
        )


def check_depmap_drug_sensitivity(
    evidence: "Evidence", ctx: "GapDetectionContext"
) -> None:
    """Check for DepMap drug sensitivity data.

    Only adds if tumor-matched cell lines exist (or no tumor type specified).
    """
    if not evidence.depmap_evidence or not evidence.depmap_evidence.drug_sensitivities:
        return

    has_tumor_matched = False
    if ctx.tumor_type and evidence.depmap_evidence.cell_line_models:
        mutant_models = [
            cl for cl in evidence.depmap_evidence.cell_line_models if cl.has_mutation
        ]
        # Args: source_disease (from model), queried_tumor (user's query)
        tumor_models = [
            m
            for m in mutant_models
            if m.primary_disease
            and tumor_types_match(m.primary_disease, ctx.tumor_type)
        ]
        has_tumor_matched = bool(tumor_models)
    elif not ctx.tumor_type:
        # No tumor type specified, so all DepMap data is valid
        has_tumor_matched = True

    if has_tumor_matched:
        n_drugs = len(evidence.depmap_evidence.drug_sensitivities)
        # DepMap drug sensitivity is variant-level (cell lines have the exact mutation)
        ctx.add_well_characterized(
            "preclinical drug sensitivity (DepMap)",
            f"{n_drugs} drugs tested",
            category=GapCategory.PRECLINICAL,
            matches_on=f"{n_drugs} variant",
        )


# =============================================================================
# RESISTANCE MECHANISMS
# =============================================================================


def check_resistance_mechanisms(
    evidence: "Evidence", ctx: "GapDetectionContext"
) -> None:
    """Check for resistance mechanism data.

    Checks multiple sources for resistance signals:
    - PubMed articles flagged as resistance evidence
    - CGI biomarkers with resistance association
    - CIViC assertions with is_resistance=True
    - VICC evidence with resistance response types
    - LLM-extracted literature knowledge mentioning resistance
    """
    # Collect resistance signals from all sources
    resistance_sources: list[str] = []
    counts = MatchCounts()

    # 1. PubMed articles with resistance evidence
    # PubMed articles don't have locus_match (count as gene-level) or tumor type info
    resistance_articles = [
        a for a in evidence.pubmed_articles if a.is_resistance_evidence()
    ]
    if resistance_articles:
        resistance_sources.append(
            f"{len(resistance_articles)} PubMed article{'s' if len(resistance_articles) != 1 else ''}"
        )
        counts.total += len(resistance_articles)
        counts.gene += len(resistance_articles)

    # 2. CGI biomarkers with resistance association
    cgi_resistance = [
        b
        for b in evidence.cgi_biomarkers
        if b.association and "RESIST" in b.association.upper()
    ]
    if cgi_resistance:
        resistance_sources.append(
            f"{len(cgi_resistance)} CGI biomarker{'s' if len(cgi_resistance) != 1 else ''}"
        )
        counts.add(count_with_levels(cgi_resistance, ctx.tumor_type))

    # 3. CIViC assertions with is_resistance=True
    civic_resistance = [a for a in evidence.civic_assertions if a.is_resistance]
    if civic_resistance:
        resistance_sources.append(
            f"{len(civic_resistance)} CIViC assertion{'s' if len(civic_resistance) != 1 else ''}"
        )
        counts.add(count_with_levels(civic_resistance, ctx.tumor_type))

    # 4. CIViC evidence items with is_resistance=True
    # Uses computed property that checks both clinical_significance AND evidence_direction
    civic_evidence_resistance = [e for e in evidence.civic_evidence if e.is_resistance]
    if civic_evidence_resistance:
        resistance_sources.append(
            f"{len(civic_evidence_resistance)} CIViC evidence item{'s' if len(civic_evidence_resistance) != 1 else ''}"
        )
        counts.add(count_with_levels(civic_evidence_resistance, ctx.tumor_type))

    # 5. VICC evidence with resistance response types
    vicc_resistance = [
        v
        for v in evidence.vicc_evidence
        if v.response_type
        and (
            "RESIST" in v.response_type.upper()
            or "REDUCED SENSITIVITY" in v.response_type.upper()
        )
    ]
    if vicc_resistance:
        resistance_sources.append(f"{len(vicc_resistance)} VICC evidence")
        counts.add(count_with_levels(vicc_resistance, ctx.tumor_type))

    # 6. LLM-extracted literature knowledge with resistance signals
    # Literature entries don't have tumor type info
    if evidence.literature_knowledge and evidence.literature_knowledge.resistant_to:
        drugs = evidence.literature_knowledge.get_resistance_drugs(predictive_only=True)
        if drugs:
            resistance_sources.append(
                f"LLM literature ({len(drugs)} drug{'s' if len(drugs) != 1 else ''})"
            )
            # Count without tumor matching (literature entries don't have tumor type)
            lit_counts = count_with_levels(evidence.literature_knowledge.resistant_to)
            counts.total += lit_counts.total
            counts.variant += lit_counts.variant
            counts.codon += lit_counts.codon
            counts.gene += lit_counts.gene

    has_resistance_data = bool(resistance_sources)

    if has_resistance_data:
        ctx.add_well_characterized(
            "resistance mechanisms",
            " + ".join(resistance_sources),
            category=GapCategory.RESISTANCE,
            matches_on=counts.matches_on_str,
            tumor_match=counts.tumor_breakdown_str,
        )

        # Add gap if evidence exists only in other tumors (not tumor-matched)
        if ctx.tumor_type and counts.tumor == 0 and counts.other_cancers:
            other_cancers_str = ", ".join(sorted(counts.other_cancers)[:3])
            ctx.add_gap(
                category=GapCategory.RESISTANCE,
                severity=GapSeverity.MODERATE,  # Data exists, just in different tumor context
                description=f"Resistance data for {ctx.gene} {ctx.variant} exists only in other cancers ({other_cancers_str}), not {ctx.tumor_type}",
                suggested_studies=[
                    "Tumor-specific resistance screen",
                    "Serial biopsy study",
                ],
                addressable_with=["Literature search", "CIViC"],
            )
    elif ctx.has_clinical or ctx.has_drug_data:
        ctx.add_gap(
            category=GapCategory.RESISTANCE,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Resistance mechanisms for {ctx.gene} {ctx.variant} not well characterized",
            suggested_studies=[
                "Serial biopsy study",
                "ctDNA monitoring",
                "Resistance screen",
            ],
            addressable_with=["Literature search", "CIViC"],
        )
        ctx.add_poorly_characterized("resistance mechanisms")


# =============================================================================
# PREVALENCE
# =============================================================================


def check_prevalence(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for prevalence/epidemiology data.

    Note: We only flag missing prevalence as a gap if cBioPortal has data for the gene
    but not the specific variant. If there's no cBioPortal data at all, we skip this
    check entirely - this handles cases like acquired resistance mutations (e.g., T790M)
    that won't appear in baseline sequencing datasets.
    """
    cbio = evidence.cbioportal_evidence

    # If no cBioPortal data at all, skip prevalence check entirely
    # This handles acquired resistance mutations that won't be in baseline sequencing
    if cbio is None or not cbio.has_data():
        return

    # Only consider it "observed" if the variant was actually found in samples
    has_variant_in_samples = cbio.samples_with_exact_variant > 0

    if has_variant_in_samples:
        study = cbio.study_name
        pct = cbio.variant_prevalence_pct
        ctx.add_well_characterized(
            f"observed in sample from study '{study}'",
            f"{pct:.1f}%",
            category=GapCategory.PREVALENCE,
            matches_on="variant",
            tumor_match="Yes",  # cBioPortal data is already tumor-specific
        )
    else:
        # We have gene-level data but not variant-specific data
        # This IS a gap worth flagging (gene is in the dataset but this variant wasn't seen)
        if ctx.is_cancer_gene and ctx.has_clinical:
            severity = GapSeverity.SIGNIFICANT
        elif ctx.is_cancer_gene:
            severity = GapSeverity.MINOR
        else:
            severity = (
                GapSeverity.INFORMATIONAL
            )  # Minor limitation for non-cancer genes
        ctx.add_gap(
            category=GapCategory.PREVALENCE,
            severity=severity,
            description=f"Prevalence of {ctx.gene} {ctx.variant} in {ctx.tumor_type or 'cancer'} unknown (gene seen in {cbio.samples_with_gene_mutation}/{cbio.total_samples} samples, but not this variant)",
            suggested_studies=["Epidemiological study", "Registry analysis"],
            addressable_with=["cBioPortal", "COSMIC", "TCGA"],
        )
        ctx.add_poorly_characterized("prevalence data")


# =============================================================================
# CLINICAL TRIALS
# =============================================================================


def check_clinical_trials(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for active clinical trials."""
    has_trials = bool(evidence.clinical_trials)

    if has_trials:
        n_trials = len(evidence.clinical_trials)

        # Count by locus match and tumor match using properties
        match_counts: dict[str, int] = {"variant": 0, "codon": 0, "gene": 0}
        tumor_match_counts: dict[str, int] = {"tumor": 0, "other": 0}

        for trial in evidence.clinical_trials:
            # Count locus match using locus_match property
            level = trial.locus_match
            if level in match_counts:
                match_counts[level] += 1

            # Count tumor match using tumor_match property
            if trial.tumor_match is True:
                tumor_match_counts["tumor"] += 1
            else:
                tumor_match_counts["other"] += 1

        # Build matches_on string (e.g., "2 variant, 1 gene")
        matches_on_parts = []
        for level in ["variant", "codon", "gene"]:
            if match_counts[level] > 0:
                matches_on_parts.append(f"{match_counts[level]} {level}")
        matches_on_str = ", ".join(matches_on_parts) if matches_on_parts else None

        # Build tumor match string
        tumor_match_parts = []
        if tumor_match_counts["tumor"] > 0:
            tumor_match_parts.append(f"{tumor_match_counts['tumor']} tumor")
        if tumor_match_counts["other"] > 0:
            tumor_match_parts.append(f"{tumor_match_counts['other']} other")
        tumor_match_str = ", ".join(tumor_match_parts) if tumor_match_parts else None

        ctx.add_well_characterized(
            "clinical trial options",
            f"{n_trials} active trial{'s' if n_trials != 1 else ''}",
            category=GapCategory.CLINICAL,
            matches_on=matches_on_str,
            tumor_match=tumor_match_str,
        )
    elif ctx.has_clinical or ctx.has_drug_data:
        ctx.add_gap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.INFORMATIONAL,  # Not a true gap, just noting no active trials
            description=f"No active clinical trials for {ctx.gene} {ctx.variant}",
            suggested_studies=["Clinical trial design", "Basket trial proposal"],
            addressable_with=["ClinicalTrials.gov"],
        )


# =============================================================================
# PRECLINICAL MODELS
# =============================================================================


def check_preclinical_models(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for cell line models with the mutation."""
    has_cell_line_models = evidence.depmap_evidence is not None and bool(
        evidence.depmap_evidence.cell_line_models
    )

    if has_cell_line_models:
        n_models = len(evidence.depmap_evidence.cell_line_models)
        mutant_models = [
            cl for cl in evidence.depmap_evidence.cell_line_models if cl.has_mutation
        ]

        if mutant_models:
            # Check for tumor-type-specific models if tumor type is specified
            if ctx.tumor_type:
                # Args: source_disease (from model), queried_tumor (user's query)
                tumor_models = [
                    m
                    for m in mutant_models
                    if m.primary_disease
                    and tumor_types_match(m.primary_disease, ctx.tumor_type)
                ]
                if tumor_models:
                    # Only add tumor-specific entry (not the general one)
                    # Cell lines have the exact variant (has_mutation=True), so match is "variant"
                    # tumor_match shows how many are tumor-matched vs other histologies
                    n_other = len(mutant_models) - len(tumor_models)
                    tumor_match_str = f"{len(tumor_models)} tumor"
                    if n_other > 0:
                        tumor_match_str += f", {n_other} other"
                    ctx.add_well_characterized(
                        f"{ctx.tumor_type} cell line models ({len(tumor_models)} available)",
                        "DepMap CCLE",
                        category=GapCategory.PRECLINICAL,
                        matches_on=f"{len(tumor_models)} variant",
                        tumor_match=tumor_match_str,
                    )
                else:
                    # No tumor-specific models found - only add gap (not well_characterized)
                    ctx.add_gap(
                        category=GapCategory.PRECLINICAL,
                        severity=GapSeverity.SIGNIFICANT,
                        description=f"Cell line models with {ctx.variant} exist in DepMap but none in {ctx.tumor_type} — cross-histology testing possible",
                        suggested_studies=[
                            f"Test in {ctx.tumor_type}-derived organoids",
                            "Compare drug response vs other histologies",
                            f"Generate isogenic model in {ctx.tumor_type} background",
                        ],
                        addressable_with=[
                            "DepMap CCLE",
                            "Patient-derived organoids",
                            "CRISPR knock-in",
                        ],
                    )
                    ctx.add_poorly_characterized(
                        f"{ctx.tumor_type}-specific preclinical models"
                    )
            else:
                # No tumor type specified - add general entry
                # Cell lines have the exact variant (has_mutation=True), so match is "variant"
                ctx.add_well_characterized(
                    f"model cell lines ({len(mutant_models)} with mutation)",
                    "DepMap CCLE",
                    category=GapCategory.PRECLINICAL,
                    matches_on=f"{len(mutant_models)} variant",
                )
        else:
            # No cell lines have the mutation - this is a gap, not well-characterized
            if ctx.has_drug_data or ctx.has_clinical or evidence.context.gene_role:
                ctx.add_gap(
                    category=GapCategory.PRECLINICAL,
                    severity=GapSeverity.MINOR,
                    description=f"DepMap has {n_models} cell lines for {ctx.gene} but none with {ctx.variant}",
                    suggested_studies=[
                        "Generate isogenic model with mutation",
                        "CRISPR knock-in",
                    ],
                    addressable_with=["DepMap CCLE", "CRISPR screens"],
                )
                ctx.add_poorly_characterized("variant-specific cell models")
    else:
        if ctx.has_drug_data or ctx.has_clinical or evidence.context.gene_role:
            ctx.add_gap(
                category=GapCategory.PRECLINICAL,
                severity=GapSeverity.MINOR,
                description=f"No cell line models identified for {ctx.gene} {ctx.variant}",
                suggested_studies=[
                    "Identify cell lines with mutation",
                    "Generate isogenic model",
                ],
                addressable_with=["DepMap CCLE", "Cellosaurus"],
            )
            ctx.add_poorly_characterized("preclinical model systems")


# =============================================================================
# LITERATURE
# =============================================================================


def check_literature_depth(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for published literature coverage."""
    if not evidence.literature_searched:
        return  # Don't report gap if user chose not to search

    pub_count = len(evidence.pubmed_articles)

    if pub_count == 0:
        severity = (
            GapSeverity.CRITICAL if ctx.is_cancer_gene else GapSeverity.SIGNIFICANT
        )
        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=severity,
            description=f"No published literature found for {ctx.gene} {ctx.variant}",
            suggested_studies=["Case report", "Functional characterization study"],
            addressable_with=["PubMed", "Semantic Scholar", "bioRxiv"],
        )
        ctx.add_poorly_characterized("published literature")
    elif pub_count < LITERATURE_WELL_CHARACTERIZED_THRESHOLD:
        ctx.add_poorly_characterized("literature depth (limited publications)")
    else:
        # Show ">= N" when we hit the limit (indicates more articles may exist)
        count_str = (
            f">= {pub_count}" if pub_count >= MAX_LITERATURE_RESULTS else str(pub_count)
        )
        ctx.add_well_characterized(
            "published literature",
            f"{count_str} PubMed articles",
            category=GapCategory.FUNCTIONAL,
        )


def check_literature_database_integration(
    evidence: "Evidence", ctx: "GapDetectionContext"
) -> None:
    """Check for literature evidence not integrated into curated databases.

    Identifies when PubMed articles mention drug sensitivity/resistance signals
    that aren't represented in CIViC, CGI, VICC, or FDA databases.
    This is an "integration gap" - evidence exists but hasn't been curated.
    """
    if not evidence.literature_searched or not evidence.pubmed_articles:
        return

    # Collect drugs mentioned in literature with sensitivity/resistance signals
    literature_drugs: set[str] = set()
    for article in evidence.pubmed_articles:
        if article.drugs_mentioned:
            if article.is_sensitivity_evidence() or article.is_resistance_evidence():
                for drug in article.drugs_mentioned:
                    literature_drugs.add(drug.lower())

    if not literature_drugs:
        return

    # Collect drugs in curated databases
    curated_drugs: set[str] = set()

    # FDA biomarker evidence (only REQUIRED_POSITIVE)
    from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement

    for ev in evidence.fda_biomarker_evidence:
        if ev.requirement != BiomarkerRequirement.REQUIRED_NEGATIVE:
            if ev.drug_name:
                curated_drugs.add(ev.drug_name.lower())

    # CIViC assertions
    for assertion in evidence.civic_assertions:
        if assertion.therapies:
            for therapy in assertion.therapies:
                curated_drugs.add(therapy.lower())

    # CIViC evidence
    for civic_ev in evidence.civic_evidence:
        if civic_ev.drugs:
            for drug in civic_ev.drugs:
                curated_drugs.add(drug.lower())

    # CGI biomarkers (all tiers)
    for biomarker in evidence.cgi_biomarkers:
        if biomarker.drug:
            curated_drugs.add(biomarker.drug.lower())
    for biomarker in evidence.preclinical_biomarkers:
        if biomarker.drug:
            curated_drugs.add(biomarker.drug.lower())
    for biomarker in evidence.early_phase_biomarkers:
        if biomarker.drug:
            curated_drugs.add(biomarker.drug.lower())

    # VICC evidence
    for vicc in evidence.vicc_evidence:
        if vicc.drugs:
            for drug in vicc.drugs:
                curated_drugs.add(drug.lower())

    # Find drugs in literature but NOT in curated databases
    uncurated_drugs = literature_drugs - curated_drugs

    if uncurated_drugs:
        # This is an integration gap - evidence exists but not in databases
        drug_list = ", ".join(sorted(uncurated_drugs)[:3])
        more_count = len(uncurated_drugs) - 3 if len(uncurated_drugs) > 3 else 0
        more_str = f" (+{more_count} more)" if more_count > 0 else ""

        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.INFORMATIONAL,  # Curation opportunity, not research gap
            description=f"Literature mentions drug signals not in curated databases: {drug_list}{more_str}",
            suggested_studies=[
                "Submit to CIViC for curation",
                "Validate findings in independent cohort",
                "Cross-reference with clinical trial results",
            ],
            addressable_with=["CIViC submission", "Literature review"],
        )


# =============================================================================
# VALIDATION GAP
# =============================================================================


def check_validation_gap(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for strong oncogenic signal with limited therapeutic validation."""
    has_strong_oncogenic_signal = (
        ctx.has_pathogenic_signal
        and evidence.depmap_evidence is not None
        and evidence.depmap_evidence.is_essential()
    )

    has_therapeutic_validation = bool(evidence.civic_assertions) or bool(
        evidence.fda_biomarker_evidence
    )

    if has_strong_oncogenic_signal:
        ctx.add_well_characterized(
            "biological driver potential",
            "Pathogenic prediction + gene essentiality",
            category=GapCategory.VALIDATION,
        )

        if not has_therapeutic_validation:
            ctx.add_gap(
                category=GapCategory.VALIDATION,
                severity=(
                    GapSeverity.CRITICAL
                    if ctx.is_cancer_gene
                    else GapSeverity.SIGNIFICANT
                ),
                description="Strong oncogenic signal but limited therapeutic validation",
                suggested_studies=[
                    "Functional validation in isogenic models",
                    "Drug sensitivity screening (PRISM/GDSC)",
                    "Patient-derived organoid testing",
                    "Structural modeling of activation mechanism",
                ],
                addressable_with=[
                    "AlphaFold",
                    "DepMap",
                    "CRISPR screens",
                    "Literature",
                ],
            )
            ctx.add_poorly_characterized("therapeutic validation")

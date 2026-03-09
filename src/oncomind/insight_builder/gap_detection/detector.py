"""Main entry point for evidence gap detection.

Orchestrates all gap detection checks and produces the final EvidenceGaps result.
"""

from typing import TYPE_CHECKING

from oncomind.models.extracted.evidence_gaps import EvidenceGaps, GapCategory
from oncomind.models.evidence.base import tumor_types_match
from oncomind.config.constants import COOCCURRENCE_STRONG_THRESHOLD_PCT

from .context import GapDetectionContext
from .helpers import (
    has_pathogenic_signal,
    sort_characterized_by_category,
    compute_overall_quality,
    compute_research_priority,
    get_primary_drug,
    get_top_sensitive_drugs,
    get_top_cooccurring_gene,
    has_strong_cooccurrence,
    get_top_codependencies,
)
from .checks import (
    check_hotspot_context,
    check_functional_predictions,
    check_gene_mechanism,
    check_clinical_evidence,
    check_tumor_type_evidence,
    check_drug_response,
    check_preclinical_biomarkers,
    check_depmap_drug_sensitivity,
    check_resistance_mechanisms,
    check_prevalence,
    check_clinical_trials,
    check_preclinical_models,
    check_literature_depth,
    check_literature_database_integration,
    check_validation_gap,
)
from .discordance import check_discordant_evidence

if TYPE_CHECKING:
    from oncomind.models.evidence import Evidence


def detect_evidence_gaps(evidence: "Evidence") -> EvidenceGaps:
    """Detect evidence gaps from aggregated evidence.

    Args:
        evidence: Aggregated evidence from all sources

    Returns:
        EvidenceGaps with identified gaps and assessment
    """
    # Initialize context
    ctx = GapDetectionContext(
        gene=evidence.identifiers.gene,
        variant=evidence.identifiers.variant,
        tumor_type=evidence.context.tumor_type,
        is_cancer_gene=evidence.context.gene_role
        in ("oncogene", "TSG", "tumor_suppressor", "ddr", "tsg_pathway_actionable"),
        has_pathogenic_signal=has_pathogenic_signal(evidence),
    )

    # Run all gap detection checks
    check_hotspot_context(evidence, ctx)
    check_functional_predictions(evidence, ctx)
    check_gene_mechanism(evidence, ctx)
    check_clinical_evidence(evidence, ctx)
    check_tumor_type_evidence(evidence, ctx)
    check_drug_response(evidence, ctx)
    check_preclinical_biomarkers(evidence, ctx)
    check_depmap_drug_sensitivity(evidence, ctx)
    check_resistance_mechanisms(evidence, ctx)
    check_discordant_evidence(evidence, ctx)
    check_prevalence(evidence, ctx)
    check_clinical_trials(evidence, ctx)
    check_preclinical_models(evidence, ctx)
    check_literature_depth(evidence, ctx)
    check_literature_database_integration(evidence, ctx)
    check_validation_gap(evidence, ctx)

    # Enrich gaps with dynamic, context-aware suggestions based on actual evidence
    _enrich_gaps_with_context(evidence, ctx)

    # Compute overall assessments
    overall_quality = compute_overall_quality(ctx.gaps, len(ctx.well_characterized))
    research_priority = compute_research_priority(
        evidence,
        ctx.gaps,
        overall_quality,
        ctx.is_cancer_gene,
        ctx.has_pathogenic_signal,
    )

    # Sort well_characterized_detailed by category for grouped display
    sorted_well_characterized = sort_characterized_by_category(
        ctx.well_characterized_detailed
    )

    return EvidenceGaps(
        gaps=ctx.gaps,
        overall_evidence_quality=overall_quality,
        well_characterized=ctx.well_characterized,
        well_characterized_detailed=sorted_well_characterized,
        poorly_characterized=ctx.poorly_characterized,
        research_priority=research_priority,
    )


def _enrich_gaps_with_context(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Enrich gap suggestions with dynamic, context-aware recommendations.

    Analyzes actual evidence to generate specific, actionable study suggestions
    based on gene, variant, tumor type, and available data.
    """
    gene = ctx.gene
    variant = ctx.variant
    tumor_type = ctx.tumor_type

    # Get primary approved drug if available
    primary_drug = get_primary_drug(evidence)

    # Get top sensitive drugs from DepMap if available (only if tumor-matched cell lines exist)
    top_sensitive_drugs = get_top_sensitive_drugs(evidence, tumor_type)

    # Get top co-occurring gene if available
    top_cooc_gene = get_top_cooccurring_gene(evidence)

    # Check for strong co-occurrence signal
    has_strong_cooc = has_strong_cooccurrence(evidence)

    # Iterate through gaps and add context-specific suggestions
    for gap in ctx.gaps:
        new_suggestions = []

        # RESISTANCE gaps: suggest bypass mechanism testing for the primary drug
        if gap.category == GapCategory.RESISTANCE and primary_drug:
            new_suggestions.append(
                f"Test bypass mechanisms for {primary_drug} resistance in {gene}-mutant models"
            )
            new_suggestions.append(
                f"ctDNA monitoring for {gene} {variant} emergence under {primary_drug} treatment"
            )

        # PRECLINICAL gaps: suggest testing with DepMap-identified sensitive drugs
        if gap.category == GapCategory.PRECLINICAL and top_sensitive_drugs:
            drugs_str = ", ".join(top_sensitive_drugs[:3])
            new_suggestions.append(
                f"Validate sensitivity to {drugs_str} in isogenic {gene} {variant} models"
            )

        # VALIDATION gaps: suggest synthetic lethality with strong co-occurring gene
        if gap.category == GapCategory.VALIDATION and has_strong_cooc and top_cooc_gene:
            new_suggestions.append(
                f"Investigate synthetic lethality with {top_cooc_gene} co-mutation"
            )
            new_suggestions.append(
                f"CRISPR screen in {gene}/{top_cooc_gene} double-mutant background"
            )

        # DRUG_RESPONSE gaps: suggest testing FDA-approved drugs in this context
        if gap.category == GapCategory.DRUG_RESPONSE:
            if primary_drug and tumor_type:
                new_suggestions.append(
                    f"Evaluate {primary_drug} efficacy in {tumor_type} with {gene} {variant}"
                )
            if top_sensitive_drugs:
                new_suggestions.append(
                    f"Confirm DepMap-predicted sensitivities ({', '.join(top_sensitive_drugs[:2])}) in PDX"
                )

        # TUMOR_TYPE gaps: suggest basket trial or cross-histology comparison
        if gap.category == GapCategory.TUMOR_TYPE and tumor_type:
            # Check if we have evidence in other tumor types
            if evidence.civic_assertions or evidence.fda_biomarker_evidence:
                new_suggestions.append(
                    f"Retrospective analysis of {gene} {variant} response in {tumor_type} vs other histologies"
                )
            if primary_drug:
                new_suggestions.append(
                    f"Basket trial cohort for {primary_drug} in {tumor_type} with {gene} {variant}"
                )

        # FUNCTIONAL gaps: suggest structural modeling with specific focus
        if gap.category == GapCategory.FUNCTIONAL:
            if evidence.context.gene_role == "oncogene":
                new_suggestions.append(
                    f"Kinase activity assay for {gene} {variant} vs wild-type"
                )
            elif evidence.context.gene_role in (
                "TSG",
                "tumor_suppressor",
                "tsg_pathway_actionable",
            ):
                new_suggestions.append(
                    f"LOF assay: assess {gene} {variant} impact on tumor suppressor function"
                )

        # Add co-dependency suggestions if relevant
        if gap.category == GapCategory.VALIDATION:
            co_deps = get_top_codependencies(evidence)
            if co_deps:
                deps_str = ", ".join(co_deps[:2])
                new_suggestions.append(
                    f"Test synthetic lethality with co-essential genes: {deps_str}"
                )

        # Append new suggestions to existing ones (avoid duplicates)
        for suggestion in new_suggestions:
            if suggestion not in gap.suggested_studies:
                gap.suggested_studies.append(suggestion)

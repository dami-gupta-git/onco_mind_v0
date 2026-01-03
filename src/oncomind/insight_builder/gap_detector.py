"""Evidence gap detection from aggregated evidence.

Analyzes what's missing from the evidence to guide research priorities.
"""

from dataclasses import dataclass, field
from oncomind.models.evidence.evidence_gaps import (
    EvidenceGaps, EvidenceGap, GapCategory, GapSeverity, CharacterizedAspect
)
from oncomind.models.gene_context import is_hotspot_variant, is_hotspot_adjacent, _extract_codon_position
from oncomind.config.constants import TUMOR_TYPE_MAPPINGS

# Import Evidence with TYPE_CHECKING to avoid circular imports
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from oncomind.models.evidence import Evidence


def _tumor_type_matches(tumor_type: str, tissue: str) -> bool:
    """Check if a tissue/disease matches the expected tumor type.

    Uses TUMOR_TYPE_MAPPINGS to handle aliases like:
    - "Melanoma" matching "SKIN"
    - "NSCLC" matching "LUNG"

    Args:
        tumor_type: The tumor type we're looking for (e.g., "Melanoma")
        tissue: The tissue/disease from the cell line (e.g., "SKIN")

    Returns:
        True if there's a match
    """
    if not tumor_type or not tissue:
        return False

    tumor_lower = tumor_type.lower()
    tissue_lower = tissue.lower()

    # Direct match
    if tumor_lower in tissue_lower or tissue_lower in tumor_lower:
        return True

    # Check via tumor type mappings
    for abbrev, aliases in TUMOR_TYPE_MAPPINGS.items():
        # Check if tumor_type matches any alias
        tumor_matches = (
            tumor_lower == abbrev or
            any(tumor_lower in alias or alias in tumor_lower for alias in aliases)
        )
        # Check if tissue matches any alias
        tissue_matches = (
            tissue_lower == abbrev or
            any(tissue_lower in alias or alias in tissue_lower for alias in aliases)
        )
        # If both match the same category, they're related
        if tumor_matches and tissue_matches:
            return True

    return False


# =============================================================================
# GAP DETECTION CONTEXT
# =============================================================================

@dataclass
class GapDetectionContext:
    """Context accumulated during gap detection.

    Tracks gaps found, well-characterized aspects, and flags for use
    across multiple detection functions.
    """
    gene: str
    variant: str
    tumor_type: str | None
    is_cancer_gene: bool
    has_pathogenic_signal: bool

    # Accumulated results
    gaps: list[EvidenceGap] = field(default_factory=list)
    well_characterized: list[str] = field(default_factory=list)
    well_characterized_detailed: list[CharacterizedAspect] = field(default_factory=list)
    poorly_characterized: list[str] = field(default_factory=list)

    # Flags set during detection (used by later checks)
    has_clinical: bool = False
    has_drug_data: bool = False

    def add_well_characterized(
        self,
        aspect: str,
        basis: str,
        category: GapCategory | None = None,
        matches_on: str | None = None,
        cancer_mismatch: str | None = None
    ) -> None:
        """Add a well-characterized aspect with its basis and category."""
        # Use title case for aspect
        aspect_title = aspect.title()
        self.well_characterized.append(aspect_title)
        self.well_characterized_detailed.append(
            CharacterizedAspect(
                aspect=aspect_title,
                basis=basis,
                category=category,
                matches_on=matches_on,
                cancer_mismatch=cancer_mismatch
            )
        )

    def add_gap(
        self,
        category: GapCategory,
        severity: GapSeverity,
        description: str,
        suggested_studies: list[str] | None = None,
        addressable_with: list[str] | None = None,
    ) -> None:
        """Add an evidence gap."""
        self.gaps.append(EvidenceGap(
            category=category,
            severity=severity,
            description=description,
            suggested_studies=suggested_studies or [],
            addressable_with=addressable_with or [],
        ))

    def add_poorly_characterized(self, aspect: str) -> None:
        """Add a poorly-characterized aspect."""
        self.poorly_characterized.append(aspect)


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

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
        is_cancer_gene=evidence.context.gene_role in (
            "oncogene", "TSG", "tumor_suppressor", "ddr", "tsg_pathway_actionable"
        ),
        has_pathogenic_signal=_has_pathogenic_signal(evidence),
    )

    # Run all gap detection checks
    _check_hotspot_context(evidence, ctx)
    _check_functional_predictions(evidence, ctx)
    _check_gene_mechanism(evidence, ctx)
    _check_clinical_evidence(evidence, ctx)
    _check_tumor_type_evidence(evidence, ctx)
    _check_drug_response(evidence, ctx)
    _check_resistance_mechanisms(evidence, ctx)
    _check_discordant_evidence(evidence, ctx)
    _check_prevalence(evidence, ctx)
    _check_clinical_trials(evidence, ctx)
    _check_preclinical_models(evidence, ctx)
    _check_literature_depth(evidence, ctx)
    _check_validation_gap(evidence, ctx)

    # Enrich gaps with dynamic, context-aware suggestions based on actual evidence
    _enrich_gaps_with_context(evidence, ctx)

    # Compute overall assessments
    overall_quality = _compute_overall_quality(ctx.gaps, len(ctx.well_characterized))
    research_priority = _compute_research_priority(evidence, ctx.gaps, overall_quality)

    # Sort well_characterized_detailed by category for grouped display
    sorted_well_characterized = _sort_characterized_by_category(ctx.well_characterized_detailed)

    return EvidenceGaps(
        gaps=ctx.gaps,
        overall_evidence_quality=overall_quality,
        well_characterized=ctx.well_characterized,
        well_characterized_detailed=sorted_well_characterized,
        poorly_characterized=ctx.poorly_characterized,
        research_priority=research_priority,
    )


# =============================================================================
# INDIVIDUAL GAP DETECTION FUNCTIONS
# =============================================================================

def _check_hotspot_context(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check if variant is at or near a known cancer hotspot."""
    is_hotspot = is_hotspot_variant(ctx.gene, ctx.variant)
    is_adjacent, nearest_hotspot = is_hotspot_adjacent(ctx.gene, ctx.variant, window=5)

    if is_hotspot:
        ctx.add_well_characterized(
            "known cancer hotspot",
            f"Codon {_extract_codon_position(ctx.variant)} is in cancerhotspots.org",
            category=GapCategory.FUNCTIONAL,
            matches_on="codon"
        )
    elif is_adjacent and nearest_hotspot:
        # Rare variant near a hotspot - high research value
        ctx.add_well_characterized(
            f"near hotspot codon {nearest_hotspot}",
            f"Within 5 codons of known hotspot — structural similarity likely",
            category=GapCategory.FUNCTIONAL
        )
        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Rare variant near known hotspot (codon {nearest_hotspot}) — functional characterization needed",
            suggested_studies=[
                f"Compare to nearby hotspot {ctx.gene} codon {nearest_hotspot}",
                "Structural modeling to assess activation mechanism",
                "Functional assay (transformation, signaling)"
            ],
            addressable_with=["AlphaFold", "Literature on nearby hotspot", "Isogenic models"]
        )
        ctx.add_poorly_characterized("rare-near-hotspot variant function")


def _check_functional_predictions(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for computational pathogenicity predictions."""
    has_functional = (
        evidence.functional.alphamissense_score is not None or
        evidence.functional.cadd_score is not None or
        evidence.functional.polyphen2_prediction is not None or
        evidence.functional.sift_prediction is not None
    )

    if has_functional:
        func_sources = []
        if evidence.functional.alphamissense_score is not None:
            func_sources.append(f"AlphaMissense={evidence.functional.alphamissense_score:.2f}")
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
            matches_on="variant"
        )
    else:
        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"No computational pathogenicity predictions for {ctx.gene} {ctx.variant}",
            suggested_studies=["Run AlphaMissense, CADD, PolyPhen2, SIFT"],
            addressable_with=["MyVariant.info", "VEP"]
        )
        ctx.add_poorly_characterized("pathogenicity predictions")

    # Check for gnomAD population frequency data (informational, not a penalty)
    # If gnomAD AF > 0.01% (0.0001), note that the variant is observed in the general population
    gnomad_af = evidence.functional.gnomad_exome_af or evidence.functional.gnomad_genome_af
    if gnomad_af is not None and gnomad_af > 0.0001:  # > 0.01%
        af_pct = gnomad_af * 100
        ctx.add_well_characterized(
            "population frequency",
            f"Observed in general population (gnomAD AF: {af_pct:.3f}%)",
            category=GapCategory.PREVALENCE
        )


def _check_gene_mechanism(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for gene function and essentiality data."""
    has_mechanism = bool(evidence.context.gene_role) or bool(evidence.context.pathway)
    has_depmap_essentiality = (
        evidence.depmap_evidence is not None and
        evidence.depmap_evidence.gene_dependency is not None
    )

    if has_mechanism:
        role = evidence.context.gene_role or "unknown"
        pathway = evidence.context.pathway or ""
        basis = f"Role: {role}" + (f", Pathway: {pathway}" if pathway else "")
        ctx.add_well_characterized("gene function", basis, category=GapCategory.FUNCTIONAL, matches_on="gene")

    # Only show gene essentiality if the gene IS essential (CERES < -0.5)
    # Gene essentiality is pan-cancer data, so don't show non-essential scores
    if has_depmap_essentiality and evidence.depmap_evidence.is_essential():
        dep = evidence.depmap_evidence.gene_dependency
        score = dep.mean_dependency_score if dep else None
        pct = dep.dependency_pct if dep else 0
        score_str = f"CERES={score:.2f}, {pct:.0f}% of cell lines depend on it"
        ctx.add_well_characterized("gene essentiality", score_str, category=GapCategory.FUNCTIONAL, matches_on="gene")

    if not has_mechanism and not has_depmap_essentiality:
        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Functional impact of {ctx.variant} on {ctx.gene} protein is unknown",
            suggested_studies=["Functional assay", "Structural modeling", "Cell-based reporter"],
            addressable_with=["UniProt", "Literature search", "DepMap"]
        )
        ctx.add_poorly_characterized("functional mechanism")


def _check_clinical_evidence(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for FDA approvals, CIViC assertions, and CIViC evidence items."""
    ctx.has_clinical = bool(evidence.civic_assertions) or bool(evidence.civic_evidence) or bool(evidence.fda_approvals)

    if ctx.has_clinical:
        n_assertions = len(evidence.civic_assertions)
        n_fda = len(evidence.fda_approvals)
        parts = []

        # Check if FDA approvals match the queried tumor type
        fda_matching = 0
        fda_other_cancers: list[str] = []  # List of cancer types that DON'T match

        if n_fda:
            tumor_type = evidence.context.tumor_type
            for approval in evidence.fda_approvals:
                if tumor_type and approval.indication:
                    parsed = approval.parse_indication_for_tumor(tumor_type)
                    if parsed.get('tumor_match'):
                        fda_matching += 1
                    else:
                        # Extract what cancer it IS approved for
                        indication_cancer = approval.extract_indication_cancer_type()
                        if indication_cancer and "pan-cancer" not in indication_cancer.lower():
                            if indication_cancer not in fda_other_cancers:
                                fda_other_cancers.append(indication_cancer)
                        else:
                            fda_matching += 1  # pan-cancer counts as matching
                else:
                    fda_matching += 1  # No tumor type specified, count as matching

            parts.append(f"{n_fda} FDA approval{'s' if n_fda > 1 else ''}")

        if n_assertions:
            parts.append(f"{n_assertions} CIViC assertion{'s' if n_assertions > 1 else ''}")

        # Determine cancer_mismatch value
        cancer_mismatch = None
        if fda_other_cancers and fda_matching == 0:
            # ALL FDA approvals are for different cancers
            cancer_mismatch = ", ".join(fda_other_cancers[:2])  # Show first 2

        ctx.add_well_characterized(
            "clinical actionability",
            " + ".join(parts) if parts else "Clinical evidence exists",
            category=GapCategory.CLINICAL,
            cancer_mismatch=cancer_mismatch
        )
    else:
        ctx.add_gap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.CRITICAL,
            description=f"No curated clinical evidence for {ctx.gene} {ctx.variant}",
            suggested_studies=["Case series", "Retrospective cohort", "Basket trial inclusion"],
            addressable_with=["CIViC submission", "Literature curation"]
        )
        ctx.add_poorly_characterized("clinical evidence")


def _check_tumor_type_evidence(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for tumor-type-specific evidence."""
    if not ctx.tumor_type:
        return

    tumor_specific = _check_tumor_specific_evidence(evidence, ctx.tumor_type)

    if tumor_specific:
        ctx.add_well_characterized(
            f"evidence in {ctx.tumor_type}",
            "Tumor-specific CIViC/FDA/VICC/CGI data",
            category=GapCategory.TUMOR_TYPE
        )
    else:
        # Severity depends on gene importance and pathogenic signal
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
                "Basket trial with histology-specific cohort"
            ],
            addressable_with=["ClinicalTrials.gov", "Literature search"]
        )
        ctx.add_poorly_characterized(f"{ctx.tumor_type}-specific data")


def _check_drug_response(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for drug sensitivity/resistance data."""
    # Check if DepMap has drug data AND tumor-matched cell lines
    has_depmap_drug_data = False
    has_tumor_matched_depmap = False
    if evidence.depmap_evidence is not None and evidence.depmap_evidence.drug_sensitivities:
        has_depmap_drug_data = True
        # Only count as tumor-matched if we have cell lines matching the tumor type
        if ctx.tumor_type and evidence.depmap_evidence.cell_line_models:
            mutant_models = [cl for cl in evidence.depmap_evidence.cell_line_models if cl.has_mutation]
            tumor_models = [
                m for m in mutant_models
                if m.primary_disease and _tumor_type_matches(ctx.tumor_type, m.primary_disease)
            ]
            has_tumor_matched_depmap = bool(tumor_models)
        elif not ctx.tumor_type:
            # No tumor type specified, so all DepMap data is valid
            has_tumor_matched_depmap = True

    ctx.has_drug_data = (
        bool(evidence.cgi_biomarkers) or
        bool(evidence.vicc_evidence) or
        bool(evidence.preclinical_biomarkers) or
        has_tumor_matched_depmap  # Only count tumor-matched DepMap data
    )

    if ctx.has_drug_data:
        n_cgi = len(evidence.cgi_biomarkers)
        n_vicc = len(evidence.vicc_evidence)
        n_preclin = len(evidence.preclinical_biomarkers)
        drug_sources = []
        if n_cgi:
            drug_sources.append(f"{n_cgi} CGI")
        if n_vicc:
            drug_sources.append(f"{n_vicc} VICC")
        if n_preclin:
            drug_sources.append(f"{n_preclin} preclinical")

        # Compute match level breakdown for drug response data
        variant_count = 0
        codon_count = 0
        gene_count = 0
        for b in evidence.cgi_biomarkers:
            level = getattr(b, 'match_level', 'gene') or 'gene'
            if level == 'variant':
                variant_count += 1
            elif level == 'codon':
                codon_count += 1
            else:
                gene_count += 1
        for v in evidence.vicc_evidence:
            level = getattr(v, 'match_level', 'gene') or 'gene'
            if level == 'variant':
                variant_count += 1
            elif level == 'codon':
                codon_count += 1
            else:
                gene_count += 1

        # Build matches_on string
        match_parts = []
        if variant_count > 0:
            match_parts.append(f"{variant_count} variant")
        if codon_count > 0:
            match_parts.append(f"{codon_count} codon")
        if gene_count > 0:
            match_parts.append(f"{gene_count} gene")
        matches_on = ", ".join(match_parts) if match_parts else None

        ctx.add_well_characterized(
            "drug response",
            " + ".join(drug_sources) if drug_sources else "Drug data available",
            category=GapCategory.DRUG_RESPONSE,
            matches_on=matches_on
        )
        # Only add DepMap drug sensitivity if tumor-matched
        if has_tumor_matched_depmap:
            n_drugs = len(evidence.depmap_evidence.drug_sensitivities)
            ctx.add_well_characterized(
                "preclinical drug sensitivity (DepMap)",
                f"{n_drugs} drugs tested",
                category=GapCategory.DRUG_RESPONSE
            )
    else:
        ctx.add_gap(
            category=GapCategory.DRUG_RESPONSE,
            severity=GapSeverity.SIGNIFICANT,
            description=f"No drug sensitivity/resistance data for {ctx.gene} {ctx.variant}",
            suggested_studies=["Cell line drug screen", "PDX drug testing", "Clinical correlative study"],
            addressable_with=["GDSC", "CTRP", "DepMap"]
        )
        ctx.add_poorly_characterized("drug response")


def _check_resistance_mechanisms(evidence: "Evidence", ctx: GapDetectionContext) -> None:
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

    # Track match levels for all resistance data
    variant_count = 0
    codon_count = 0
    gene_count = 0

    # 1. PubMed articles with resistance evidence
    resistance_articles = [a for a in evidence.pubmed_articles if a.is_resistance_evidence()]
    if resistance_articles:
        resistance_sources.append(f"{len(resistance_articles)} PubMed article{'s' if len(resistance_articles) != 1 else ''}")
        # PubMed articles don't have match_level, count as gene-level
        gene_count += len(resistance_articles)

    # 2. CGI biomarkers with resistance association
    cgi_resistance = [
        b for b in evidence.cgi_biomarkers
        if b.association and "RESIST" in b.association.upper()
    ]
    if cgi_resistance:
        resistance_sources.append(f"{len(cgi_resistance)} CGI biomarker{'s' if len(cgi_resistance) != 1 else ''}")
        for b in cgi_resistance:
            level = getattr(b, 'match_level', 'gene') or 'gene'
            if level == 'variant':
                variant_count += 1
            elif level == 'codon':
                codon_count += 1
            else:
                gene_count += 1

    # 3. CIViC assertions with is_resistance=True
    civic_resistance = [a for a in evidence.civic_assertions if a.is_resistance]
    if civic_resistance:
        resistance_sources.append(f"{len(civic_resistance)} CIViC assertion{'s' if len(civic_resistance) != 1 else ''}")
        for a in civic_resistance:
            level = getattr(a, 'match_level', 'gene') or 'gene'
            if level == 'variant':
                variant_count += 1
            elif level == 'codon':
                codon_count += 1
            else:
                gene_count += 1

    # 4. VICC evidence with resistance response types
    vicc_resistance = [
        v for v in evidence.vicc_evidence
        if v.response_type and ("RESIST" in v.response_type.upper() or "REDUCED SENSITIVITY" in v.response_type.upper())
    ]
    if vicc_resistance:
        resistance_sources.append(f"{len(vicc_resistance)} VICC evidence")
        for v in vicc_resistance:
            level = getattr(v, 'match_level', 'gene') or 'gene'
            if level == 'variant':
                variant_count += 1
            elif level == 'codon':
                codon_count += 1
            else:
                gene_count += 1

    # 5. LLM-extracted literature knowledge with resistance signals
    if evidence.literature_knowledge and evidence.literature_knowledge.resistant_to:
        drugs = evidence.literature_knowledge.get_resistance_drugs(predictive_only=True)
        if drugs:
            resistance_sources.append(f"LLM literature ({len(drugs)} drug{'s' if len(drugs) != 1 else ''})")
            # Count match levels from literature resistance entries
            for entry in evidence.literature_knowledge.resistant_to:
                level = getattr(entry, 'match_level', 'gene') or 'gene'
                if level == 'variant':
                    variant_count += 1
                elif level == 'codon':
                    codon_count += 1
                else:
                    gene_count += 1

    has_resistance_data = bool(resistance_sources)

    if has_resistance_data:
        # Build matches_on string
        match_parts = []
        if variant_count > 0:
            match_parts.append(f"{variant_count} variant")
        if codon_count > 0:
            match_parts.append(f"{codon_count} codon")
        if gene_count > 0:
            match_parts.append(f"{gene_count} gene")
        matches_on = ", ".join(match_parts) if match_parts else None

        ctx.add_well_characterized(
            "resistance mechanisms",
            " + ".join(resistance_sources),
            category=GapCategory.RESISTANCE,
            matches_on=matches_on
        )
    elif ctx.has_clinical or ctx.has_drug_data:
        ctx.add_gap(
            category=GapCategory.RESISTANCE,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Resistance mechanisms for {ctx.gene} {ctx.variant} not well characterized",
            suggested_studies=["Serial biopsy study", "ctDNA monitoring", "Resistance screen"],
            addressable_with=["Literature search", "CIViC"]
        )
        ctx.add_poorly_characterized("resistance mechanisms")


def _check_discordant_evidence(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for conflicting evidence between sources."""
    discordant_findings = _detect_discordant_evidence_internal(evidence)

    for finding in discordant_findings:
        ctx.add_gap(
            category=GapCategory.DISCORDANT,
            severity=GapSeverity.SIGNIFICANT,
            description=finding,
            suggested_studies=["Meta-analysis", "Prospective validation study"],
            addressable_with=["Literature review", "Expert consensus"]
        )
        ctx.add_poorly_characterized("conflicting evidence")


def _check_prevalence(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for prevalence/epidemiology data."""
    cbio = evidence.cbioportal_evidence

    # Only consider it "observed" if the variant was actually found in samples
    has_variant_in_samples = (
        cbio is not None and
        cbio.has_data() and
        cbio.samples_with_exact_variant > 0
    )

    if has_variant_in_samples:
        study = cbio.study_name if cbio else "cBioPortal"
        pct = cbio.variant_prevalence_pct if cbio else 0
        ctx.add_well_characterized(
            "observed in samples",
            f"{pct:.1f}% in {study}",
            category=GapCategory.PREVALENCE,
            matches_on="variant"
        )
    else:
        severity = GapSeverity.SIGNIFICANT if (ctx.is_cancer_gene and ctx.has_clinical) else GapSeverity.MINOR
        ctx.add_gap(
            category=GapCategory.PREVALENCE,
            severity=severity,
            description=f"Prevalence of {ctx.gene} {ctx.variant} in {ctx.tumor_type or 'cancer'} unknown",
            suggested_studies=["Epidemiological study", "Registry analysis"],
            addressable_with=["cBioPortal", "COSMIC", "TCGA"]
        )
        ctx.add_poorly_characterized("prevalence data")


def _check_clinical_trials(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for active clinical trials."""
    has_trials = bool(evidence.clinical_trials)

    if has_trials:
        n_trials = len(evidence.clinical_trials)
        ctx.add_well_characterized(
            "clinical trial options",
            f"{n_trials} active trial{'s' if n_trials != 1 else ''}",
            category=GapCategory.CLINICAL
        )
    elif ctx.has_clinical or ctx.has_drug_data:
        ctx.add_gap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.MINOR,
            description=f"No active clinical trials for {ctx.gene} {ctx.variant}",
            suggested_studies=["Clinical trial design", "Basket trial proposal"],
            addressable_with=["ClinicalTrials.gov"]
        )


def _check_preclinical_models(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for cell line models with the mutation."""
    has_cell_line_models = (
        evidence.depmap_evidence is not None and
        bool(evidence.depmap_evidence.cell_line_models)
    )

    if has_cell_line_models:
        n_models = len(evidence.depmap_evidence.cell_line_models)
        mutant_models = [
            cl for cl in evidence.depmap_evidence.cell_line_models if cl.has_mutation
        ]

        if mutant_models:
            # Check for tumor-type-specific models if tumor type is specified
            if ctx.tumor_type:
                tumor_models = [
                    m for m in mutant_models
                    if m.primary_disease and _tumor_type_matches(ctx.tumor_type, m.primary_disease)
                ]
                if tumor_models:
                    # Only add tumor-specific entry (not the general one)
                    # Cell lines have the exact variant (has_mutation=True), so match is "variant"
                    ctx.add_well_characterized(
                        f"{ctx.tumor_type} cell line models ({len(tumor_models)} available)",
                        "DepMap CCLE",
                        category=GapCategory.PRECLINICAL,
                        matches_on=f"{len(tumor_models)} variant"
                    )
                else:
                    # No tumor-specific models found - only add gap (not well_characterized)
                    ctx.add_gap(
                        category=GapCategory.PRECLINICAL,
                        severity=GapSeverity.SIGNIFICANT,
                        description=f"Models with {ctx.variant} exist but none in {ctx.tumor_type} — cross-histology testing possible",
                        suggested_studies=[
                            f"Test in {ctx.tumor_type}-derived organoids",
                            "Compare drug response vs other histologies",
                            f"Generate isogenic model in {ctx.tumor_type} background"
                        ],
                        addressable_with=["DepMap CCLE", "Patient-derived organoids", "CRISPR knock-in"]
                    )
                    ctx.add_poorly_characterized(f"{ctx.tumor_type}-specific preclinical models")
            else:
                # No tumor type specified - add general entry
                # Cell lines have the exact variant (has_mutation=True), so match is "variant"
                ctx.add_well_characterized(
                    f"model cell lines ({len(mutant_models)} with mutation)",
                    "DepMap CCLE",
                    category=GapCategory.PRECLINICAL,
                    matches_on=f"{len(mutant_models)} variant"
                )
        else:
            ctx.add_well_characterized(
                f"model cell lines ({n_models} available)",
                "DepMap CCLE",
                category=GapCategory.PRECLINICAL
            )
    else:
        if ctx.has_drug_data or ctx.has_clinical or evidence.context.gene_role:
            ctx.add_gap(
                category=GapCategory.PRECLINICAL,
                severity=GapSeverity.MINOR,
                description=f"No cell line models identified for {ctx.gene} {ctx.variant}",
                suggested_studies=["Identify cell lines with mutation", "Generate isogenic model"],
                addressable_with=["DepMap CCLE", "Cellosaurus"]
            )
            ctx.add_poorly_characterized("preclinical model systems")


def _check_literature_depth(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for published literature coverage."""
    if not evidence.literature_searched:
        return  # Don't report gap if user chose not to search

    pub_count = len(evidence.pubmed_articles)

    if pub_count == 0:
        severity = GapSeverity.CRITICAL if ctx.is_cancer_gene else GapSeverity.SIGNIFICANT
        ctx.add_gap(
            category=GapCategory.FUNCTIONAL,
            severity=severity,
            description=f"No published literature found for {ctx.gene} {ctx.variant}",
            suggested_studies=["Case report", "Functional characterization study"],
            addressable_with=["PubMed", "Semantic Scholar", "bioRxiv"]
        )
        ctx.add_poorly_characterized("published literature")
    elif pub_count < 5:
        ctx.add_poorly_characterized("literature depth (limited publications)")
    else:
        ctx.add_well_characterized(
            "published literature",
            f"{pub_count} PubMed articles",
            category=GapCategory.FUNCTIONAL
        )


def _check_validation_gap(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Check for strong oncogenic signal with limited therapeutic validation."""
    has_strong_oncogenic_signal = (
        ctx.has_pathogenic_signal and
        evidence.depmap_evidence is not None and
        evidence.depmap_evidence.is_essential()
    )

    has_therapeutic_validation = (
        bool(evidence.civic_assertions) or
        bool(evidence.fda_approvals) or
        bool(evidence.vicc_evidence)
    )

    if has_strong_oncogenic_signal:
        ctx.add_well_characterized(
            "biological driver potential",
            "Pathogenic prediction + gene essentiality",
            category=GapCategory.VALIDATION
        )

        if not has_therapeutic_validation:
            ctx.add_gap(
                category=GapCategory.VALIDATION,
                severity=GapSeverity.CRITICAL if ctx.is_cancer_gene else GapSeverity.SIGNIFICANT,
                description="Strong oncogenic signal but limited therapeutic validation",
                suggested_studies=[
                    "Functional validation in isogenic models",
                    "Drug sensitivity screening (PRISM/GDSC)",
                    "Patient-derived organoid testing",
                    "Structural modeling of activation mechanism"
                ],
                addressable_with=["AlphaFold", "DepMap", "CRISPR screens", "Literature"]
            )
            ctx.add_poorly_characterized("therapeutic validation")


# =============================================================================
# CONTEXT-AWARE SUGGESTION ENRICHMENT
# =============================================================================

def _enrich_gaps_with_context(evidence: "Evidence", ctx: GapDetectionContext) -> None:
    """Enrich gap suggestions with dynamic, context-aware recommendations.

    Analyzes actual evidence to generate specific, actionable study suggestions
    based on gene, variant, tumor type, and available data.
    """
    gene = ctx.gene
    variant = ctx.variant
    tumor_type = ctx.tumor_type

    # Get primary approved drug if available
    primary_drug = _get_primary_drug(evidence)

    # Get top sensitive drugs from DepMap if available (only if tumor-matched cell lines exist)
    top_sensitive_drugs = _get_top_sensitive_drugs(evidence, tumor_type)

    # Get top co-occurring gene if available
    top_cooc_gene = _get_top_cooccurring_gene(evidence)

    # Check for strong co-occurrence signal
    has_strong_cooc = _has_strong_cooccurrence(evidence)

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
            if evidence.civic_assertions or evidence.fda_approvals:
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
            elif evidence.context.gene_role in ("TSG", "tumor_suppressor"):
                new_suggestions.append(
                    f"LOF assay: assess {gene} {variant} impact on tumor suppressor function"
                )

        # Add co-dependency suggestions if relevant
        if gap.category == GapCategory.VALIDATION:
            co_deps = _get_top_codependencies(evidence)
            if co_deps:
                deps_str = ", ".join(co_deps[:2])
                new_suggestions.append(
                    f"Test synthetic lethality with co-essential genes: {deps_str}"
                )

        # Append new suggestions to existing ones (avoid duplicates)
        for suggestion in new_suggestions:
            if suggestion not in gap.suggested_studies:
                gap.suggested_studies.append(suggestion)


def _get_primary_drug(evidence: "Evidence") -> str | None:
    """Get the primary approved drug for this variant."""
    if evidence.fda_approvals:
        approval = evidence.fda_approvals[0]
        return approval.brand_name or approval.generic_name or approval.drug_name
    # Fall back to CIViC assertion therapies
    for assertion in evidence.civic_assertions:
        if assertion.therapies:
            return assertion.therapies[0] if isinstance(assertion.therapies, list) else assertion.therapies
    return None


def _get_top_sensitive_drugs(evidence: "Evidence", tumor_type: str | None = None) -> list[str]:
    """Get top sensitive drugs from DepMap data.

    Only returns drugs if there are tumor-matched cell lines with the mutation,
    since sensitivity data from unrelated tumor types is not applicable.
    """
    if not evidence.depmap_evidence or not evidence.depmap_evidence.drug_sensitivities:
        return []

    # Check for tumor-matched cell lines
    if tumor_type and evidence.depmap_evidence.cell_line_models:
        mutant_models = [cl for cl in evidence.depmap_evidence.cell_line_models if cl.has_mutation]
        tumor_matched = [
            m for m in mutant_models
            if m.primary_disease and _tumor_type_matches(tumor_type, m.primary_disease)
        ]
        if not tumor_matched:
            # No tumor-matched cell lines - don't suggest drugs from unrelated tumor types
            return []

    top_drugs = evidence.depmap_evidence.get_top_sensitive_drugs(5)
    return [ds.drug_name for ds in top_drugs]


def _get_top_cooccurring_gene(evidence: "Evidence") -> str | None:
    """Get the top co-occurring gene from cBioPortal data."""
    if evidence.cbioportal_evidence and evidence.cbioportal_evidence.co_occurring:
        return evidence.cbioportal_evidence.co_occurring[0].gene
    return None


def _has_strong_cooccurrence(evidence: "Evidence", threshold_pct: float = 20.0) -> bool:
    """Check if there's a strong co-occurrence signal (>threshold% co-mutation rate)."""
    if evidence.cbioportal_evidence and evidence.cbioportal_evidence.co_occurring:
        top_cooc = evidence.cbioportal_evidence.co_occurring[0]
        return top_cooc.pct >= threshold_pct
    return False


def _get_top_codependencies(evidence: "Evidence") -> list[str]:
    """Get top co-dependent genes from DepMap data."""
    if evidence.depmap_evidence and evidence.depmap_evidence.co_dependencies:
        return [cd.gene for cd in evidence.depmap_evidence.co_dependencies[:3]]
    return []


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Define category ordering for sorting well-characterized aspects
_CATEGORY_ORDER: dict[GapCategory | None, int] = {
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


def _sort_characterized_by_category(
    items: list[CharacterizedAspect]
) -> list[CharacterizedAspect]:
    """Sort well-characterized aspects by category for grouped display.

    Args:
        items: List of CharacterizedAspect items

    Returns:
        Sorted list with items grouped by category
    """
    return sorted(items, key=lambda x: _CATEGORY_ORDER.get(x.category, 99))


def _check_tumor_specific_evidence(evidence: "Evidence", tumor_type: str) -> bool:
    """Check if any evidence is specific to this tumor type."""
    tumor_lower = tumor_type.lower()

    # Check CIViC assertions
    for assertion in evidence.civic_assertions:
        if assertion.disease and tumor_lower in assertion.disease.lower():
            return True

    # Check CIViC evidence items
    for civic in evidence.civic_evidence:
        if civic.disease and civic.disease_match:
            # disease_match is True when disease matches the queried tumor type
            return True

    # Check FDA approvals
    for approval in evidence.fda_approvals:
        if approval.indication:
            parsed = approval.parse_indication_for_tumor(tumor_type)
            if parsed.get('tumor_match'):
                return True

    # Check VICC evidence
    for vicc in evidence.vicc_evidence:
        if vicc.disease and tumor_lower in vicc.disease.lower():
            return True

    # Check CGI biomarkers
    for cgi in evidence.cgi_biomarkers:
        if cgi.tumor_type and tumor_lower in cgi.tumor_type.lower():
            return True

    return False


def _has_pathogenic_signal(evidence: "Evidence") -> bool:
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
        "p", "pathogenic", "likely_pathogenic", "lp"
    ):
        return True

    # CADD score >= 20 suggests deleteriousness
    if func.cadd_score is not None and func.cadd_score >= 20:
        return True

    # PolyPhen2 damaging prediction
    if func.polyphen2_prediction and func.polyphen2_prediction.lower() in (
        "d", "damaging", "probably_damaging", "possibly_damaging"
    ):
        return True

    # SIFT deleterious prediction (D = deleterious, T = tolerated)
    if func.sift_prediction and func.sift_prediction.lower() in (
        "d", "deleterious"
    ):
        return True

    # Has clinical evidence (strongest signal)
    if evidence.civic_assertions or evidence.fda_approvals:
        return True

    # ClinVar pathogenic entries
    for entry in evidence.clinvar_entries:
        if entry.clinical_significance and "pathogenic" in entry.clinical_significance.lower():
            return True

    # Check if overall ClinVar significance is pathogenic
    if evidence.clinvar_significance and "pathogenic" in evidence.clinvar_significance.lower():
        return True

    # Truncating variants (nonsense, frameshift) are generally pathogenic
    if func.snpeff_effect:
        effect_lower = func.snpeff_effect.lower()
        if any(term in effect_lower for term in [
            "stop_gained", "frameshift", "splice_donor", "splice_acceptor",
            "start_lost", "nonsense"
        ]):
            return True

    return False


def _normalize_source(source: str) -> str:
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


def _detect_discordant_evidence_internal(evidence: "Evidence") -> list[str]:
    """Detect conflicting evidence between different sources.

    Only flags TRUE cross-source conflicts (e.g., CGI says sensitive, CIViC says resistant).
    Intra-source conflicts (multiple CIViC entries disagreeing) are not flagged as they
    often represent context-dependent responses (different tumor types, combinations, etc.)
    rather than genuine discordance.

    Returns list of human-readable conflict descriptions.
    """
    conflicts: list[str] = []

    # Collect drug response signals from different sources
    sensitive_drugs: dict[str, set[str]] = {}
    resistant_drugs: dict[str, set[str]] = {}

    # Check FDA approvals (always sensitivity)
    for approval in evidence.fda_approvals:
        drug = approval.generic_name or approval.brand_name or approval.drug_name
        if drug:
            sensitive_drugs.setdefault(drug.lower(), set()).add("FDA")

    # Check CGI biomarkers (FDA-approved)
    for cgi in evidence.cgi_biomarkers:
        if not cgi.drug:
            continue
        drug = cgi.drug.lower()
        if cgi.association:
            assoc_upper = cgi.association.upper()
            if "RESIST" in assoc_upper:
                resistant_drugs.setdefault(drug, set()).add("CGI")
            elif "SENS" in assoc_upper or "RESPON" in assoc_upper:
                sensitive_drugs.setdefault(drug, set()).add("CGI")

    # Check CGI preclinical biomarkers
    for cgi in evidence.preclinical_biomarkers:
        if not cgi.drug:
            continue
        drug = cgi.drug.lower()
        if cgi.association:
            assoc_upper = cgi.association.upper()
            if "RESIST" in assoc_upper:
                resistant_drugs.setdefault(drug, set()).add("CGI (preclinical)")
            elif "SENS" in assoc_upper or "RESPON" in assoc_upper:
                sensitive_drugs.setdefault(drug, set()).add("CGI (preclinical)")

    # Check CGI early-phase biomarkers
    for cgi in evidence.early_phase_biomarkers:
        if not cgi.drug:
            continue
        drug = cgi.drug.lower()
        if cgi.association:
            assoc_upper = cgi.association.upper()
            if "RESIST" in assoc_upper:
                resistant_drugs.setdefault(drug, set()).add("CGI (early phase)")
            elif "SENS" in assoc_upper or "RESPON" in assoc_upper:
                sensitive_drugs.setdefault(drug, set()).add("CGI (early phase)")

    # Check VICC evidence
    for vicc in evidence.vicc_evidence:
        if not vicc.drugs or len(vicc.drugs) > 1:  # Skip combinations
            continue
        drug_lower = vicc.drugs[0].lower()
        if vicc.response_type:
            resp_upper = vicc.response_type.upper()
            source_name = _normalize_source(f"VICC/{vicc.source}" if vicc.source else "VICC")
            if "RESIST" in resp_upper:
                resistant_drugs.setdefault(drug_lower, set()).add(source_name)
            elif "SENS" in resp_upper or "RESPON" in resp_upper:
                sensitive_drugs.setdefault(drug_lower, set()).add(source_name)

    # Check CIViC assertions
    for assertion in evidence.civic_assertions:
        if not assertion.therapies:
            continue
        for therapy in assertion.therapies:
            drug_lower = therapy.lower()
            if assertion.is_resistance:
                resistant_drugs.setdefault(drug_lower, set()).add("CIViC")
            elif assertion.is_sensitivity:
                sensitive_drugs.setdefault(drug_lower, set()).add("CIViC")

    # Check CIViC evidence items
    for civic in evidence.civic_evidence:
        if not civic.drugs or len(civic.drugs) > 1:  # Skip combinations
            continue
        drug_lower = civic.drugs[0].lower()
        if civic.clinical_significance:
            sig_upper = civic.clinical_significance.upper()
            if "RESIST" in sig_upper:
                resistant_drugs.setdefault(drug_lower, set()).add("CIViC")
            elif "SENS" in sig_upper or "RESPON" in sig_upper:
                sensitive_drugs.setdefault(drug_lower, set()).add("CIViC")

    # Find drugs with TRUE CROSS-SOURCE conflicts only
    for drug in set(sensitive_drugs.keys()) & set(resistant_drugs.keys()):
        sens_sources = sensitive_drugs[drug]
        resist_sources = resistant_drugs[drug]

        # Only flag if sources are truly different (cross-source conflict)
        if sens_sources == resist_sources:
            continue

        sens_only = sens_sources - resist_sources
        resist_only = resist_sources - sens_sources

        if sens_only and resist_only:
            conflicts.append(
                f"Conflicting drug response for {drug}: "
                f"sensitive ({', '.join(sorted(sens_sources))}) vs "
                f"resistant ({', '.join(sorted(resist_sources))})"
            )

    # Check ClinVar significance conflicts
    clinvar_sigs = set()
    for entry in evidence.clinvar_entries:
        if entry.clinical_significance:
            sig = entry.clinical_significance.lower()
            if "pathogenic" in sig and "benign" not in sig:
                clinvar_sigs.add("pathogenic")
            elif "benign" in sig and "pathogenic" not in sig:
                clinvar_sigs.add("benign")

    if "pathogenic" in clinvar_sigs and "benign" in clinvar_sigs:
        conflicts.append(
            "ClinVar has conflicting interpretations: both pathogenic and benign submissions"
        )

    return conflicts


# =============================================================================
# SCORING FUNCTIONS
# =============================================================================

# Research-oriented gap weights: biological unknowns weighted higher than clinical gaps
GAP_CATEGORY_WEIGHTS: dict[GapCategory, float] = {
    GapCategory.VALIDATION: 3.5,      # Strong signal + no validation = high research value
    GapCategory.FUNCTIONAL: 3.0,      # Mechanism unknown
    GapCategory.PRECLINICAL: 2.5,     # No models to test hypotheses
    GapCategory.RESISTANCE: 2.0,      # Resistance mechanisms unknown
    GapCategory.DISCORDANT: 2.0,      # Conflicting evidence needs resolution
    GapCategory.DRUG_RESPONSE: 1.5,   # Drug sensitivity unknown
    GapCategory.TUMOR_TYPE: 1.5,      # Not studied in this tumor
    GapCategory.PREVALENCE: 1.0,      # Epidemiology unknown
    GapCategory.CLINICAL: 1.0,        # Lower weight for research context
    GapCategory.PROGNOSTIC: 1.0,      # Prognostic impact unknown
}

SEVERITY_MULTIPLIERS: dict[GapSeverity, float] = {
    GapSeverity.CRITICAL: 3.0,
    GapSeverity.SIGNIFICANT: 2.0,
    GapSeverity.MINOR: 1.0,
}


def _compute_overall_quality(gaps: list[EvidenceGap], well_characterized_count: int) -> str:
    """Compute overall evidence quality using net scoring (gaps vs well-characterized).

    A variant with many well-characterized aspects and few gaps scores better than
    one with few well-characterized aspects and the same gaps.

    Args:
        gaps: List of evidence gaps found
        well_characterized_count: Number of well-characterized aspects

    Returns:
        Quality rating: "comprehensive" | "moderate" | "limited" | "minimal"
    """
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

    # Apply thresholds to net score
    if net_score >= 12.0:
        return "minimal"
    elif net_score >= 6.0:
        return "limited"
    elif net_score >= 0.0:
        return "moderate"
    else:
        return "comprehensive"


def _compute_research_priority(
    evidence: "Evidence",
    gaps: list[EvidenceGap],
    overall_quality: str,
) -> str:
    """Compute research priority based on gene importance and gap profile.

    Args:
        evidence: The aggregated evidence
        gaps: List of identified gaps
        overall_quality: The computed overall quality ("comprehensive", "moderate", etc.)

    Returns: "very_high" | "high" | "medium" | "low"
    """

    # Count gaps first
    critical_count = sum(1 for g in gaps if g.severity == GapSeverity.CRITICAL)
    significant_count = sum(1 for g in gaps if g.severity == GapSeverity.SIGNIFICANT)

    # Only return low if comprehensive AND no significant/critical gaps
    if overall_quality == "comprehensive" and critical_count == 0 and significant_count == 0:
        return "low"

    gene = evidence.identifiers.gene
    variant = evidence.identifiers.variant

    is_cancer_gene = evidence.context.gene_role in (
        "oncogene", "TSG", "tumor_suppressor", "ddr", "tsg_pathway_actionable"
    )

    has_strong_oncogenic_signal = (
        _has_pathogenic_signal(evidence) and
        evidence.depmap_evidence is not None and
        evidence.depmap_evidence.is_essential()
    )

    has_biological_gaps = any(
        g.category in (GapCategory.VALIDATION, GapCategory.FUNCTIONAL, GapCategory.PRECLINICAL)
        for g in gaps
    )

    is_adjacent, _ = is_hotspot_adjacent(gene, variant, window=5)

    # Very high: strong oncogenic signal + biological gaps = prime research target
    if has_strong_oncogenic_signal and has_biological_gaps:
        return "very_high"

    # Very high: hotspot-adjacent variant with pathogenic signal
    if is_adjacent and _has_pathogenic_signal(evidence) and has_biological_gaps:
        return "very_high"

    # High: cancer gene with critical gaps
    if is_cancer_gene and critical_count > 0:
        return "high"

    # High: hotspot-adjacent variant in cancer gene
    if is_adjacent and is_cancer_gene:
        return "high"

    # Medium: any critical gaps OR cancer gene with significant gaps
    if critical_count > 0 or (is_cancer_gene and significant_count > 0):
        return "medium"

    return "low"



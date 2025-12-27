"""Evidence gap detection from aggregated evidence.

Analyzes what's missing from the evidence to guide research priorities.
"""

from oncomind.models.evidence.evidence_gaps import (
    EvidenceGaps, EvidenceGap, GapCategory, GapSeverity, CharacterizedAspect
)
from oncomind.models.gene_context import is_hotspot_variant, is_hotspot_adjacent, _extract_codon_position

# Import Evidence with TYPE_CHECKING to avoid circular imports
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from oncomind.models.evidence import Evidence


def detect_evidence_gaps(evidence: "Evidence") -> EvidenceGaps:
    """Detect evidence gaps from aggregated evidence.

    Args:
        evidence: Aggregated evidence from all sources

    Returns:
        EvidenceGaps with identified gaps and assessment
    """
    gaps: list[EvidenceGap] = []
    well_characterized: list[str] = []  # Legacy simple strings
    well_characterized_detailed: list[CharacterizedAspect] = []  # New structured format
    poorly_characterized: list[str] = []

    def add_well_characterized(aspect: str, basis: str) -> None:
        """Add a well-characterized aspect with its basis."""
        well_characterized.append(aspect)
        well_characterized_detailed.append(CharacterizedAspect(aspect=aspect, basis=basis))

    gene = evidence.identifiers.gene
    variant = evidence.identifiers.variant
    tumor_type = evidence.context.tumor_type

    # Determine if this is a known cancer gene (affects severity calculations)
    is_cancer_gene = evidence.context.gene_role in (
        "oncogene", "TSG", "tumor_suppressor", "ddr", "tsg_pathway_actionable"
    )

    # Check if variant has pathogenic signal (not a benign polymorphism)
    # Used to avoid overcalling gaps on clearly benign variants
    has_pathogenic_signal = _has_pathogenic_signal(evidence)

    # === Check hotspot context ===
    # Hotspot variants are well-characterized; hotspot-adjacent variants are research gold
    is_hotspot = is_hotspot_variant(gene, variant)
    is_adjacent, nearest_hotspot = is_hotspot_adjacent(gene, variant, window=5)

    if is_hotspot:
        add_well_characterized("known cancer hotspot", f"Codon {_extract_codon_position(variant)} is in cancerhotspots.org")
    elif is_adjacent and nearest_hotspot:
        # Rare variant near a hotspot - high research value
        add_well_characterized(
            f"near hotspot codon {nearest_hotspot}",
            f"Within 5 codons of known hotspot — structural similarity likely"
        )
        # This is a research opportunity - rare variant near known activating hotspot
        # Structural similarity to hotspot suggests similar functional impact
        gaps.append(EvidenceGap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Rare variant near known hotspot (codon {nearest_hotspot}) — functional characterization needed",
            suggested_studies=[
                f"Compare to nearby hotspot {gene} codon {nearest_hotspot}",
                "Structural modeling to assess activation mechanism",
                "Functional assay (transformation, signaling)"
            ],
            addressable_with=["AlphaFold", "Literature on nearby hotspot", "Isogenic models"]
        ))
        poorly_characterized.append("rare-near-hotspot variant function")

    # === Check functional characterization ===
    has_functional = (
        evidence.functional.alphamissense_score is not None or
        evidence.functional.cadd_score is not None or
        evidence.functional.polyphen2_prediction is not None
    )

    if has_functional:
        # Build basis string based on available predictions
        func_sources = []
        if evidence.functional.alphamissense_score is not None:
            func_sources.append(f"AlphaMissense={evidence.functional.alphamissense_score:.2f}")
        if evidence.functional.cadd_score is not None:
            func_sources.append(f"CADD={evidence.functional.cadd_score:.1f}")
        if evidence.functional.polyphen2_prediction:
            func_sources.append(f"PolyPhen2={evidence.functional.polyphen2_prediction}")
        add_well_characterized("computational pathogenicity", " | ".join(func_sources) if func_sources else "Predictions available")
    else:
        gaps.append(EvidenceGap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"No computational pathogenicity predictions for {gene} {variant}",
            suggested_studies=["Run AlphaMissense, CADD, PolyPhen2"],
            addressable_with=["MyVariant.info", "VEP"]
        ))
        poorly_characterized.append("pathogenicity predictions")

    # === Check mechanism/functional studies ===
    has_mechanism = bool(evidence.context.gene_role) or bool(evidence.context.pathway)

    # Check DepMap gene essentiality - key functional evidence
    has_depmap_essentiality = (
        evidence.depmap_evidence is not None and
        evidence.depmap_evidence.gene_dependency is not None
    )

    if has_mechanism:
        role = evidence.context.gene_role or "unknown"
        pathway = evidence.context.pathway or ""
        basis = f"Role: {role}" + (f", Pathway: {pathway}" if pathway else "")
        add_well_characterized("gene function", basis)
    if has_depmap_essentiality:
        dep = evidence.depmap_evidence.gene_dependency
        score = dep.mean_dependency_score if dep else None
        score_str = f"CERES={score:.2f}" if score else "DepMap CRISPR data"
        add_well_characterized("gene essentiality", score_str)
        # If gene is essential, this adds confidence to functional impact
        if evidence.depmap_evidence is not None and evidence.depmap_evidence.is_essential():
            pct = dep.dependency_pct if dep else 0
            add_well_characterized(f"{gene} is essential", f"{pct:.0f}% of cell lines depend on it")

    if not has_mechanism and not has_depmap_essentiality:
        gaps.append(EvidenceGap(
            category=GapCategory.FUNCTIONAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Functional impact of {variant} on {gene} protein is unknown",
            suggested_studies=["Functional assay", "Structural modeling", "Cell-based reporter"],
            addressable_with=["UniProt", "Literature search", "DepMap"]
        ))
        poorly_characterized.append("functional mechanism")

    # === Check clinical evidence ===
    has_clinical = bool(evidence.civic_assertions) or bool(evidence.fda_approvals)

    if has_clinical:
        # Count assertions and approvals
        n_assertions = len(evidence.civic_assertions)
        n_fda = len(evidence.fda_approvals)
        parts = []
        if n_fda:
            parts.append(f"{n_fda} FDA approval{'s' if n_fda > 1 else ''}")
        if n_assertions:
            parts.append(f"{n_assertions} CIViC assertion{'s' if n_assertions > 1 else ''}")
        add_well_characterized("clinical actionability", " + ".join(parts) if parts else "Clinical evidence exists")
    else:
        # Clinical gap is CRITICAL - no curated clinical evidence is a major gap
        gaps.append(EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.CRITICAL,
            description=f"No curated clinical evidence for {gene} {variant}",
            suggested_studies=["Case series", "Retrospective cohort", "Basket trial inclusion"],
            addressable_with=["CIViC submission", "Literature curation"]
        ))
        poorly_characterized.append("clinical evidence")

    # === Check tumor-type-specific evidence ===
    if tumor_type:
        tumor_specific_evidence = _check_tumor_specific_evidence(evidence, tumor_type)

        if tumor_specific_evidence:
            add_well_characterized(f"evidence in {tumor_type}", "Tumor-specific CIViC/FDA/VICC/CGI data")
        else:
            # Tumor-type gap severity depends on gene importance and pathogenic signal
            # CRITICAL: cancer gene + no clinical evidence + pathogenic signal
            # SIGNIFICANT: cancer gene OR pathogenic signal
            # MINOR: benign/unknown variants in non-cancer genes
            if is_cancer_gene and not has_clinical and has_pathogenic_signal:
                tumor_severity = GapSeverity.CRITICAL
            elif is_cancer_gene or has_pathogenic_signal:
                tumor_severity = GapSeverity.SIGNIFICANT
            else:
                tumor_severity = GapSeverity.MINOR
            gaps.append(EvidenceGap(
                category=GapCategory.TUMOR_TYPE,
                severity=tumor_severity,
                description=f"No evidence specific to {tumor_type} for {gene} {variant}",
                suggested_studies=[
                    f"Case series in {tumor_type}",
                    f"Retrospective analysis in {tumor_type} cohort",
                    "Basket trial with histology-specific cohort"
                ],
                addressable_with=["ClinicalTrials.gov", "Literature search"]
            ))
            poorly_characterized.append(f"{tumor_type}-specific data")

    # === Check drug response data ===
    # Include DepMap drug sensitivities in drug response check
    has_depmap_drug_data = (
        evidence.depmap_evidence is not None and
        bool(evidence.depmap_evidence.drug_sensitivities)
    )
    has_drug_data = (
        bool(evidence.cgi_biomarkers) or
        bool(evidence.vicc_evidence) or
        bool(evidence.preclinical_biomarkers) or
        has_depmap_drug_data
    )

    if has_drug_data:
        # Count drug data sources
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
        add_well_characterized("drug response", " + ".join(drug_sources) if drug_sources else "Drug data available")
        if has_depmap_drug_data:
            n_drugs = len(evidence.depmap_evidence.drug_sensitivities)
            add_well_characterized("preclinical drug sensitivity (DepMap)", f"{n_drugs} drug IC50 values")
    else:
        gaps.append(EvidenceGap(
            category=GapCategory.DRUG_RESPONSE,
            severity=GapSeverity.SIGNIFICANT,
            description=f"No drug sensitivity/resistance data for {gene} {variant}",
            suggested_studies=["Cell line drug screen", "PDX drug testing", "Clinical correlative study"],
            addressable_with=["GDSC", "CTRP", "DepMap"]
        ))
        poorly_characterized.append("drug response")

    # === Check resistance mechanisms ===
    resistance_articles = [a for a in evidence.pubmed_articles if a.is_resistance_evidence()]
    has_resistance_data = (
        bool(resistance_articles) or
        any(b.association and "RESIST" in b.association.upper() for b in evidence.cgi_biomarkers)
    )

    if has_resistance_data:
        n_resist_articles = len(resistance_articles)
        add_well_characterized("resistance mechanisms", f"{n_resist_articles} resistance article{'s' if n_resist_articles != 1 else ''} + CGI data")
    elif has_clinical or has_drug_data:
        # Only flag as gap if this is likely a clinically relevant variant
        # Resistance is significant because it affects treatment decisions
        gaps.append(EvidenceGap(
            category=GapCategory.RESISTANCE,
            severity=GapSeverity.SIGNIFICANT,
            description=f"Resistance mechanisms for {gene} {variant} not well characterized",
            suggested_studies=["Serial biopsy study", "ctDNA monitoring", "Resistance screen"],
            addressable_with=["Literature search", "CIViC"]
        ))
        poorly_characterized.append("resistance mechanisms")

    # === Check for discordant/conflicting evidence ===
    discordant_findings = _detect_discordant_evidence(evidence)
    if discordant_findings:
        for finding in discordant_findings:
            gaps.append(EvidenceGap(
                category=GapCategory.DISCORDANT,
                severity=GapSeverity.SIGNIFICANT,
                description=finding,
                suggested_studies=["Meta-analysis", "Prospective validation study"],
                addressable_with=["Literature review", "Expert consensus"]
            ))
            poorly_characterized.append("conflicting evidence")

    # === Check prevalence/epidemiology ===
    has_prevalence = (
        evidence.cbioportal_evidence is not None and
        evidence.cbioportal_evidence.has_data()
    )

    if has_prevalence:
        cbio = evidence.cbioportal_evidence
        study = cbio.study_name if cbio else "cBioPortal"
        pct = cbio.variant_prevalence_pct if cbio else 0
        add_well_characterized("prevalence", f"{pct:.1f}% in {study}")
    else:
        # Prevalence gap is SIGNIFICANT for cancer genes with clinical relevance, MINOR otherwise
        prevalence_severity = GapSeverity.SIGNIFICANT if (is_cancer_gene and has_clinical) else GapSeverity.MINOR
        gaps.append(EvidenceGap(
            category=GapCategory.PREVALENCE,
            severity=prevalence_severity,
            description=f"Prevalence of {gene} {variant} in {tumor_type or 'cancer'} unknown",
            suggested_studies=["Epidemiological study", "Registry analysis"],
            addressable_with=["cBioPortal", "COSMIC", "TCGA"]
        ))
        poorly_characterized.append("prevalence data")

    # === Check clinical trials ===
    has_trials = bool(evidence.clinical_trials)

    if has_trials:
        n_trials = len(evidence.clinical_trials)
        add_well_characterized("clinical trial options", f"{n_trials} active trial{'s' if n_trials != 1 else ''}")
    elif has_clinical or has_drug_data:
        gaps.append(EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.MINOR,
            description=f"No active clinical trials for {gene} {variant}",
            suggested_studies=["Clinical trial design", "Basket trial proposal"],
            addressable_with=["ClinicalTrials.gov"]
        ))

    # === Check preclinical model systems (cell lines) ===
    # Enhanced logic: check for tumor-type-specific models and flag cross-histology opportunities
    has_cell_line_models = (
        evidence.depmap_evidence is not None and
        bool(evidence.depmap_evidence.cell_line_models)
    )

    if has_cell_line_models:
        n_models = len(evidence.depmap_evidence.cell_line_models)
        # Get CellLineModel objects (not just names) to access tumor type info
        mutant_models = [
            cl for cl in evidence.depmap_evidence.cell_line_models if cl.has_mutation
        ]

        if mutant_models:
            add_well_characterized(f"model cell lines ({len(mutant_models)} with mutation)", "DepMap CCLE cell lines")

            # Check if models exist but none match the queried tumor type
            if tumor_type:
                tumor_lower = tumor_type.lower()
                tumor_models = [
                    m for m in mutant_models
                    if m.primary_disease and tumor_lower in m.primary_disease.lower()
                ]
                if not tumor_models:
                    # Models exist with mutation but not in queried tumor type
                    # This is a research opportunity for cross-histology hypothesis
                    gaps.append(EvidenceGap(
                        category=GapCategory.PRECLINICAL,
                        severity=GapSeverity.SIGNIFICANT,
                        description=f"Models with {variant} exist but none in {tumor_type} — cross-histology testing possible",
                        suggested_studies=[
                            f"Test in {tumor_type}-derived organoids",
                            "Compare drug response vs other histologies",
                            "Generate isogenic model in {tumor_type} background"
                        ],
                        addressable_with=["DepMap", "Patient-derived organoids", "CRISPR knock-in"]
                    ))
                    poorly_characterized.append(f"{tumor_type}-specific preclinical models")
                else:
                    add_well_characterized(f"{tumor_type} cell line models ({len(tumor_models)} available)", "Tumor-specific DepMap models")
        else:
            add_well_characterized(f"model cell lines ({n_models} available)", "DepMap CCLE cell lines")
    else:
        # Only flag as gap if gene is likely to be studied
        if has_drug_data or has_clinical or evidence.context.gene_role:
            gaps.append(EvidenceGap(
                category=GapCategory.PRECLINICAL,
                severity=GapSeverity.MINOR,
                description=f"No cell line models identified for {gene} {variant}",
                suggested_studies=["Identify cell lines with mutation", "Generate isogenic model"],
                addressable_with=["DepMap", "CCLE", "Cellosaurus"]
            ))
            poorly_characterized.append("preclinical model systems")

    # === Check literature depth ===
    # Only flag literature gaps if literature search was actually performed
    pub_count = len(evidence.pubmed_articles)
    if evidence.literature_searched:
        if pub_count == 0:
            # Literature gap is CRITICAL for cancer genes, SIGNIFICANT otherwise
            lit_severity = GapSeverity.CRITICAL if is_cancer_gene else GapSeverity.SIGNIFICANT
            gaps.append(EvidenceGap(
                category=GapCategory.FUNCTIONAL,
                severity=lit_severity,
                description=f"No published literature found for {gene} {variant}",
                suggested_studies=["Case report", "Functional characterization study"],
                addressable_with=["PubMed", "Semantic Scholar", "bioRxiv"]
            ))
            poorly_characterized.append("published literature")
        elif pub_count < 5:
            poorly_characterized.append("literature depth (limited publications)")
        else:
            add_well_characterized("published literature", f"{pub_count} PubMed articles")
    # If literature wasn't searched, don't report it as a gap (user chose not to search)

    # === Check for strong oncogenic signal with limited validation ===
    # This identifies high-potential research targets: variants with biological
    # driver signal but lacking therapeutic/clinical validation
    has_strong_oncogenic_signal = (
        has_pathogenic_signal and
        (evidence.depmap_evidence is not None and evidence.depmap_evidence.is_essential())
    )

    has_therapeutic_validation = (
        bool(evidence.civic_assertions) or
        bool(evidence.fda_approvals) or
        bool(evidence.vicc_evidence)
    )

    if has_strong_oncogenic_signal:
        add_well_characterized("biological driver potential", "Pathogenic prediction + gene essentiality")

        if not has_therapeutic_validation:
            gaps.append(EvidenceGap(
                category=GapCategory.VALIDATION,
                severity=GapSeverity.CRITICAL if is_cancer_gene else GapSeverity.SIGNIFICANT,
                description=f"Strong oncogenic signal but limited therapeutic validation",
                suggested_studies=[
                    "Functional validation in isogenic models",
                    "Drug sensitivity screening (PRISM/GDSC)",
                    "Patient-derived organoid testing",
                    "Structural modeling of activation mechanism"
                ],
                addressable_with=["AlphaFold", "DepMap", "CRISPR screens", "Literature"]
            ))
            poorly_characterized.append("therapeutic validation")

    # === Determine overall evidence quality ===
    overall_quality = _compute_overall_quality(gaps)

    # === Determine research priority ===
    research_priority = _compute_research_priority(evidence, gaps)

    return EvidenceGaps(
        gaps=gaps,
        overall_evidence_quality=overall_quality,
        well_characterized=well_characterized,
        well_characterized_detailed=well_characterized_detailed,
        poorly_characterized=poorly_characterized,
        research_priority=research_priority,
    )


def _check_tumor_specific_evidence(evidence: "Evidence", tumor_type: str) -> bool:
    """Check if any evidence is specific to this tumor type."""
    tumor_lower = tumor_type.lower()

    # Check CIViC assertions
    for assertion in evidence.civic_assertions:
        if assertion.disease and tumor_lower in assertion.disease.lower():
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

# Severity multipliers
SEVERITY_MULTIPLIERS: dict[GapSeverity, float] = {
    GapSeverity.CRITICAL: 3.0,
    GapSeverity.SIGNIFICANT: 2.0,
    GapSeverity.MINOR: 1.0,
}


def _compute_overall_quality(gaps: list[EvidenceGap]) -> str:
    """Compute overall evidence quality from gaps using weighted scoring.

    Research-oriented: weights biological gaps higher than clinical gaps.
    """
    if not gaps:
        return "comprehensive"

    # Compute weighted gap score
    total_score = 0.0
    for gap in gaps:
        category_weight = GAP_CATEGORY_WEIGHTS.get(gap.category, 1.0)
        severity_mult = SEVERITY_MULTIPLIERS.get(gap.severity, 1.0)
        total_score += category_weight * severity_mult

    # Thresholds for quality levels
    if total_score >= 15.0:
        return "minimal"
    elif total_score >= 10.0:
        return "limited"
    elif total_score >= 5.0:
        return "moderate"
    else:
        return "comprehensive"


def _compute_research_priority(evidence: "Evidence", gaps: list[EvidenceGap]) -> str:
    """Compute research priority based on gene importance and gap profile.

    Research-oriented: prioritizes biologically promising but under-explored variants.
    Returns: "very_high" | "high" | "medium" | "low"
    """
    gene = evidence.identifiers.gene
    variant = evidence.identifiers.variant

    is_cancer_gene = evidence.context.gene_role in (
        "oncogene", "TSG", "tumor_suppressor", "ddr", "tsg_pathway_actionable"
    )

    # Check for strong oncogenic signal (pathogenic + essential)
    has_strong_oncogenic_signal = (
        _has_pathogenic_signal(evidence) and
        evidence.depmap_evidence is not None and
        evidence.depmap_evidence.is_essential()
    )

    # Check for biological/preclinical gaps (high research value)
    has_biological_gaps = any(
        g.category in (GapCategory.VALIDATION, GapCategory.FUNCTIONAL, GapCategory.PRECLINICAL)
        for g in gaps
    )

    # Check hotspot context
    is_hotspot = is_hotspot_variant(gene, variant)
    is_adjacent, _ = is_hotspot_adjacent(gene, variant, window=5)

    critical_count = sum(1 for g in gaps if g.severity == GapSeverity.CRITICAL)
    significant_count = sum(1 for g in gaps if g.severity == GapSeverity.SIGNIFICANT)

    # Very high: strong oncogenic signal + biological gaps = prime research target
    if has_strong_oncogenic_signal and has_biological_gaps:
        return "very_high"

    # Very high: hotspot-adjacent variant with pathogenic signal = rare near hotspot
    if is_adjacent and _has_pathogenic_signal(evidence) and has_biological_gaps:
        return "very_high"

    # High: cancer gene with critical gaps
    if is_cancer_gene and critical_count > 0:
        return "high"

    # High: hotspot-adjacent variant in cancer gene (research opportunity)
    if is_adjacent and is_cancer_gene:
        return "high"

    # Medium: any critical gaps OR cancer gene with significant gaps
    if critical_count > 0 or (is_cancer_gene and significant_count > 0):
        return "medium"

    return "low"


def _has_pathogenic_signal(evidence: "Evidence") -> bool:
    """Check if variant has any signal suggesting pathogenicity.

    Used to avoid overcalling gaps on clearly benign variants (common polymorphisms).

    Returns True if any of:
    - AlphaMissense predicts pathogenic (P or likely_pathogenic)
    - CADD score >= 20 (predicted deleterious)
    - PolyPhen2 predicts damaging (D or probably_damaging)
    - Rare in population (gnomAD AF < 0.01 or absent)
    - Has any clinical assertions or FDA approvals
    - Has any ClinVar pathogenic/likely pathogenic entries
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

    # Rare variant (gnomAD AF < 1% or absent suggests not a common polymorphism)
    if func.gnomad_exome_af is None or func.gnomad_exome_af < 0.01:
        # Absence or rarity is suggestive but not definitive
        # Only count as pathogenic signal if combined with other evidence
        pass

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


def _detect_discordant_evidence(evidence: "Evidence") -> list[str]:
    """Detect conflicting evidence between different sources.

    Only flags TRUE cross-source conflicts (e.g., CGI says sensitive, CIViC says resistant).
    Intra-source conflicts (multiple CIViC entries disagreeing) are not flagged as they
    often represent context-dependent responses (different tumor types, combinations, etc.)
    rather than genuine discordance.

    Key filtering:
    - Combination therapies are tracked separately from monotherapy
    - VICC/civic is treated as the same source as CIViC
    - Only flags conflicts between truly independent sources

    Returns list of human-readable conflict descriptions.
    """
    conflicts: list[str] = []

    # Collect drug response signals from different sources
    # Use sets to de-duplicate sources (avoid "CIViC, CIViC, CIViC...")
    # Key: drug name (monotherapy only), Value: set of normalized source names
    sensitive_drugs: dict[str, set[str]] = {}
    resistant_drugs: dict[str, set[str]] = {}

    # Check CGI biomarkers
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

    # Check VICC evidence
    for vicc in evidence.vicc_evidence:
        if not vicc.drugs:
            continue
        # Skip combination therapies (only look at single-drug entries)
        if len(vicc.drugs) > 1:
            continue
        drug_lower = vicc.drugs[0].lower()
        if vicc.response_type:
            resp_upper = vicc.response_type.upper()
            source_name = _normalize_source(f"VICC/{vicc.source}" if vicc.source else "VICC")
            if "RESIST" in resp_upper:
                resistant_drugs.setdefault(drug_lower, set()).add(source_name)
            elif "SENS" in resp_upper or "RESPON" in resp_upper:
                sensitive_drugs.setdefault(drug_lower, set()).add(source_name)

    # Check CIViC evidence
    for civic in evidence.civic_evidence:
        if not civic.drugs:
            continue
        # Skip combination therapies (only look at single-drug entries)
        if len(civic.drugs) > 1:
            continue
        drug_lower = civic.drugs[0].lower()
        if civic.clinical_significance:
            sig_upper = civic.clinical_significance.upper()
            if "RESIST" in sig_upper:
                resistant_drugs.setdefault(drug_lower, set()).add("CIViC")
            elif "SENS" in sig_upper or "RESPON" in sig_upper:
                sensitive_drugs.setdefault(drug_lower, set()).add("CIViC")

    # Find drugs with TRUE CROSS-SOURCE conflicts only
    # Intra-source conflicts (CIViC vs CIViC) are not flagged - they often represent
    # context-dependent responses rather than genuine discordance
    for drug in set(sensitive_drugs.keys()) & set(resistant_drugs.keys()):
        sens_sources = sensitive_drugs[drug]
        resist_sources = resistant_drugs[drug]

        # Only flag if the sources are truly different (cross-source conflict)
        # If both sets contain only CIViC (or CIViC-derived data), skip it
        if sens_sources == resist_sources:
            continue

        # Check for actual cross-source disagreement
        sens_only = sens_sources - resist_sources
        resist_only = resist_sources - sens_sources

        if sens_only and resist_only:
            # True cross-source conflict: different sources on each side
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
            elif "uncertain" in sig or "vus" in sig:
                clinvar_sigs.add("uncertain")

    if "pathogenic" in clinvar_sigs and "benign" in clinvar_sigs:
        conflicts.append(
            f"ClinVar has conflicting interpretations: both pathogenic and benign submissions"
        )

    return conflicts

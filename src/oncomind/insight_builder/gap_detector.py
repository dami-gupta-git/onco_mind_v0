"""Evidence gap detection from aggregated evidence.

Analyzes what's missing from the evidence to guide research priorities.
"""

from oncomind.models.evidence.evidence_gaps import (
    EvidenceGaps, EvidenceGap, GapCategory, GapSeverity
)

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
    well_characterized: list[str] = []
    poorly_characterized: list[str] = []

    gene = evidence.identifiers.gene
    variant = evidence.identifiers.variant
    tumor_type = evidence.context.tumor_type

    # === Check functional characterization ===
    has_functional = (
        evidence.functional.alphamissense_score is not None or
        evidence.functional.cadd_score is not None or
        evidence.functional.polyphen2_prediction is not None
    )

    if has_functional:
        well_characterized.append("computational pathogenicity")
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
        well_characterized.append("gene function")
    if has_depmap_essentiality:
        well_characterized.append("gene essentiality (DepMap CRISPR)")
        # If gene is essential, this adds confidence to functional impact
        if evidence.depmap_evidence.is_essential():
            well_characterized.append(f"{gene} is essential in cancer cells")

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
        well_characterized.append("clinical actionability")
    else:
        gaps.append(EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.SIGNIFICANT,
            description=f"No curated clinical evidence for {gene} {variant}",
            suggested_studies=["Case series", "Retrospective cohort", "Basket trial inclusion"],
            addressable_with=["CIViC submission", "Literature curation"]
        ))
        poorly_characterized.append("clinical evidence")

    # === Check tumor-type-specific evidence ===
    if tumor_type:
        tumor_specific_evidence = _check_tumor_specific_evidence(evidence, tumor_type)

        if tumor_specific_evidence:
            well_characterized.append(f"evidence in {tumor_type}")
        else:
            gaps.append(EvidenceGap(
                category=GapCategory.TUMOR_TYPE,
                severity=GapSeverity.CRITICAL,
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
        well_characterized.append("drug response")
        if has_depmap_drug_data:
            well_characterized.append("preclinical drug sensitivity (DepMap)")
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
        well_characterized.append("resistance mechanisms")
    elif has_clinical or has_drug_data:
        # Only flag as gap if this is likely a clinically relevant variant
        gaps.append(EvidenceGap(
            category=GapCategory.RESISTANCE,
            severity=GapSeverity.MINOR,
            description=f"Resistance mechanisms for {gene} {variant} not well characterized",
            suggested_studies=["Serial biopsy study", "ctDNA monitoring", "Resistance screen"],
            addressable_with=["Literature search", "CIViC"]
        ))
        poorly_characterized.append("resistance mechanisms")

    # === Check prevalence/epidemiology ===
    has_prevalence = (
        evidence.cbioportal_evidence is not None and
        evidence.cbioportal_evidence.has_data()
    )

    if has_prevalence:
        well_characterized.append("prevalence")
    else:
        gaps.append(EvidenceGap(
            category=GapCategory.PREVALENCE,
            severity=GapSeverity.MINOR,
            description=f"Prevalence of {gene} {variant} in {tumor_type or 'cancer'} unknown",
            suggested_studies=["Epidemiological study", "Registry analysis"],
            addressable_with=["cBioPortal", "COSMIC", "TCGA"]
        ))
        poorly_characterized.append("prevalence data")

    # === Check clinical trials ===
    has_trials = bool(evidence.clinical_trials)

    if has_trials:
        well_characterized.append("clinical trial options")
    elif has_clinical or has_drug_data:
        gaps.append(EvidenceGap(
            category=GapCategory.CLINICAL,
            severity=GapSeverity.MINOR,
            description=f"No active clinical trials for {gene} {variant}",
            suggested_studies=["Clinical trial design", "Basket trial proposal"],
            addressable_with=["ClinicalTrials.gov"]
        ))

    # === Check preclinical model systems (cell lines) ===
    has_cell_line_models = (
        evidence.depmap_evidence is not None and
        bool(evidence.depmap_evidence.cell_line_models)
    )

    if has_cell_line_models:
        n_models = len(evidence.depmap_evidence.cell_line_models)
        mutant_models = evidence.depmap_evidence.get_model_cell_lines(with_mutation_only=True)
        if mutant_models:
            well_characterized.append(f"model cell lines ({len(mutant_models)} with mutation)")
        else:
            well_characterized.append(f"model cell lines ({n_models} available)")
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
            gaps.append(EvidenceGap(
                category=GapCategory.FUNCTIONAL,
                severity=GapSeverity.CRITICAL,
                description=f"No published literature found for {gene} {variant}",
                suggested_studies=["Case report", "Functional characterization study"],
                addressable_with=["PubMed", "Semantic Scholar", "bioRxiv"]
            ))
            poorly_characterized.append("published literature")
        elif pub_count < 5:
            poorly_characterized.append("literature depth (limited publications)")
        else:
            well_characterized.append("published literature")
    # If literature wasn't searched, don't report it as a gap (user chose not to search)

    # === Determine overall evidence quality ===
    overall_quality = _compute_overall_quality(gaps)

    # === Determine research priority ===
    research_priority = _compute_research_priority(evidence, gaps)

    return EvidenceGaps(
        gaps=gaps,
        overall_evidence_quality=overall_quality,
        well_characterized=well_characterized,
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


def _compute_overall_quality(gaps: list[EvidenceGap]) -> str:
    """Compute overall evidence quality from gaps."""
    critical_count = sum(1 for g in gaps if g.severity == GapSeverity.CRITICAL)
    significant_count = sum(1 for g in gaps if g.severity == GapSeverity.SIGNIFICANT)

    if critical_count >= 2:
        return "minimal"
    elif critical_count == 1:
        return "limited"
    elif significant_count >= 2:
        return "limited"
    elif significant_count == 1:
        return "moderate"
    elif len(gaps) == 0:
        return "comprehensive"
    else:
        return "moderate"


def _compute_research_priority(evidence: "Evidence", gaps: list[EvidenceGap]) -> str:
    """Compute research priority based on gene importance and gaps."""
    critical_count = sum(1 for g in gaps if g.severity == GapSeverity.CRITICAL)
    significant_count = sum(1 for g in gaps if g.severity == GapSeverity.SIGNIFICANT)

    # High priority: clinically important gene with critical gaps
    is_cancer_gene = evidence.context.gene_role in (
        "oncogene", "TSG", "tumor_suppressor", "ddr", "tsg_pathway_actionable"
    )

    if is_cancer_gene and critical_count > 0:
        return "high"
    elif critical_count > 0 or (is_cancer_gene and significant_count > 0):
        return "medium"
    else:
        return "low"

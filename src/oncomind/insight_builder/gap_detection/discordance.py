"""Discordant evidence detection for gap analysis.

Detects conflicting evidence between different sources (e.g., CIViC says sensitive,
CGI says resistant). This is separated from checks.py due to its complexity and
distinct purpose.
"""

from typing import TYPE_CHECKING

from oncomind.models.extracted.evidence_gaps import GapCategory, GapSeverity
from .helpers import normalize_source

if TYPE_CHECKING:
    from oncomind.models.evidence import Evidence
    from .context import GapDetectionContext


def check_discordant_evidence(evidence: "Evidence", ctx: "GapDetectionContext") -> None:
    """Check for conflicting evidence between sources."""
    high_conflicts, significant_conflicts = detect_discordant_evidence_internal(evidence)

    # HIGH severity conflicts (cross-source drug response, ClinVar conflicts)
    for finding in high_conflicts:
        ctx.add_gap(
            category=GapCategory.DISCORDANT,
            severity=GapSeverity.HIGH,  # Conflicting data deserves high urgency
            description=finding,
            suggested_studies=["Meta-analysis", "Prospective validation study"],
            addressable_with=["Literature review", "Expert consensus"]
        )
        ctx.add_poorly_characterized("conflicting evidence")

    # SIGNIFICANT severity conflicts (FDA vs VICC at variant level)
    for finding in significant_conflicts:
        ctx.add_gap(
            category=GapCategory.DISCORDANT,
            severity=GapSeverity.SIGNIFICANT,  # VICC is less reliable than FDA
            description=finding,
            suggested_studies=["Independent validation study", "Case series analysis"],
            addressable_with=["Literature review", "Clinical correlative data"]
        )
        ctx.add_poorly_characterized("conflicting evidence")


def detect_discordant_evidence_internal(evidence: "Evidence") -> tuple[list[str], list[str]]:
    """Detect conflicting evidence between different sources.

    Only flags TRUE cross-source conflicts (e.g., CGI says sensitive, CIViC says resistant).
    Intra-source conflicts (multiple CIViC entries disagreeing) are not flagged as they
    often represent context-dependent responses (different tumor types, combinations, etc.)
    rather than genuine discordance.

    Returns:
        Tuple of (high_severity_conflicts, significant_severity_conflicts)
        - high: Cross-source drug conflicts, ClinVar internal/external conflicts
        - significant: FDA sensitive vs VICC resistant at variant level (VICC less reliable)
    """
    high_conflicts: list[str] = []
    significant_conflicts: list[str] = []

    # Collect drug response signals from different sources
    # We track BOTH all evidence AND variant-level-only evidence
    # Cross-source conflicts are only flagged when BOTH sides have variant-level evidence
    sensitive_drugs: dict[str, set[str]] = {}  # All evidence (any locus level)
    resistant_drugs: dict[str, set[str]] = {}  # All evidence (any locus level)
    sensitive_variant_level: dict[str, set[str]] = {}  # Variant-level only
    resistant_variant_level: dict[str, set[str]] = {}  # Variant-level only

    # Track FDA evidence at variant level for FDA vs VICC conflict detection
    fda_sensitive_variant_level: dict[str, bool] = {}

    # Check FDA biomarker evidence (only REQUIRED_POSITIVE, always sensitivity)
    from oncomind.models.evidence.fda_biomarker import BiomarkerRequirement
    for ev in evidence.fda_biomarker_evidence:
        if ev.requirement == BiomarkerRequirement.REQUIRED_NEGATIVE:
            continue
        if not ev.drug_name or ev.drug_name == "Unknown":
            continue
        drug_lower = ev.drug_name.lower()
        sensitive_drugs.setdefault(drug_lower, set()).add("FDA")
        # Check match level using locus_match property (consistent with other evidence types)
        if ev.locus_match == "variant":
            fda_sensitive_variant_level[drug_lower] = True
            sensitive_variant_level.setdefault(drug_lower, set()).add("FDA")

    # Check CGI biomarkers (FDA-approved)
    for cgi in evidence.cgi_biomarkers:
        if not cgi.drug:
            continue
        drug = cgi.drug.lower()
        if cgi.association:
            assoc_upper = cgi.association.upper()
            if "RESIST" in assoc_upper:
                resistant_drugs.setdefault(drug, set()).add("CGI")
                if cgi.locus_match == "variant":
                    resistant_variant_level.setdefault(drug, set()).add("CGI")
            elif "SENS" in assoc_upper or "RESPON" in assoc_upper:
                sensitive_drugs.setdefault(drug, set()).add("CGI")
                if cgi.locus_match == "variant":
                    sensitive_variant_level.setdefault(drug, set()).add("CGI")

    # Check CGI preclinical biomarkers
    for cgi in evidence.preclinical_biomarkers:
        if not cgi.drug:
            continue
        drug = cgi.drug.lower()
        if cgi.association:
            assoc_upper = cgi.association.upper()
            if "RESIST" in assoc_upper:
                resistant_drugs.setdefault(drug, set()).add("CGI (preclinical)")
                if cgi.locus_match == "variant":
                    resistant_variant_level.setdefault(drug, set()).add("CGI (preclinical)")
            elif "SENS" in assoc_upper or "RESPON" in assoc_upper:
                sensitive_drugs.setdefault(drug, set()).add("CGI (preclinical)")
                if cgi.locus_match == "variant":
                    sensitive_variant_level.setdefault(drug, set()).add("CGI (preclinical)")

    # Check CGI early-phase biomarkers
    for cgi in evidence.early_phase_biomarkers:
        if not cgi.drug:
            continue
        drug = cgi.drug.lower()
        if cgi.association:
            assoc_upper = cgi.association.upper()
            if "RESIST" in assoc_upper:
                resistant_drugs.setdefault(drug, set()).add("CGI (early phase)")
                if cgi.locus_match == "variant":
                    resistant_variant_level.setdefault(drug, set()).add("CGI (early phase)")
            elif "SENS" in assoc_upper or "RESPON" in assoc_upper:
                sensitive_drugs.setdefault(drug, set()).add("CGI (early phase)")
                if cgi.locus_match == "variant":
                    sensitive_variant_level.setdefault(drug, set()).add("CGI (early phase)")

    # Track VICC resistance at variant level for FDA vs VICC conflict detection
    vicc_resistant_variant_level: dict[str, bool] = {}

    # Check VICC evidence
    for vicc in evidence.vicc_evidence:
        if not vicc.drugs or len(vicc.drugs) > 1:  # Skip combinations
            continue
        drug_lower = vicc.drugs[0].lower()
        if vicc.response_type:
            resp_upper = vicc.response_type.upper()
            source_name = normalize_source(f"VICC/{vicc.source}" if vicc.source else "VICC")
            if "RESIST" in resp_upper:
                resistant_drugs.setdefault(drug_lower, set()).add(source_name)
                # Track if variant-level
                if vicc.locus_match == "variant":
                    vicc_resistant_variant_level[drug_lower] = True
                    resistant_variant_level.setdefault(drug_lower, set()).add(source_name)
            elif "SENS" in resp_upper or "RESPON" in resp_upper:
                sensitive_drugs.setdefault(drug_lower, set()).add(source_name)
                if vicc.locus_match == "variant":
                    sensitive_variant_level.setdefault(drug_lower, set()).add(source_name)

    # Check CIViC assertions
    for assertion in evidence.civic_assertions:
        if not assertion.therapies:
            continue
        for therapy in assertion.therapies:
            drug_lower = therapy.lower()
            if assertion.is_resistance:
                resistant_drugs.setdefault(drug_lower, set()).add("CIViC")
                if assertion.locus_match == "variant":
                    resistant_variant_level.setdefault(drug_lower, set()).add("CIViC")
            elif assertion.is_sensitivity:
                sensitive_drugs.setdefault(drug_lower, set()).add("CIViC")
                if assertion.locus_match == "variant":
                    sensitive_variant_level.setdefault(drug_lower, set()).add("CIViC")

    # Check CIViC evidence items
    # Use computed properties that check both clinical_significance AND evidence_direction
    for civic in evidence.civic_evidence:
        if not civic.drugs or len(civic.drugs) > 1:  # Skip combinations
            continue
        drug_lower = civic.drugs[0].lower()
        if civic.is_resistance:
            resistant_drugs.setdefault(drug_lower, set()).add("CIViC")
            if civic.locus_match == "variant":
                resistant_variant_level.setdefault(drug_lower, set()).add("CIViC")
        elif civic.is_sensitivity:
            sensitive_drugs.setdefault(drug_lower, set()).add("CIViC")
            if civic.locus_match == "variant":
                sensitive_variant_level.setdefault(drug_lower, set()).add("CIViC")

    # Find drugs with TRUE CROSS-SOURCE conflicts at VARIANT LEVEL only
    # Gene/codon-level conflicts are not flagged as different variants can behave differently
    for drug in set(sensitive_variant_level.keys()) & set(resistant_variant_level.keys()):
        sens_sources = sensitive_variant_level[drug]
        resist_sources = resistant_variant_level[drug]

        # Only flag if sources are truly different (cross-source conflict)
        if sens_sources == resist_sources:
            continue

        sens_only = sens_sources - resist_sources
        resist_only = resist_sources - sens_sources

        if sens_only and resist_only:
            # FDA vs VICC conflicts are handled separately (SIGNIFICANT, not HIGH)
            # because VICC data is less reliable than FDA
            if sens_only == {"FDA"} and resist_only == {"VICC"}:
                continue  # Will be caught in FDA vs VICC check below

            high_conflicts.append(
                f"Conflicting drug response for {drug} at variant level: "
                f"sensitive ({', '.join(sorted(sens_sources))}) vs "
                f"resistant ({', '.join(sorted(resist_sources))})"
            )

    # Check FDA (sensitive) vs VICC (resistant) at variant level
    # SKIP conflicts for known sensitizing variants where acquired resistance is expected.
    # For example, BRAF V600E + dabrafenib: FDA approval proves initial sensitivity,
    # VICC resistance reports represent acquired resistance (expected clinical behavior).
    from oncomind.config.constants import is_sensitizing_variant_for_drug

    query_gene = evidence.identifiers.gene if evidence.identifiers else None
    query_variant = evidence.identifiers.variant if evidence.identifiers else None

    for drug in fda_sensitive_variant_level.keys():
        if drug in vicc_resistant_variant_level:
            # Check if this is a known sensitizing variant for this drug
            # If so, VICC resistance is acquired (expected), not intrinsic (conflict)
            if query_gene and query_variant and is_sensitizing_variant_for_drug(
                query_gene, query_variant, drug
            ):
                # Skip - this is expected acquired resistance, not a true conflict
                continue

            # True conflict: FDA approves but VICC reports resistance for unknown variant
            significant_conflicts.append(
                f"FDA approves {drug} (sensitive) at variant level but VICC reports resistance at variant level — requires validation"
            )

    # Check ClinVar significance conflicts (internal) - differentiate by locus level
    # Only flag as a true conflict if both pathogenic and benign exist at the SAME locus level
    # (different variants in the same gene can legitimately have different pathogenicity)
    clinvar_by_locus: dict[str, set[str]] = {
        "variant": set(),
        "codon": set(),
        "gene": set(),
    }
    for entry in evidence.clinvar_entries:
        if entry.clinical_significance:
            sig = entry.clinical_significance.lower()
            locus = entry.locus_match  # "variant", "codon", or "gene"
            if "pathogenic" in sig and "benign" not in sig:
                clinvar_by_locus[locus].add("pathogenic")
            elif "benign" in sig and "pathogenic" not in sig:
                clinvar_by_locus[locus].add("benign")

    # Check for conflicts at each locus level (most specific first)
    # Variant-level conflict is most concerning
    if "pathogenic" in clinvar_by_locus["variant"] and "benign" in clinvar_by_locus["variant"]:
        high_conflicts.append(
            "ClinVar has conflicting interpretations at variant level: both pathogenic and benign submissions for this exact variant"
        )
    elif "pathogenic" in clinvar_by_locus["codon"] and "benign" in clinvar_by_locus["codon"]:
        high_conflicts.append(
            "ClinVar has conflicting interpretations at codon level: both pathogenic and benign submissions for variants at this position"
        )
    # Note: gene-level "conflicts" are not flagged — different variants legitimately differ

    # Check ClinVar vs CIViC pathogenicity conflicts (cross-source)
    # Only flag if BOTH are at variant level — gene-level ClinVar benign doesn't conflict with
    # variant-level CIViC actionable (different variants can have different significance)
    clinvar_benign_at_variant = "benign" in clinvar_by_locus["variant"] and "pathogenic" not in clinvar_by_locus["variant"]

    if clinvar_benign_at_variant:
        civic_actionable_at_variant: list[str] = []

        # Check CIViC assertions at variant level for ONCOGENIC type or actionable evidence
        for assertion in evidence.civic_assertions:
            if assertion.locus_match != "variant":
                continue
            if assertion.assertion_type and assertion.assertion_type.upper() == "ONCOGENIC":
                civic_actionable_at_variant.append("CIViC assertion (oncogenic)")
                break
            if assertion.significance and "ONCOGENIC" in assertion.significance.upper():
                civic_actionable_at_variant.append("CIViC assertion (oncogenic)")
                break
            # Predictive assertions with therapies suggest actionability
            if assertion.assertion_type and assertion.assertion_type.upper() == "PREDICTIVE":
                if assertion.therapies:
                    civic_actionable_at_variant.append("CIViC assertion (predictive)")
                    break

        # Check CIViC evidence items at variant level for PREDISPOSING/ONCOGENIC types
        for evi in evidence.civic_evidence:
            if evi.locus_match != "variant":
                continue
            if evi.evidence_type:
                etype = evi.evidence_type.upper()
                if etype in ("PREDISPOSING", "ONCOGENIC"):
                    civic_actionable_at_variant.append(f"CIViC evidence ({etype.lower()})")
                    break
            # Predictive evidence with drugs suggests actionability
            if evi.evidence_type and evi.evidence_type.upper() == "PREDICTIVE":
                if evi.drugs:
                    civic_actionable_at_variant.append("CIViC evidence (predictive)")
                    break

        if civic_actionable_at_variant:
            source_str = civic_actionable_at_variant[0]
            high_conflicts.append(
                f"ClinVar classifies as benign at variant level but {source_str} suggests clinical relevance for this exact variant"
            )

    return high_conflicts, significant_conflicts

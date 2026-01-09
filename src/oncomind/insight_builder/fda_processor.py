"""FDA label processing and evidence conversion.

This module handles FDA label data transformation, including:
- Biomarker specificity parsing from indication text
- Variant coverage determination
- Label to evidence model conversion
- Combination therapy partner linking
- Tumor type matching
"""

import re
from typing import Any

from oncomind.api.fda_label_service import get_fda_labels_for_drugs
from oncomind.models.evidence import (
    FDAApproval,
    FDALabelEvidence,
)
from oncomind.models.evidence.base import EvidenceLevel, is_pan_cancer_term, extract_variant_codon
from oncomind.models.evidence.fda import extract_combination_partners
from oncomind.config.debug import get_logger

logger = get_logger(__name__)


def parse_biomarker_specificity(text: str, gene: str | None = None) -> dict | None:
    """Parse biomarker specificity from FDA indication text.

    Determines what level of biomarker specificity the FDA approval requires:
    - variant: specific variant like "KRAS G12C" (also captures codon)
    - codon: specific codon like "BRAF V600"
    - exon: specific exon like "EGFR exon 19"
    - gene: any alteration in the gene like "AKT1 alteration"
    - phenotype: biomarker phenotype like "MSI-H", "dMMR", "HRD", "TMB-H"

    Args:
        text: FDA indication text to parse
        gene: Gene symbol to look for (optional for phenotype-only parsing)

    Returns:
        Dict with parsing results:
        - For variant: {"level": "variant", "codon": "G12", "specified_variant": "G12C",
                        "target_genes": ["KRAS"]}
        - For codon: {"level": "codon", "codon": "V600", "target_genes": ["BRAF"]}
        - For exon: {"level": "exon", "exon": "19", "target_genes": ["EGFR"]}
        - For gene: {"level": "gene", "target_genes": ["BRCA1", "BRCA2"]}
        - For phenotype: {"level": "phenotype", "phenotypes": ["MSI-H", "dMMR"]}
        Or None if no match
    """
    if not text:
        return None

    text_upper = text.upper()

    # Check for phenotype-based approvals first (these are gene-agnostic)
    phenotype_patterns = [
        (r'MSI-H|MICROSATELLITE INSTABILITY[- ]HIGH', 'MSI-H'),
        (r'DMMR|MISMATCH REPAIR DEFICIEN', 'dMMR'),
        (r'HRD|HOMOLOGOUS RECOMBINATION DEFICIEN', 'HRD'),
        (r'TMB-H|TUMOR MUTATIONAL BURDEN[- ]HIGH', 'TMB-H'),
    ]

    phenotypes_found = []
    for pattern, phenotype in phenotype_patterns:
        if re.search(pattern, text_upper):
            phenotypes_found.append(phenotype)

    if phenotypes_found:
        return {"level": "phenotype", "phenotypes": phenotypes_found}

    # If no gene specified, we can't do gene-specific parsing
    if not gene:
        return None

    # Variant-specific: "KRAS G12C" - captures ref, position, alt
    variant_match = re.search(rf"{gene}\s+([A-Z])(\d+)([A-Z])", text, re.IGNORECASE)
    if variant_match:
        ref, codon_num, alt = variant_match.groups()
        ref = ref.upper()
        alt = alt.upper()
        return {
            "level": "variant",
            "codon": f"{ref}{codon_num}",  # "G12"
            "specified_variant": f"{ref}{codon_num}{alt}",  # "G12C"
            "target_genes": [gene.upper()],
        }

    # Codon-specific: "BRAF V600"
    codon_match = re.search(rf"{gene}\s+([A-Z]\d+)\s", text, re.IGNORECASE)
    if codon_match:
        return {
            "level": "codon",
            "codon": codon_match.group(1).upper(),
            "target_genes": [gene.upper()],
        }

    # Exon-specific: "EGFR exon 19"
    exon_match = re.search(rf"{gene}\s+exon\s+(\d+)", text, re.IGNORECASE)
    if exon_match:
        return {
            "level": "exon",
            "exon": exon_match.group(1),
            "target_genes": [gene.upper()],
        }

    # Gene-level: "AKT1 alteration" or "BRCA-mutated"
    gene_level_patterns = [
        rf"{gene}.*(?:alter|mutat)",
        rf"{gene}[- ]POSITIVE",
        rf"{gene}\s+REARRANGEMENT",
        rf"{gene}\s+FUSION",
    ]
    for pattern in gene_level_patterns:
        if re.search(pattern, text, re.IGNORECASE):
            return {"level": "gene", "target_genes": [gene.upper()]}

    return None


def get_variant_match_level(
    query_variant: str,
    specificity: dict,
    allow_codon_fallback: bool = True,
) -> str | None:
    """Get the match level between a query variant and FDA approval specificity.

    Determines what level of match exists between the queried variant and
    the biomarker specificity parsed from FDA indication text.

    Args:
        query_variant: The variant being queried (e.g., "E17K", "V600E")
        specificity: Dict from parse_biomarker_specificity
        allow_codon_fallback: If True, a variant-level approval can match at
            codon level (e.g., G12D matches G12C approval at codon level).
            If False, variant-level approvals require exact match.

    Returns:
        Match level: "variant", "codon", "gene", or None if no match
    """
    if specificity["level"] == "gene":
        return "gene"

    if specificity["level"] == "variant":
        # Check for exact variant match
        if query_variant.upper() == specificity["specified_variant"].upper():
            return "variant"
        # Optionally check for codon match (same position, different AA)
        if allow_codon_fallback:
            query_codon = extract_variant_codon(query_variant)
            if query_codon and query_codon == specificity["codon"]:
                return "codon"
        return None

    if specificity["level"] == "codon":
        # Extract codon from variant: E17K → E17
        query_codon = extract_variant_codon(query_variant)
        if query_codon and query_codon == specificity["codon"]:
            return "codon"
        return None

    if specificity["level"] == "exon":
        # For exon-level, we'd need to know what exon the variant is in
        # This requires additional mapping - return None for now
        return None

    return None


def is_variant_covered(query_variant: str, specificity: dict) -> bool:
    """Check if a query variant is covered by the FDA approval specificity.

    Determines whether a specific variant (e.g., "E17K") is covered by
    the biomarker specificity parsed from FDA indication text.

    Note: For variant-level approvals, only exact matches count as "covered".
    A G12D variant is NOT covered by a G12C-specific approval, even though
    they share the same codon.

    Args:
        query_variant: The variant being queried (e.g., "E17K", "V600E")
        specificity: Dict from parse_biomarker_specificity

    Returns:
        True if the variant is covered by this approval, False otherwise
    """
    # Don't allow codon fallback for is_variant_covered - must be exact match
    # for variant-level approvals
    match_level = get_variant_match_level(
        query_variant, specificity, allow_codon_fallback=False
    )
    return match_level is not None


def populate_locus_variant_match(
    fda_labels: list[FDALabelEvidence], query_variant: str | None = None
) -> None:
    """Populate the locus_variant_match field on FDALabelEvidence objects.

    Determines the relationship between the queried variant and the approved
    variant (variant/codon/gene match) and sets the locus_variant_match field.

    Modifies fda_labels in place.

    Args:
        fda_labels: List of FDALabelEvidence to populate
        query_variant: The variant being queried (e.g., "E17K", "G12C")
    """
    for label in fda_labels:
        # Parse biomarker specificity to determine locus_variant_match
        # locus_variant_match reflects the relationship between the queried
        # variant and the approved variant (variant/codon/gene match)
        if label.indications_and_usage and label.gene:
            biomarker_spec = parse_biomarker_specificity(
                label.indications_and_usage, label.gene
            )
            if biomarker_spec:
                # Determine locus_level based on queried variant vs approved variant
                if query_variant:
                    match_level = get_variant_match_level(query_variant, biomarker_spec)
                    if match_level:
                        label.locus_variant_match = EvidenceLevel(
                            level=match_level,
                            scope="specific",
                            origin="kb"
                        )
                elif biomarker_spec["level"] == "gene":
                    # No query variant but gene-level approval
                    label.locus_variant_match = EvidenceLevel(
                        level="gene",
                        scope="unspecified",
                        origin="kb"
                    )


def convert_fda_labels_to_approvals(
    fda_label_evidence: list[FDALabelEvidence]
) -> list[FDAApproval]:
    """Convert FDALabelEvidence to FDAApproval objects for LLM consumption.

    Takes the already-processed FDALabelEvidence (which has locus_variant_match
    already computed) and converts to FDAApproval for LLM prompt.

    Note: FDA labels do not contain sensitivity/resistance associations.
    Those signals come from CGI/CIViC/VICC therapeutic evidence separately.

    Args:
        fda_label_evidence: List of FDALabelEvidence (already processed)

    Returns:
        List of FDAApproval objects for LLM prompt
    """
    fda_approvals = []
    for label in fda_label_evidence:
        try:
            fda_approvals.append(label.to_approval())
        except Exception as e:
            logger.warning(f"Failed to convert FDA label to FDA approval: {e}")
    return fda_approvals


def enrich_fda_with_tumor_match(
    approvals: list[FDAApproval],
    tumor_type: str | None,
) -> list[FDAApproval]:
    """Enrich FDA approvals with cancer_type_match based on tumor type.

    Uses parse_indication_for_tumor() to determine if the approval
    matches the queried tumor type, then sets cancer_type_match
    for consistency with ClinicalTrialEvidence.

    Args:
        approvals: List of FDAApproval objects
        tumor_type: The queried tumor type (may be None)

    Returns:
        Same list with cancer_type_match populated
    """
    if not tumor_type:
        # No tumor type to match against - leave cancer_type_match as None
        return approvals

    for approval in approvals:
        # Use existing parse_indication_for_tumor method
        parsed = approval.parse_indication_for_tumor(tumor_type)
        tumor_match = parsed.get('tumor_match', False)

        # Check for pan-cancer / tumor-agnostic approvals (MSI-H, NTRK, etc.)
        # Uses centralized is_pan_cancer_term() which handles substring matching
        is_pan_cancer = is_pan_cancer_term(approval.indication)

        if tumor_match:
            if is_pan_cancer:
                level = "pan_cancer"
            else:
                level = "cancer_specific"
        else:
            # Not a match - extract what cancer it IS for
            extracted_cancer = approval.extract_indication_cancer_type()
            level = extracted_cancer if extracted_cancer else "pan_cancer"

        approval.cancer_type_match = EvidenceLevel(
            level=level,
            scope="specific" if tumor_match else "unspecified",
            origin="kb",
        )

    return approvals


def sort_fda_by_association(approvals: list[FDAApproval]) -> list[FDAApproval]:
    """Sort FDA approvals by association: Responsive/Sensitive first, Resistant last.

    Sort order:
    1. Responsive/Sensitivity (drug works for this variant)
    2. None/Unknown (from FDA labels without CGI association data)
    3. Resistant (drug does NOT work for this variant)

    This ensures that when displaying FDA approvals, the drugs most likely
    to be effective appear first, with resistance mutations clearly shown last.

    Args:
        approvals: List of FDAApproval objects

    Returns:
        Sorted list with Responsive first, Resistant last
    """
    def sort_key(approval: FDAApproval) -> int:
        assoc = (approval.association or "").lower()
        if assoc in ("responsive", "sensitivity", "sensitive"):
            return 0  # Best - drug works
        elif assoc in ("resistant", "resistance"):
            return 2  # Worst - drug doesn't work
        else:
            return 1  # Unknown - from FDA labels

    return sorted(approvals, key=sort_key)



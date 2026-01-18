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

from oncomind.models.evidence import (
    BiomarkerMatch,
    FDAApproval,
    FDALabelEvidence,
)
from oncomind.models.evidence.base import EvidenceLevel, is_pan_cancer_term, extract_variant_codon, tumor_types_match
from oncomind.models.evidence.fda import extract_combination_partners
from oncomind.config.debug import get_logger

logger = get_logger(__name__)


# Common cancer/tumor type patterns in FDA indications
# These are used to extract the disease context from indication text
FDA_TUMOR_PATTERNS = [
    # Specific cancers (order matters - more specific first)
    (r"acute myeloid leukemia|AML", "AML"),
    (r"non-?small cell lung cancer|NSCLC", "NSCLC"),
    (r"small cell lung cancer|SCLC", "SCLC"),
    (r"lung cancer|lung adenocarcinoma|lung carcinoma", "Lung Cancer"),
    (r"colorectal cancer|colon cancer|rectal cancer|CRC", "Colorectal Cancer"),
    (r"breast cancer|breast carcinoma", "Breast Cancer"),
    (r"melanoma", "Melanoma"),
    (r"glioma|astrocytoma|oligodendroglioma|glioblastoma|GBM", "Glioma"),
    (r"cholangiocarcinoma|bile duct cancer", "Cholangiocarcinoma"),
    (r"thyroid cancer|thyroid carcinoma", "Thyroid Cancer"),
    (r"prostate cancer|prostate carcinoma", "Prostate Cancer"),
    (r"ovarian cancer|ovarian carcinoma", "Ovarian Cancer"),
    (r"pancreatic cancer|pancreatic carcinoma", "Pancreatic Cancer"),
    (r"gastric cancer|stomach cancer", "Gastric Cancer"),
    (r"esophageal cancer", "Esophageal Cancer"),
    (r"bladder cancer|urothelial cancer|urothelial carcinoma", "Bladder Cancer"),
    (r"renal cell carcinoma|kidney cancer|RCC", "Kidney Cancer"),
    (r"hepatocellular carcinoma|liver cancer|HCC", "Liver Cancer"),
    (r"head and neck.*cancer|HNSCC", "Head and Neck Cancer"),
    (r"cervical cancer", "Cervical Cancer"),
    (r"endometrial cancer|uterine cancer", "Endometrial Cancer"),
    (r"gastrointestinal stromal tumor|GIST", "GIST"),
    (r"chronic myeloid leukemia|CML", "CML"),
    (r"chronic lymphocytic leukemia|CLL", "CLL"),
    (r"multiple myeloma", "Multiple Myeloma"),
    (r"lymphoma", "Lymphoma"),
    (r"leukemia", "Leukemia"),
    # Pan-cancer/tumor-agnostic biomarker approvals
    # These are tumor-agnostic approvals based on molecular biomarkers
    # IMPORTANT: MSI-H/dMMR/TMB-H only count as pan-cancer if combined with "solid tumor"
    # "MSI-H colorectal cancer" is NOT pan-cancer - it's specific to CRC
    (r"(MSI-?H|microsatellite instability.high).*solid tumou?r", "Solid Tumor"),
    (r"(dMMR|mismatch repair deficien).*solid tumou?r", "Solid Tumor"),
    (r"(TMB-?H|tumor mutational burden.high).*solid tumou?r", "Solid Tumor"),
    (r"NTRK.*fusion.*solid tumou?r", "Solid Tumor"),
    (r"RET.*fusion.*solid tumou?r", "Solid Tumor"),
    # Note: "advanced [specific] cancer" is stage, not pan-cancer
    # Only "solid tumor(s)" without specific type is truly tumor-agnostic
    (r"solid tumou?rs?(?!\s+(?:breast|lung|colon|prostate|pancreatic))", "Solid Tumor"),
]


def extract_tumor_from_indication(indication_text: str) -> list[str]:
    """Extract tumor/disease types from FDA indication text.

    Args:
        indication_text: FDA indications and usage text

    Returns:
        List of normalized tumor types found (e.g., ["AML", "Cholangiocarcinoma"])
    """
    if not indication_text:
        return []

    tumors_found = []
    text_lower = indication_text.lower()

    for pattern, tumor_name in FDA_TUMOR_PATTERNS:
        if re.search(pattern, text_lower, re.IGNORECASE):
            if tumor_name not in tumors_found:
                tumors_found.append(tumor_name)

    return tumors_found


def match_fda_tumor(
    indication_text: str,
    query_tumor: str | None,
) -> tuple[bool, str | None]:
    """Match query tumor against FDA indication text.

    Args:
        indication_text: FDA indications and usage text
        query_tumor: User's tumor type query (e.g., "Glioma", "NSCLC")

    Returns:
        Tuple of (matched: bool, match_type: str | None)
        match_type is "cancer_specific" if tumor matches, "pan_cancer" if solid tumor, None if no match
    """
    if not query_tumor or not indication_text:
        return False, None

    # Extract tumors from indication
    indication_tumors = extract_tumor_from_indication(indication_text)

    if not indication_tumors:
        return False, None

    # Check for pan-cancer approvals (e.g., "solid tumors" without specific type)
    if "Solid Tumor" in indication_tumors:
        return True, "pan_cancer"

    # Check if query tumor matches any indication tumor
    for ind_tumor in indication_tumors:
        if tumor_types_match(ind_tumor, query_tumor):
            return True, "cancer_specific"

    return False, None


def parse_biomarker_specificity(text: str, gene: str | None = None) -> dict | None:
    """Parse biomarker specificity from FDA indication text.

    Determines what level of biomarker specificity the FDA approval requires:
    - variant: specific variant like "KRAS G12C" (also captures codon)
    - variant_list: multiple specific variants like "BRAF V600E or V600K"
    - codon: specific codon like "BRAF V600"
    - exon: specific exon like "EGFR exon 19"
    - gene: any alteration in the gene like "AKT1 alteration"
    - phenotype: biomarker phenotype like "MSI-H", "dMMR", "HRD", "TMB-H"
    - contraindication: drug NOT indicated for this biomarker
    - wild_type_required: requires wild-type (unmutated) gene

    Args:
        text: FDA indication text to parse
        gene: Gene symbol to look for (optional for phenotype-only parsing)

    Returns:
        Dict with parsing results:
        - For variant: {"level": "variant", "codon": "G12", "specified_variant": "G12C",
                        "target_genes": ["KRAS"]}
        - For variant_list: {"level": "variant_list", "specified_variants": ["V600E", "V600K"],
                             "target_genes": ["BRAF"]}
        - For codon: {"level": "codon", "codon": "V600", "target_genes": ["BRAF"]}
        - For exon: {"level": "exon", "exon": "19", "target_genes": ["EGFR"]}
        - For gene: {"level": "gene", "target_genes": ["BRCA1", "BRCA2"]}
        - For phenotype: {"level": "phenotype", "phenotypes": ["MSI-H", "dMMR"]}
        - For contraindication: {"level": "contraindication", "target_genes": ["KRAS"]}
        - For wild_type_required: {"level": "wild_type_required", "target_genes": ["KRAS"]}
        Or None if no match
    """
    if not text:
        return None

    text_upper = text.upper()

    # 1. Phenotype-based approvals (gene-agnostic)
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

    # 2. Multi-gene patterns: "PIK3CA/AKT1/PTEN-alteration" or "PIK3CA/AKT1/PTEN -alteration"
    # Handle both with and without space before hyphen
    multi_gene_match = re.search(
        r'([A-Z0-9]+(?:/[A-Z0-9]+)+)\s*-?\s*(?:ALTER|MUTAT)',
        text_upper
    )
    if multi_gene_match:
        genes = multi_gene_match.group(1).split('/')
        return {"level": "gene", "target_genes": genes}

    # If no gene specified, can't continue with gene-specific parsing
    if not gene:
        return None

    gene_upper = gene.upper()

    # 3. Contraindication: "not indicated for KRAS-mutant"
    if re.search(rf"not indicated.*{gene}.*mutant", text, re.IGNORECASE):
        return {"level": "contraindication", "target_genes": [gene_upper]}

    # 4. Wild-type required: "KRAS wild-type"
    if re.search(rf"{gene}\s+wild[- ]?type", text, re.IGNORECASE):
        return {"level": "wild_type_required", "target_genes": [gene_upper]}

    # 5. EGFR "exon 19 deletions or exon 21 L858R" pattern (MUST come before multi-variant)
    # This is extremely common in EGFR TKI labels and needs special handling
    # because it combines two different specificity types in one approval
    # Match both "EGFR exon 19" and "(EGFR) exon 19" formats
    egfr_combined_match = re.search(
        rf"(?:{gene}|\({gene}\))\s+exon\s+19\s+deletion.*?(?:or|and).*?(?:exon\s+21\s+)?L858R",
        text, re.IGNORECASE
    )
    if egfr_combined_match:
        return {
            "level": "variant_list",
            "specified_variants": ["L858R"],
            "exons": ["19"],  # Also covers exon 19 deletions
            "target_genes": [gene_upper],
        }

    # 6. Exon-specific: "EGFR exon 19" or "(EGFR) exon 19" (before multi-variant)
    exon_match = re.search(rf"(?:{gene}|\({gene}\))\s+exon\s+(\d+)", text, re.IGNORECASE)
    if exon_match:
        return {
            "level": "exon",
            "exon": exon_match.group(1),
            "target_genes": [gene_upper],
        }

    # 7. Multiple variants: "BRAF V600E or V600K"
    multi_variant_match = re.findall(
        rf"{gene}\s+([A-Z]\d+[A-Z])", text, re.IGNORECASE
    )
    if len(multi_variant_match) > 1:
        return {
            "level": "variant_list",
            "specified_variants": [v.upper() for v in multi_variant_match],
            "target_genes": [gene_upper],
        }

    # 8. Variant-specific: "KRAS G12C" - captures ref, position, alt
    variant_match = re.search(
        rf"{gene}\s+([A-Z])(\d+)([A-Z])", text, re.IGNORECASE
    )
    if variant_match:
        ref, codon_num, alt = variant_match.groups()
        return {
            "level": "variant",
            "codon": f"{ref.upper()}{codon_num}",
            "specified_variant": f"{ref.upper()}{codon_num}{alt.upper()}",
            "target_genes": [gene_upper],
        }

    # 9. Codon-specific: "BRAF V600" (no trailing AA)
    codon_match = re.search(
        rf"{gene}\s+([A-Z]\d+)(?:\s|[-]|$|[^A-Z0-9])", text, re.IGNORECASE
    )
    if codon_match:
        return {
            "level": "codon",
            "codon": codon_match.group(1).upper(),
            "target_genes": [gene_upper],
        }

    # 9. Gene-level: "AKT1 alteration", "BRCA-mutated", "ALK-positive"
    # Also handles "(IDH1) or ... (IDH2) mutation" format from FDA labels
    # NOTE: The "(GENE)...mutation" pattern should NOT match if followed by exon/variant info
    gene_level_patterns = [
        rf"{gene}\s*-?\s*(?:alter|mutat)",
        rf"{gene}[- ]?positive",
        rf"{gene}\s+rearrangement",
        rf"{gene}\s+fusion",
        rf"\({gene}\)\s+(?!exon).*?mutation",  # "(IDH1) mutation" but NOT "(EGFR) exon 19...mutation"
    ]
    for pattern in gene_level_patterns:
        if re.search(pattern, text, re.IGNORECASE):
            return {"level": "gene", "target_genes": [gene_upper]}

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
        Match level: "variant", "codon", "gene", "contraindicated", or None if no match
    """
    level = specificity.get("level")

    # Contraindication/wild-type: explicitly NOT covered
    if level in ("contraindication", "wild_type_required"):
        return "contraindicated"

    # Phenotype: not variant-based, can't determine match
    if level == "phenotype":
        return None

    if level == "gene":
        return "gene"

    if level == "variant":
        # Check for exact variant match
        if query_variant.upper() == specificity["specified_variant"].upper():
            return "variant"
        # Optionally check for codon match (same position, different AA)
        if allow_codon_fallback:
            query_codon = extract_variant_codon(query_variant)
            if query_codon and query_codon == specificity.get("codon"):
                return "codon"
        return None

    if level == "variant_list":
        # Check if query variant is in the list
        specified = specificity.get("specified_variants", [])
        if query_variant.upper() in [v.upper() for v in specified]:
            return "variant"
        # Check for codon match with any variant in the list
        if allow_codon_fallback:
            query_codon = extract_variant_codon(query_variant)
            for v in specified:
                if query_codon and query_codon == extract_variant_codon(v):
                    return "codon"
        return None

    if level == "codon":
        # Extract codon from variant: E17K → E17
        query_codon = extract_variant_codon(query_variant)
        if query_codon and query_codon == specificity.get("codon"):
            return "codon"
        return None

    if level == "exon":
        # For exon-level, we'd need to know what exon the variant is in
        # This requires additional mapping - return None for now
        return None

    return None

def is_variant_covered(query_variant: str, specificity: dict) -> tuple[bool, str | None]:
    if not specificity:
        return (False, None)

    level = specificity.get("level")

    if level in ("wild_type_required", "contraindication"):
        return (False, "contraindicated")

    if level == "phenotype":
        return (False, None)

    if level == "gene":
        return (True, "gene")

    if level == "variant":
        # Exact match
        if query_variant.upper() == specificity.get("specified_variant", "").upper():
            return (True, "variant")
        # Same codon, different variant
        query_codon = extract_variant_codon(query_variant)
        if query_codon and query_codon == specificity.get("codon"):
            return (False, "codon")  # <-- THIS IS THE KEY LINE
        return (False, None)

    if level == "codon":
        query_codon = extract_variant_codon(query_variant)
        if query_codon and query_codon == specificity.get("codon"):
            return (True, "codon")
        return (False, None)

    if level == "variant_list":
        if query_variant.upper() in [v.upper() for v in specificity.get("specified_variants", [])]:
            return (True, "variant")
        query_codon = extract_variant_codon(query_variant)
        for v in specificity.get("specified_variants", []):
            if query_codon == extract_variant_codon(v):
                return (False, "codon")
        return (False, None)

    return (False, None)

def match_fda_approval(
    query_gene: str,
    query_variant: str,
    query_tumor: str,
    fda_indication: str,
) -> dict:
    """Match a query against FDA indication text."""

    specificity = parse_biomarker_specificity(fda_indication, query_gene)

    if not specificity:
        return {"matched": False, "reason": "gene_not_found"}

    # Check if query gene is in target genes
    target_genes = specificity.get("target_genes", [])
    if target_genes and query_gene.upper() not in [g.upper() for g in target_genes]:
        return {"matched": False, "reason": "gene_not_in_targets"}

    covered, match_level = is_variant_covered(query_variant, specificity)
    partners = extract_combination_partners(fda_indication)

    return {
        "matched": covered,
        "match_level": match_level,  # "gene", "codon", "variant", "contraindicated"
        "specificity": specificity,
        "combination_partners": partners,
        # TODO: add tumor matching
    }


def populate_locus_variant_match(
    fda_labels: list[FDALabelEvidence],
    query_variant: str | None = None,
    query_tumor: str | None = None,
) -> None:
    """Populate locus_variant_match, cancer_type_match, and biomarker_match on FDALabelEvidence.

    Determines the relationship between the queried variant and the approved
    variant (variant/codon/gene match) and sets:
    - locus_variant_match: EvidenceLevel for variant/codon/gene match level
    - cancer_type_match: EvidenceLevel for tumor type match (cancer_specific/pan_cancer/other)
    - biomarker_match: BiomarkerMatch with matched bool, match_level, and tumor_matched

    Modifies fda_labels in place.

    Args:
        fda_labels: List of FDALabelEvidence to populate
        query_variant: The variant being queried (e.g., "E17K", "G12C")
        query_tumor: The tumor type being queried (e.g., "colorectal cancer", "NSCLC")
    """
    for label in fda_labels:
        if not label.indications_and_usage or not label.gene:
            continue

        # Check tumor match first (applies to all cases)
        tumor_matched, tumor_match_type = match_fda_tumor(
            label.indications_and_usage, query_tumor
        )

        # Use match_fda_approval for biomarker matching logic
        if query_variant:
            result = match_fda_approval(
                query_gene=label.gene,
                query_variant=query_variant,
                query_tumor=query_tumor or "",
                fda_indication=label.indications_and_usage,
            )

            biomarker_covered = result["matched"]
            match_level = result.get("match_level")
            partners = result.get("combination_partners", [])

            # Final "matched" requires BOTH biomarker AND tumor match
            # A drug approved for AML shouldn't match for a Glioma patient
            fully_matched = biomarker_covered and tumor_matched

            # Set biomarker_match
            label.biomarker_match = BiomarkerMatch(
                matched=fully_matched,
                match_level=match_level if fully_matched else None,
                tumor_matched=tumor_matched,
                tumor_match_type=tumor_match_type,
                combination_partners=partners,
            )

            # Also set locus_variant_match for backward compatibility
            if match_level and fully_matched:
                label.locus_variant_match = EvidenceLevel(
                    level=match_level,
                    scope="specific",
                    origin="kb"
                )

            # Set cancer_type_match - always set this so we know what cancer the drug is approved for
            if tumor_matched:
                label.cancer_type_match = EvidenceLevel(
                    level=tumor_match_type or "cancer_specific",
                    scope="specific",
                    origin="kb"
                )
            else:
                # Not matched - extract what cancer the drug IS approved for
                indication_tumors = extract_tumor_from_indication(label.indications_and_usage)
                if indication_tumors:
                    # Use the first extracted tumor as the cancer type
                    label.cancer_type_match = EvidenceLevel(
                        level=indication_tumors[0].lower(),  # e.g., "breast cancer"
                        scope="unspecified",
                        origin="kb"
                    )
                # If no tumors extracted, cancer_type_match remains None (not pan_cancer)
        else:
            # No query variant - check if gene-level approval
            biomarker_spec = parse_biomarker_specificity(
                label.indications_and_usage, label.gene
            )
            if biomarker_spec and biomarker_spec.get("level") == "gene":
                partners = extract_combination_partners(label.indications_and_usage)

                # Gene-level match still requires tumor match
                fully_matched = tumor_matched

                label.locus_variant_match = EvidenceLevel(
                    level="gene",
                    scope="unspecified",
                    origin="kb"
                )
                label.biomarker_match = BiomarkerMatch(
                    matched=fully_matched,
                    match_level="gene" if fully_matched else None,
                    tumor_matched=tumor_matched,
                    tumor_match_type=tumor_match_type,
                    combination_partners=partners,
                )

                if tumor_matched:
                    label.cancer_type_match = EvidenceLevel(
                        level=tumor_match_type or "cancer_specific",
                        scope="specific",
                        origin="kb"
                    )
                else:
                    # Not matched - extract what cancer the drug IS approved for
                    indication_tumors = extract_tumor_from_indication(label.indications_and_usage)
                    if indication_tumors:
                        label.cancer_type_match = EvidenceLevel(
                            level=indication_tumors[0].lower(),
                            scope="unspecified",
                            origin="kb"
                        )


def sort_fda_by_association(approvals: list) -> list:
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



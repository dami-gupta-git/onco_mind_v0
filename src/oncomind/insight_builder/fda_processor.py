"""FDA label processing utilities.

This module handles FDA label data transformation, including:
- Biomarker specificity parsing from indication text
- Variant coverage determination (is_variant_covered)
"""

import re

from oncomind.models.evidence.base import extract_variant_codon
from oncomind.config.debug import get_logger

logger = get_logger(__name__)


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
        (r"MSI-H|MICROSATELLITE INSTABILITY[- ]HIGH", "MSI-H"),
        (r"DMMR|MISMATCH REPAIR DEFICIEN", "dMMR"),
        (r"HRD|HOMOLOGOUS RECOMBINATION DEFICIEN", "HRD"),
        (r"TMB-H|TUMOR MUTATIONAL BURDEN[- ]HIGH", "TMB-H"),
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
        r"([A-Z0-9]+(?:/[A-Z0-9]+)+)\s*-?\s*(?:ALTER|MUTAT)", text_upper
    )
    if multi_gene_match:
        genes = multi_gene_match.group(1).split("/")
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
        text,
        re.IGNORECASE,
    )
    if egfr_combined_match:
        return {
            "level": "variant_list",
            "specified_variants": ["L858R"],
            "exons": ["19"],  # Also covers exon 19 deletions
            "target_genes": [gene_upper],
        }

    # 6. Exon-specific: "EGFR exon 19" or "(EGFR) exon 19" (before multi-variant)
    exon_match = re.search(
        rf"(?:{gene}|\({gene}\))\s+exon\s+(\d+)", text, re.IGNORECASE
    )
    if exon_match:
        return {
            "level": "exon",
            "exon": exon_match.group(1),
            "target_genes": [gene_upper],
        }

    # 7. Multiple variants: "BRAF V600E or V600K"
    multi_variant_match = re.findall(rf"{gene}\s+([A-Z]\d+[A-Z])", text, re.IGNORECASE)
    if len(multi_variant_match) > 1:
        return {
            "level": "variant_list",
            "specified_variants": [v.upper() for v in multi_variant_match],
            "target_genes": [gene_upper],
        }

    # 8. Variant-specific: "KRAS G12C" - captures ref, position, alt
    variant_match = re.search(rf"{gene}\s+([A-Z])(\d+)([A-Z])", text, re.IGNORECASE)
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

    # 10. Gene-level: "AKT1 alteration", "BRCA-mutated", "ALK-positive"
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


def is_variant_covered(
    query_variant: str, specificity: dict
) -> tuple[bool, str | None]:
    """Check if a query variant is covered by the parsed biomarker specificity.

    Args:
        query_variant: The variant being queried (e.g., "E17K", "V600E")
        specificity: Dict from parse_biomarker_specificity

    Returns:
        Tuple of (covered: bool, match_level: str | None)
        match_level is "variant", "codon", "gene", "contraindicated", or None
    """
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
            return (False, "codon")
        return (False, None)

    if level == "codon":
        query_codon = extract_variant_codon(query_variant)
        if query_codon and query_codon == specificity.get("codon"):
            return (True, "codon")
        return (False, None)

    if level == "variant_list":
        if query_variant.upper() in [
            v.upper() for v in specificity.get("specified_variants", [])
        ]:
            return (True, "variant")
        query_codon = extract_variant_codon(query_variant)
        for v in specificity.get("specified_variants", []):
            if query_codon == extract_variant_codon(v):
                return (False, "codon")
        return (False, None)

    return (False, None)

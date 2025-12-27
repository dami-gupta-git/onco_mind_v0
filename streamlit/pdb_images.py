"""PDB structure image URL generation for cancer genes."""

# Gene to PDB ID mapping for common cancer genes
# Selected structures show the most relevant domain (kinase, DNA-binding, etc.)
GENE_PDB_MAPPING: dict[str, dict] = {
    "BRAF": {"pdb_id": "6u2g", "description": "Kinase domain"},
    "EGFR": {"pdb_id": "1m17", "description": "Kinase domain with erlotinib"},
    "KRAS": {"pdb_id": "4obe", "description": "GTPase"},
    "NRAS": {"pdb_id": "5uhv", "description": "GTPase"},
    "HRAS": {"pdb_id": "4q21", "description": "GTPase"},
    "PIK3CA": {"pdb_id": "4ovv", "description": "PI3K catalytic subunit"},
    "TP53": {"pdb_id": "2ocj", "description": "DNA binding domain"},
    "ALK": {"pdb_id": "4mkc", "description": "Kinase domain"},
    "MET": {"pdb_id": "3dkc", "description": "Kinase domain"},
    "HER2": {"pdb_id": "3pp0", "description": "Kinase domain"},
    "ERBB2": {"pdb_id": "3pp0", "description": "Kinase domain"},
    "KIT": {"pdb_id": "1t46", "description": "Kinase domain"},
    "BRCA1": {"pdb_id": "3wrz", "description": "BRCT domain"},
    "BRCA2": {"pdb_id": "1mje", "description": "DNA binding domain"},
    "PDGFRA": {"pdb_id": "5k5x", "description": "Kinase domain"},
    "FGFR1": {"pdb_id": "4v05", "description": "Kinase domain"},
    "FGFR2": {"pdb_id": "2psq", "description": "Kinase domain"},
    "FGFR3": {"pdb_id": "4k33", "description": "Kinase domain"},
    "RET": {"pdb_id": "2ivt", "description": "Kinase domain"},
    "ROS1": {"pdb_id": "3zbf", "description": "Kinase domain"},
    "JAK2": {"pdb_id": "4ytc", "description": "Kinase domain"},
    "ABL1": {"pdb_id": "1iep", "description": "Kinase domain with imatinib"},
    "BTK": {"pdb_id": "3gen", "description": "Kinase domain"},
    "AKT1": {"pdb_id": "3o96", "description": "Kinase domain"},
    "MTOR": {"pdb_id": "4dri", "description": "FRB domain"},
    "CDK4": {"pdb_id": "2w96", "description": "Kinase with cyclin D"},
    "CDK6": {"pdb_id": "1bi7", "description": "Kinase domain"},
    "ERBB3": {"pdb_id": "4leo", "description": "Kinase domain"},
    "MAP2K1": {"pdb_id": "3eqc", "description": "MEK1 kinase"},
    "MAP2K2": {"pdb_id": "1s9j", "description": "MEK2 kinase"},
    "IDH1": {"pdb_id": "3inm", "description": "Dehydrogenase"},
    "IDH2": {"pdb_id": "4ja8", "description": "Dehydrogenase"},
    "PTEN": {"pdb_id": "1d5r", "description": "Phosphatase domain"},
    "VHL": {"pdb_id": "1lm8", "description": "VHL-ElonginC-ElonginB complex"},
    "SMAD4": {"pdb_id": "1dd1", "description": "MH2 domain"},
    "APC": {"pdb_id": "3nmz", "description": "Armadillo repeat"},
    "ATM": {"pdb_id": "5np0", "description": "FAT domain"},
    "ARID1A": {"pdb_id": "6lth", "description": "SWI/SNF complex"},
    "FBXW7": {"pdb_id": "2ovq", "description": "WD40 domain"},
    "NOTCH1": {"pdb_id": "3eto", "description": "Ankyrin domain"},
    "CTNNB1": {"pdb_id": "1jdh", "description": "Armadillo repeat"},
    "STK11": {"pdb_id": "2wtk", "description": "Kinase domain"},
    "NF1": {"pdb_id": "1nf1", "description": "GAP domain"},
    "RB1": {"pdb_id": "1ad6", "description": "Pocket domain"},
    "KEAP1": {"pdb_id": "2flu", "description": "Kelch domain"},
    "NFE2L2": {"pdb_id": "2flu", "description": "With KEAP1"},
    "ESR1": {"pdb_id": "1err", "description": "Ligand binding domain"},
    "AR": {"pdb_id": "1e3g", "description": "Ligand binding domain"},
}


def get_pdb_image_url(gene: str) -> dict | None:
    """Get PDB structure image URL for a gene.

    Args:
        gene: Gene symbol (e.g., "BRAF", "EGFR")

    Returns:
        Dict with keys: url, pdb_id, description
        None if gene not in mapping
    """
    gene_upper = gene.upper()

    if gene_upper not in GENE_PDB_MAPPING:
        return None

    info = GENE_PDB_MAPPING[gene_upper]
    pdb_id = info["pdb_id"].lower()

    # URL format: middle 2 characters of PDB ID determine subdirectory
    middle2 = pdb_id[1:3]
    url = f"https://cdn.rcsb.org/images/structures/{middle2}/{pdb_id}/{pdb_id}_assembly-1.jpeg"

    return {
        "url": url,
        "pdb_id": pdb_id.upper(),
        "description": info["description"],
    }


def get_pdb_page_url(gene: str) -> str | None:
    """Get RCSB PDB page URL for a gene.

    Args:
        gene: Gene symbol

    Returns:
        URL to RCSB structure page, or None if gene not mapped
    """
    gene_upper = gene.upper()

    if gene_upper not in GENE_PDB_MAPPING:
        return None

    pdb_id = GENE_PDB_MAPPING[gene_upper]["pdb_id"].upper()
    return f"https://www.rcsb.org/structure/{pdb_id}"

"""Centralized constants and mappings for OncoMind.

This module consolidates all hardcoded mappings used across the codebase:
- Gene aliases (different nomenclatures for the same gene)
- Tumor type mappings (abbreviations to full names)
- Amino acid codes
- Variant type classifications

Centralizing these makes maintenance easier and ensures consistency.
"""

# =============================================================================
# GENE ALIASES
# =============================================================================
# Maps gene symbols to their alternative names used in different databases
# FDA labels, CIViC, CGI, etc. may use different nomenclature

GENE_ALIASES: dict[str, list[str]] = {
    "ERBB2": ["HER2", "NEU"],
    "HER2": ["ERBB2", "NEU"],
    "EGFR": ["HER1", "ERBB1"],
    "ERBB1": ["EGFR", "HER1"],
    "HER1": ["EGFR", "ERBB1"],
    "MET": ["HGFR"],
    "KIT": ["CD117", "C-KIT"],
    "PDGFRA": ["CD140A"],
    "PDGFRB": ["CD140B"],
    # BRCA1/BRCA2: FDA labels often use generic "BRCA" or "gBRCA" (germline BRCA)
    # without specifying BRCA1 vs BRCA2. PARP inhibitors like Olaparib (Lynparza)
    # are approved for "BRCA-mutated" cancers.
    "BRCA1": ["BRCA", "gBRCA"],
    "BRCA2": ["BRCA", "gBRCA"],
}


# =============================================================================
# TUMOR TYPE MAPPINGS
# =============================================================================
# Maps tumor type abbreviations to their full names and synonyms
# Used for matching user input to database entries (CGI, FDA, CIViC)

TUMOR_TYPE_MAPPINGS: dict[str, list[str]] = {
    # Lung
    "nsclc": ["non-small cell lung", "nsclc", "lung non-small cell", "lung adenocarcinoma", "lung squamous", "lung non-small cell carcinoma"],
    "l": ["lung", "non-small cell lung", "nsclc", "small cell lung", "lung carcinoma"],
    "sclc": ["small cell lung", "sclc"],
    "luad": ["lung adenocarcinoma", "luad"],
    "lusc": ["lung squamous", "lusc", "squamous cell lung"],

    # Colorectal
    "crc": ["colorectal", "colon", "crc", "rectal", "rectum"],
    "coad": ["colon adenocarcinoma", "coad"],
    "read": ["rectal adenocarcinoma", "read", "rectum"],
    "coread": ["colorectal", "colon", "rectal", "colorectal adenocarcinoma", "crc"],  # CGI uses COREAD

    # Melanoma
    "mel": ["melanoma", "mel", "cutaneous melanoma", "skin melanoma"],
    "skcm": ["skin cutaneous melanoma", "skcm", "melanoma"],

    # Breast
    "bc": ["breast", "bc", "breast cancer", "breast carcinoma"],
    "brca": ["breast carcinoma", "brca", "breast cancer"],
    "idc": ["invasive ductal carcinoma", "idc", "breast ductal"],
    "ilc": ["invasive lobular carcinoma", "ilc", "breast lobular"],

    # Thyroid
    "atc": ["anaplastic thyroid", "atc", "thyroid anaplastic"],
    "thca": ["thyroid carcinoma", "thca", "thyroid cancer"],
    "ptc": ["papillary thyroid", "ptc"],

    # Gastrointestinal
    "gist": ["gastrointestinal stromal", "gist"],
    "stad": ["stomach adenocarcinoma", "stad", "gastric"],
    "esca": ["esophageal carcinoma", "esca", "esophageal"],
    "paad": ["pancreatic adenocarcinoma", "paad", "pancreatic", "pancreas"],
    "lihc": ["liver hepatocellular", "lihc", "hepatocellular", "hcc"],
    "chol": ["cholangiocarcinoma", "chol", "bile duct"],

    # Genitourinary
    "prad": ["prostate adenocarcinoma", "prad", "prostate"],
    "blca": ["bladder carcinoma", "blca", "bladder", "urothelial"],
    "rcc": ["renal cell carcinoma", "rcc", "kidney"],
    "ccrcc": ["clear cell renal", "ccrcc", "kidney clear cell"],

    # Gynecologic
    "ov": ["ovarian", "ov", "ovary"],
    "ucec": ["uterine corpus endometrial", "ucec", "endometrial", "uterine"],
    "cesc": ["cervical squamous", "cesc", "cervical"],

    # Head and Neck
    "hnsc": ["head and neck squamous", "hnsc", "head neck"],

    # Brain
    "gbm": ["glioblastoma", "gbm", "glioblastoma multiforme"],
    "lgg": ["low grade glioma", "lgg", "glioma"],

    # Hematologic
    "aml": ["acute myeloid leukemia", "aml"],
    "cml": ["chronic myeloid leukemia", "cml"],
    "all": ["acute lymphoblastic leukemia", "all"],
    "cll": ["chronic lymphocytic leukemia", "cll"],
    "dlbcl": ["diffuse large b-cell lymphoma", "dlbcl"],
    "mm": ["multiple myeloma", "mm", "myeloma"],
}


# =============================================================================
# PRIORITY TUMOR TYPES FOR UI
# =============================================================================
# OncoTree codes for commonly assessed tumor types
# These appear first in UI dropdowns/autocomplete

PRIORITY_TUMOR_CODES: list[str] = [
    # Lung
    "NSCLC", "LUAD", "LUSC", "SCLC",
    # Breast
    "BRCA", "IDC", "ILC",
    # Colorectal
    "CRC", "COAD", "READ",
    # Melanoma
    "MEL", "SKCM",
    # Pancreatic
    "PAAD",
    # Prostate
    "PRAD",
    # Ovarian
    "OV",
    # Glioblastoma
    "GBM",
    # Bladder
    "BLCA",
    # Kidney
    "RCC", "CCRCC",
    # Head and Neck
    "HNSC",
    # Gastric
    "STAD",
    # Liver
    "LIHC",
    # Thyroid
    "THCA",
    # Endometrial
    "UCEC",
    # Hematologic
    "AML", "CML",
]


# =============================================================================
# AMINO ACID CODES
# =============================================================================
# Standard amino acid code conversions (3-letter to 1-letter and vice versa)

AMINO_ACID_3TO1: dict[str, str] = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'TER': '*', 'STP': '*', 'STOP': '*',  # Stop codons
    'SEC': 'U',  # Selenocysteine (rare)
    'PYL': 'O',  # Pyrrolysine (rare)
}

AMINO_ACID_1TO3: dict[str, str] = {v: k for k, v in AMINO_ACID_3TO1.items() if v not in ('*',)}
AMINO_ACID_1TO3['*'] = 'TER'  # Prefer TER for stop codon


# =============================================================================
# VARIANT TYPE CLASSIFICATIONS
# =============================================================================
# Variant types allowed for SNP/small indel analysis

ALLOWED_VARIANT_TYPES: set[str] = {
    'missense',
    'nonsense',
    'insertion',
    'deletion',
    'frameshift',
}

# Variant types that are structural (not supported in current system)
STRUCTURAL_VARIANT_TYPES: set[str] = {
    'fusion',
    'amplification',
    'deletion_large',
    'rearrangement',
    'copy_number',
    'truncating',
    'splice',
}


# =============================================================================
# GENE CHROMOSOMES
# =============================================================================
# Gene to chromosome mapping for common cancer genes
# Used for genomic coordinate lookups

GENE_CHROMOSOMES: dict[str, str] = {
    "BRAF": "7",
    "KRAS": "12",
    "NRAS": "1",
    "EGFR": "7",
    "TP53": "17",
    "PIK3CA": "3",
    "ALK": "2",
    "ROS1": "6",
    "RET": "10",
    "MET": "7",
    "ERBB2": "17",
    "HER2": "17",
    "BRCA1": "17",
    "BRCA2": "13",
    "APC": "5",
    "PTEN": "10",
    "KIT": "4",
    "PDGFRA": "4",
    "IDH1": "2",
    "IDH2": "15",
    "FGFR1": "8",
    "FGFR2": "10",
    "FGFR3": "4",
    "ATM": "11",
    "BRIP1": "17",
    "PALB2": "16",
    "CHEK2": "22",
    "CDH1": "16",
    "STK11": "19",
    "SMAD4": "18",
    "VHL": "3",
    "NF1": "17",
    "NF2": "22",
    "RB1": "13",
    "CDKN2A": "9",
    "MLH1": "3",
    "MSH2": "2",
    "MSH6": "2",
    "PMS2": "7",
    "ARID1A": "1",
    "NOTCH1": "9",
    "FBXW7": "4",
    "CTNNB1": "3",
    "AKT1": "14",
    "MTOR": "1",
    "ESR1": "6",
    "AR": "X",
    "JAK2": "9",
    "MPL": "1",
    "CALR": "19",
    "NPM1": "5",
    "FLT3": "13",
    "DNMT3A": "2",
    "TET2": "4",
    "ASXL1": "20",
    "SF3B1": "2",
    "U2AF1": "21",
    "SRSF2": "17",
    "RUNX1": "21",
    "CEBPA": "19",
    "WT1": "11",
    "GATA2": "3",
    "MYD88": "3",
    "BTK": "X",
    "BCL2": "18",
    "BCL6": "3",
    "MYC": "8",
    "CCND1": "11",
    "CDK4": "12",
    "CDK6": "7",
    "MDM2": "12",
    "TERT": "5",
}



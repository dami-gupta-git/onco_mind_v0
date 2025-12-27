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
    "mel": ["melanoma", "mel", "cutaneous melanoma", "skin melanoma", "skin"],
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

# =============================================================================
# CBIOPORTAL STUDY MAPPINGS
# =============================================================================
# Maps tumor types to cBioPortal study IDs
# First study in list is the primary study used for queries
# TCGA PanCancer Atlas studies are preferred for broad coverage
# MSK-IMPACT and institutional studies provide additional depth

CBIOPORTAL_STUDY_MAPPINGS: dict[str, list[str]] = {
    # Lung
    "nsclc": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018", "nsclc_mskcc_2018"],
    "lung": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"],
    "luad": ["luad_tcga_pan_can_atlas_2018"],
    "lusc": ["lusc_tcga_pan_can_atlas_2018"],
    "sclc": ["sclc_ucologne_2015"],

    # Melanoma
    "melanoma": ["skcm_tcga_pan_can_atlas_2018", "mel_ucla_2016", "skcm_mskcc_2014"],
    "skcm": ["skcm_tcga_pan_can_atlas_2018"],

    # Colorectal
    "colorectal": ["coadread_tcga_pan_can_atlas_2018", "crc_msk_2017"],
    "crc": ["coadread_tcga_pan_can_atlas_2018", "crc_msk_2017"],
    "colon": ["coadread_tcga_pan_can_atlas_2018"],
    "coad": ["coadread_tcga_pan_can_atlas_2018"],
    "rectal": ["coadread_tcga_pan_can_atlas_2018"],

    # Breast
    "breast": ["brca_tcga_pan_can_atlas_2018", "breast_msk_2018"],
    "brca": ["brca_tcga_pan_can_atlas_2018"],

    # Pancreatic
    "pancreatic": ["paad_tcga_pan_can_atlas_2018", "paad_qcmg_uq_2016"],
    "paad": ["paad_tcga_pan_can_atlas_2018"],
    "pancreas": ["paad_tcga_pan_can_atlas_2018"],

    # Brain
    "glioblastoma": ["gbm_tcga_pan_can_atlas_2018"],
    "gbm": ["gbm_tcga_pan_can_atlas_2018"],
    "glioma": ["lgg_tcga_pan_can_atlas_2018", "gbm_tcga_pan_can_atlas_2018"],
    "lgg": ["lgg_tcga_pan_can_atlas_2018"],

    # Ovarian
    "ovarian": ["ov_tcga_pan_can_atlas_2018"],
    "ov": ["ov_tcga_pan_can_atlas_2018"],

    # Prostate
    "prostate": ["prad_tcga_pan_can_atlas_2018", "prad_mskcc_2017"],
    "prad": ["prad_tcga_pan_can_atlas_2018"],

    # Bladder
    "bladder": ["blca_tcga_pan_can_atlas_2018"],
    "blca": ["blca_tcga_pan_can_atlas_2018"],
    "urothelial": ["blca_tcga_pan_can_atlas_2018"],

    # Kidney
    "kidney": ["kirc_tcga_pan_can_atlas_2018", "kirp_tcga_pan_can_atlas_2018"],
    "renal": ["kirc_tcga_pan_can_atlas_2018"],
    "rcc": ["kirc_tcga_pan_can_atlas_2018"],
    "kirc": ["kirc_tcga_pan_can_atlas_2018"],
    "kirp": ["kirp_tcga_pan_can_atlas_2018"],

    # Thyroid
    "thyroid": ["thca_tcga_pan_can_atlas_2018"],
    "thca": ["thca_tcga_pan_can_atlas_2018"],

    # Head and Neck
    "head and neck": ["hnsc_tcga_pan_can_atlas_2018"],
    "hnsc": ["hnsc_tcga_pan_can_atlas_2018"],

    # Liver
    "liver": ["lihc_tcga_pan_can_atlas_2018"],
    "lihc": ["lihc_tcga_pan_can_atlas_2018"],
    "hepatocellular": ["lihc_tcga_pan_can_atlas_2018"],
    "hcc": ["lihc_tcga_pan_can_atlas_2018"],

    # Gastric/Stomach
    "gastric": ["stad_tcga_pan_can_atlas_2018"],
    "stomach": ["stad_tcga_pan_can_atlas_2018"],
    "stad": ["stad_tcga_pan_can_atlas_2018"],

    # Esophageal
    "esophageal": ["esca_tcga_pan_can_atlas_2018"],
    "esca": ["esca_tcga_pan_can_atlas_2018"],

    # Uterine/Endometrial
    "uterine": ["ucec_tcga_pan_can_atlas_2018"],
    "endometrial": ["ucec_tcga_pan_can_atlas_2018"],
    "ucec": ["ucec_tcga_pan_can_atlas_2018"],

    # Cervical
    "cervical": ["cesc_tcga_pan_can_atlas_2018"],
    "cesc": ["cesc_tcga_pan_can_atlas_2018"],

    # Cholangiocarcinoma
    "cholangiocarcinoma": ["chol_tcga_pan_can_atlas_2018"],
    "chol": ["chol_tcga_pan_can_atlas_2018"],
    "bile duct": ["chol_tcga_pan_can_atlas_2018"],

    # Sarcoma
    "sarcoma": ["sarc_tcga_pan_can_atlas_2018"],
    "sarc": ["sarc_tcga_pan_can_atlas_2018"],

    # GIST
    "gist": ["gist_mskcc"],
    "gastrointestinal stromal": ["gist_mskcc"],

    # Mesothelioma
    "mesothelioma": ["meso_tcga_pan_can_atlas_2018"],
    "meso": ["meso_tcga_pan_can_atlas_2018"],

    # Adrenocortical
    "adrenocortical": ["acc_tcga_pan_can_atlas_2018"],
    "acc": ["acc_tcga_pan_can_atlas_2018"],

    # Pheochromocytoma
    "pheochromocytoma": ["pcpg_tcga_pan_can_atlas_2018"],
    "pcpg": ["pcpg_tcga_pan_can_atlas_2018"],

    # Testicular
    "testicular": ["tgct_tcga_pan_can_atlas_2018"],
    "tgct": ["tgct_tcga_pan_can_atlas_2018"],

    # Thymoma
    "thymoma": ["thym_tcga_pan_can_atlas_2018"],
    "thym": ["thym_tcga_pan_can_atlas_2018"],

    # Uveal Melanoma
    "uveal melanoma": ["uvm_tcga_pan_can_atlas_2018"],
    "uvm": ["uvm_tcga_pan_can_atlas_2018"],

    # AML
    "aml": ["laml_tcga_pan_can_atlas_2018"],
    "laml": ["laml_tcga_pan_can_atlas_2018"],
    "acute myeloid leukemia": ["laml_tcga_pan_can_atlas_2018"],

    # Diffuse Large B-Cell Lymphoma
    "dlbcl": ["dlbc_tcga_pan_can_atlas_2018"],
    "dlbc": ["dlbc_tcga_pan_can_atlas_2018"],
}

# Default pan-cancer study when no tumor-specific study is found
CBIOPORTAL_DEFAULT_STUDY = "msk_impact_2017"


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


# =============================================================================
# CANCER GENES FOR CO-OCCURRENCE ANALYSIS
# =============================================================================
# Curated list of cancer genes used for co-mutation analysis in cBioPortal.
# These genes are selected for biologically meaningful co-occurrence patterns.
# Used by: cbioportal.py for co-mutation queries

CANCER_GENES_CO_OCCURRENCE: list[str] = [
    # Top tumor suppressors (frequently co-mutated)
    "TP53",      # Most frequently mutated TSG, co-mutates with many genes
    "PTEN",      # Often co-mutates with PIK3CA
    "CDKN2A",    # Cell cycle regulator, co-mutates with TP53
    "RB1",       # Cell cycle TSG, co-mutates with TP53
    "SMAD4",     # Co-mutates with KRAS in pancreatic/CRC
    "STK11",     # Co-mutates with KRAS in lung adenocarcinoma
    "NF1",       # RAS pathway, co-mutates with TP53
    "APC",       # Colorectal driver
    "ARID1A",    # Chromatin remodeling, co-mutates with TP53
    "KEAP1",     # Oxidative stress pathway (lung)
    "ATM",       # DDR gene, co-mutates across tumor types
    "BRCA1", "BRCA2",  # DDR genes
    # Key oncogenes
    "KRAS",      # RAS pathway, often mutually exclusive with BRAF/NRAS
    "NRAS",      # RAS pathway, mutually exclusive with KRAS/BRAF
    "BRAF",      # MAPK pathway
    "PIK3CA",    # PI3K pathway, co-mutates with PTEN
    "EGFR",      # Receptor tyrosine kinase
    "ERBB2",     # HER2
    "MET",       # Often co-mutates as resistance mechanism
    "IDH1", "IDH2",  # Glioma/AML drivers
    "CTNNB1",    # Wnt pathway
    # Oxidative stress
    "NFE2L2",    # Co-mutates with KEAP1 in lung
    # Chromatin/epigenetic
    "KMT2D",     # Frequently mutated in solid tumors
    "ARID2",     # Chromatin remodeling
    # Receptor tyrosine kinases
    "FGFR1", "FGFR2", "FGFR3",
    "KIT", "PDGFRA",
    "ALK", "ROS1", "RET",
]


# =============================================================================
# DEPMAP FALLBACK DATA
# =============================================================================
# Pre-computed fallback data for DepMap queries when API is unavailable.
# WARNING: This is static data that may become outdated. It serves as a
# fallback cache only - real-time API queries are preferred.
#
# CANCER_GENE_DEPENDENCIES: CERES scores from CRISPR screens
#   - score: Mean dependency score (more negative = more essential)
#   - dependent_pct: % of cell lines where gene is essential
#
# GENE_DRUG_SENSITIVITIES: IC50 values from PRISM drug screens
#   - ic50_mutant: IC50 in nM for mutant cell lines
#   - ic50_wt: IC50 in nM for wild-type cell lines
#   - target: Drug target annotation

DEPMAP_GENE_DEPENDENCIES_FALLBACK: dict[str, dict] = {
    "BRAF": {"score": -0.8, "dependent_pct": 45},
    "KRAS": {"score": -1.2, "dependent_pct": 65},
    "EGFR": {"score": -0.6, "dependent_pct": 35},
    "TP53": {"score": -0.1, "dependent_pct": 5},  # TSG - not essential when lost
    "PIK3CA": {"score": -0.5, "dependent_pct": 30},
    "NRAS": {"score": -0.7, "dependent_pct": 40},
    "MYC": {"score": -1.5, "dependent_pct": 80},
    "ERBB2": {"score": -0.9, "dependent_pct": 50},
    "ALK": {"score": -0.4, "dependent_pct": 20},
    "RET": {"score": -0.3, "dependent_pct": 15},
    "MET": {"score": -0.5, "dependent_pct": 25},
    "FGFR1": {"score": -0.4, "dependent_pct": 20},
    "FGFR2": {"score": -0.5, "dependent_pct": 25},
    "FGFR3": {"score": -0.4, "dependent_pct": 20},
    "KIT": {"score": -0.6, "dependent_pct": 35},
    "PDGFRA": {"score": -0.3, "dependent_pct": 15},
}

DEPMAP_DRUG_SENSITIVITIES_FALLBACK: dict[str, list[dict]] = {
    "BRAF": [
        {"drug": "vemurafenib", "ic50_mutant": 50, "ic50_wt": 2000, "target": "BRAF V600"},
        {"drug": "dabrafenib", "ic50_mutant": 30, "ic50_wt": 1500, "target": "BRAF V600"},
        {"drug": "encorafenib", "ic50_mutant": 20, "ic50_wt": 1200, "target": "BRAF V600"},
        {"drug": "trametinib", "ic50_mutant": 8, "ic50_wt": 200, "target": "MEK1/2"},
        {"drug": "cobimetinib", "ic50_mutant": 15, "ic50_wt": 250, "target": "MEK1/2"},
    ],
    "KRAS": [
        {"drug": "sotorasib", "ic50_mutant": 100, "ic50_wt": None, "target": "KRAS G12C"},
        {"drug": "adagrasib", "ic50_mutant": 80, "ic50_wt": None, "target": "KRAS G12C"},
        {"drug": "trametinib", "ic50_mutant": 25, "ic50_wt": 300, "target": "MEK1/2"},
    ],
    "EGFR": [
        {"drug": "erlotinib", "ic50_mutant": 40, "ic50_wt": 3000, "target": "EGFR"},
        {"drug": "gefitinib", "ic50_mutant": 35, "ic50_wt": 2500, "target": "EGFR"},
        {"drug": "osimertinib", "ic50_mutant": 15, "ic50_wt": 500, "target": "EGFR T790M"},
        {"drug": "afatinib", "ic50_mutant": 20, "ic50_wt": 1000, "target": "EGFR/HER2"},
    ],
    "PIK3CA": [
        {"drug": "alpelisib", "ic50_mutant": 50, "ic50_wt": 500, "target": "PI3K alpha"},
        {"drug": "everolimus", "ic50_mutant": 5, "ic50_wt": 50, "target": "mTOR"},
        {"drug": "capivasertib", "ic50_mutant": 100, "ic50_wt": 800, "target": "AKT"},
    ],
    "ERBB2": [
        {"drug": "lapatinib", "ic50_mutant": 30, "ic50_wt": 2000, "target": "EGFR/HER2"},
        {"drug": "neratinib", "ic50_mutant": 10, "ic50_wt": 500, "target": "pan-HER"},
        {"drug": "tucatinib", "ic50_mutant": 8, "ic50_wt": 1000, "target": "HER2"},
    ],
}

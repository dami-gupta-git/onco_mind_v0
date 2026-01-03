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
    "nsclc": [
        "non-small cell lung", "nsclc", "lung non-small cell",
        "non small cell lung cancer", "non-small cell lung carcinoma",
        "lung adenocarcinoma", "lung squamous", "lung adenosquamous"
    ],
    "lung": ["lung", "lung cancer", "lung carcinoma"],
    "sclc": ["small cell lung", "sclc", "small cell lung cancer", "small cell lung carcinoma"],
    "luad": ["lung adenocarcinoma", "luad", "adenocarcinoma of lung"],
    "lusc": ["lung squamous", "lusc", "squamous cell lung", "squamous cell carcinoma of lung"],

    # Colorectal
    "crc": ["colorectal", "colon", "rectal", "crc", "colorectal cancer", "colorectal carcinoma"],
    "coad": ["colon adenocarcinoma", "coad", "colon cancer"],
    "read": ["rectal adenocarcinoma", "read", "rectum adenocarcinoma", "rectal cancer"],
    "coread": ["colorectal", "colon", "rectal", "colorectal adenocarcinoma", "crc", "coadread"],

    # Melanoma
    "mel": ["melanoma", "mel", "cutaneous melanoma", "skin melanoma"],
    "cm": ["cutaneous melanoma", "cm", "melanoma", "skin melanoma"],  # CGI uses CM
    "skcm": ["skin cutaneous melanoma", "skcm", "cutaneous melanoma"],
    "uvm": ["uveal melanoma", "uvm", "ocular melanoma"],

    # Breast
    "bc": ["breast", "bc", "breast cancer", "breast carcinoma"],
    "brca": ["breast invasive carcinoma", "brca", "breast cancer", "breast carcinoma"],
    "idc": ["invasive ductal carcinoma", "idc", "ductal carcinoma", "breast ductal"],
    "ilc": ["invasive lobular carcinoma", "ilc", "lobular carcinoma", "breast lobular"],
    "tnbc": ["triple negative breast cancer", "tnbc", "triple-negative breast"],

    # Thyroid
    "thca": ["thyroid carcinoma", "thca", "thyroid cancer"],
    "ptc": ["papillary thyroid carcinoma", "ptc", "papillary thyroid"],
    "ftc": ["follicular thyroid carcinoma", "ftc"],
    "atc": ["anaplastic thyroid carcinoma", "atc", "anaplastic thyroid cancer"],

    # Gastrointestinal
    "gist": ["gastrointestinal stromal tumor", "gist", "gastrointestinal stromal"],
    "stad": ["stomach adenocarcinoma", "stad", "gastric adenocarcinoma", "gastric cancer"],
    "esca": ["esophageal carcinoma", "esca", "esophageal cancer", "oesophageal"],
    "escc": ["esophageal squamous cell carcinoma", "escc"],
    "eadc": ["esophageal adenocarcinoma", "eadc", "barrett"],
    "paad": ["pancreatic adenocarcinoma", "paad", "pancreatic cancer", "pancreas adenocarcinoma"],
    "pdac": ["pancreatic ductal adenocarcinoma", "pdac", "pancreatic ductal"],
    "lihc": ["liver hepatocellular carcinoma", "lihc", "hepatocellular carcinoma", "hcc"],
    "hcc": ["hepatocellular carcinoma", "hcc"],
    "chol": ["cholangiocarcinoma", "chol", "bile duct cancer", "intrahepatic cholangiocarcinoma"],

    # Genitourinary
    "prad": ["prostate adenocarcinoma", "prad", "prostate cancer", "prostate carcinoma"],
    "prostate": ["prostate", "prostate cancer"],
    "blca": ["bladder urothelial carcinoma", "blca", "bladder cancer", "urothelial carcinoma"],
    "rcc": ["renal cell carcinoma", "rcc", "kidney cancer"],
    "ccrcc": ["clear cell renal cell carcinoma", "ccrcc", "clear cell kidney"],
    "prcc": ["papillary renal cell carcinoma", "prcc"],
    "chrcc": ["chromophobe renal cell carcinoma", "chrcc"],

    # Gynecologic
    "ov": ["ovarian carcinoma", "ov", "ovarian cancer", "ovary"],
    "hgsc": ["high grade serous ovarian", "hgsc", "serous ovarian"],
    "ucec": ["uterine corpus endometrial carcinoma", "ucec", "endometrial cancer", "endometrial carcinoma"],
    "ucs": ["uterine carcinosarcoma", "ucs"],
    "cesc": ["cervical squamous cell carcinoma", "cesc", "cervical cancer"],

    # Head and Neck
    "hnsc": [
        "head and neck squamous cell carcinoma", "hnsc",
        "head neck", "head and neck cancer", "oral cavity", "oropharynx", "larynx", "hypopharynx"
    ],

    # Brain / CNS
    "gbm": ["glioblastoma multiforme", "gbm", "glioblastoma"],
    "lgg": ["low grade glioma", "lgg", "diffuse glioma"],
    "astrocytoma": ["astrocytoma", "idhm"],
    "oligodendroglioma": ["oligodendroglioma", "oligo"],

    # Sarcoma
    "sarc": ["sarcoma", "sarc", "soft tissue sarcoma"],
    "ups": ["undifferentiated pleomorphic sarcoma", "ups"],
    "leiomyosarcoma": ["leiomyosarcoma", "lms"],
    "liposarcoma": ["liposarcoma"],

    # Hematologic
    "aml": ["acute myeloid leukemia", "aml"],
    "all": ["acute lymphoblastic leukemia", "all", "acute lymphocytic leukemia"],
    "cll": ["chronic lymphocytic leukemia", "cll", "sll"],
    "cml": ["chronic myeloid leukemia", "cml"],
    "dlbcl": ["diffuse large b-cell lymphoma", "dlbcl"],
    "fl": ["follicular lymphoma", "fl"],
    "mm": ["multiple myeloma", "mm", "plasma cell myeloma", "myeloma"],

    # Other / Rare
    "acc": ["adrenocortical carcinoma", "acc"],
    "meso": ["mesothelioma", "meso"],
    "tgct": ["testicular germ cell tumor", "tgct"],
    "thym": ["thymoma", "thym"],
    "pnet": ["pancreatic neuroendocrine tumor", "pnet", "pancreatic net"],
    "nets": ["neuroendocrine tumor", "nets", "net"],
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
    # Lung / NSCLC - Use TCGA PanCancer Atlas studies (verified working)
    "nsclc": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"],
    "lung": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"],
    "luad": ["luad_tcga_pan_can_atlas_2018"],
    "lusc": ["lusc_tcga_pan_can_atlas_2018"],
    "sclc": ["sclc_ucologne_2015"],

    # Melanoma
    "melanoma": ["mel_mskimpact_2020", "skcm_tcga_pan_can_atlas_2018", "mel_ucla_2016"],
    "skcm": ["skcm_tcga_pan_can_atlas_2018"],

    # Colorectal
    "colorectal": ["crc_msk_2017", "crc_eo_2020", "coadread_tcga_pan_can_atlas_2018"],
    "crc": ["crc_msk_2017", "crc_eo_2020", "coadread_tcga_pan_can_atlas_2018"],
    "colon": ["coadread_tcga_pan_can_atlas_2018"],
    "coad": ["coadread_tcga_pan_can_atlas_2018"],
    "rectal": ["coadread_tcga_pan_can_atlas_2018"],

    # Breast - Use latest MSK if available; METABRIC often used in aggregates
    "breast": ["brca_metabric", "breast_msk_2018", "brca_tcga_pan_can_atlas_2018"],
    "brca": ["brca_metabric", "brca_tcga_pan_can_atlas_2018"],

    # Pancreatic - Large 2024 MSK cohort first
    "pancreatic": ["pdac_msk_2024", "pancreas_msk_2024", "paad_qcmg_uq_2016"],
    "paad": ["pdac_msk_2024", "pancreas_msk_2024"],
    "pancreas": ["pdac_msk_2024", "pancreas_msk_2024"],

    # Brain
    "glioblastoma": ["gbm_tcga_pan_can_atlas_2018"],
    "gbm": ["gbm_tcga_pan_can_atlas_2018"],
    "glioma": ["lgg_tcga_pan_can_atlas_2018", "gbm_tcga_pan_can_atlas_2018"],
    "lgg": ["lgg_tcga_pan_can_atlas_2018"],

    # Ovarian
    "ovarian": ["ov_tcga_pan_can_atlas_2018"],
    "ov": ["ov_tcga_pan_can_atlas_2018"],

    # Prostate
    "prostate": ["prostate_msk_2024", "prad_mskcc_2017", "prad_tcga_pan_can_atlas_2018"],
    "prad": ["prostate_msk_2024", "prad_tcga_pan_can_atlas_2018"],

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
    "gist": ["gist_msk_2023"],
    "gastrointestinal stromal": ["gist_msk_2023"],

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

# CBIOPORTAL_STUDY_MAPPINGS: dict[str, list[str]] = {
#     # Lung / NSCLC
#     "nsclc": ["genie_bpc_nsclc", "nsclc_ctdx_msk_2022", "luad_mskcc_2023_or_later_if_available"],
#     "lung": ["genie_bpc_nsclc", "luad_tcga_pan_can_atlas_2018"],
#     "luad": ["genie_bpc_nsclc_luad_subset_or_multi_study_luad"],
#     "lusc": ["lusc_tcga_pan_can_atlas_2018"],
#     "sclc": ["sclc_ucologne_2015"],
#
#     # Melanoma
#     "melanoma": ["mel_mskimpact_2020", "skcm_tcga_pan_can_atlas_2018", "mel_ucla_2016", "skcm_mskcc_2014"],
#     "skcm": ["skcm_tcga_pan_can_atlas_2018"],
#
#     # Colorectal
#     "colorectal": ["crc_eo_2020", "crc_msk_2017", "coadread_tcga_pan_can_atlas_2018"],
#     "crc": ["crc_eo_2020", "crc_msk_2017", "coadread_tcga_pan_can_atlas_2018"],
#     "colon": ["coadread_tcga_pan_can_atlas_2018"],
#     "coad": ["coadread_tcga_pan_can_atlas_2018"],
#     "rectal": ["coadread_tcga_pan_can_atlas_2018"],
#
#     # Breast
#     "breast": ["breast_msk_2025", "breast_msk_2018", "brca_tcga_pan_can_atlas_2018"],
#     "brca": ["brca_tcga_pan_can_atlas_2018"],
#
#     # Pancreatic - Updated with latest large cohorts (as of late 2025)
#     "pancreatic": ["pdac_msk_2024", "pancreas_msk_2024", "paad_qcmg_uq_2016"],
#     "paad": ["pdac_msk_2024", "pancreas_msk_2024"],
#     "pancreas": ["pdac_msk_2024", "pancreas_msk_2024"],
#
#     # Brain
#     "glioblastoma": ["gbm_tcga_pan_can_atlas_2018"],
#     "gbm": ["gbm_tcga_pan_can_atlas_2018"],
#     "glioma": ["lgg_tcga_pan_can_atlas_2018", "gbm_tcga_pan_can_atlas_2018"],
#     "lgg": ["lgg_tcga_pan_can_atlas_2018"],
#
#     # Ovarian
#     "ovarian": ["ov_tcga_pan_can_atlas_2018"],
#     "ov": ["ov_tcga_pan_can_atlas_2018"],
#
#     # Prostate
#     "prostate": ["prostate_msk_2024", "prad_mskcc_2017", "prad_tcga_pan_can_atlas_2018"],
#     "prad": ["prostate_msk_2024", "prad_tcga_pan_can_atlas_2018"],
#
#     # Bladder
#     "bladder": ["blca_tcga_pan_can_atlas_2018"],
#     "blca": ["blca_tcga_pan_can_atlas_2018"],
#     "urothelial": ["blca_tcga_pan_can_atlas_2018"],
#
#     # Kidney
#     "kidney": ["kirc_tcga_pan_can_atlas_2018", "kirp_tcga_pan_can_atlas_2018"],
#     "renal": ["kirc_tcga_pan_can_atlas_2018"],
#     "rcc": ["kirc_tcga_pan_can_atlas_2018"],
#     "kirc": ["kirc_tcga_pan_can_atlas_2018"],
#     "kirp": ["kirp_tcga_pan_can_atlas_2018"],
#
#     # Thyroid
#     "thyroid": ["thca_tcga_pan_can_atlas_2018"],
#     "thca": ["thca_tcga_pan_can_atlas_2018"],
#
#     # Head and Neck
#     "head and neck": ["hnsc_tcga_pan_can_atlas_2018"],
#     "hnsc": ["hnsc_tcga_pan_can_atlas_2018"],
#
#     # Liver
#     "liver": ["lihc_tcga_pan_can_atlas_2018"],
#     "lihc": ["lihc_tcga_pan_can_atlas_2018"],
#     "hepatocellular": ["lihc_tcga_pan_can_atlas_2018"],
#     "hcc": ["lihc_tcga_pan_can_atlas_2018"],
#
#     # Gastric/Stomach
#     "gastric": ["stad_tcga_pan_can_atlas_2018"],
#     "stomach": ["stad_tcga_pan_can_atlas_2018"],
#     "stad": ["stad_tcga_pan_can_atlas_2018"],
#
#     # Esophageal
#     "esophageal": ["esca_tcga_pan_can_atlas_2018"],
#     "esca": ["esca_tcga_pan_can_atlas_2018"],
#
#     # Uterine/Endometrial
#     "uterine": ["ucec_tcga_pan_can_atlas_2018"],
#     "endometrial": ["ucec_tcga_pan_can_atlas_2018"],
#     "ucec": ["ucec_tcga_pan_can_atlas_2018"],
#
#     # Cervical
#     "cervical": ["cesc_tcga_pan_can_atlas_2018"],
#     "cesc": ["cesc_tcga_pan_can_atlas_2018"],
#
#     # Cholangiocarcinoma
#     "cholangiocarcinoma": ["chol_tcga_pan_can_atlas_2018"],
#     "chol": ["chol_tcga_pan_can_atlas_2018"],
#     "bile duct": ["chol_tcga_pan_can_atlas_2018"],
#
#     # Sarcoma
#     "sarcoma": ["sarc_tcga_pan_can_atlas_2018"],
#     "sarc": ["sarc_tcga_pan_can_atlas_2018"],
#
#     # GIST
#     "gist": ["gist_mskcc"],
#     "gastrointestinal stromal": ["gist_mskcc"],
#
#     # Mesothelioma
#     "mesothelioma": ["meso_tcga_pan_can_atlas_2018"],
#     "meso": ["meso_tcga_pan_can_atlas_2018"],
#
#     # Adrenocortical
#     "adrenocortical": ["acc_tcga_pan_can_atlas_2018"],
#     "acc": ["acc_tcga_pan_can_atlas_2018"],
#
#     # Pheochromocytoma
#     "pheochromocytoma": ["pcpg_tcga_pan_can_atlas_2018"],
#     "pcpg": ["pcpg_tcga_pan_can_atlas_2018"],
#
#     # Testicular
#     "testicular": ["tgct_tcga_pan_can_atlas_2018"],
#     "tgct": ["tgct_tcga_pan_can_atlas_2018"],
#
#     # Thymoma
#     "thymoma": ["thym_tcga_pan_can_atlas_2018"],
#     "thym": ["thym_tcga_pan_can_atlas_2018"],
#
#     # Uveal Melanoma
#     "uveal melanoma": ["uvm_tcga_pan_can_atlas_2018"],
#     "uvm": ["uvm_tcga_pan_can_atlas_2018"],
#
#     # AML
#     "aml": ["laml_tcga_pan_can_atlas_2018"],
#     "laml": ["laml_tcga_pan_can_atlas_2018"],
#     "acute myeloid leukemia": ["laml_tcga_pan_can_atlas_2018"],
#
#     # Diffuse Large B-Cell Lymphoma
#     "dlbcl": ["dlbc_tcga_pan_can_atlas_2018"],
#     "dlbc": ["dlbc_tcga_pan_can_atlas_2018"],
# }

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


LITERATURE_SENSITIVITY_TERMS = [
    'sensitivity', 'sensitive', 'response',
    'efficacy', 'effective', 'benefit',
]


# =============================================================================
# AMBIGUOUS VARIANTS
# =============================================================================
# Variants that are ambiguous - same notation exists in multiple genes or
# the variant notation is too broad to be considered specific.
# These should have scope="ambiguous" when variant_level.level="variant"
# Format: (gene, variant) tuples with exact match required

BROAD_VARIANTS: set[tuple[str, str]] = {
    ("BRAF", "V600"), ("IDH1", "R132"), ("IDH2", "R140"),
    ("IDH2", "R172"), ("EGFR", "exon 19 deletion"), ("PIK3CA", "E542"),
    ("PIK3CA", "E545"), ("PIK3CA", "H1047"),
    ("KRAS", "G12"), ("NRAS", "G12"), ("HRAS", "G12"),
    ("KRAS", "G13"), ("NRAS", "G13"), ("HRAS", "G13"),
    ("KRAS", "Q61"), ("NRAS", "Q61"), ("HRAS", "Q61"),
    ("PIK3CA", "hotspot"), ("MET", "exon 14 skipping"),
}
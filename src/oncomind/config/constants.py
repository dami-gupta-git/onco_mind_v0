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
    "l": ["lung", "lung cancer", "lung carcinoma"],  # CGI uses single-letter "L"
    "sclc": ["small cell lung", "sclc", "small cell lung cancer", "small cell lung carcinoma"],
    "luad": ["lung adenocarcinoma", "luad", "adenocarcinoma of lung"],
    "lusc": ["lung squamous", "lusc", "squamous cell lung", "squamous cell carcinoma of lung"],

    # Colorectal
    "crc": ["colorectal", "colon", "rectal", "crc", "mcrc", "colorectal cancer", "colorectal carcinoma"],
    "coad": ["colon adenocarcinoma", "coad", "colon cancer"],
    "read": ["rectal adenocarcinoma", "read", "rectum adenocarcinoma", "rectal cancer"],
    "coread": ["colorectal", "colon", "rectal", "colorectal adenocarcinoma", "crc", "coadread"],

    # Melanoma
    "mel": ["melanoma", "mel", "cutaneous melanoma", "skin melanoma"],
    "cm": ["cutaneous melanoma", "cm", "melanoma", "skin melanoma"],  # CGI uses CM
    "skcm": ["skin cutaneous melanoma", "skcm", "cutaneous melanoma"],
    "uvm": ["uveal melanoma", "uvm", "ocular melanoma", "eye cancer", "eye"],

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
    "chol": ["cholangiocarcinoma", "chol", "bile duct cancer", "biliary", "intrahepatic cholangiocarcinoma"],

    # Genitourinary
    "prad": ["prostate adenocarcinoma", "prad", "prostate cancer", "prostate carcinoma"],
    "prostate": ["prostate", "prostate cancer"],
    "blca": ["bladder urothelial carcinoma", "blca", "bladder cancer", "urothelial carcinoma", "transitional cell"],
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
    # Myeloproliferative neoplasms
    "mpn": ["myeloproliferative neoplasm", "mpn", "myeloproliferative"],
    "mf": ["myelofibrosis", "mf", "primary myelofibrosis"],
    "pv": ["polycythemia vera", "pv"],

    # Other / Rare
    "acc": ["adrenocortical carcinoma", "acc"],
    "meso": ["mesothelioma", "meso"],
    "tgct": ["testicular germ cell tumor", "tgct"],
    "thym": ["thymoma", "thym"],
    "pnet": ["pancreatic neuroendocrine tumor", "pnet", "pancreatic net"],
    "nets": ["neuroendocrine tumor", "nets", "net"],
}


# =============================================================================
# PAN-CANCER TERMS
# =============================================================================
# Generic tumor-agnostic disease terms that should NOT be considered a specific
# tumor type match. When evidence is labeled with these terms, it applies broadly
# across tumor types and should be marked as "pan_cancer" rather than "cancer_specific".

# Terms that require EXACT match (case-insensitive)
# These are generic terms that are too short/common for substring matching
PAN_CANCER_TERMS: set[str] = {
    "cancer",
    "solid tumor",
    "solid tumors",
    "solid tumour",
    "solid tumours",
    "advanced solid tumor",
    "advanced solid tumors",
    "malignant neoplasm",
    "malignant neoplasms",
    "neoplasm",
    "tumor",
    "tumour",
    "all tumor types",
    "any tumor",
    "any cancer",
}

# Terms that use SUBSTRING match (for FDA indications like "MSI-H Solid Tumors")
# These are specific enough to safely use substring matching
PAN_CANCER_BIOMARKERS: set[str] = {
    "msi-h",
    "msi-high",
    "dmmr",
    "tumor-agnostic",
    "tumour-agnostic",
    "ntrk",
    "tmb-h",
    "tmb-high",
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
# ACQUIRED RESISTANCE MUTATIONS
# =============================================================================
# Mutations that arise under treatment pressure (not present at baseline).
# These won't appear in cBioPortal/TCGA datasets (which sequence treatment-naive tumors).
# Used to skip prevalence gap flagging and provide appropriate context to LLM.
#
# Format: {gene: [list of acquired resistance variants]}
# References:
#   - EGFR T790M: 50-60% of EGFR-TKI resistant NSCLC (PMID: 15737014, 16014882)
#   - EGFR C797S: Osimertinib resistance (PMID: 26286086)
#   - ALK resistance: Crizotinib/alectinib resistance (PMID: 24675041, 28425994)
#   - ROS1 resistance: Crizotinib resistance (PMID: 26698910)
#   - ESR1: Aromatase inhibitor resistance in breast cancer (PMID: 26698910)

ACQUIRED_RESISTANCE_MUTATIONS: dict[str, list[str]] = {
    # EGFR - TKI resistance mutations
    "EGFR": [
        "T790M",   # 1st/2nd-gen TKI resistance (gefitinib, erlotinib, afatinib)
        "C797S",   # 3rd-gen TKI resistance (osimertinib)
        "G724S",   # Osimertinib resistance
        "L718Q",   # Osimertinib resistance
        "L718V",   # Osimertinib resistance
        "L792F",   # Osimertinib resistance
        "G796D",   # Osimertinib resistance
    ],
    # ALK - ALK inhibitor resistance mutations
    "ALK": [
        "G1202R",  # Pan-ALK inhibitor resistance
        "L1196M",  # Crizotinib resistance (gatekeeper)
        "I1171T",  # Alectinib resistance
        "I1171N",  # Alectinib resistance
        "F1174L",  # Crizotinib resistance
        "F1174C",  # Crizotinib resistance
        "C1156Y",  # Crizotinib resistance
        "L1152R",  # Crizotinib resistance
        "G1269A",  # Crizotinib resistance
    ],
    # ROS1 - ROS1 inhibitor resistance mutations
    "ROS1": [
        "G2032R",  # Crizotinib resistance (solvent front)
        "D2033N",  # Crizotinib resistance
        "L2026M",  # Crizotinib resistance (gatekeeper)
        "S1986F",  # Entrectinib/lorlatinib resistance
        "S1986Y",  # Entrectinib/lorlatinib resistance
    ],
    # BRAF - Acquired splice variants and resistance mutations
    "BRAF": [
        # Splice variants causing vemurafenib/dabrafenib resistance
        # Note: These are typically detected as RNA alterations
    ],
    # ESR1 - Aromatase inhibitor resistance in HR+ breast cancer
    "ESR1": [
        "D538G",   # Most common AI resistance mutation
        "Y537S",   # AI resistance
        "Y537N",   # AI resistance
        "Y537C",   # AI resistance
        "E380Q",   # AI resistance
        "L536Q",   # AI resistance
    ],
    # KIT - Imatinib resistance in GIST (secondary mutations in exon 17)
    "KIT": [
        "D816V",   # Primary imatinib resistance (also germline in mastocytosis)
        "D816H",   # Imatinib resistance
        "D816Y",   # Imatinib resistance
        "D820Y",   # Sunitinib resistance
        "N822K",   # Imatinib resistance
        "Y823D",   # Secondary imatinib resistance
    ],
    # FGFR2 - Acquired resistance to FGFR inhibitors
    "FGFR2": [
        "N549H",   # Pemigatinib/infigratinib resistance
        "N549K",   # Pemigatinib/infigratinib resistance
        "V564F",   # Gatekeeper - FGFR inhibitor resistance
        "V564I",   # Gatekeeper
        "E565A",   # FGFR inhibitor resistance
        "L617V",   # FGFR inhibitor resistance
        "K659M",   # FGFR inhibitor resistance
    ],
    # BTK - Ibrutinib resistance in CLL/lymphoma
    "BTK": [
        "C481S",   # Ibrutinib resistance (covalent binding site)
        "C481R",   # Ibrutinib resistance
        "C481F",   # Ibrutinib resistance
    ],
    # BCR-ABL1 - TKI resistance in CML
    "ABL1": [
        "T315I",   # Pan-TKI resistance (gatekeeper) - except ponatinib
        "F317L",   # Dasatinib resistance
        "V299L",   # Dasatinib resistance
        "Y253H",   # Imatinib resistance
        "E255K",   # Imatinib resistance
        "E255V",   # Imatinib resistance
    ],
}

# Flattened set for quick lookup: {(gene, variant), ...}
ACQUIRED_RESISTANCE_MUTATIONS_SET: set[tuple[str, str]] = {
    (gene.upper(), variant.upper())
    for gene, variants in ACQUIRED_RESISTANCE_MUTATIONS.items()
    for variant in variants
}


def is_acquired_resistance_mutation(gene: str, variant: str) -> bool:
    """Check if a gene/variant is a known acquired resistance mutation.

    Args:
        gene: Gene symbol (e.g., "EGFR")
        variant: Variant notation (e.g., "T790M")

    Returns:
        True if this is a known acquired resistance mutation
    """
    return (gene.upper(), variant.upper()) in ACQUIRED_RESISTANCE_MUTATIONS_SET


# =============================================================================
# CBIOPORTAL STUDY MAPPINGS
# =============================================================================
CBIOPORTAL_STUDY_MAPPINGS: dict[str, list[str]] = {
    # Lung / NSCLC – prioritize recent MSK + GENIE BPC for clinical depth
    "nsclc": [
        "nsclc_msk_2023",                    # Large recent MSK-IMPACT NSCLC cohort
        "lung_msk_2024",                     # Broader recent MSK lung
        "nsclc_public_genie_bpc",            # GENIE BPC NSCLC v2.0-public (rich clinical annotations)
        "genie_public",                      # Latest GENIE public release (18.0-public as of mid-2025)
        "luad_tcga_pan_can_atlas_2018",       # Reliable harmonized fallback
        "lusc_tcga_pan_can_atlas_2018"
    ],
    "lung": [
        "lung_msk_2024",
        "nsclc_public_genie_bpc",
        "genie_public",
        "luad_tcga_pan_can_atlas_2018",
        "lusc_tcga_pan_can_atlas_2018"
    ],
    "luad": [
        "lung_msk_2024",                     # Covers LUAD well in recent MSK
        "nsclc_public_genie_bpc",
        "luad_tcga_pan_can_atlas_2018"
    ],
    "lusc": [
        "lusc_msk_2023",                     # Recent MSK squamous if available
        "lusc_tcga_pan_can_atlas_2018"
    ],
    "sclc": [
        "sclc_msk_2024",                     # Recent MSK SCLC if available
        "sclc_ucologne_2015"                 # Classic small reference
    ],

    # Melanoma
    "melanoma": ["mel_mskimpact_2020", "skcm_tcga_pan_can_atlas_2018"],
    "skcm": ["skcm_tcga_pan_can_atlas_2018"],

    # Colorectal – replace old 2017 with TCGA (no major newer public MSK-specific)
    "colorectal": ["coadread_tcga_pan_can_atlas_2018"],
    "crc": ["coadread_tcga_pan_can_atlas_2018"],
    "colon": ["coadread_tcga_pan_can_atlas_2018"],
    "coad": ["coadread_tcga_pan_can_atlas_2018"],
    "rectal": ["coadread_tcga_pan_can_atlas_2018"],

    # Breast
    "breast": ["breast_msk_2025", "brca_metabric", "brca_tcga_pan_can_atlas_2018"],
    "brca": ["breast_msk_2025", "brca_tcga_pan_can_atlas_2018"],

    # Pancreatic
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
    "prostate": ["prostate_msk_2024", "prad_tcga_pan_can_atlas_2018"],
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
#
# CBIOPORTAL_STUDY_MAPPINGS: dict[str, list[str]] = {
#     # Lung / NSCLC – keep TCGA PanCan; consider adding msk_impact_2017 for pan‑solid prevalence if you want
#     "nsclc": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"],
#     "lung": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"],
#     "luad": ["luad_tcga_pan_can_atlas_2018"],
#     "lusc": ["lusc_tcga_pan_can_atlas_2018"],
#     "sclc": ["sclc_ucologne_2015"],
#
#     # Melanoma
#     "melanoma": ["mel_mskimpact_2020", "skcm_tcga_pan_can_atlas_2018"],
#     "skcm": ["skcm_tcga_pan_can_atlas_2018"],
#
#     # Colorectal
#     "colorectal": ["crc_msk_2017", "coadread_tcga_pan_can_atlas_2018"],
#     "crc": ["crc_msk_2017", "coadread_tcga_pan_can_atlas_2018"],
#     "colon": ["coadread_tcga_pan_can_atlas_2018"],
#     "coad": ["coadread_tcga_pan_can_atlas_2018"],
#     "rectal": ["coadread_tcga_pan_can_atlas_2018"],
#
#     # Breast – updated
#     "breast": ["breast_msk_2025", "brca_metabric", "brca_tcga_pan_can_atlas_2018"],
#     "brca": ["breast_msk_2025", "brca_tcga_pan_can_atlas_2018"],
#
#     # Pancreatic – keep 2024 MSK cohorts
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
#     # Prostate – keep 2024 MSK cohort
#     "prostate": ["prostate_msk_2024", "prad_tcga_pan_can_atlas_2018"],
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
#     # GIST – keep MSK GIST cohort
#     "gist": ["gist_msk_2023"],
#     "gastrointestinal stromal": ["gist_msk_2023"],
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
# GAP DETECTION THRESHOLDS
# =============================================================================
# Magic numbers used in gap_detector.py for evidence quality assessment

# gnomAD allele frequency threshold (0.01%) - variants above this are
# considered "observed in general population"
GNOMAD_COMMON_AF_THRESHOLD = 0.0001

# Minimum number of publications to be considered "well characterized"
# 0 pubs = gap, 1-4 pubs = poorly characterized, 5+ = well characterized
LITERATURE_WELL_CHARACTERIZED_THRESHOLD = 5

# Co-occurrence percentage threshold for "strong" co-mutation signal
COOCCURRENCE_STRONG_THRESHOLD_PCT = 20.0

# Hotspot adjacency window (codons) - variants within this distance
# of a known hotspot are flagged as "near hotspot"
HOTSPOT_ADJACENCY_WINDOW = 5


# =============================================================================
# AMBIGUOUS VARIANTS
# =============================================================================
# Variants that are ambiguous - same notation exists in multiple genes or
# the variant notation is too broad to be considered specific.
# These should have scope="ambiguous" when locus_match.level="variant"
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


# =============================================================================
# LLM SERVICE CONSTANTS
# =============================================================================
# Configuration for LLM-based evidence synthesis and paper scoring

# Default model for LLM calls
LLM_DEFAULT_MODEL = "claude-sonnet-4-20250514"

# Fast/cheap model for high-volume operations (paper scoring, knowledge extraction)
LLM_FAST_MODEL = "gpt-4o-mini"

# Timeout settings (seconds)
LLM_TIMEOUT_DEFAULT = 60
LLM_TIMEOUT_CLAUDE = 120

# Token limits for different LLM tasks
LLM_MAX_TOKENS_SYNTHESIS = 1500
LLM_MAX_TOKENS_HYPOTHESIS = 800
LLM_MAX_TOKENS_PAPER_SCORING = 500
LLM_MAX_TOKENS_KNOWLEDGE_EXTRACTION = 1500

# Content truncation limits (characters)
LLM_PAPER_CONTENT_TRUNCATION = 1500
LLM_PAPER_ABSTRACT_TRUNCATION = 1000

# Paper processing limits
LLM_MAX_PAPERS_FOR_EXTRACTION = 5

# Relevance score threshold for paper filtering
LLM_PAPER_RELEVANCE_THRESHOLD = 0.6

# Default temperature for LLM calls
LLM_DEFAULT_TEMPERATURE = 0.1

# =============================================================================
# EVIDENCE AGGREGATOR CONSTANTS
# =============================================================================
# Result limits for evidence fetching
MAX_VICC_RESULTS = 500

# VICC sources that have dedicated tabs/counts elsewhere
# These should be filtered out when displaying VICC to avoid double-counting
# VICC aggregates from: CIViC, CGI, JAX-CKB, OncoKB, PMKB, MolecularMatch
VICC_DUPLICATE_SOURCES: set[str] = {"civic", "cgi"}
MAX_CIVIC_ASSERTIONS = 500
MAX_CLINICAL_TRIALS = 500
MAX_LITERATURE_RESULTS = 30

# =============================================================================
# UI DISPLAY CONSTANTS
# =============================================================================
# Row limits for evidence display tables
UI_MAX_CIVIC_EVIDENCE_ROWS = 50
UI_MAX_RESEARCH_HYPOTHESES = 3
UI_MAX_REFERENCES = 5

# Functional prediction thresholds
CADD_DELETERIOUS_THRESHOLD = 20  # CADD score > 20 considered deleterious
GNOMAD_RARE_THRESHOLD = 0.01  # gnomAD AF < 0.01 considered rare


# =============================================================================
# DEPMAP CONSTANTS
# =============================================================================
# Configuration for DepMap API and data access

# API endpoints
DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/download/api/download"
DEPMAP_CUSTOM_DOWNLOAD_API = "https://depmap.org/portal/api/download/custom"
DEPMAP_TASK_STATUS_API = "https://depmap.org/portal/api/task"

# PRISM drug sensitivity data URLs (24Q2 release)
DEPMAP_PRISM_SENSITIVITY_URL = (
    "https://depmap.org/portal/download/api/download?"
    "file_name=Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv"
)
DEPMAP_PRISM_DRUG_INFO_URL = (
    "https://depmap.org/portal/download/api/download?"
    "file_name=Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"
)

# Data file names for local cache
DEPMAP_MUTATIONS_FILE = "OmicsSomaticMutations.csv"
DEPMAP_CRISPR_FILE = "CRISPRGeneEffect.csv"

# DepMap release version
DEPMAP_RELEASE = "DepMap+Public+24Q4"
DEPMAP_DATA_VERSION = "DepMap Public 24Q4"

# Timeout for DepMap requests (seconds)
DEPMAP_DEFAULT_TIMEOUT = 60.0
DEPMAP_TASK_POLL_TIMEOUT = 10.0
DEPMAP_DOWNLOAD_TIMEOUT = 30.0

# Task polling configuration
DEPMAP_MAX_TASK_POLL_ATTEMPTS = 10
DEPMAP_TASK_POLL_INTERVAL = 1  # seconds

# Drug sensitivity thresholds
DEPMAP_SENSITIVITY_THRESHOLD = -1.7  # Log2FC threshold for sensitivity
DEPMAP_MIN_CELL_LINES = 3  # Minimum cell lines required for drug inclusion
DEPMAP_TOP_DRUGS_LIMIT = 10  # Maximum number of drugs to return
DEPMAP_SENSITIVE_LINES_DISPLAY = 5  # Max sensitive cell lines to display per drug

# Gene dependency threshold
DEPMAP_DEPENDENCY_THRESHOLD = -0.5  # Score below this = dependent


# =============================================================================
# FDA ONCOLOGY FILTER TERMS
# =============================================================================
# Terms used to identify oncology drugs and filter out false positives
# (e.g., "eGFR" kidney function vs "EGFR" cancer gene)

FDA_ONCOLOGY_TERMS = [
    # Cancer types
    "cancer", "tumor", "tumour", "carcinoma", "sarcoma", "lymphoma",
    "leukemia", "leukaemia", "melanoma", "neoplasm", "malignant",
    "myeloma", "glioma", "glioblastoma", "blastoma", "mesothelioma",
    # Clinical terms
    "metastatic", "oncology", "chemotherapy", "antineoplastic",
    "locally advanced", "unresectable", "refractory",
    # Biomarker terms
    "mutation-positive", "mutation", "mutated", "biomarker", "pd-l1",
    "msi-h", "microsatellite", "tmb", "tumor mutational",
    # Drug class terms
    "kinase inhibitor", "tyrosine kinase", "checkpoint inhibitor",
    "immunotherapy", "targeted therapy", "antibody-drug conjugate",
    "monoclonal antibody", "bispecific",
    # Specific cancer contexts
    "egfr-positive", "egfr mutation", "non-small cell lung",
    "nsclc", "sclc", "gist", "gastrointestinal stromal",
    "colorectal", "breast cancer", "prostate cancer", "pancreatic",
    "hepatocellular", "cholangiocarcinoma", "ovarian", "endometrial",
]

# Terms indicating non-oncology context for EGFR (kidney function "eGFR")
FDA_EGFR_EXCLUSION_TERMS = [
    "glomerular filtration", "egfr ml/min", "egfr <", "egfr >",
    "renal function", "kidney", "diabetes", "type 2 diabetes",
    "glycemic", "glucose", "insulin", "sglt2", "metformin",
    "dapagliflozin", "empagliflozin", "canagliflozin",
]

# Terms indicating non-oncology context for KIT gene
# "kit" is a common English word that appears in many drug labels
# (e.g., "Kit for the preparation of...", "test kit", "diagnostic kit")
# Only exclude if these patterns appear WITHOUT oncology context
FDA_KIT_EXCLUSION_PATTERNS = [
    "kit for the preparation",
    "diagnostic kit",
    "test kit",
    "reagent kit",
    "assay kit",
    "collection kit",
    "constipation",  # Lubiprostone
    "irritable bowel",
    "opioid-induced",
]


# =============================================================================
# BIOMARKER SELECTION DRUGS
# =============================================================================
# Drugs where the queried gene is a PATIENT SELECTION BIOMARKER, not the drug target.
# These drugs are approved for patients with mutations in certain genes, but the drug
# mechanism targets a different protein.
#
# Example: DATROWAY (datopotamab deruxtecan) is a TROP2 ADC approved for EGFR-mutant
# NSCLC patients who progressed on EGFR therapy. The drug targets TROP2, not EGFR.
# Without this mapping, it would incorrectly appear as an "EGFR drug".
#
# Format: {drug_name_lowercase: {"target": actual_target, "biomarker_genes": [genes]}}
#
# TODO: Replace with DGIdb integration for comprehensive drug-target relationships.
# See TODO.md for details.

BIOMARKER_SELECTION_DRUGS: dict[str, dict] = {
    # TROP2 ADCs - indicated for EGFR-mutant NSCLC but target TROP2
    "datopotamab deruxtecan": {"target": "TROP2", "biomarker_genes": ["EGFR"]},
    "datroway": {"target": "TROP2", "biomarker_genes": ["EGFR"]},
    # Chemotherapies often indicated alongside biomarker-selected patients
    "pemetrexed": {"target": "TYMS", "biomarker_genes": ["EGFR", "ALK"]},
    "pemetrexed disodium": {"target": "TYMS", "biomarker_genes": ["EGFR", "ALK"]},
    # VEGFR2 inhibitors - ramucirumab (Cyramza) approved with erlotinib for EGFR-mutant NSCLC
    # but targets VEGFR2, not EGFR
    "ramucirumab": {"target": "VEGFR2", "biomarker_genes": ["EGFR"]},
    "cyramza": {"target": "VEGFR2", "biomarker_genes": ["EGFR"]},
    # PD-1 inhibitors - approved based on PD-L1 expression or TMB-high, not specific mutations
    # These should not appear as "EGFR drugs" when querying EGFR variants
    "pembrolizumab": {"target": "PD-1", "biomarker_genes": ["EGFR", "KRAS", "BRAF", "ALK", "ROS1", "RET", "MET", "NTRK1", "NTRK2", "NTRK3"]},
    "keytruda": {"target": "PD-1", "biomarker_genes": ["EGFR", "KRAS", "BRAF", "ALK", "ROS1", "RET", "MET", "NTRK1", "NTRK2", "NTRK3"]},
    "keytruda qlex": {"target": "PD-1", "biomarker_genes": ["EGFR", "KRAS", "BRAF", "ALK", "ROS1", "RET", "MET", "NTRK1", "NTRK2", "NTRK3"]},
    "pembrolizumab and hyaluronidase": {"target": "PD-1", "biomarker_genes": ["EGFR", "KRAS", "BRAF", "ALK", "ROS1", "RET", "MET", "NTRK1", "NTRK2", "NTRK3"]},
    # VEGFR inhibitors approved for CRC - FDA label mentions EGFR therapy for patient selection
    # (RAS wild-type CRC patients should have received anti-EGFR therapy before fruquintinib)
    # This is NOT an EGFR-targeting drug - it targets VEGFR1/2/3
    "fruquintinib": {"target": "VEGFR", "biomarker_genes": ["EGFR", "KRAS", "RAS"]},
    "fruzaqla": {"target": "VEGFR", "biomarker_genes": ["EGFR", "KRAS", "RAS"]},
}


def is_biomarker_selection_drug(drug_name: str, gene: str) -> bool:
    """Check if a drug uses the gene as a biomarker for patient selection, not as its target.

    For example, datopotamab deruxtecan (DATROWAY) is approved for EGFR-mutant NSCLC,
    but it targets TROP2, not EGFR. For EGFR queries, we should filter it out because
    it's not an "EGFR drug" - it's a drug for patients who happen to have EGFR mutations.

    Args:
        drug_name: Drug name (case-insensitive)
        gene: Gene symbol to check (e.g., "EGFR")

    Returns:
        True if this drug uses the gene as a biomarker (should be filtered for that gene)
    """
    if not drug_name or not gene:
        return False

    drug_lower = drug_name.lower()
    gene_upper = gene.upper()

    for drug_key, info in BIOMARKER_SELECTION_DRUGS.items():
        # Check if drug name matches (substring match in both directions)
        if drug_key in drug_lower or drug_lower in drug_key:
            if gene_upper in [g.upper() for g in info.get("biomarker_genes", [])]:
                return True
    return False


# =============================================================================
# LOCUS MATCH LEVEL DESCRIPTIONS
# =============================================================================
# Human-readable descriptions for locus match levels (variant, codon, gene)
# Used for hover text/tooltips in the UI to explain what each level means
#
# These descriptions make clear that codon/gene level matches are EXTRAPOLATIONS
# from different variants, not direct evidence for the queried variant.

def get_locus_match_description(
    level: str,
    gene: str | None = None,
    variant: str | None = None,
    codon: str | int | None = None,
) -> str:
    """Get human-readable description for a locus match level.

    Args:
        level: The locus match level ("variant", "codon", or "gene")
        gene: Gene symbol for substitution in description
        variant: Variant notation for substitution in description
        codon: Codon number for substitution in codon-level description

    Returns:
        Human-readable description suitable for hover text/tooltip
    """
    if level == "variant":
        if variant:
            return f"Direct evidence for {variant}"
        return "Direct evidence for this exact variant"

    if level == "codon":
        if codon:
            return f"Evidence for other variants at codon {codon} — extrapolation from similar position"
        return "Evidence for other variants at this codon — extrapolation from similar position"

    if level == "gene":
        if gene and variant:
            return f"Evidence for other {gene} variants — may not apply to {variant}"
        if gene:
            return f"Evidence for other {gene} variants — may not apply to this variant"
        return "Evidence for other variants in this gene — may not apply to this variant"

    # Fallback for unknown levels
    return f"Evidence at {level} level"


# Static descriptions for use in documentation/UI without dynamic substitution
LOCUS_MATCH_DESCRIPTIONS: dict[str, str] = {
    "variant": "Direct evidence for this exact variant",
    "codon": "Evidence for other variants at this codon — extrapolation from similar position",
    "gene": "Evidence for other variants in this gene — may not apply to this variant",
}

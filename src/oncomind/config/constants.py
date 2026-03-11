"""Centralized constants and mappings for OncoMind.

This module consolidates all hardcoded mappings used across the codebase:
- Gene aliases (different nomenclatures for the same gene)
- Tumor type mappings (abbreviations to full names)
- Amino acid codes
- Variant type classifications

Centralizing these makes maintenance easier and ensures consistency.
"""

from pathlib import Path

# =============================================================================
# DATA PATHS
# =============================================================================
# Centralized paths for prefetched/cached data files.
# Each data source has its own subdirectory under data/.

# Project root (onco_mind_v0/)
_PROJECT_ROOT = Path(__file__).parent.parent.parent.parent

# Data directories - each source has its own subdirectory
CGI_DATA_DIR = _PROJECT_ROOT / "data" / "cgi"
FDA_DATA_DIR = _PROJECT_ROOT / "data" / "fda"
PROCESSING_DATA_DIR = _PROJECT_ROOT / "data" / "processing"

# Individual data files
CGI_BIOMARKERS_FILE = CGI_DATA_DIR / "cgi_biomarkers.tsv"
FDA_LABELS_FILE = PROCESSING_DATA_DIR / "labels" / "fda_labels.json"
DRUG_NAMES_CACHE_FILE = PROCESSING_DATA_DIR / "drug_names.json"

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
# GENE FULL NAMES
# =============================================================================
# Maps gene symbols to their full protein/receptor names as used in FDA labels
# FDA labels often use full names like "epidermal growth factor receptor" instead of "EGFR"

GENE_FULL_NAMES: dict[str, list[str]] = {
    "EGFR": ["epidermal growth factor receptor", "egfr"],
    "ERBB2": ["human epidermal growth factor receptor 2", "her2", "her-2", "erbb2"],
    "HER2": ["human epidermal growth factor receptor 2", "her2", "her-2", "erbb2"],
    "ALK": ["anaplastic lymphoma kinase", "alk"],
    "ROS1": ["ros proto-oncogene 1", "ros1", "c-ros"],
    "BRAF": ["b-raf", "braf", "b-raf proto-oncogene"],
    "KRAS": ["kirsten rat sarcoma", "kras", "k-ras"],
    "NRAS": ["neuroblastoma ras", "nras", "n-ras"],
    "MET": ["met proto-oncogene", "hepatocyte growth factor receptor", "hgfr", "c-met"],
    "RET": ["ret proto-oncogene", "rearranged during transfection"],
    "NTRK1": ["neurotrophic tyrosine receptor kinase 1", "trka", "ntrk1"],
    "NTRK2": ["neurotrophic tyrosine receptor kinase 2", "trkb", "ntrk2"],
    "NTRK3": ["neurotrophic tyrosine receptor kinase 3", "trkc", "ntrk3"],
    "FGFR1": ["fibroblast growth factor receptor 1", "fgfr1"],
    "FGFR2": ["fibroblast growth factor receptor 2", "fgfr2"],
    "FGFR3": ["fibroblast growth factor receptor 3", "fgfr3"],
    "KIT": ["kit proto-oncogene", "c-kit", "cd117", "stem cell factor receptor"],
    "PDGFRA": [
        "platelet-derived growth factor receptor alpha",
        "pdgfra",
        "cd140a",
        "PDGFR",
    ],
    "PDGFRB": ["platelet-derived growth factor receptor beta", "pdgfrb", "cd140b"],
    "PIK3CA": [
        "phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha",
        "pi3k",
        "pik3ca",
    ],
    "AKT1": ["akt serine/threonine kinase 1", "protein kinase b", "pkb", "akt1"],
    "BRCA1": ["brca1", "breast cancer 1", "brca1 dna repair associated"],
    "BRCA2": ["brca2", "breast cancer 2", "brca2 dna repair associated"],
    "IDH1": ["isocitrate dehydrogenase 1", "idh1"],
    "IDH2": ["isocitrate dehydrogenase 2", "idh2"],
    "FLT3": ["fms-related tyrosine kinase 3", "flt3", "fms-like tyrosine kinase 3"],
    "JAK2": ["janus kinase 2", "jak2"],
    "MPL": ["mpl proto-oncogene", "thrombopoietin receptor", "tpo receptor"],
    "CALR": ["calreticulin", "calr"],
    "BCR-ABL1": ["bcr-abl1", "bcr-abl", "philadelphia chromosome"],
    "ABL1": ["abl proto-oncogene 1", "abl1", "c-abl"],
    "TP53": ["tumor protein p53", "p53", "tp53"],
    "PTEN": ["phosphatase and tensin homolog", "pten"],
    "ATM": ["ataxia telangiectasia mutated", "atm"],
    "PALB2": ["partner and localizer of brca2", "palb2"],
    "RAD51C": ["rad51 paralog c", "rad51c"],
    "RAD51D": ["rad51 paralog d", "rad51d"],
    "CHEK2": ["checkpoint kinase 2", "chek2", "chk2"],
    "MSH2": ["muts homolog 2", "msh2"],
    "MSH6": ["muts homolog 6", "msh6"],
    "MLH1": ["mutl homolog 1", "mlh1"],
    "PMS2": ["pms1 homolog 2", "pms2"],
}


# =============================================================================
# TUMOR TYPE MAPPINGS
# =============================================================================
# Maps tumor type abbreviations to their full names and synonyms
# Used for matching user input to database entries (CGI, FDA, CIViC)

TUMOR_TYPE_MAPPINGS: dict[str, list[str]] = {
    # Lung
    "nsclc": [
        "non-small cell lung",
        "nsclc",
        "lung non-small cell",
        "non small cell lung cancer",
        "non-small cell lung carcinoma",
        "lung adenocarcinoma",
        "lung squamous",
        "lung adenosquamous",
    ],
    "lung": ["lung", "lung cancer", "lung carcinoma"],
    "l": [
        "lung",
        "lung cancer",
        "lung carcinoma",
        "nsclc",
        "non-small cell lung",
        "sclc",
        "small cell lung",
    ],  # CGI uses single-letter "L" for all lung cancers
    "sclc": [
        "small cell lung",
        "sclc",
        "small cell lung cancer",
        "small cell lung carcinoma",
    ],
    "luad": ["lung adenocarcinoma", "luad", "adenocarcinoma of lung"],
    "lusc": [
        "lung squamous",
        "lusc",
        "squamous cell lung",
        "squamous cell carcinoma of lung",
    ],
    # Colorectal
    "crc": [
        "colorectal",
        "colon",
        "rectal",
        "crc",
        "mcrc",
        "colorectal cancer",
        "colorectal carcinoma",
    ],
    "coad": ["colon adenocarcinoma", "coad", "colon cancer"],
    "read": ["rectal adenocarcinoma", "read", "rectum adenocarcinoma", "rectal cancer"],
    "coread": [
        "colorectal",
        "colon",
        "rectal",
        "colorectal adenocarcinoma",
        "crc",
        "coadread",
    ],
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
    "mtc": [
        "medullary thyroid carcinoma",
        "mtc",
        "medullary thyroid cancer",
        "medullary thyroid",
        "thyroid cancer",
    ],
    # Gastrointestinal
    "gist": ["gastrointestinal stromal tumor", "gist", "gastrointestinal stromal"],
    "stad": [
        "stomach adenocarcinoma",
        "stad",
        "gastric adenocarcinoma",
        "gastric cancer",
    ],
    "esca": ["esophageal carcinoma", "esca", "esophageal cancer", "oesophageal"],
    "escc": ["esophageal squamous cell carcinoma", "escc"],
    "eadc": ["esophageal adenocarcinoma", "eadc", "barrett"],
    "paad": [
        "pancreatic adenocarcinoma",
        "paad",
        "pancreatic cancer",
        "pancreas adenocarcinoma",
    ],
    "pdac": ["pancreatic ductal adenocarcinoma", "pdac", "pancreatic ductal"],
    "lihc": [
        "liver hepatocellular carcinoma",
        "lihc",
        "hepatocellular carcinoma",
        "hcc",
    ],
    "hcc": ["hepatocellular carcinoma", "hcc"],
    "chol": [
        "cholangiocarcinoma",
        "chol",
        "bile duct cancer",
        "biliary",
        "intrahepatic cholangiocarcinoma",
    ],
    # Genitourinary
    "prad": [
        "prostate adenocarcinoma",
        "prad",
        "prostate cancer",
        "prostate carcinoma",
    ],
    "prostate": ["prostate", "prostate cancer"],
    "blca": [
        "bladder urothelial carcinoma",
        "blca",
        "bladder cancer",
        "urothelial carcinoma",
        "transitional cell",
    ],
    "rcc": ["renal cell carcinoma", "rcc", "kidney cancer"],
    "ccrcc": ["clear cell renal cell carcinoma", "ccrcc", "clear cell kidney"],
    "prcc": ["papillary renal cell carcinoma", "prcc"],
    "chrcc": ["chromophobe renal cell carcinoma", "chrcc"],
    # Gynecologic
    "ov": ["ovarian carcinoma", "ov", "ovarian cancer", "ovary"],
    "hgsc": ["high grade serous ovarian", "hgsc", "serous ovarian"],
    "ucec": [
        "uterine corpus endometrial carcinoma",
        "ucec",
        "endometrial cancer",
        "endometrial carcinoma",
    ],
    "ucs": ["uterine carcinosarcoma", "ucs"],
    "cesc": ["cervical squamous cell carcinoma", "cesc", "cervical cancer"],
    # Head and Neck
    "hnsc": [
        "head and neck squamous cell carcinoma",
        "hnsc",
        "head neck",
        "head and neck cancer",
        "oral cavity",
        "oropharynx",
        "larynx",
        "hypopharynx",
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
    "NSCLC",
    "LUAD",
    "LUSC",
    "SCLC",
    # Breast
    "BRCA",
    "IDC",
    "ILC",
    # Colorectal
    "CRC",
    "COAD",
    "READ",
    # Melanoma
    "MEL",
    "SKCM",
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
    "RCC",
    "CCRCC",
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
    "AML",
    "CML",
]


# =============================================================================
# AMINO ACID CODES
# =============================================================================
# Standard amino acid code conversions (3-letter to 1-letter and vice versa)

AMINO_ACID_3TO1: dict[str, str] = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "TER": "*",
    "STP": "*",
    "STOP": "*",  # Stop codons
    "SEC": "U",  # Selenocysteine (rare)
    "PYL": "O",  # Pyrrolysine (rare)
}

AMINO_ACID_1TO3: dict[str, str] = {
    v: k for k, v in AMINO_ACID_3TO1.items() if v not in ("*",)
}
AMINO_ACID_1TO3["*"] = "TER"  # Prefer TER for stop codon


# =============================================================================
# VARIANT TYPE CLASSIFICATIONS
# =============================================================================
# Variant types allowed for SNP/small indel analysis

ALLOWED_VARIANT_TYPES: set[str] = {
    "missense",
    "nonsense",
    "insertion",
    "deletion",
    "frameshift",
}

# Variant types that are structural (not supported in current system)
STRUCTURAL_VARIANT_TYPES: set[str] = {
    "fusion",
    "amplification",
    "deletion_large",
    "rearrangement",
    "copy_number",
    "truncating",
    "splice",
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

# Each entry is either a plain variant string (unrestricted — applies to any tumor type)
# or a dict {"variant": str, "tumor_types": list[str]} (restricted — only applies when
# the tumor type matches one of the listed substrings, case-insensitive).
ACQUIRED_RESISTANCE_MUTATIONS: dict[str, list] = {
    # EGFR - TKI resistance mutations
    "EGFR": [
        "T790M",  # 1st/2nd-gen TKI resistance (gefitinib, erlotinib, afatinib)
        "C797S",  # 3rd-gen TKI resistance (osimertinib)
        "G724S",  # Osimertinib resistance
        "L718Q",  # Osimertinib resistance
        "L718V",  # Osimertinib resistance
        "L792F",  # Osimertinib resistance
        "G796D",  # Osimertinib resistance
    ],
    # ALK - ALK inhibitor resistance mutations
    "ALK": [
        "G1202R",  # Pan-ALK inhibitor resistance
        "L1196M",  # Crizotinib resistance (gatekeeper)
        "I1171T",  # Alectinib resistance
        "I1171N",  # Alectinib resistance
        # F1174L/C are acquired resistance in lung cancers but primary oncogenic
        # drivers in neuroblastoma — restrict to lung tumor types.
        {"variant": "F1174L", "tumor_types": ["lung", "nsclc"]},
        {"variant": "F1174C", "tumor_types": ["lung", "nsclc"]},
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
        "D538G",  # Most common AI resistance mutation
        "Y537S",  # AI resistance
        "Y537N",  # AI resistance
        "Y537C",  # AI resistance
        "E380Q",  # AI resistance
        "L536Q",  # AI resistance
    ],
    # KIT - Imatinib resistance in GIST (secondary mutations in exon 17)
    "KIT": [
        "D816V",  # Primary imatinib resistance (also germline in mastocytosis)
        "D816H",  # Imatinib resistance
        "D816Y",  # Imatinib resistance
        "D820Y",  # Sunitinib resistance
        "N822K",  # Imatinib resistance
        "Y823D",  # Secondary imatinib resistance
    ],
    # FGFR2 - Acquired resistance to FGFR inhibitors
    "FGFR2": [
        "N549H",  # Pemigatinib/infigratinib resistance
        "N549K",  # Pemigatinib/infigratinib resistance
        "V564F",  # Gatekeeper - FGFR inhibitor resistance
        "V564I",  # Gatekeeper
        "E565A",  # FGFR inhibitor resistance
        "L617V",  # FGFR inhibitor resistance
        "K659M",  # FGFR inhibitor resistance
    ],
    # BTK - Ibrutinib resistance in CLL/lymphoma
    "BTK": [
        "C481S",  # Ibrutinib resistance (covalent binding site)
        "C481R",  # Ibrutinib resistance
        "C481F",  # Ibrutinib resistance
    ],
    # BCR-ABL1 - TKI resistance in CML
    "ABL1": [
        "T315I",  # Pan-TKI resistance (gatekeeper) - except ponatinib
        "F317L",  # Dasatinib resistance
        "V299L",  # Dasatinib resistance
        "Y253H",  # Imatinib resistance
        "E255K",  # Imatinib resistance
        "E255V",  # Imatinib resistance
    ],
}

# Unrestricted lookup: variants that apply regardless of tumor type.
_ACQUIRED_RESISTANCE_UNRESTRICTED: set[tuple[str, str]] = {
    (gene.upper(), entry.upper())
    for gene, entries in ACQUIRED_RESISTANCE_MUTATIONS.items()
    for entry in entries
    if isinstance(entry, str)
}

# Restricted lookup: variants that only apply to specific tumor types.
# Maps (GENE, VARIANT) -> list of allowed tumor-type substrings (lowercase).
_ACQUIRED_RESISTANCE_RESTRICTED: dict[tuple[str, str], list[str]] = {
    (gene.upper(), entry["variant"].upper()): [t.lower() for t in entry["tumor_types"]]
    for gene, entries in ACQUIRED_RESISTANCE_MUTATIONS.items()
    for entry in entries
    if isinstance(entry, dict)
}


def is_acquired_resistance_mutation(
    gene: str, variant: str, tumor_type: str | None = None
) -> bool:
    """Check if a gene/variant is a known acquired resistance mutation.

    Args:
        gene: Gene symbol (e.g., "EGFR")
        variant: Variant notation (e.g., "T790M")
        tumor_type: Optional tumor type string. Required for context-dependent
            variants (e.g., ALK F1174L). When absent and the variant is
            restricted, returns False (safe default).

    Returns:
        True if this is a known acquired resistance mutation in the given context.
    """
    key = (gene.upper(), variant.upper())

    if key in _ACQUIRED_RESISTANCE_UNRESTRICTED:
        return True

    if key in _ACQUIRED_RESISTANCE_RESTRICTED:
        if tumor_type is None:
            return False
        tumor_lower = tumor_type.lower()
        return any(
            allowed in tumor_lower for allowed in _ACQUIRED_RESISTANCE_RESTRICTED[key]
        )

    return False


# =============================================================================
# TARGETABLE SENSITIZING VARIANTS
# =============================================================================
# Variants where FDA-approved drugs exist and acquired resistance is expected.
# When VICC reports resistance for these drug+variant pairs, it represents
# acquired resistance (develops under treatment) NOT intrinsic resistance
# (drug never works). This distinction is critical for discordant evidence detection.
#
# Format: {drug_name_lower: [(gene, variant), ...]}
# Drug names should be lowercase for case-insensitive matching.
#
# Used by: gap_detector._detect_discordant_evidence_internal()
# to avoid flagging FDA approval + VICC resistance as a conflict.

TARGETABLE_SENSITIZING_VARIANTS: dict[str, list[tuple[str, str]]] = {
    # BRAF V600 inhibitors - melanoma, NSCLC, CRC
    "vemurafenib": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "zelboraf": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "dabrafenib": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "tafinlar": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "encorafenib": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "braftovi": [("BRAF", "V600E"), ("BRAF", "V600K")],
    # MEK inhibitors (used with BRAF inhibitors)
    "trametinib": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "mekinist": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "cobimetinib": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "cotellic": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "binimetinib": [("BRAF", "V600E"), ("BRAF", "V600K")],
    "mektovi": [("BRAF", "V600E"), ("BRAF", "V600K")],
    # EGFR TKIs - NSCLC
    "erlotinib": [("EGFR", "L858R"), ("EGFR", "exon 19 deletion"), ("EGFR", "del19")],
    "tarceva": [("EGFR", "L858R"), ("EGFR", "exon 19 deletion"), ("EGFR", "del19")],
    "gefitinib": [("EGFR", "L858R"), ("EGFR", "exon 19 deletion"), ("EGFR", "del19")],
    "iressa": [("EGFR", "L858R"), ("EGFR", "exon 19 deletion"), ("EGFR", "del19")],
    "afatinib": [("EGFR", "L858R"), ("EGFR", "exon 19 deletion"), ("EGFR", "del19")],
    "gilotrif": [("EGFR", "L858R"), ("EGFR", "exon 19 deletion"), ("EGFR", "del19")],
    "osimertinib": [
        ("EGFR", "L858R"),
        ("EGFR", "exon 19 deletion"),
        ("EGFR", "del19"),
        ("EGFR", "T790M"),
    ],
    "tagrisso": [
        ("EGFR", "L858R"),
        ("EGFR", "exon 19 deletion"),
        ("EGFR", "del19"),
        ("EGFR", "T790M"),
    ],
    "dacomitinib": [("EGFR", "L858R"), ("EGFR", "exon 19 deletion"), ("EGFR", "del19")],
    "vizimpro": [("EGFR", "L858R"), ("EGFR", "exon 19 deletion"), ("EGFR", "del19")],
    # EGFR exon 20 insertion inhibitors
    "amivantamab": [("EGFR", "exon 20 insertion")],
    "rybrevant": [("EGFR", "exon 20 insertion")],
    "mobocertinib": [("EGFR", "exon 20 insertion")],
    "exkivity": [("EGFR", "exon 20 insertion")],
    # KRAS G12C inhibitors - NSCLC, CRC
    "sotorasib": [("KRAS", "G12C")],
    "lumakras": [("KRAS", "G12C")],
    "adagrasib": [("KRAS", "G12C")],
    "krazati": [("KRAS", "G12C")],
    # ALK inhibitors - NSCLC
    "crizotinib": [("ALK", "fusion"), ("ALK", "rearrangement"), ("ROS1", "fusion")],
    "xalkori": [("ALK", "fusion"), ("ALK", "rearrangement"), ("ROS1", "fusion")],
    "ceritinib": [("ALK", "fusion"), ("ALK", "rearrangement")],
    "zykadia": [("ALK", "fusion"), ("ALK", "rearrangement")],
    "alectinib": [("ALK", "fusion"), ("ALK", "rearrangement")],
    "alecensa": [("ALK", "fusion"), ("ALK", "rearrangement")],
    "brigatinib": [("ALK", "fusion"), ("ALK", "rearrangement")],
    "alunbrig": [("ALK", "fusion"), ("ALK", "rearrangement")],
    "lorlatinib": [("ALK", "fusion"), ("ALK", "rearrangement")],
    "lorbrena": [("ALK", "fusion"), ("ALK", "rearrangement")],
    # ROS1 inhibitors - NSCLC
    "entrectinib": [
        ("ROS1", "fusion"),
        ("NTRK1", "fusion"),
        ("NTRK2", "fusion"),
        ("NTRK3", "fusion"),
    ],
    "rozlytrek": [
        ("ROS1", "fusion"),
        ("NTRK1", "fusion"),
        ("NTRK2", "fusion"),
        ("NTRK3", "fusion"),
    ],
    "repotrectinib": [("ROS1", "fusion")],
    "augtyro": [("ROS1", "fusion")],
    # RET inhibitors - NSCLC, thyroid
    "selpercatinib": [("RET", "fusion"), ("RET", "mutation")],
    "retevmo": [("RET", "fusion"), ("RET", "mutation")],
    "pralsetinib": [("RET", "fusion"), ("RET", "mutation")],
    "gavreto": [("RET", "fusion"), ("RET", "mutation")],
    # MET inhibitors - NSCLC
    "capmatinib": [("MET", "exon 14 skipping"), ("MET", "amplification")],
    "tabrecta": [("MET", "exon 14 skipping"), ("MET", "amplification")],
    "tepotinib": [("MET", "exon 14 skipping")],
    "tepmetko": [("MET", "exon 14 skipping")],
    # NTRK inhibitors - tumor agnostic
    "larotrectinib": [("NTRK1", "fusion"), ("NTRK2", "fusion"), ("NTRK3", "fusion")],
    "vitrakvi": [("NTRK1", "fusion"), ("NTRK2", "fusion"), ("NTRK3", "fusion")],
    # HER2 inhibitors - breast, gastric
    "trastuzumab": [("ERBB2", "amplification"), ("HER2", "amplification")],
    "herceptin": [("ERBB2", "amplification"), ("HER2", "amplification")],
    "pertuzumab": [("ERBB2", "amplification"), ("HER2", "amplification")],
    "perjeta": [("ERBB2", "amplification"), ("HER2", "amplification")],
    "trastuzumab deruxtecan": [
        ("ERBB2", "amplification"),
        ("HER2", "amplification"),
        ("ERBB2", "mutation"),
    ],
    "enhertu": [
        ("ERBB2", "amplification"),
        ("HER2", "amplification"),
        ("ERBB2", "mutation"),
    ],
    "neratinib": [("ERBB2", "mutation"), ("HER2", "mutation")],
    "nerlynx": [("ERBB2", "mutation"), ("HER2", "mutation")],
    "tucatinib": [("ERBB2", "amplification"), ("HER2", "amplification")],
    "tukysa": [("ERBB2", "amplification"), ("HER2", "amplification")],
    # BCR-ABL inhibitors - CML
    "imatinib": [("ABL1", "fusion"), ("BCR-ABL1", "fusion")],
    "gleevec": [("ABL1", "fusion"), ("BCR-ABL1", "fusion")],
    "dasatinib": [("ABL1", "fusion"), ("BCR-ABL1", "fusion")],
    "sprycel": [("ABL1", "fusion"), ("BCR-ABL1", "fusion")],
    "nilotinib": [("ABL1", "fusion"), ("BCR-ABL1", "fusion")],
    "tasigna": [("ABL1", "fusion"), ("BCR-ABL1", "fusion")],
    "bosutinib": [("ABL1", "fusion"), ("BCR-ABL1", "fusion")],
    "bosulif": [("ABL1", "fusion"), ("BCR-ABL1", "fusion")],
    "ponatinib": [("ABL1", "fusion"), ("BCR-ABL1", "fusion"), ("ABL1", "T315I")],
    "iclusig": [("ABL1", "fusion"), ("BCR-ABL1", "fusion"), ("ABL1", "T315I")],
    # KIT/PDGFRA inhibitors - GIST
    "avapritinib": [("KIT", "D816V"), ("PDGFRA", "D842V")],
    "ayvakit": [("KIT", "D816V"), ("PDGFRA", "D842V")],
    "ripretinib": [("KIT", "mutation")],
    "qinlock": [("KIT", "mutation")],
    # IDH inhibitors - AML, cholangiocarcinoma
    "ivosidenib": [("IDH1", "R132H"), ("IDH1", "R132C"), ("IDH1", "R132")],
    "tibsovo": [("IDH1", "R132H"), ("IDH1", "R132C"), ("IDH1", "R132")],
    "enasidenib": [
        ("IDH2", "R140Q"),
        ("IDH2", "R172K"),
        ("IDH2", "R140"),
        ("IDH2", "R172"),
    ],
    "idhifa": [
        ("IDH2", "R140Q"),
        ("IDH2", "R172K"),
        ("IDH2", "R140"),
        ("IDH2", "R172"),
    ],
    # FGFR inhibitors - cholangiocarcinoma, bladder
    "pemigatinib": [("FGFR2", "fusion"), ("FGFR2", "rearrangement")],
    "pemazyre": [("FGFR2", "fusion"), ("FGFR2", "rearrangement")],
    "infigratinib": [("FGFR2", "fusion"), ("FGFR2", "rearrangement")],
    "truseltiq": [("FGFR2", "fusion"), ("FGFR2", "rearrangement")],
    "futibatinib": [("FGFR2", "fusion"), ("FGFR2", "rearrangement")],
    "lytgobi": [("FGFR2", "fusion"), ("FGFR2", "rearrangement")],
    "erdafitinib": [("FGFR3", "mutation"), ("FGFR2", "mutation"), ("FGFR3", "fusion")],
    "balversa": [("FGFR3", "mutation"), ("FGFR2", "mutation"), ("FGFR3", "fusion")],
    # PIK3CA inhibitors - breast
    "alpelisib": [
        ("PIK3CA", "H1047R"),
        ("PIK3CA", "E545K"),
        ("PIK3CA", "E542K"),
        ("PIK3CA", "mutation"),
    ],
    "piqray": [
        ("PIK3CA", "H1047R"),
        ("PIK3CA", "E545K"),
        ("PIK3CA", "E542K"),
        ("PIK3CA", "mutation"),
    ],
    # AKT inhibitors - breast
    "capivasertib": [("AKT1", "E17K"), ("AKT1", "mutation"), ("PTEN", "loss")],
    "truqap": [("AKT1", "E17K"), ("AKT1", "mutation"), ("PTEN", "loss")],
}

# Build reverse lookup: {(gene_upper, variant_upper): [drug_names]}
# For efficient checking if a gene+variant has known sensitizing drugs
SENSITIZING_VARIANT_TO_DRUGS: dict[tuple[str, str], list[str]] = {}
for drug, variants in TARGETABLE_SENSITIZING_VARIANTS.items():
    for gene, variant in variants:
        key = (gene.upper(), variant.upper())
        if key not in SENSITIZING_VARIANT_TO_DRUGS:
            SENSITIZING_VARIANT_TO_DRUGS[key] = []
        SENSITIZING_VARIANT_TO_DRUGS[key].append(drug)


def is_sensitizing_variant_for_drug(gene: str, variant: str, drug: str) -> bool:
    """Check if a gene+variant is a known sensitizing target for a drug.

    When True, VICC resistance reports for this drug+variant pair represent
    acquired resistance (expected clinical behavior), not intrinsic resistance.

    Args:
        gene: Gene symbol (e.g., "BRAF")
        variant: Variant notation (e.g., "V600E")
        drug: Drug name (case-insensitive)

    Returns:
        True if this variant is a known sensitizing target for this drug
    """
    drug_lower = drug.lower()
    if drug_lower not in TARGETABLE_SENSITIZING_VARIANTS:
        return False

    gene_upper = gene.upper()
    variant_upper = variant.upper()

    for target_gene, target_variant in TARGETABLE_SENSITIZING_VARIANTS[drug_lower]:
        if (
            target_gene.upper() == gene_upper
            and target_variant.upper() == variant_upper
        ):
            return True
    return False


# =============================================================================
# CBIOPORTAL STUDY MAPPINGS
# =============================================================================
CBIOPORTAL_STUDY_MAPPINGS: dict[str, list[str]] = {
    # Lung / NSCLC – prioritize recent MSK + GENIE BPC for clinical depth
    "nsclc": [
        "nsclc_msk_2023",  # Large recent MSK-IMPACT NSCLC cohort
        "lung_msk_2024",  # Broader recent MSK lung
        "nsclc_public_genie_bpc",  # GENIE BPC NSCLC v2.0-public (rich clinical annotations)
        "genie_public",  # Latest GENIE public release (18.0-public as of mid-2025)
        "luad_tcga_pan_can_atlas_2018",  # Reliable harmonized fallback
        "lusc_tcga_pan_can_atlas_2018",
    ],
    "lung": [
        "lung_msk_2024",
        "nsclc_public_genie_bpc",
        "genie_public",
        "luad_tcga_pan_can_atlas_2018",
        "lusc_tcga_pan_can_atlas_2018",
    ],
    "luad": [
        "lung_msk_2024",  # Covers LUAD well in recent MSK
        "nsclc_public_genie_bpc",
        "luad_tcga_pan_can_atlas_2018",
    ],
    "lusc": [
        "lusc_msk_2023",  # Recent MSK squamous if available
        "lusc_tcga_pan_can_atlas_2018",
    ],
    "sclc": [
        "sclc_msk_2024",  # Recent MSK SCLC if available
        "sclc_ucologne_2015",  # Classic small reference
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
    "TP53",  # Most frequently mutated TSG, co-mutates with many genes
    "PTEN",  # Often co-mutates with PIK3CA
    "CDKN2A",  # Cell cycle regulator, co-mutates with TP53
    "RB1",  # Cell cycle TSG, co-mutates with TP53
    "SMAD4",  # Co-mutates with KRAS in pancreatic/CRC
    "STK11",  # Co-mutates with KRAS in lung adenocarcinoma
    "NF1",  # RAS pathway, co-mutates with TP53
    "APC",  # Colorectal driver
    "ARID1A",  # Chromatin remodeling, co-mutates with TP53
    "KEAP1",  # Oxidative stress pathway (lung)
    "ATM",  # DDR gene, co-mutates across tumor types
    "BRCA1",
    "BRCA2",  # DDR genes
    # Key oncogenes
    "KRAS",  # RAS pathway, often mutually exclusive with BRAF/NRAS
    "NRAS",  # RAS pathway, mutually exclusive with KRAS/BRAF
    "BRAF",  # MAPK pathway
    "PIK3CA",  # PI3K pathway, co-mutates with PTEN
    "EGFR",  # Receptor tyrosine kinase
    "ERBB2",  # HER2
    "MET",  # Often co-mutates as resistance mechanism
    "IDH1",
    "IDH2",  # Glioma/AML drivers
    "CTNNB1",  # Wnt pathway
    # Oxidative stress
    "NFE2L2",  # Co-mutates with KEAP1 in lung
    # Chromatin/epigenetic
    "KMT2D",  # Frequently mutated in solid tumors
    "ARID2",  # Chromatin remodeling
    # Receptor tyrosine kinases
    "FGFR1",
    "FGFR2",
    "FGFR3",
    "KIT",
    "PDGFRA",
    "ALK",
    "ROS1",
    "RET",
]


LITERATURE_SENSITIVITY_TERMS = [
    "sensitivity",
    "sensitive",
    "response",
    "efficacy",
    "effective",
    "benefit",
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
    ("BRAF", "V600"),
    ("IDH1", "R132"),
    ("IDH2", "R140"),
    ("IDH2", "R172"),
    ("EGFR", "exon 19 deletion"),
    ("PIK3CA", "E542"),
    ("PIK3CA", "E545"),
    ("PIK3CA", "H1047"),
    ("KRAS", "G12"),
    ("NRAS", "G12"),
    ("HRAS", "G12"),
    ("KRAS", "G13"),
    ("NRAS", "G13"),
    ("HRAS", "G13"),
    ("KRAS", "Q61"),
    ("NRAS", "Q61"),
    ("HRAS", "Q61"),
    ("PIK3CA", "hotspot"),
    ("MET", "exon 14 skipping"),
}


# =============================================================================
# LLM SERVICE CONSTANTS
# =============================================================================
# Configuration for LLM-based evidence synthesis and paper scoring

# Default model for LLM calls
LLM_DEFAULT_MODEL = "gemini/gemini-2.0-flash"

# Fast/cheap model for high-volume operations (paper scoring, knowledge extraction)
LLM_FAST_MODEL = "gpt-4o-mini"

# Timeout settings (seconds)
LLM_TIMEOUT_DEFAULT = 60
LLM_TIMEOUT_CLAUDE = 120

# Token limits for different LLM tasks
LLM_MAX_TOKENS_SYNTHESIS = 3000
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
MAX_CIVIC_ASSERTIONS = 50  # Per evidence level (A, B, C, D, E)
MAX_CLINICAL_TRIALS = 30
MAX_LITERATURE_RESULTS = 30

# =============================================================================
# UI DISPLAY CONSTANTS
# =============================================================================
# Row limits for evidence display tables
UI_MAX_CIVIC_EVIDENCE_ROWS = 50
UI_MAX_RESEARCH_HYPOTHESES = 4
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
    "cancer",
    "tumor",
    "tumour",
    "carcinoma",
    "sarcoma",
    "lymphoma",
    "leukemia",
    "leukaemia",
    "melanoma",
    "neoplasm",
    "malignant",
    "myeloma",
    "glioma",
    "glioblastoma",
    "blastoma",
    "mesothelioma",
    # Clinical terms
    "metastatic",
    "oncology",
    "chemotherapy",
    "antineoplastic",
    "locally advanced",
    "unresectable",
    "refractory",
    # Biomarker terms
    "mutation-positive",
    "mutation",
    "mutated",
    "biomarker",
    "pd-l1",
    "msi-h",
    "microsatellite",
    "tmb",
    "tumor mutational",
    # Drug class terms
    "kinase inhibitor",
    "tyrosine kinase",
    "checkpoint inhibitor",
    "immunotherapy",
    "targeted therapy",
    "antibody-drug conjugate",
    "monoclonal antibody",
    "bispecific",
    # Specific cancer contexts
    "egfr-positive",
    "egfr mutation",
    "non-small cell lung",
    "nsclc",
    "sclc",
    "gist",
    "gastrointestinal stromal",
    "colorectal",
    "breast cancer",
    "prostate cancer",
    "pancreatic",
    "hepatocellular",
    "cholangiocarcinoma",
    "ovarian",
    "endometrial",
]

# Terms indicating non-oncology context for EGFR (kidney function "eGFR")
FDA_EGFR_EXCLUSION_TERMS = [
    "glomerular filtration",
    "egfr ml/min",
    "egfr <",
    "egfr >",
    "renal function",
    "kidney",
    "diabetes",
    "type 2 diabetes",
    "glycemic",
    "glucose",
    "insulin",
    "sglt2",
    "metformin",
    "dapagliflozin",
    "empagliflozin",
    "canagliflozin",
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
    "pembrolizumab": {
        "target": "PD-1",
        "biomarker_genes": [
            "EGFR",
            "KRAS",
            "BRAF",
            "ALK",
            "ROS1",
            "RET",
            "MET",
            "NTRK1",
            "NTRK2",
            "NTRK3",
        ],
    },
    "keytruda": {
        "target": "PD-1",
        "biomarker_genes": [
            "EGFR",
            "KRAS",
            "BRAF",
            "ALK",
            "ROS1",
            "RET",
            "MET",
            "NTRK1",
            "NTRK2",
            "NTRK3",
        ],
    },
    "keytruda qlex": {
        "target": "PD-1",
        "biomarker_genes": [
            "EGFR",
            "KRAS",
            "BRAF",
            "ALK",
            "ROS1",
            "RET",
            "MET",
            "NTRK1",
            "NTRK2",
            "NTRK3",
        ],
    },
    "pembrolizumab and hyaluronidase": {
        "target": "PD-1",
        "biomarker_genes": [
            "EGFR",
            "KRAS",
            "BRAF",
            "ALK",
            "ROS1",
            "RET",
            "MET",
            "NTRK1",
            "NTRK2",
            "NTRK3",
        ],
    },
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

    # For combination drugs like "Osimertinib + pemetrexed", only check the PRIMARY drug
    # (before the first '+'), not combination partners
    primary_drug = drug_lower.split("+")[0].strip()

    for drug_key, info in BIOMARKER_SELECTION_DRUGS.items():
        # Check if PRIMARY drug name matches (substring match in both directions)
        if drug_key in primary_drug or primary_drug in drug_key:
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
        return (
            "Evidence for other variants in this gene — may not apply to this variant"
        )

    # Fallback for unknown levels
    return f"Evidence at {level} level"


# Static descriptions for use in documentation/UI without dynamic substitution
LOCUS_MATCH_DESCRIPTIONS: dict[str, str] = {
    "variant": "Direct evidence for this exact variant",
    "codon": "Evidence for other variants at this codon — extrapolation from similar position",
    "gene": "Evidence for other variants in this gene — may not apply to this variant",
}

SENSITIVITY_KEYWORDS = [
    "inhibit",
    "potent",
    "antiproliferative",
    "proapoptotic",
    "regress",
    "nanomolar activity",
    "selective",
    "selectively",
    "significant impact",
    "effic",
    "sensitiv",
    "respons",
]

# =============================================================================
# HOTSPOT TUMOR TYPE MAPPING
# =============================================================================
# Maps clinical tumor types to hotspot database organ terms
# Hotspot database uses organ-based terms (brain, lung, skin, bowel, etc.)
# Clinical queries use disease names (Glioma, NSCLC, Melanoma, CRC, etc.)
#
# Format: {clinical_term_lowercase: [list of hotspot organ terms]}

# =============================================================================
# VICC DRUG NAME CORRECTIONS
# =============================================================================
# VICC MetaKB returns truncated/malformed drug names from some sources (especially JAX).
# This mapping corrects known bad drug names to their proper forms.
#
# Example: JAX data returns "AgI" instead of "AGI-5198", "881" instead of "AG-881"
# These truncations appear to be upstream data quality issues in VICC/JAX.
#
# Format: {malformed_name: corrected_name}
# Matching is case-insensitive.

VICC_DRUG_NAME_CORRECTIONS: dict[str, str] = {
    # IDH inhibitors (JAX truncations)
    "AgI": "AGI-5198",
    "AGI": "AGI-5198",
    "881": "AG-881",
    "AG-881": "AG-881",  # Ensure proper form stays proper
    # GSK compounds (JAX truncations) - GSK2118436 is dabrafenib
    "GSK": "GSK2118436",
    # BGB compounds (JAX truncations)
    "BGB": "BGB-3290",
    # VX compounds (JAX truncations) - VX-11e is an ERK inhibitor
    "11E": "VX-11E",
    "11e": "VX-11E",
    # InChIKey to drug name mappings (JAX sometimes uses InChIKeys instead of names)
    "DPMYVVGAYAPQNS-UHFFFAOYSA-N": "CCT241161",  # ERK5 inhibitor
    # Add more corrections as discovered
}

# Build lowercase lookup for case-insensitive matching
_VICC_DRUG_CORRECTIONS_LOWER: dict[str, str] = {
    k.lower(): v for k, v in VICC_DRUG_NAME_CORRECTIONS.items()
}

# Regex pattern for InChIKey format (e.g., DPMYVVGAYAPQNS-UHFFFAOYSA-N)
# Format: 14 uppercase letters - 10 uppercase letters - 1 uppercase letter
import re

_INCHIKEY_PATTERN = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

# Patterns that indicate a compound is NOT a drug (sugars, alcohols, acids, etc.)
# These are case-insensitive patterns
_NON_DRUG_PATTERNS = [
    r"-?lactose$",  # Sugars: lactose, 4-deoxylactose
    r"-?glucose$",  # Sugars
    r"-?fructose$",  # Sugars
    r"-?galactose$",  # Sugars
    r"-?sucrose$",  # Sugars
    r"-?maltose$",  # Sugars
    r"^d-",  # D-glucose, D-mannose etc (stereoisomer prefix for sugars)
    r"^l-",  # L-glucose etc
    r"sulfone$",  # Industrial solvents: methyl phenyl sulfone, etc.
    r"^ethephon$",  # Plant growth regulator (agricultural chemical, not a drug)
]
_NON_DRUG_REGEX = re.compile("|".join(_NON_DRUG_PATTERNS), re.IGNORECASE)


def is_valid_drug_name(drug_name: str) -> bool:
    """Check if a drug name appears to be a valid therapeutic compound.

    Filters out:
    - Sugars (lactose, glucose, etc.)
    - Common metabolites and non-therapeutic compounds

    Args:
        drug_name: Drug name to validate

    Returns:
        True if the name appears to be a valid drug, False otherwise
    """
    if not drug_name:
        return False

    # Filter out non-drug compounds based on patterns
    if _NON_DRUG_REGEX.search(drug_name):
        return False

    return True


def correct_vicc_drug_name(drug_name: str) -> str:
    """Correct known malformed drug names from VICC MetaKB and normalize to uppercase.

    Handles:
    - Known truncated drug names (e.g., "AgI" -> "AGI-5198")
    - InChIKey identifiers (e.g., "DPMYVVGAYAPQNS-UHFFFAOYSA-N" -> "InChIKey: DPMYVVGA...")
    - Normalizes all drug names to uppercase for consistency

    Args:
        drug_name: Drug name as returned by VICC API

    Returns:
        Corrected and uppercase drug name, or original uppercase if no correction needed
    """
    if not drug_name:
        return drug_name

    # Check for known corrections first
    corrected = _VICC_DRUG_CORRECTIONS_LOWER.get(drug_name.lower())
    if corrected:
        return corrected.upper()

    # Format InChIKey identifiers for readability
    if _INCHIKEY_PATTERN.match(drug_name):
        # Show abbreviated InChIKey (first 8 chars of the identifier)
        return f"InChIKey: {drug_name[:8]}..."

    # Normalize to uppercase
    return drug_name.upper()


# =============================================================================
# CIVIC CLINICAL SIGNIFICANCE FORMATTING
# =============================================================================
# CIViC API returns significance values like "SENSITIVITYRESPONSE" (no separator).
# This mapping formats them for human-readable display.
#
# Format: {raw_value_upper: formatted_display}

CIVIC_SIGNIFICANCE_DISPLAY: dict[str, str] = {
    "SENSITIVITYRESPONSE": "Sensitivity/Response",
    "SENSITIVITY/RESPONSE": "Sensitivity/Response",
    "RESISTANCE": "Resistance",
    "REDUCED SENSITIVITY": "Reduced Sensitivity",
    "ADVERSE RESPONSE": "Adverse Response",
    "ONCOGENIC": "Oncogenic",
    "BENIGN": "Benign",
    "LIKELY BENIGN": "Likely Benign",
    "PATHOGENIC": "Pathogenic",
    "LIKELY PATHOGENIC": "Likely Pathogenic",
    "UNCERTAIN SIGNIFICANCE": "Uncertain Significance",
    "N/A": "N/A",
}


def format_civic_significance(raw_significance: str | None) -> str:
    """Format CIViC clinical significance for display.

    Args:
        raw_significance: Raw significance value from CIViC API (e.g., "SENSITIVITYRESPONSE")

    Returns:
        Human-readable formatted string (e.g., "Sensitivity/Response")
    """
    if not raw_significance:
        return "Unknown"
    upper = raw_significance.upper()
    return CIVIC_SIGNIFICANCE_DISPLAY.get(upper, raw_significance.title())


HOTSPOT_TUMOR_MAPPING: dict[str, list[str]] = {
    # Brain/CNS tumors
    "glioma": ["brain"],
    "gbm": ["brain"],
    "glioblastoma": ["brain"],
    "astrocytoma": ["brain"],
    "oligodendroglioma": ["brain"],
    "lgg": ["brain"],
    "brain": ["brain"],
    "cns": ["brain"],
    # Lung cancers
    "nsclc": ["lung"],
    "lung": ["lung"],
    "luad": ["lung"],
    "lusc": ["lung"],
    "sclc": ["lung"],
    "lung adenocarcinoma": ["lung"],
    "non-small cell lung": ["lung"],
    "non small cell lung": ["lung"],
    # Skin/Melanoma
    "melanoma": ["skin"],
    "mel": ["skin"],
    "skcm": ["skin"],
    "cutaneous melanoma": ["skin"],
    "skin": ["skin"],
    # Colorectal/Bowel
    "crc": ["bowel"],
    "colorectal": ["bowel"],
    "colorectal cancer": ["bowel"],
    "colon": ["bowel"],
    "colon cancer": ["bowel"],
    "rectal": ["bowel"],
    "rectal cancer": ["bowel"],
    "coad": ["bowel"],
    "read": ["bowel"],
    "bowel": ["bowel"],
    # Thyroid
    "thyroid": ["thyroid"],
    "thca": ["thyroid"],
    "ptc": ["thyroid"],
    "ftc": ["thyroid"],
    "atc": ["thyroid"],
    # Breast
    "breast": ["breast"],
    "brca": ["breast"],
    "bc": ["breast"],
    "tnbc": ["breast"],
    # Pancreatic
    "pancreatic": ["pancreas"],
    "paad": ["pancreas"],
    "pdac": ["pancreas"],
    "pancreas": ["pancreas"],
    # Stomach/Gastric
    "gastric": ["stomach"],
    "stad": ["stomach"],
    "stomach": ["stomach"],
    # Bladder
    "bladder": ["bladder"],
    "blca": ["bladder"],
    "urothelial": ["bladder"],
    # Kidney/Renal
    "renal": ["kidney"],
    "rcc": ["kidney"],
    "kidney": ["kidney"],
    "ccrcc": ["kidney"],
    # Ovarian
    "ovarian": ["ovary"],
    "ov": ["ovary"],
    "hgsc": ["ovary"],
    "ovary": ["ovary"],
    # Liver
    "liver": ["liver"],
    "hcc": ["liver"],
    "hepatocellular": ["liver"],
    "lihc": ["liver"],
    # Biliary/Cholangiocarcinoma
    "cholangiocarcinoma": ["biliary_tract"],
    "chol": ["biliary_tract"],
    "biliary": ["biliary_tract"],
    "bile duct": ["biliary_tract"],
    # Head and Neck
    "head and neck": ["head_neck"],
    "hnsc": ["head_neck"],
    "oral": ["head_neck"],
    "oropharynx": ["head_neck"],
    # Uterine/Endometrial
    "endometrial": ["uterus"],
    "uterine": ["uterus"],
    "ucec": ["uterus"],
    "uterus": ["uterus"],
    # Cervical
    "cervical": ["cervix"],
    "cesc": ["cervix"],
    "cervix": ["cervix"],
    # Prostate
    "prostate": ["prostate"],
    "prad": ["prostate"],
    # Testicular
    "testicular": ["testis"],
    "tgct": ["testis"],
    "testis": ["testis"],
    # Soft tissue/Sarcoma
    "sarcoma": ["soft_tissue"],
    "sarc": ["soft_tissue"],
    "soft tissue": ["soft_tissue"],
    # Hematologic (hotspot uses "blood")
    "aml": ["blood"],
    "all": ["blood"],
    "cml": ["blood"],
    "cll": ["blood"],
    "leukemia": ["blood"],
    "blood": ["blood"],
    # Bone
    "bone": ["bone"],
    "osteosarcoma": ["bone"],
    # Adrenal
    "adrenal": ["adrenal_gland"],
    "acc": ["adrenal_gland"],
    # Thymus
    "thymoma": ["thymus"],
    "thym": ["thymus"],
    "thymus": ["thymus"],
    # Lymphoma (hotspot uses "lymph")
    "lymphoma": ["lymph"],
    "dlbcl": ["lymph"],
    "fl": ["lymph"],
}

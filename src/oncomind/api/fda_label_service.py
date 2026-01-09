"""
FDA Label Service - Unified interface for FDA drug label lookups.

This module provides functions to:
1. Collect drug names from multiple sources (CGI, CIViC, VICC, FDA oncology biomarkers)
2. Look up FDA label info from fda_labels.json cache
3. Fetch missing drugs from OpenFDA API and update cache

Usage:
    from oncomind.api.fda_label_service import get_fda_labels_for_gene

    # Get FDA label info for all drugs associated with a gene
    labels = await get_fda_labels_for_gene(
        gene="BRAF",
        cgi_drugs=["vemurafenib", "dabrafenib"],
        civic_drugs=["encorafenib"],
        vicc_drugs=["trametinib"],
    )
"""

import json
import logging
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from oncomind.config.constants import (
    PROCESSING_DATA_DIR,
    FDA_LABELS_FILE,
    CGI_BIOMARKERS_FILE,
    DRUG_NAMES_CACHE_FILE,
)

from oncomind.models.evidence.fda import (
    FDALabelEvidence,
    ClinicalStudyEvidence,
    MechanismEvidence,
    AdverseReactionsEvidence,
    extract_combination_partners,
)

logger = logging.getLogger(__name__)

# File paths
FDA_ONCOLOGY_BIOMARKERS_FILE = PROCESSING_DATA_DIR / "fda_oncology_biomarkers.xlsx"

# OpenFDA API
OPENFDA_BASE = "https://api.fda.gov/drug/label.json"


def fetch_latest_labels(drug_name: str) -> list[dict[str, Any]]:
    """Fetch the latest FDA labels from OpenFDA for a drug.

    A drug may have multiple FDA labels (different set_ids), each with
    multiple versions. This returns one label per unique set_id, each
    at its latest (highest) version.

    Args:
        drug_name: Drug name to search for

    Returns:
        List of processed label data dicts (one per unique set_id).
        Empty list if no labels found.
    """
    # Group results by set_id, keeping track of latest version for each
    # Key: set_id, Value: (version_num, raw_result)
    latest_by_set_id: dict[str, tuple[int, dict[str, Any]]] = {}

    for field in ["openfda.brand_name", "openfda.generic_name"]:
        params = {
            "search": f'{field}:"{drug_name}"',
            "limit": 10,
        }
        try:
            resp = requests.get(OPENFDA_BASE, params=params, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                for result in data.get("results", []):
                    set_id = result.get("set_id")
                    version = result.get("version")
                    if set_id and version:
                        try:
                            version_num = int(version)
                            # Check if this set_id is new or has a higher version
                            if set_id not in latest_by_set_id or version_num > latest_by_set_id[set_id][0]:
                                latest_by_set_id[set_id] = (version_num, result)
                        except ValueError:
                            continue
        except requests.RequestException as e:
            logger.debug(f"Error fetching label for {drug_name}: {e}")
            continue

    if not latest_by_set_id:
        return []

    # Process each latest result
    return [_process_raw_label(result, drug_name) for _, result in latest_by_set_id.values()]


# In-memory cache for drug names (loaded from disk on first access)
_drug_names_cache: dict[str, dict[str, str | None]] | None = None


def _load_drug_names_cache() -> dict[str, dict[str, str | None]]:
    """Load drug names cache from disk.

    Returns:
        Dict mapping drug name (lowercase) to {generic, brand} dict
    """
    global _drug_names_cache

    if _drug_names_cache is not None:
        return _drug_names_cache

    if DRUG_NAMES_CACHE_FILE.exists():
        try:
            with open(DRUG_NAMES_CACHE_FILE, "r") as f:
                data = json.load(f)
            _drug_names_cache = data.get("drugs", {})
            logger.debug(f"Loaded {len(_drug_names_cache)} drugs from cache")
        except Exception as e:
            logger.warning(f"Error loading drug names cache: {e}")
            _drug_names_cache = {}
    else:
        _drug_names_cache = {}

    return _drug_names_cache


def _save_drug_names_cache() -> None:
    """Save drug names cache to disk."""
    global _drug_names_cache

    if _drug_names_cache is None:
        return

    from datetime import datetime

    DRUG_NAMES_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)

    data = {
        "last_updated": datetime.now().isoformat(),
        "drugs": _drug_names_cache,
    }

    with open(DRUG_NAMES_CACHE_FILE, "w") as f:
        json.dump(data, f, indent=2)

    logger.debug(f"Saved {len(_drug_names_cache)} drugs to cache")


def _add_to_drug_names_cache(drug_name: str, info: dict[str, str | None]) -> None:
    """Add a drug to the cache and save.

    Args:
        drug_name: Drug name (will be lowercased for lookup key)
        info: Dict with 'generic' and 'brand' keys
    """
    cache = _load_drug_names_cache()

    # Add by input name
    cache[drug_name.lower()] = info

    # Also add by generic and brand names for reverse lookups
    if info.get("generic"):
        cache[info["generic"].lower()] = info
    if info.get("brand"):
        cache[info["brand"].lower()] = info

    _save_drug_names_cache()


def get_drug_info(drug_name: str) -> dict[str, str | None] | None:
    """Look up brand/generic name for a drug.

    Checks local cache first, then queries OpenFDA API.
    Results are cached to avoid repeated API calls.

    Args:
        drug_name: Drug name to look up (generic or brand)

    Returns:
        Dict with 'generic' and 'brand' keys, or None if not found
    """
    # Check cache first
    cache = _load_drug_names_cache()
    cached = cache.get(drug_name.lower())
    if cached:
        logger.debug(f"Cache hit for drug: {drug_name}")
        return cached

    # Query OpenFDA API
    search = f'openfda.generic_name:"{drug_name}" OR openfda.brand_name:"{drug_name}"'
    params = {"search": search, "limit": 1}

    try:
        resp = requests.get(OPENFDA_BASE, params=params, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            if data.get("results"):
                openfda = data["results"][0].get("openfda", {})
                info = {
                    "generic": openfda.get("generic_name", [None])[0],
                    "brand": openfda.get("brand_name", [None])[0],
                }
                # Cache the result
                _add_to_drug_names_cache(drug_name, info)
                logger.debug(f"Fetched and cached drug info: {drug_name} -> {info}")
                return info
    except requests.RequestException as e:
        logger.debug(f"Error fetching drug info for {drug_name}: {e}")

    return None


# =============================================================================
# Data Models
# =============================================================================

@dataclass
class ClinicalStudyData:
    """Structured clinical study data extracted from FDA label."""
    trial_name: str | None = None
    nct_id: str | None = None
    patients_n: int | None = None
    pfs_months_treatment: float | None = None
    pfs_months_control: float | None = None
    hazard_ratio: float | None = None
    hazard_ratio_ci: tuple[float, float] | None = None
    orr_treatment: float | None = None  # Percentage
    orr_control: float | None = None  # Percentage
    biomarker_breakdown: dict[str, float] | None = None  # e.g., {"PIK3CA": 76, "AKT1": 13}


@dataclass
class MechanismOfActionData:
    """Structured mechanism of action data."""
    targets: list[str] = field(default_factory=list)
    mechanism: str | None = None
    preclinical: str | None = None


@dataclass
class AdverseReactionsData:
    """Structured adverse reactions data."""
    common_toxicities: list[tuple[str, float]] = field(default_factory=list)  # (name, percentage)
    serious_rate: float | None = None
    discontinuation_rate: float | None = None


@dataclass
class RecentMajorChangesData:
    """Structured recent major changes data."""
    last_label_update: str | None = None
    update_reason: str | None = None


# FDALabelData removed - use FDALabelEvidence from models/evidence/fda.py directly


# =============================================================================
# Drug Collection Functions
# =============================================================================

def get_drugs_from_fda_biomarkers(gene: str) -> set[str]:
    """Get all drugs associated with a gene from fda_oncology_biomarkers.xlsx."""
    drugs: set[str] = set()

    if not FDA_ONCOLOGY_BIOMARKERS_FILE.exists():
        logger.warning(f"FDA biomarkers file not found: {FDA_ONCOLOGY_BIOMARKERS_FILE}")
        return drugs

    try:
        df = pd.read_excel(FDA_ONCOLOGY_BIOMARKERS_FILE)

        # Find the gene column (might be "Biomarker†", "Biomarker", or "Gene")
        gene_col = None
        for col in ["Biomarker†", "Biomarker", "Gene", "gene", "biomarker"]:
            if col in df.columns:
                gene_col = col
                break

        drug_col = None
        for col in ["Drug", "drug"]:
            if col in df.columns:
                drug_col = col
                break

        if gene_col is None or drug_col is None:
            logger.warning(f"Could not find gene/drug columns in {FDA_ONCOLOGY_BIOMARKERS_FILE}")
            return drugs

        # Filter by gene (case-insensitive partial match)
        gene_upper = gene.upper()
        mask = df[gene_col].str.upper().str.contains(gene_upper, na=False)
        gene_drugs = df.loc[mask, drug_col].dropna().unique()

        for drug in gene_drugs:
            drugs.add(str(drug).strip())

        logger.debug(f"Found {len(drugs)} drugs for {gene} in FDA biomarkers")

    except Exception as e:
        logger.error(f"Error reading FDA biomarkers: {e}")

    return drugs


def collect_all_drugs(
    gene: str,
    cgi_drugs: list[str] | None = None,
    civic_drugs: list[str] | None = None,
    vicc_drugs: list[str] | None = None,
) -> set[str]:
    """Collect all unique drug names from all sources for a gene.

    Args:
        gene: Gene symbol (e.g., "BRAF")
        cgi_drugs: Drug names from CGI biomarkers
        civic_drugs: Drug names from CIViC evidence
        vicc_drugs: Drug names from VICC evidence

    Returns:
        Set of unique drug names
    """
    all_drugs: set[str] = set()

    # Add drugs from each source
    if cgi_drugs:
        all_drugs.update(d.strip() for d in cgi_drugs if d)

    if civic_drugs:
        all_drugs.update(d.strip() for d in civic_drugs if d)

    if vicc_drugs:
        all_drugs.update(d.strip() for d in vicc_drugs if d)

    # Add drugs from FDA oncology biomarkers for this gene
    fda_drugs = get_drugs_from_fda_biomarkers(gene)
    all_drugs.update(fda_drugs)

    logger.info(f"Collected {len(all_drugs)} unique drugs for gene {gene}")
    return all_drugs


# =============================================================================
# Cache Loading Functions
# =============================================================================

def load_fda_labels_json() -> dict[str, Any]:
    """Load fda_labels.json cache.

    Returns the 'drugs' dictionary from the JSON file structure:
    {"last_updated": "...", "drugs": {...}}
    """
    if not FDA_LABELS_FILE.exists():
        logger.warning(f"fda_labels.json not found: {FDA_LABELS_FILE}")
        return {}

    try:
        with open(FDA_LABELS_FILE, "r") as f:
            data = json.load(f)

        # Handle the nested structure: {"last_updated": ..., "drugs": {...}}
        drugs = data.get("drugs", {})
        logger.debug(f"Loaded {len(drugs)} entries from fda_labels.json")
        return drugs
    except Exception as e:
        logger.error(f"Error loading fda_labels.json: {e}")
        return {}


# =============================================================================
# Lookup Functions
# =============================================================================

def lookup_by_set_id_version(set_id: str, version: str, labels: dict[str, Any]) -> dict[str, Any] | None:
    """Look up a drug in the cache by set_id and version.

    Args:
        set_id: FDA label set_id (UUID)
        version: FDA label version
        labels: Dictionary from fda_labels.json

    Returns:
        Label data dict if found, None otherwise
    """
    cache_key = f"{set_id}_{version}"
    return labels.get(cache_key)


def lookup_by_drug_name(drug: str, labels: dict[str, Any]) -> dict[str, Any] | None:
    """Look up a drug in the cache by drug name (for cache-only mode).

    Searches all cache entries for matching drug_name field.

    Args:
        drug: Drug name to search for
        labels: Dictionary from fda_labels.json

    Returns:
        Label data dict if found, None otherwise
    """
    drug_lower = drug.lower().strip()
    for entry in labels.values():
        entry_drug_name = entry.get("drug_name", "").lower().strip()
        if entry_drug_name == drug_lower:
            return entry
    return None


# =============================================================================
# Parsing Functions for Structured Data Extraction
# =============================================================================

def parse_clinical_studies(text: str) -> ClinicalStudyData | None:
    """Parse clinical studies text to extract structured data.

    Extracts trial name, NCT ID, patient counts, PFS, HR, ORR, biomarker breakdown.
    """
    if not text:
        return None

    data = ClinicalStudyData()

    # Extract NCT ID (e.g., NCT04305496)
    nct_match = re.search(r"(NCT\d{8})", text)
    if nct_match:
        data.nct_id = nct_match.group(1)

    # Extract trial name - look for patterns like "in TRIAL-NAME (NCT...)"
    # Common patterns: CAPItello-291, KEYNOTE-123, CheckMate-123
    trial_patterns = [
        r"in\s+([A-Z][A-Za-z]+(?:ello|NCER|Mate)?-?\d+)\s*\(",
        r"evaluated in\s+([A-Z][A-Za-z]+-?\d+)",
        r"study\s+([A-Z][A-Za-z]+-?\d+)",
    ]
    for pattern in trial_patterns:
        trial_match = re.search(pattern, text)
        if trial_match:
            data.trial_name = trial_match.group(1)
            break

    # Extract patient count for biomarker-altered subset
    # Pattern: "289 patients had tumors with eligible PIK3CA/AKT1/PTEN"
    patient_patterns = [
        r"(\d+)\s+patients\s+had\s+tumors?\s+with\s+(?:eligible\s+)?(?:PIK3CA|AKT1|BRAF|EGFR)",
        r"of\s+which\s+(\d+)\s+patients\s+had",
        r"(\d+)\s+patients\s+with\s+(?:PIK3CA|AKT1|BRAF|biomarker).?alter",
    ]
    for pattern in patient_patterns:
        pt_match = re.search(pattern, text, re.IGNORECASE)
        if pt_match:
            data.patients_n = int(pt_match.group(1))
            break

    # Extract PFS median months
    # Pattern: "7.3 (5.5, 9.0)" vs "3.1 (1.9, 3.9)"
    # Or: "Median, months (95%CI)" followed by values
    pfs_pattern = r"(?:Median|PFS)[^0-9]*?(\d+\.?\d*)\s*(?:\([^)]+\))?\s*(?:vs\.?|versus)\s*(\d+\.?\d*)"
    pfs_match = re.search(pfs_pattern, text, re.IGNORECASE)
    if pfs_match:
        data.pfs_months_treatment = float(pfs_match.group(1))
        data.pfs_months_control = float(pfs_match.group(2))

    # Extract Hazard Ratio
    # Pattern: "Hazard Ratio ... 0.50 (0.38, 0.65)" - may have text in between
    # Look for a decimal number followed by CI in parentheses after "Hazard ratio"
    hr_patterns = [
        # Standard pattern: HR followed immediately by value and CI
        r"[Hh]azard\s+[Rr]atio[^0-9]*?(\d+\.?\d*)\s*\((\d+\.?\d*)[,\s]+(\d+\.?\d*)\)",
        # Pattern with stratification text in between: find HR header then look for pattern
        r"(\d+\.\d+)\s*\((\d+\.\d+)[,\s]+(\d+\.\d+)\)\s*p-value",
    ]
    for hr_pattern in hr_patterns:
        hr_match = re.search(hr_pattern, text)
        if hr_match:
            data.hazard_ratio = float(hr_match.group(1))
            data.hazard_ratio_ci = (float(hr_match.group(2)), float(hr_match.group(3)))
            break

    # Extract ORR (Objective Response Rate)
    # Pattern: "ORR 26% vs 8%" or "ORR (95% CI) 26% (19, 34) 8% (4, 14)"
    # Must skip "95% CI" which is confidence interval notation, not ORR value
    orr_patterns = [
        # Standard vs pattern: "ORR 26% vs 8%"
        r"(?:ORR|response\s+rate)[^0-9]*?(\d+)%\s*(?:vs\.?|versus)\s*(\d+)%",
        # FDA label table format: "ORR (95% CI) 26% (19, 34) 8% (4, 14)"
        # Skip the "95% CI" part and capture the actual ORR values
        r"ORR\s*\(95%\s*CI\)\s*(\d+)%\s*\([^)]+\)\s*(\d+)%",
        # Simpler pattern for "ORR ... X% ... Y%" but excluding 95 from CI
        r"ORR[^0-9]*?(?:95%\s*CI[^0-9]*)?(\d+)%\s*\([^)]+\)\s*(\d+)%",
    ]
    for pattern in orr_patterns:
        orr_match = re.search(pattern, text, re.IGNORECASE)
        if orr_match:
            val1 = float(orr_match.group(1))
            val2 = float(orr_match.group(2))
            # Sanity check: ORR should typically be < 90%
            # Skip if we accidentally captured "95" from "95% CI"
            if val1 == 95 or val2 == 95:
                continue
            data.orr_treatment = val1
            data.orr_control = val2
            break

    # Extract biomarker breakdown
    # Look for percentages associated with gene names
    # Patterns can be "PIK3CA: 76%" or "13% had an alteration in AKT1"
    biomarker_breakdown = {}
    genes = ["PIK3CA", "AKT1", "PTEN", "BRAF", "EGFR", "KRAS", "ALK", "HER2", "NTRK", "ROS1", "RET", "MET"]

    # Map written numbers to digits
    written_numbers = {
        "seventy-six": 76, "seventy six": 76,
        "thirteen": 13, "seventeen": 17,
        "seventy": 70, "sixty": 60, "fifty": 50, "forty": 40, "thirty": 30, "twenty": 20,
    }

    for gene in genes:
        # Pattern 1: "gene: XX%" or "gene XX%"
        pattern1 = rf"{gene}[:\s]+(\d+)%"
        # Pattern 2: "XX% had ... in gene" - tight match within same clause
        pattern2 = rf"(\d+)%\s+(?:of\s+patients\s+)?had\s+(?:an\s+)?(?:alteration|mutation)\s+in\s+{gene}\b"
        # Pattern 3: "Written-out percent ... in gene" (e.g., "Seventy-six percent ... in PIK3CA")
        pattern3 = rf"([A-Za-z-]+)\s+percent\s+(?:of\s+patients\s+)?had\s+(?:an\s+)?(?:alteration|mutation)\s+in\s+{gene}\b"

        match = re.search(pattern1, text, re.IGNORECASE)
        if match:
            biomarker_breakdown[gene] = float(match.group(1))
        else:
            match = re.search(pattern2, text, re.IGNORECASE)
            if match:
                biomarker_breakdown[gene] = float(match.group(1))
            else:
                match = re.search(pattern3, text, re.IGNORECASE)
                if match:
                    written = match.group(1).lower()
                    if written in written_numbers:
                        biomarker_breakdown[gene] = float(written_numbers[written])

    if biomarker_breakdown:
        data.biomarker_breakdown = biomarker_breakdown

    # Return None if we didn't find any meaningful data
    if (data.nct_id or data.trial_name or data.pfs_months_treatment or
            data.hazard_ratio or data.orr_treatment or data.patients_n or
            data.biomarker_breakdown):
        return data
    return None


def parse_mechanism_of_action(text: str) -> MechanismOfActionData | None:
    """Parse mechanism of action text to extract structured data.

    Extracts targets, mechanism description, and preclinical findings.
    """
    if not text:
        return None

    data = MechanismOfActionData()

    # Extract targets - look for AKT isoforms, kinases, receptors
    targets = []

    # Pattern 1: "inhibitor of... AKT (AKT1, AKT2 and AKT3)"
    akt_pattern = r"\(AKT1,?\s*AKT2\s*(?:and\s+)?AKT3\)"
    akt_match = re.search(akt_pattern, text)
    if akt_match:
        targets = ["AKT1", "AKT2", "AKT3"]
    else:
        # Pattern 2: Generic kinase/gene targets
        target_patterns = [
            r"inhibitor of[^(]*\(([A-Z0-9]+(?:\s*,\s*[A-Z0-9]+)*(?:\s*and\s+[A-Z0-9]+)?)\)",
            r"inhibits?\s+([A-Z0-9]+(?:\s*,\s*[A-Z0-9]+)*)",
            r"targets?\s*(?:include)?\s*([A-Z0-9]+(?:\s*,\s*[A-Z0-9]+)*)",
        ]
        for pattern in target_patterns:
            target_match = re.search(pattern, text)
            if target_match:
                target_str = target_match.group(1)
                # Clean up and split
                target_str = re.sub(r"\s+and\s+", ", ", target_str)
                for t in re.split(r"[,/]", target_str):
                    t = t.strip()
                    if t and len(t) >= 2 and t.upper() == t:  # Uppercase gene names
                        targets.append(t)
                break

    if targets:
        data.targets = targets

    # Extract mechanism - look for "inhibits..." statements
    mech_patterns = [
        r"(?:and\s+)?inhibits\s+phosphorylation\s+of\s+downstream\s+[^.]+",
        r"inhibits?\s+[^.]*?(?:signaling|pathway|phosphorylation)[^.]*",
    ]
    for pattern in mech_patterns:
        mech_match = re.search(pattern, text, re.IGNORECASE)
        if mech_match:
            data.mechanism = mech_match.group(0).strip()
            break

    # Extract preclinical findings - look for "In vitro" or "In vivo" statements
    preclinical_patterns = [
        r"In vitro[^.]+\.",
        r"In vivo[^.]+\.",
        r"reduced growth of[^.]+\.",
    ]
    preclinical_findings = []
    for pattern in preclinical_patterns:
        for match in re.finditer(pattern, text, re.IGNORECASE):
            preclinical_findings.append(match.group(0).strip())

    if preclinical_findings:
        data.preclinical = " ".join(preclinical_findings)

    if data.targets or data.mechanism or data.preclinical:
        return data
    return None


def parse_adverse_reactions(text: str) -> AdverseReactionsData | None:
    """Parse adverse reactions text to extract structured data.

    Extracts common toxicities with percentages, serious rate, discontinuation rate.
    """
    if not text:
        return None

    data = AdverseReactionsData()

    # First extract serious/discontinuation rates before toxicities to avoid false matches
    serious_rate = None
    discont_rate = None

    serious_match = re.search(r"[Ss]erious\s+adverse\s+reactions?\s+(?:occurred\s+in\s+)?(\d+)%", text)
    if serious_match:
        serious_rate = float(serious_match.group(1))

    discont_match = re.search(r"(?:[Pp]ermanent\s+)?[Dd]iscontinuation[^0-9]*?(\d+)%", text)
    if discont_match:
        discont_rate = float(discont_match.group(1))

    # Look for "were" followed by toxicity names - pattern like "were diarrhea, nausea"
    # Also look for toxicities with percentages in a list format
    toxicities = []

    # Pattern for toxicity list: "were diarrhea 77%, nausea 50%, fatigue 40%"
    # or "diarrhea (77%)"
    known_toxicities = [
        "diarrhea", "nausea", "vomiting", "fatigue", "rash", "stomatitis",
        "hyperglycemia", "decreased appetite", "mucositis", "constipation",
        "abdominal pain", "headache", "arthralgia", "myalgia", "pruritus",
        "alopecia", "peripheral neuropathy", "anemia", "neutropenia",
        "thrombocytopenia", "leukopenia", "lymphopenia",
        "cutaneous adverse reactions", "cutaneous reactions",
    ]

    text_lower = text.lower()
    for tox in known_toxicities:
        # Pattern: toxicity name followed by percentage
        pattern = rf"\b{re.escape(tox)}\b[^0-9]*?(\d+)%"
        match = re.search(pattern, text_lower)
        if match:
            pct = float(match.group(1))
            if pct >= 10:  # Only include common ones
                toxicities.append((tox.title(), pct))

    # Sort by percentage descending and take top 10
    toxicities.sort(key=lambda x: x[1], reverse=True)
    data.common_toxicities = toxicities[:10]

    data.serious_rate = serious_rate
    data.discontinuation_rate = discont_rate

    if data.common_toxicities or data.serious_rate is not None or data.discontinuation_rate is not None:
        return data
    return None


def parse_recent_major_changes(text: str) -> RecentMajorChangesData | None:
    """Parse recent major changes text to extract structured data.

    Extracts last label update date and reason.
    """
    if not text:
        return None

    data = RecentMajorChangesData()

    # Extract date - pattern: "02/2025" or "2/2025"
    date_pattern = r"(\d{1,2}/\d{4})"
    date_match = re.search(date_pattern, text)
    if date_match:
        data.last_label_update = date_match.group(1)

    # Extract reason - typically the section name before the date
    # Pattern: "Warnings and Precautions, Hyperglycemia (5.1) 02/2025"
    reason_patterns = [
        r"([^,]+),\s*([^(]+)\s*\(",
        r"^([^,\d]+)",
    ]
    for pattern in reason_patterns:
        reason_match = re.search(pattern, text)
        if reason_match:
            reason = reason_match.group(1).strip()
            if len(reason_match.groups()) > 1:
                reason += ": " + reason_match.group(2).strip()
            data.update_reason = reason
            break

    # Only return data if we found a date - that's the meaningful indicator
    # update_reason alone (without date) is just raw text, not structured data
    if data.last_label_update:
        return data
    return None


# =============================================================================
# OpenFDA Fetch Functions
# =============================================================================

def _join_list_to_text(value: list | str | None) -> str:
    """Convert a list or string to a single text string."""
    if value is None:
        return ""
    if isinstance(value, list):
        return " ".join(str(v) for v in value)
    return str(value)


def extract_approved_indications(indications_text: str) -> list[str]:
    """Extract list of approved diseases/conditions from FDA indication text.

    Parses the indications_and_usage text to find all approved cancer types
    and other diseases the drug is approved for.

    Args:
        indications_text: Raw indications and usage text from FDA label

    Returns:
        List of disease/condition names (e.g., ["breast cancer", "GIST", "RCC"])
    """
    if not indications_text:
        return []

    text_lower = indications_text.lower()
    found = []

    # Cancer types and other conditions to look for
    # Order matters - more specific patterns first
    indication_patterns = [
        # Specific cancer types
        (r"gastrointestinal stromal tumor", "GIST"),
        (r"\bgist\b", "GIST"),
        (r"renal cell carcinoma", "Renal Cell Carcinoma"),
        (r"\brcc\b", "Renal Cell Carcinoma"),
        (r"non-small cell lung cancer", "NSCLC"),
        (r"\bnsclc\b", "NSCLC"),
        (r"small cell lung cancer", "Small Cell Lung Cancer"),
        (r"breast cancer", "Breast Cancer"),
        (r"colorectal cancer", "Colorectal Cancer"),
        (r"colon cancer", "Colorectal Cancer"),
        (r"pancreatic neuroendocrine tumor", "Pancreatic NET"),
        (r"\bpnet\b", "Pancreatic NET"),
        (r"pancreatic cancer", "Pancreatic Cancer"),
        (r"hepatocellular carcinoma", "Hepatocellular Carcinoma"),
        (r"thyroid cancer", "Thyroid Cancer"),
        (r"anaplastic thyroid", "Anaplastic Thyroid Cancer"),
        (r"\bmelanoma\b", "Melanoma"),
        (r"ovarian cancer", "Ovarian Cancer"),
        (r"prostate cancer", "Prostate Cancer"),
        (r"bladder cancer", "Bladder Cancer"),
        (r"urothelial carcinoma", "Urothelial Carcinoma"),
        (r"cholangiocarcinoma", "Cholangiocarcinoma"),
        (r"bile duct cancer", "Cholangiocarcinoma"),
        (r"gastric cancer", "Gastric Cancer"),
        (r"esophageal cancer", "Esophageal Cancer"),
        (r"head and neck", "Head and Neck Cancer"),
        (r"acute myeloid leukemia", "AML"),
        (r"\baml\b", "AML"),
        (r"chronic myeloid leukemia", "CML"),
        (r"\bcml\b", "CML"),
        (r"multiple myeloma", "Multiple Myeloma"),
        (r"\blymphoma\b", "Lymphoma"),
        (r"myelofibrosis", "Myelofibrosis"),
        (r"polycythemia vera", "Polycythemia Vera"),
        (r"systemic mastocytosis", "Systemic Mastocytosis"),
        (r"endometrial cancer", "Endometrial Cancer"),
        (r"cervical cancer", "Cervical Cancer"),
        (r"glioblastoma", "Glioblastoma"),
        # Pan-cancer/biomarker-driven
        (r"msi-h", "MSI-H Tumors"),
        (r"microsatellite instability-high", "MSI-H Tumors"),
        (r"dmmr", "dMMR Tumors"),
        (r"mismatch repair deficient", "dMMR Tumors"),
        (r"ntrk.+fusion", "NTRK Fusion+ Tumors"),
        (r"tumor.agnostic", "Tumor-Agnostic"),
        (r"solid tumor", "Solid Tumors"),
    ]

    seen = set()
    for pattern, display_name in indication_patterns:
        if re.search(pattern, text_lower):
            if display_name not in seen:
                found.append(display_name)
                seen.add(display_name)

    return found


def _fetch_initial_approval_date(drug_name: str, brand_names: list[str] | None = None) -> str | None:
    """Fetch initial FDA approval date from drugsfda endpoint.

    The drugsfda endpoint contains the full approval history with submissions.
    The initial approval is the ORIG submission with status AP.

    Args:
        drug_name: Drug name to search for
        brand_names: Optional list of brand names to try

    Returns:
        Initial approval date in YYYY-MM-DD format, or None if not found
    """
    drugsfda_url = "https://api.fda.gov/drug/drugsfda.json"

    # Try brand names first (more reliable), then drug name
    search_terms = []
    if brand_names:
        search_terms.extend(brand_names)
    search_terms.append(drug_name)

    for term in search_terms:
        params = {
            "search": f'openfda.brand_name:"{term}"',
            "limit": 1
        }
        try:
            resp = requests.get(drugsfda_url, params=params, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                if data.get("results"):
                    result = data["results"][0]
                    submissions = result.get("submissions", [])

                    # Find the original approval (ORIG with status AP)
                    for sub in submissions:
                        if sub.get("submission_type") == "ORIG" and sub.get("submission_status") == "AP":
                            status_date = sub.get("submission_status_date")
                            if status_date and len(status_date) == 8:
                                # Convert YYYYMMDD to YYYY-MM-DD
                                return f"{status_date[:4]}-{status_date[4:6]}-{status_date[6:8]}"

        except requests.RequestException as e:
            logger.debug(f"Error fetching approval date for {term}: {e}")
            continue

    return None


def _process_raw_label(result: dict[str, Any], drug_name: str) -> dict[str, Any]:
    """Process raw OpenFDA label result into structured label data.

    Args:
        result: Raw result dict from OpenFDA API
        drug_name: Drug name (for fetching approval date)

    Returns:
        Processed label data dict with all extracted fields
    """
    # Extract indications_and_usage as a string
    indications = result.get("indications_and_usage", [])
    indications_str = _join_list_to_text(indications)

    # Extract brand/generic names
    openfda = result.get("openfda", {})
    brand_names = openfda.get("brand_name", [])
    generic_names = openfda.get("generic_name", [])
    manufacturer_names = openfda.get("manufacturer_name", [])

    # Get initial approval date from drugsfda endpoint
    initial_approval_date = _fetch_initial_approval_date(drug_name, brand_names)

    # Extract effective_time (label revision date) - format is YYYYMMDD
    effective_time_raw = result.get("effective_time")
    effective_time = None
    if effective_time_raw and len(effective_time_raw) == 8:
        effective_time = f"{effective_time_raw[:4]}-{effective_time_raw[4:6]}-{effective_time_raw[6:8]}"

    # Extract set_id and version for cache key
    set_id = result.get("set_id")
    version = result.get("version")

    # Extract additional sections as text
    clinical_studies_text = _join_list_to_text(result.get("clinical_studies", []))
    mechanism_of_action_text = _join_list_to_text(result.get("mechanism_of_action", []))
    adverse_reactions_text = _join_list_to_text(result.get("adverse_reactions", []))
    recent_major_changes_text = _join_list_to_text(result.get("recent_major_changes", []))

    # Parse structured data from text
    clinical_studies_parsed = parse_clinical_studies(clinical_studies_text)
    mechanism_parsed = parse_mechanism_of_action(mechanism_of_action_text)
    adverse_parsed = parse_adverse_reactions(adverse_reactions_text)
    recent_changes_parsed = parse_recent_major_changes(recent_major_changes_text)

    # Build result dict
    label_data: dict[str, Any] = {
        "indications_and_usage": indications,
        "indications_and_usage_text": indications_str,
        "brand_name": brand_names,
        "generic_name": generic_names,
        "manufacturer_name": manufacturer_names,
        "initial_approval_date": initial_approval_date,
        "effective_time": effective_time,
        "set_id": set_id,
        "version": version,
        "clinical_studies_text": clinical_studies_text,
        "mechanism_of_action_text": mechanism_of_action_text,
        "adverse_reactions_text": adverse_reactions_text,
        "recent_major_changes_text": recent_major_changes_text,
    }

    # Add parsed structured data if available
    if clinical_studies_parsed:
        label_data["clinical_studies"] = {
            "trial_name": clinical_studies_parsed.trial_name,
            "nct_id": clinical_studies_parsed.nct_id,
            "patients_n": clinical_studies_parsed.patients_n,
            "pfs_months_treatment": clinical_studies_parsed.pfs_months_treatment,
            "pfs_months_control": clinical_studies_parsed.pfs_months_control,
            "hazard_ratio": clinical_studies_parsed.hazard_ratio,
            "hazard_ratio_ci": clinical_studies_parsed.hazard_ratio_ci,
            "orr_treatment": clinical_studies_parsed.orr_treatment,
            "orr_control": clinical_studies_parsed.orr_control,
            "biomarker_breakdown": clinical_studies_parsed.biomarker_breakdown,
        }

    if mechanism_parsed:
        label_data["mechanism_of_action"] = {
            "targets": mechanism_parsed.targets,
            "mechanism": mechanism_parsed.mechanism,
            "preclinical": mechanism_parsed.preclinical,
        }

    if adverse_parsed:
        label_data["adverse_reactions"] = {
            "common_toxicities": adverse_parsed.common_toxicities,
            "serious_rate": adverse_parsed.serious_rate,
            "discontinuation_rate": adverse_parsed.discontinuation_rate,
        }

    if recent_changes_parsed:
        label_data["recent_major_changes"] = {
            "last_label_update": recent_changes_parsed.last_label_update,
            "update_reason": recent_changes_parsed.update_reason,
        }

    return label_data


# =============================================================================
# Cache Update Functions
# =============================================================================

def update_fda_labels_json(drug: str, label_data: dict[str, Any], genes: list[str]) -> bool:
    """Update fda_labels.json with a new drug entry.

    Maintains the file structure: {"last_updated": "...", "drugs": {...}}

    Uses composite key format: "set_id_version" (e.g., "d698c106-2322-401e-b738-cbd83c843ecf_8")
    where set_id is the FDA label's unique identifier and version is the label version.
    This allows tracking multiple versions of the same label and different labels
    for different approvals.

    Args:
        drug: Drug name
        label_data: Label data from OpenFDA
        genes: List of genes associated with this drug

    Returns:
        True if the entry was saved, False if set_id or version was missing
    """
    from datetime import datetime

    # Require set_id and version for cache key - do not fall back to drug name
    set_id = label_data.get("set_id")
    version = label_data.get("version")
    if not set_id or not version:
        logger.warning(f"Cannot cache {drug}: missing set_id or version")
        return False

    # Load existing file to get full structure
    if FDA_LABELS_FILE.exists():
        with open(FDA_LABELS_FILE, "r") as f:
            full_data = json.load(f)
        drugs = full_data.get("drugs", {})
    else:
        full_data = {}
        drugs = {}
    entry = {
        "drug_name": drug,  # Store original drug name in the entry
        "genes": genes,
        "indications_and_usage": label_data.get("indications_and_usage", []),
        "brand_name": label_data.get("brand_name", []),
        "generic_name": label_data.get("generic_name", []),
        "manufacturer_name": label_data.get("manufacturer_name", []),
        "initial_approval_date": label_data.get("initial_approval_date"),
        "effective_time": label_data.get("effective_time"),  # Label revision date
        "set_id": set_id,
        "version": version,
    }

    # Add extended fields if available
    if label_data.get("clinical_studies_text"):
        entry["clinical_studies_text"] = label_data["clinical_studies_text"]
    if label_data.get("clinical_studies"):
        entry["clinical_studies"] = label_data["clinical_studies"]

    if label_data.get("mechanism_of_action_text"):
        entry["mechanism_of_action_text"] = label_data["mechanism_of_action_text"]
    if label_data.get("mechanism_of_action"):
        entry["mechanism_of_action"] = label_data["mechanism_of_action"]

    if label_data.get("adverse_reactions_text"):
        entry["adverse_reactions_text"] = label_data["adverse_reactions_text"]
    if label_data.get("adverse_reactions"):
        entry["adverse_reactions"] = label_data["adverse_reactions"]

    if label_data.get("recent_major_changes_text"):
        entry["recent_major_changes_text"] = label_data["recent_major_changes_text"]
    if label_data.get("recent_major_changes"):
        entry["recent_major_changes"] = label_data["recent_major_changes"]

    # Create composite key: set_id_version (e.g., "d698c106-2322-401e-b738-cbd83c843ecf_8")
    cache_key = f"{set_id}_{version}"
    drugs[cache_key] = entry

    # Rebuild full structure with updated drugs
    full_data["last_updated"] = datetime.now().isoformat()
    full_data["drugs"] = drugs

    # Save
    FDA_LABELS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(FDA_LABELS_FILE, "w") as f:
        json.dump(full_data, f, indent=2)

    logger.info(f"Updated fda_labels.json with {cache_key}")
    return True


# =============================================================================
# Main Service Functions
# =============================================================================

def _create_fda_label_evidence(
    drug: str,
    gene: str,
    label_data: dict[str, Any],
) -> FDALabelEvidence:
    """Create FDALabelEvidence from label data dict.

    Helper function that converts raw label data to FDALabelEvidence model.

    Args:
        drug: Drug name
        gene: Gene symbol
        label_data: Label data dict from cache or OpenFDA

    Returns:
        FDALabelEvidence instance
    """
    # Extract basic fields
    brand_name = label_data.get("brand_name", [])
    if isinstance(brand_name, list):
        brand_name = brand_name[0] if brand_name else None

    generic_name = label_data.get("generic_name", [])
    if isinstance(generic_name, list):
        generic_name = generic_name[0] if generic_name else None

    manufacturer = label_data.get("manufacturer_name", [])
    if isinstance(manufacturer, list):
        manufacturer = manufacturer[0] if manufacturer else None

    indications = label_data.get("indications_and_usage_text", "")
    if not indications:
        indications_list = label_data.get("indications_and_usage", [])
        if isinstance(indications_list, list):
            indications = " ".join(indications_list)

    # Build structured data from parsed fields using Pydantic models
    clinical_studies_evidence = None
    if label_data.get("clinical_studies"):
        cs = label_data["clinical_studies"]
        clinical_studies_evidence = ClinicalStudyEvidence(
            trial_name=cs.get("trial_name"),
            nct_id=cs.get("nct_id"),
            patients_n=cs.get("patients_n"),
            pfs_months_treatment=cs.get("pfs_months_treatment"),
            pfs_months_control=cs.get("pfs_months_control"),
            hazard_ratio=cs.get("hazard_ratio"),
            hazard_ratio_ci=tuple(cs["hazard_ratio_ci"]) if cs.get("hazard_ratio_ci") else None,
            orr_treatment=cs.get("orr_treatment"),
            orr_control=cs.get("orr_control"),
            biomarker_breakdown=cs.get("biomarker_breakdown"),
        )

    mechanism_evidence = None
    if label_data.get("mechanism_of_action"):
        moa = label_data["mechanism_of_action"]
        mechanism_evidence = MechanismEvidence(
            targets=moa.get("targets", []),
            mechanism=moa.get("mechanism"),
            preclinical=moa.get("preclinical"),
        )

    adverse_evidence = None
    if label_data.get("adverse_reactions"):
        ar = label_data["adverse_reactions"]
        adverse_evidence = AdverseReactionsEvidence(
            common_toxicities=[tuple(t) for t in ar.get("common_toxicities", [])],
            serious_rate=ar.get("serious_rate"),
            discontinuation_rate=ar.get("discontinuation_rate"),
        )

    # Extract last_label_update from recent_major_changes
    last_label_update = None
    update_reason = None
    if label_data.get("recent_major_changes"):
        rmc = label_data["recent_major_changes"]
        last_label_update = rmc.get("last_label_update")
        update_reason = rmc.get("update_reason")

    # Extract approved indications from the full indication text
    approved_indications = extract_approved_indications(indications or "")

    # Build therapy name with combination partners (if any)
    # Extract partner names from indication text directly
    partner_names = extract_combination_partners(indications or "")
    if partner_names:
        therapy_name = f"{drug} + {' + '.join(partner_names)}"
        logger.debug(f"Combination therapy: {therapy_name}")
    else:
        therapy_name = drug

    # Create FDALabelEvidence directly (locus_variant_match will be None initially)
    x= FDALabelEvidence(
        drug=therapy_name,  # e.g., "capivasertib + FASLODEX"
        gene=gene,
        brand_name=brand_name,
        generic_name=generic_name,
        manufacturer=manufacturer,
        indications_and_usage=indications,
        initial_approval_date=label_data.get("initial_approval_date"),
        effective_time=label_data.get("effective_time"),
        set_id=label_data.get("set_id"),
        version=str(label_data.get("version")) if label_data.get("version") else None,
        approved_indications=approved_indications,
        clinical_studies=clinical_studies_evidence,
        mechanism_of_action=mechanism_evidence,
        adverse_reactions=adverse_evidence,
        last_label_update=last_label_update,
        update_reason=update_reason,
        clinical_studies_text=label_data.get("clinical_studies_text"),
        mechanism_of_action_text=label_data.get("mechanism_of_action_text"),
        adverse_reactions_text=label_data.get("adverse_reactions_text"),
    )
    return x


def get_fda_labels_for_drugs(
    drugs: set[str],
    gene: str,
    fetch_missing: bool = True,
) -> list[FDALabelEvidence]:
    """Get FDA label data for a set of drugs.

    When fetch_missing=True:
    1. Call fetch_label_version(drug) to get set_id and version from OpenFDA
    2. Look up cache by set_id_version key
    3. If found in cache, use it; if not, fetch full label

    When fetch_missing=False (cache-only mode):
    1. Look up cache by drug_name field (no API calls)

    Args:
        drugs: Set of drug names to look up
        gene: Gene symbol for context
        fetch_missing: Whether to fetch from OpenFDA (if False, only uses cache)

    Returns:
        List of FDALabelEvidence for found drugs (without locus_variant_match populated)
    """
    results: list[FDALabelEvidence] = []

    # Load cache
    labels_cache = load_fda_labels_json()

    for drug in drugs:
        if fetch_missing:
            # Online mode: fetch latest labels (one per set_id)
            labels_data = fetch_latest_labels(drug)

            if not labels_data:
                logger.debug(f"No FDA labels found for {drug}")
                continue

            for label_data in labels_data:
                set_id = label_data.get("set_id")
                version = label_data.get("version")
                if not set_id or not version:
                    logger.warning(f"FDA label for {drug} missing set_id or version, skipping")
                    continue

                # Check if already in cache
                cached_entry = lookup_by_set_id_version(set_id, version, labels_cache)
                if cached_entry:
                    logger.debug(f"Cache hit for {drug} (set_id={set_id}, version={version})")
                    label_evidence = _create_fda_label_evidence(drug, gene, cached_entry)
                else:
                    # Cache miss - save to cache
                    logger.debug(f"Cache miss for {drug}, saving to cache...")
                    update_fda_labels_json(drug, label_data, [gene])
                    label_evidence = _create_fda_label_evidence(drug, gene, label_data)

                results.append(label_evidence)

            # Rate limiting between API calls
            time.sleep(0.2)
        else:
            # Cache-only mode: search by drug_name field (no API calls)
            cached_entry = lookup_by_drug_name(drug, labels_cache)
            if cached_entry:
                logger.debug(f"Cache hit for {drug}")
                label_evidence = _create_fda_label_evidence(drug, gene, cached_entry)
                results.append(label_evidence)
            else:
                logger.debug(f"Cache miss for {drug} (cache-only mode)")

    return results


async def get_fda_labels_for_gene(
    gene: str,
    cgi_drugs: list[str] | None = None,
    civic_drugs: list[str] | None = None,
    vicc_drugs: list[str] | None = None,
    fetch_missing: bool = True,
) -> list[FDALabelEvidence]:
    """Get FDA label data for all drugs associated with a gene.

    Collects drug names from CGI, CIViC, VICC, and FDA oncology biomarkers,
    then looks up FDA label data for each drug.

    Args:
        gene: Gene symbol (e.g., "BRAF")
        cgi_drugs: Drug names from CGI biomarkers
        civic_drugs: Drug names from CIViC evidence
        vicc_drugs: Drug names from VICC evidence
        fetch_missing: Whether to fetch missing drugs from OpenFDA

    Returns:
        List of FDALabelEvidence for all found drugs
    """
    # Collect all drugs from all sources
    all_drugs = collect_all_drugs(
        gene=gene,
        cgi_drugs=cgi_drugs,
        civic_drugs=civic_drugs,
        vicc_drugs=vicc_drugs,
    )

    if not all_drugs:
        logger.info(f"No drugs found for gene {gene}")
        return []

    # Get FDA labels for all drugs
    return get_fda_labels_for_drugs(all_drugs, gene, fetch_missing)


# Synchronous wrapper
def get_fda_labels_for_gene_sync(
    gene: str,
    cgi_drugs: list[str] | None = None,
    civic_drugs: list[str] | None = None,
    vicc_drugs: list[str] | None = None,
    fetch_missing: bool = True,
) -> list[FDALabelEvidence]:
    """Synchronous version of get_fda_labels_for_gene."""
    import asyncio
    return asyncio.run(get_fda_labels_for_gene(
        gene=gene,
        cgi_drugs=cgi_drugs,
        civic_drugs=civic_drugs,
        vicc_drugs=vicc_drugs,
        fetch_missing=fetch_missing,
    ))

"""General utility functions."""

from typing import Any, Dict, List

from oncomind.config.constants import FDA_KIT_EXCLUSION_PATTERNS


def normalize_drug_name(drug_name: str) -> str:
    """Normalize drug name by removing salt/formulation suffixes.

    Used for deduplication of FDA evidence where the same drug may appear
    with different salt forms (e.g., ERLOTINIB vs ERLOTINIB HYDROCHLORIDE).

    NOTE: This function intentionally does NOT normalize biologic formulations
    (e.g., IV vs SC formulations like AMIVANTAMAB-VMJW vs AMIVANTAMAB AND
    HYALURONIDASE-LPUJ) because these are clinically distinct products that
    should be shown separately.

    Args:
        drug_name: Original drug name (e.g., "ERLOTINIB HYDROCHLORIDE")

    Returns:
        Normalized drug name in lowercase (e.g., "erlotinib")

    Examples:
        >>> normalize_drug_name("ERLOTINIB HYDROCHLORIDE")
        'erlotinib'
        >>> normalize_drug_name("DABRAFENIB MESYLATE")
        'dabrafenib'
        >>> normalize_drug_name("Pembrolizumab")
        'pembrolizumab'
    """
    if not drug_name:
        return ""

    name = drug_name.upper().strip()

    # Common salt/formulation suffixes to remove
    suffixes = [
        " HYDROCHLORIDE",
        " MESYLATE",
        " MALEATE",
        " FUMARATE",
        " SUCCINATE",
        " TARTRATE",
        " CITRATE",
        " SULFATE",
        " PHOSPHATE",
        " SODIUM",
        " POTASSIUM",
        " CALCIUM",
        " DIMESYLATE",
        " DIHYDROCHLORIDE",
        " TOSYLATE",
        " ACETATE",
        " BESYLATE",
    ]

    for suffix in suffixes:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    return name.lower()


def is_kit_false_positive(indication_text: str | None) -> bool:
    """Check if an FDA label is a KIT false positive (diagnostic kit, not KIT gene).

    Args:
        indication_text: The indications_and_usage text from FDA label

    Returns:
        True if this appears to be a diagnostic/preparation kit, not a KIT oncogene drug
    """
    if not indication_text:
        return False

    text_lower = indication_text.lower()

    # Check for KIT exclusion patterns (diagnostic kit, test kit, etc.)
    has_exclusion = any(pattern in text_lower for pattern in FDA_KIT_EXCLUSION_PATTERNS)

    # Check for oncology context
    oncology_terms = [
        "cancer",
        "tumor",
        "malignant",
        "neoplasm",
        "carcinoma",
        "leukemia",
        "lymphoma",
        "melanoma",
        "sarcoma",
        "gist",
        "mastocytosis",
    ]
    has_oncology = any(term in text_lower for term in oncology_terms)

    # If it has exclusion patterns and no oncology context, it's a false positive
    return has_exclusion and not has_oncology


def dedupe_civic_evidence(civic_evidence_list) -> List[Dict[str, Any]]:
    """Deduplicate CIViC evidence items by evidence_id.

    Defensive measure to ensure no duplicate EIDs appear in the UI,
    even if the API returns them.

    Args:
        civic_evidence_list: List of CIViC evidence objects with evidence_id attribute

    Returns:
        List of deduplicated evidence dicts ready for display
    """
    seen_ids = set()
    deduped = []
    for e in civic_evidence_list:
        if e.evidence_id in seen_ids:
            continue
        seen_ids.add(e.evidence_id)
        deduped.append(
            {
                "evidence_id": e.evidence_id,
                "eid": e.eid,  # Formatted ID (e.g., "EID5586")
                "civic_url": e.civic_url,  # Direct link to CIViC
                "evidence_type": e.evidence_type,
                "evidence_level": e.evidence_level,
                "clinical_significance": e.clinical_significance,
                "disease": e.disease,
                "drugs": e.drugs,
                "description": e.description,
                "pmid": e.pmid,
                "source_url": e.source_url,
                "trust_rating": e.trust_rating
                or e.rating,  # Use trust_rating if available, else rating
                "evidence_direction": e.evidence_direction,
                # Match specificity tracking
                "locus_match": e.locus_match,
                "matched_profile": e.matched_profile,
                "tumor_match": e.tumor_match,
            }
        )
    return deduped


def dedupe_vicc_evidence(vicc_evidence_list) -> List[Dict[str, Any]]:
    """Deduplicate VICC evidence items by (drugs, disease, response_type).

    VICC MetaKB entries don't have unique IDs like CIViC evidence_id,
    so we deduplicate by the combination of drugs, disease, and response type.
    This prevents the same drug-disease-response combination from appearing
    multiple times in the UI.

    Args:
        vicc_evidence_list: List of VICC evidence objects

    Returns:
        List of deduplicated evidence dicts ready for display
    """
    seen_keys = set()
    deduped = []
    for v in vicc_evidence_list:
        # Create a deduplication key from drugs (as sorted tuple), disease, and response
        drugs_key = tuple(sorted(v.drugs)) if v.drugs else ()
        disease_key = (v.disease or "").lower().strip()
        response_key = (v.response_type or "").lower().strip()
        dedup_key = (drugs_key, disease_key, response_key)

        if dedup_key in seen_keys:
            continue
        seen_keys.add(dedup_key)

        deduped.append(
            {
                "source": v.source,
                "drugs": v.drugs,
                "disease": v.disease,
                "response_type": v.response_type,
                "evidence_level": v.evidence_level,
                "molecular_profile": v.molecular_profile,
                "molecular_profile_score": v.molecular_profile_score,
                "publication_url": (
                    v.publication_url[0]
                    if isinstance(v.publication_url, list) and v.publication_url
                    else v.publication_url
                ),
                # Match specificity tracking
                "locus_match": v.locus_match,
                "matched_profile": v.matched_profile,
                "tumor_match": v.tumor_match,
            }
        )
    return deduped

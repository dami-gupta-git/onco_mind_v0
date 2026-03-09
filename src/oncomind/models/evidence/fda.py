"""FDA-related utility functions for evidence processing."""

import re


def extract_combination_partners(indication_text: str | None) -> list[str]:
    """Extract combination drug partners from indication text."""
    if not indication_text:
        return []

    patterns = [
        r"in combination with\s+([^,\.]+?)(?:\s+for|\s+in|\.|,|$)",
        r"combined with\s+([^,\.]+?)(?:\s+for|\s+in|\.|,|$)",
        r"coadministered with\s+([^,\.]+?)(?:\s+for|\s+in|\.|,|$)",
    ]

    partners = []
    for pattern in patterns:
        matches = re.findall(pattern, indication_text, re.IGNORECASE)
        partners.extend([m.strip() for m in matches])

    # Handle "X and Y" or "X or Y" in partner string
    expanded = []
    for p in partners:
        parts = re.split(r"\s+(?:and|or|\+|with)\s+", p, flags=re.IGNORECASE)
        expanded.extend([x.strip() for x in parts])

    # Clean up artifacts
    cleaned = []
    for p in expanded:
        # Remove leading "either", "without", "or without"
        p = re.sub(r"^(either|without|or without)\s+", "", p, flags=re.IGNORECASE)
        # Remove trailing "or without X"
        p = re.sub(r"\s+(or without|without).*$", "", p, flags=re.IGNORECASE)
        p = p.strip()
        if p:
            cleaned.append(p)

    # Filter to only drug-like names
    drug_suffixes = (
        "mab",
        "nib",
        "lib",
        "tib",
        "tin",
        "zumab",
        "ximab",
        "mumab",
        "platin",
        "taxel",
        "mustine",
        "rubicin",
        "mycin",
        "tinib",
        "ciclib",
        "lisib",
        "rafenib",
        "metinib",
        "sertib",
        "parin",
        "mide",
        "uracil",
        "citabine",
        "trexed",
        "fosfamide",
        "imod",
        "olimus",
        "fulvestrant",
        "pemetrexed",
        "asertib",
    )

    filtered = []
    for p in cleaned:
        p_lower = p.lower()
        # Skip generic terms
        skip_terms = (
            "chemotherapy",
            "radiotherapy",
            "radiation",
            "therapy",
            "platinum-based",
            "platinum",
            "chemoradiotherapy",
            "neoadjuvant treatment",
            "fluoropyrimidine-based",
            "platinum-containing",
            "single agent",
        )
        if any(term in p_lower for term in skip_terms):
            continue
        # Check suffix
        if any(p_lower.endswith(suffix) for suffix in drug_suffixes):
            filtered.append(p)

    # Deduplicate while preserving order
    seen = set()
    unique = []
    for p in filtered:
        p_lower = p.lower()
        if p_lower not in seen:
            seen.add(p_lower)
            unique.append(p)

    return unique

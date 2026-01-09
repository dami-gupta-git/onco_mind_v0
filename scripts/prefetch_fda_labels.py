#!/usr/bin/env python3
"""
Prefetch FDA drug labels for all drugs in the oncology biomarker table.

This script:
1. Reads unique drugs from fda_oncology_biomarkers.xlsx
2. Fetches the indications_and_usage section from OpenFDA for each drug
3. Saves results to data/processing/fda_labels.json
"""

import json
import logging
import sys
import time
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pandas as pd
import requests

from oncomind.config.constants import PROCESSING_DATA_DIR, FDA_LABELS_FILE

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

OPENFDA_BASE = "https://api.fda.gov/drug/label.json"
OUTPUT_FILE = FDA_LABELS_FILE


def fetch_drug_label(drug_name: str) -> dict | None:
    """Fetch FDA label for a drug from OpenFDA API."""
    # Try brand name first, then generic name
    for field in ["openfda.brand_name", "openfda.generic_name"]:
        params = {
            "search": f'{field}:"{drug_name}"',
            "limit": 1
        }
        try:
            resp = requests.get(OPENFDA_BASE, params=params, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                if data.get("results"):
                    result = data["results"][0]
                    return {
                        "indications_and_usage": result.get("indications_and_usage", []),
                        "brand_name": result.get("openfda", {}).get("brand_name", []),
                        "generic_name": result.get("openfda", {}).get("generic_name", []),
                        "manufacturer_name": result.get("openfda", {}).get("manufacturer_name", []),
                    }
            elif resp.status_code == 404:
                continue  # Try next field
            else:
                logger.warning(f"HTTP {resp.status_code} for {drug_name} ({field})")
        except requests.RequestException as e:
            logger.error(f"Request failed for {drug_name}: {e}")

    return None


def main():
    # Read unique drugs from biomarker file
    biomarker_file = PROCESSING_DATA_DIR / "fda_oncology_biomarkers.xlsx"
    df = pd.read_excel(biomarker_file)
    unique_drugs = sorted(df["Drug"].dropna().unique())

    logger.info(f"Found {len(unique_drugs)} unique drugs to fetch")

    # Create output directory
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Load existing labels if any (for resuming)
    existing_labels = {}
    if OUTPUT_FILE.exists():
        with open(OUTPUT_FILE) as f:
            existing_labels = json.load(f)
        logger.info(f"Loaded {len(existing_labels)} existing labels")

    # Fetch labels
    labels = existing_labels.copy()
    fetched = 0
    not_found = []

    for i, drug in enumerate(unique_drugs):
        if drug in labels:
            logger.debug(f"Skipping {drug} (already fetched)")
            continue

        logger.info(f"[{i+1}/{len(unique_drugs)}] Fetching: {drug}")
        label = fetch_drug_label(drug)

        if label:
            labels[drug] = label
            fetched += 1
        else:
            not_found.append(drug)
            labels[drug] = None  # Mark as attempted

        # Rate limiting - be nice to the API
        time.sleep(0.2)

        # Save periodically
        if fetched % 10 == 0 and fetched > 0:
            with open(OUTPUT_FILE, "w") as f:
                json.dump(labels, f, indent=2)
            logger.info(f"Saved {len(labels)} labels")

    # Final save
    with open(OUTPUT_FILE, "w") as f:
        json.dump(labels, f, indent=2)

    # Report
    found = sum(1 for v in labels.values() if v is not None)
    logger.info(f"\nDone! Fetched {fetched} new labels")
    logger.info(f"Total labels: {found}/{len(unique_drugs)}")

    if not_found:
        logger.warning(f"\nNot found ({len(not_found)} drugs):")
        for drug in not_found:
            logger.warning(f"  - {drug}")


if __name__ == "__main__":
    main()

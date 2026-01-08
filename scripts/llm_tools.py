#!/usr/bin/env python3
"""LLM debugging and inspection tools.

This script provides utilities for debugging and inspecting the data pipeline
that feeds into the LLM for variant annotation.

Commands:
    baseline  - Test LLM baseline knowledge (no evidence provided)
    data      - Show all data that would be sent to the LLM prompt

Usage:
    # Run from the project root directory

    # Test what the LLM knows from training data alone (no database evidence)
    python scripts/llm_tools.py baseline PIK3CA H1047R "Thyroid Cancer"

    # Show all evidence data that would be sent to the LLM prompt
    # This runs the full pipeline but does NOT call the LLM
    python scripts/llm_tools.py data BRAF V600E Melanoma

    # Without tumor type (pan-cancer context)
    python scripts/llm_tools.py data KRAS G12D

    # Save output to file
    python scripts/llm_tools.py data EGFR T790M NSCLC > data/egfr_t790m_nsclc.txt

The 'data' command outputs:
    - DATA AVAILABILITY FLAGS: What data sources have evidence
    - EVIDENCE ASSESSMENT: Overall quality, well-characterized aspects, gaps
    - LOCUS MATCH SUMMARY: Variant vs codon vs gene-level evidence breakdown
    - TUMOR MATCH SUMMARY: Tumor-specific vs pan-cancer evidence breakdown
    - THERAPEUTIC SIGNALS: Sensitivity and resistance summaries
    - BIOLOGICAL CONTEXT: Gene role, pathway, cBioPortal, DepMap, hotspots data
    - EVIDENCE SUMMARY: Compact therapeutic evidence from FDA, CIViC, VICC, CGI
    - LITERATURE SUMMARY: PubMed articles and extracted knowledge
"""

import asyncio
import json
import sys


async def run_baseline(gene: str, variant: str, tumor_type: str | None):
    """Test LLM baseline knowledge for a variant (no evidence provided)."""
    from oncomind.llm.service import LLMService

    service = LLMService()
    result = await service.test_llm_baseline_knowledge(gene, variant, tumor_type)

    print("\n" + "=" * 80)
    print("PARSED RESULT:")
    print("=" * 80)
    print(json.dumps(result, indent=2))


async def show_prompt_data(gene: str, variant: str, tumor_type: str | None):
    """Run the evidence pipeline and show all data that would be sent to the LLM.

    This runs the full pipeline up to prompt creation, but does NOT call the LLM.
    Useful for debugging what data the LLM sees.
    """
    from oncomind.config.constants import TUMOR_TYPE_MAPPINGS
    from oncomind.insight_builder.conductor import Conductor, ConductorConfig

    # Run with LLM disabled - we just want the evidence
    config = ConductorConfig(enable_llm=False)

    variant_str = f"{gene} {variant}"

    async with Conductor(config) as conductor:
        result = await conductor.run(variant_str, tumor_type)

    evidence = result.evidence

    # Compute gaps
    evidence_gaps = evidence.compute_evidence_gaps()

    # Build all the data that would go to the LLM
    evidence_summary = evidence.get_evidence_summary_for_llm()
    biological_context = evidence.get_biological_context_for_llm()
    literature_summary = evidence.get_literature_summary_for_llm()
    evidence_assessment = evidence_gaps.to_dict_for_llm()
    resistance_summary = evidence.get_resistance_summary()
    sensitivity_summary = evidence.get_sensitivity_summary()
    locus_match_summary = evidence.get_locus_match_summary()
    tumor_match_summary = evidence.get_tumor_match_summary()

    # Compute data availability flags
    data_availability = {
        "has_tumor_specific_cbioportal": _check_tumor_specific_cbioportal(
            evidence, tumor_type
        ),
        "has_civic_assertions": bool(evidence.civic_assertions),
        "has_fda_approvals": bool(evidence.fda_approvals),
        "has_vicc_evidence": bool(evidence.get_vicc_unique()),
    }

    # Print everything
    print("\n" + "=" * 80)
    print(f"LLM PROMPT DATA FOR: {gene} {variant} ({tumor_type or 'Pan-cancer'})")
    print("=" * 80)

    print("\n### DATA AVAILABILITY FLAGS ###")
    print(json.dumps(data_availability, indent=2))

    print("\n### EVIDENCE ASSESSMENT ###")
    print(json.dumps(evidence_assessment, indent=2))

    print("\n### LOCUS MATCH SUMMARY ###")
    print(json.dumps(locus_match_summary, indent=2))

    print("\n### TUMOR MATCH SUMMARY ###")
    print(json.dumps(tumor_match_summary, indent=2))

    print("\n### THERAPEUTIC SIGNALS ###")
    print(f"Sensitivity: {sensitivity_summary}")
    print(f"Resistance: {resistance_summary}")

    print("\n### BIOLOGICAL CONTEXT ###")
    print(biological_context or "(empty)")

    print("\n### EVIDENCE SUMMARY ###")
    print(evidence_summary or "(empty)")

    print("\n### LITERATURE SUMMARY ###")
    print(literature_summary or "(empty)")

    print("\n" + "=" * 80)


def _check_tumor_specific_cbioportal(evidence, tumor_type: str | None) -> bool:
    """Check if cBioPortal data is tumor-specific (not pan-cancer)."""
    from oncomind.config.constants import TUMOR_TYPE_MAPPINGS

    if not tumor_type:
        return False

    if not evidence.cbioportal_evidence or not evidence.cbioportal_evidence.has_data():
        return False

    study_id = evidence.cbioportal_evidence.study_id or ""
    study_lower = study_id.lower()

    pan_cancer_indicators = [
        "pan_cancer",
        "pancancer",
        "all_tcga",
        "msk_impact",
        "mixed",
        "combined",
        "multi",
    ]

    if any(indicator in study_lower for indicator in pan_cancer_indicators):
        return False

    tumor_lower = tumor_type.lower().replace(" ", "_")
    tumor_words = tumor_type.lower().split()

    for word in tumor_words:
        if word in study_lower:
            return True

    for abbrev, synonyms in TUMOR_TYPE_MAPPINGS.items():
        if tumor_lower in [s.lower() for s in synonyms] or tumor_lower == abbrev:
            if abbrev in study_lower:
                return True
        for word in tumor_words:
            if any(word in syn.lower() for syn in synonyms):
                if abbrev in study_lower:
                    return True

    return False


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        print("\nExamples:")
        print("  python scripts/llm_tools.py baseline BRAF V600E")
        print("  python scripts/llm_tools.py data EGFR T790M NSCLC")
        sys.exit(1)

    command = sys.argv[1]
    gene = sys.argv[2]
    variant = sys.argv[3]
    tumor_type = sys.argv[4] if len(sys.argv) > 4 else None

    if command == "baseline":
        asyncio.run(run_baseline(gene, variant, tumor_type))
    elif command == "data":
        asyncio.run(show_prompt_data(gene, variant, tumor_type))
    else:
        print(f"Unknown command: {command}")
        print("Available commands: baseline, data")
        sys.exit(1)


if __name__ == "__main__":
    main()

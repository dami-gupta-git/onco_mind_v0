#!/usr/bin/env python3
"""LLM debugging and inspection tools.

This script provides utilities for debugging and inspecting the data pipeline
that feeds into the LLM for variant annotation.

Commands:
    baseline  - Test LLM baseline knowledge (no evidence provided)
    data      - Show all data that would be sent to the LLM prompt
    prompt    - Generate and optionally run the full LLM prompt

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

    # Generate the full prompt (but don't run LLM)
    python scripts/llm_tools.py prompt BRAF V600E Melanoma

    # Generate and RUN the prompt through the LLM
    python scripts/llm_tools.py prompt BRAF V600E Melanoma --run

    # Run with a specific model
    python scripts/llm_tools.py prompt AKT1 E17K "Breast Cancer" --run --model gpt-4o-mini

The 'data' command outputs:
    - DATA AVAILABILITY FLAGS: What data sources have evidence
    - EVIDENCE ASSESSMENT: Overall quality, well-characterized aspects, gaps
    - LOCUS MATCH SUMMARY: Variant vs codon vs gene-level evidence breakdown
    - TUMOR MATCH SUMMARY: Tumor-specific vs pan-cancer evidence breakdown
    - THERAPEUTIC SIGNALS: Sensitivity and resistance summaries
    - BIOLOGICAL CONTEXT: Gene role, pathway, cBioPortal, DepMap, hotspots data
    - EVIDENCE SUMMARY: Compact therapeutic evidence from FDA, CIViC, VICC, CGI
    - LITERATURE SUMMARY: PubMed articles and extracted knowledge

The 'prompt' command outputs:
    - The full system prompt
    - The full user prompt (with all evidence inserted)
    - Optionally: the LLM response (with --run flag)
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


async def generate_and_run_prompt(
    gene: str,
    variant: str,
    tumor_type: str | None,
    run_llm: bool = False,
    model: str | None = None,
):
    """Generate the full LLM prompt and optionally run it.

    This is useful for iterating on prompt design - you can see exactly what
    the LLM receives and what it outputs.
    """
    from oncomind.config.constants import TUMOR_TYPE_MAPPINGS
    from oncomind.insight_builder.conductor import Conductor, ConductorConfig
    from oncomind.llm.prompts import create_synthesis_prompt

    # Run with LLM disabled initially - we build the prompt ourselves
    config = ConductorConfig(enable_llm=False)

    variant_str = f"{gene} {variant}"

    print(f"\n{'='*80}")
    print(f"GATHERING EVIDENCE FOR: {gene} {variant} ({tumor_type or 'Pan-cancer'})")
    print("="*80)

    async with Conductor(config) as conductor:
        result = await conductor.run(variant_str, tumor_type)

    evidence = result.evidence

    # Compute all the data needed for the prompt
    evidence_gaps = evidence.compute_evidence_gaps()
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

    # Build the prompt
    messages = create_synthesis_prompt(
        gene=gene,
        variant=variant,
        tumor_type=tumor_type,
        biological_context=biological_context or "",
        evidence_summary=evidence_summary or "",
        evidence_assessment=evidence_assessment,
        literature_summary=literature_summary or "",
        data_availability=data_availability,
        resistance_summary=resistance_summary or "",
        sensitivity_summary=sensitivity_summary or "",
        locus_match_summary=locus_match_summary,
        tumor_match_summary=tumor_match_summary,
    )

    # Print the prompt
    print("\n" + "="*80)
    print("SYSTEM PROMPT")
    print("="*80)
    print(messages[0]["content"])

    print("\n" + "="*80)
    print("USER PROMPT")
    print("="*80)
    print(messages[1]["content"])

    # Optionally run through LLM
    if run_llm:
        print("\n" + "="*80)
        print(f"RUNNING LLM ({model or 'default model'})...")
        print("="*80)

        from oncomind.llm.service import LLMService

        service = LLMService()
        if model:
            service.model = model

        # Call the LLM directly with our messages
        response = await service._call_llm(messages)

        print("\n" + "="*80)
        print("LLM RAW RESPONSE")
        print("="*80)
        print(response)

        # Try to parse as JSON
        try:
            import re
            # Extract JSON from response (handle markdown code blocks)
            json_match = re.search(r'```(?:json)?\s*([\s\S]*?)\s*```', response)
            if json_match:
                json_str = json_match.group(1)
            else:
                json_str = response

            parsed = json.loads(json_str)
            print("\n" + "="*80)
            print("PARSED JSON RESPONSE")
            print("="*80)
            print(json.dumps(parsed, indent=2))
        except json.JSONDecodeError as e:
            print(f"\n(Could not parse response as JSON: {e})")

    print("\n" + "="*80)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="LLM debugging and inspection tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/llm_tools.py baseline BRAF V600E
  python scripts/llm_tools.py data EGFR T790M NSCLC
  python scripts/llm_tools.py prompt BRAF V600E Melanoma
  python scripts/llm_tools.py prompt AKT1 E17K "Breast Cancer" --run
  python scripts/llm_tools.py prompt KRAS G12D --run --model gpt-4o-mini
        """
    )

    parser.add_argument(
        "command",
        choices=["baseline", "data", "prompt"],
        help="Command to run: baseline (test LLM knowledge), data (show evidence), prompt (generate/run prompt)"
    )
    parser.add_argument("gene", help="Gene symbol (e.g., BRAF)")
    parser.add_argument("variant", help="Variant notation (e.g., V600E)")
    parser.add_argument("tumor_type", nargs="?", default=None, help="Tumor type (optional)")
    parser.add_argument("--run", action="store_true", help="Run the prompt through the LLM (for 'prompt' command)")
    parser.add_argument("--model", type=str, default=None, help="LLM model to use (for 'prompt --run')")

    args = parser.parse_args()

    if args.command == "baseline":
        asyncio.run(run_baseline(args.gene, args.variant, args.tumor_type))
    elif args.command == "data":
        asyncio.run(show_prompt_data(args.gene, args.variant, args.tumor_type))
    elif args.command == "prompt":
        asyncio.run(generate_and_run_prompt(
            args.gene,
            args.variant,
            args.tumor_type,
            run_llm=args.run,
            model=args.model,
        ))


if __name__ == "__main__":
    main()

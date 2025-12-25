"""
Backend logic for OncoMind Streamlit app.

Integrates with the oncomind package for:
- Variant insight generation (single and batch)
- Evidence gathering from multiple databases

ARCHITECTURE:
    All paths use InsightBuilder → Insight
    With optional LLMService → Insight.llm for narrative generation

Modes:
    - Lite mode: InsightBuilder only (fast, no LLM)
    - Default mode: InsightBuilder + LLMService narrative
    - Full mode: InsightBuilder + Literature + LLMService narrative
"""

from typing import Any, Dict, List, Optional, Callable

from oncomind.insight_builder import InsightBuilder, InsightBuilderConfig
from oncomind.llm.service import LLMService


async def get_variant_insight(
    gene: str,
    variant: str,
    tumor_type: Optional[str] = None,
    enable_llm: bool = False,
    enable_literature: bool = False,
    model: str = "gpt-4o-mini",
    temperature: float = 0.1
) -> Dict[str, Any]:
    """
    Generate insight for a single variant.

    Args:
        gene: Gene symbol (e.g., BRAF)
        variant: Variant notation (e.g., V600E)
        tumor_type: Optional tumor type (e.g., Melanoma)
        enable_llm: Enable LLM narrative generation
        enable_literature: Enable literature search
        model: LLM model to use (if enable_llm=True)
        temperature: LLM temperature (0.0-1.0)

    Returns:
        Dict containing insight results with identifiers, evidence, etc.
    """
    try:
        # Step 1: Build insight using InsightBuilder
        config = InsightBuilderConfig(enable_literature=enable_literature)
        async with InsightBuilder(config) as builder:
            insight = await builder.build_insight(f"{gene} {variant}", tumor_type=tumor_type)

        # Step 2: Generate LLM narrative if enabled
        if enable_llm:
            # Get evidence summary for LLM
            evidence_summary = insight.get_evidence_summary_for_llm()

            # Generate LLM insight
            llm_service = LLMService(model=model, temperature=temperature)
            llm_insight = await llm_service.get_llm_insight(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                evidence_summary=evidence_summary,
                has_clinical_trials=bool(insight.clinical.clinical_trials),
            )

            # Embed LLM insight in the main insight object
            insight.llm = llm_insight

        # Build response
        return _build_response(insight)

    except Exception as e:
        return {"error": f"Insight generation failed: {str(e)}"}


async def batch_get_variant_insights(
    variants: List[Dict[str, str]],
    enable_llm: bool = False,
    enable_literature: bool = False,
    model: str = "gpt-4o-mini",
    temperature: float = 0.1,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> List[Dict[str, Any]]:
    """
    Generate insights for multiple variants.

    Args:
        variants: List of dicts with 'gene', 'variant', and optional 'tumor_type'
        enable_llm: Enable LLM narrative generation
        enable_literature: Enable literature search
        model: LLM model to use (if enable_llm=True)
        temperature: LLM temperature (0.0-1.0)
        progress_callback: Optional callback(current, total) for progress updates

    Returns:
        List of insight results
    """
    results = []
    total = len(variants)

    for i, v in enumerate(variants):
        if progress_callback:
            progress_callback(i + 1, total)

        try:
            result = await get_variant_insight(
                gene=v['gene'],
                variant=v['variant'],
                tumor_type=v.get('tumor_type'),
                enable_llm=enable_llm,
                enable_literature=enable_literature,
                model=model,
                temperature=temperature
            )
            results.append(result)
        except Exception as e:
            results.append({
                "variant": {
                    "gene": v['gene'],
                    "variant": v['variant'],
                    "tumor_type": v.get('tumor_type'),
                },
                "error": str(e)
            })

    return results


# === Private helper functions ===


def _build_response(insight) -> Dict[str, Any]:
    """Build the standard response dict from an Insight object."""
    llm = insight.llm  # May be None if LLM was not enabled

    return {
        "variant": {
            "gene": insight.identifiers.gene,
            "variant": insight.identifiers.variant,
            "tumor_type": insight.clinical.tumor_type,
        },
        "insight": {
            "summary": insight.get_summary(),  # Always 1-line summary
            "llm_narrative": llm.llm_summary if llm else None,  # LLM insight when available
            "rationale": llm.rationale if llm else None,
        },
        "identifiers": {
            "cosmic_id": insight.identifiers.cosmic_id,
            "ncbi_gene_id": insight.identifiers.ncbi_gene_id,
            "dbsnp_id": insight.identifiers.dbsnp_id,
            "clinvar_id": insight.identifiers.clinvar_id,
        },
        "hgvs": {
            "genomic": insight.identifiers.hgvs_genomic,
            "protein": insight.identifiers.hgvs_protein,
            "transcript": insight.identifiers.hgvs_transcript,
        },
        "clinvar": {
            "clinical_significance": insight.clinical.clinvar_clinical_significance,
            "accession": insight.clinical.clinvar_accession,
        },
        "annotations": {
            "snpeff_effect": insight.functional.snpeff_effect,
            "polyphen2_prediction": insight.functional.polyphen2_prediction,
            "cadd_score": insight.functional.cadd_score,
            "gnomad_exome_af": insight.functional.gnomad_exome_af,
            "alphamissense_score": insight.functional.alphamissense_score,
            "alphamissense_prediction": insight.functional.alphamissense_prediction,
        },
        "transcript": {
            "id": insight.identifiers.transcript_id,
            "consequence": insight.identifiers.transcript_consequence,
        },
        "evidence_panel": {
            "clinical": {
                "clinical_trials": [
                    {
                        "nct_id": t.nct_id,
                        "title": t.title,
                        "phase": t.phase,
                        "status": t.status,
                        "drugs": t.interventions,  # interventions contains drug names
                        "conditions": t.conditions,
                        "url": t.url,
                        "variant_specific": t.variant_specific,
                    }
                    for t in insight.clinical.clinical_trials
                ],
                "fda_approvals": [
                    {
                        "drug_name": a.drug_name,
                        "brand_name": a.brand_name,
                        "indication": a.indication,
                    }
                    for a in insight.clinical.fda_approvals
                ],
            },
        },
        "recommended_therapies": (
            [
                {
                    "drug_name": t.drug_name,
                    "evidence_level": t.evidence_level,
                    "approval_status": t.approval_status,
                    "clinical_context": t.clinical_context,
                }
                for t in llm.recommended_therapies
            ]
            if llm
            else _extract_therapies(insight)
        ),
        "insight_data": insight.model_dump(mode="json"),
    }


def _extract_therapies(insight) -> List[Dict[str, Any]]:
    """Extract therapy recommendations from an Insight object."""
    therapies = []

    # From FDA approvals
    for approval in insight.clinical.fda_approvals:
        drug_name = approval.brand_name or approval.generic_name or approval.drug_name
        if drug_name:
            therapies.append({
                "drug_name": drug_name,
                "evidence_level": "A",
                "approval_status": "FDA Approved",
                "clinical_context": approval.indication if approval.indication else "",
            })

    # From CGI biomarkers
    for biomarker in insight.kb.cgi_biomarkers:
        if biomarker.fda_approved and biomarker.drug:
            therapies.append({
                "drug_name": biomarker.drug,
                "evidence_level": biomarker.evidence_level or "A",
                "approval_status": "FDA Approved (CGI)",
                "clinical_context": biomarker.tumor_type or "",
            })

    # Deduplicate by drug name
    seen_drugs = set()
    unique_therapies = []
    for t in therapies:
        if t["drug_name"] not in seen_drugs:
            seen_drugs.add(t["drug_name"])
            unique_therapies.append(t)

    return unique_therapies[:10]

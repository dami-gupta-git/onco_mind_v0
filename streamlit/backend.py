"""
Backend logic for OncoMind Streamlit app.

Integrates with the oncomind package for:
- Variant insight generation (single and batch)
- Evidence gathering from multiple databases

ARCHITECTURE:
    All paths use Conductor → Result (evidence + optional LLM narrative)

Modes:
    - Lite mode: Evidence only (fast, no LLM)
    - Default mode: Evidence + LLM narrative
    - Full mode: Evidence + Literature + LLM narrative
"""

from typing import Any, Dict, List, Optional, Callable

from oncomind.insight_builder import Conductor, ConductorConfig


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
        # Configure and run the Conductor
        config = ConductorConfig(
            enable_literature=enable_literature,
            enable_llm=enable_llm,
            llm_model=model,
            llm_temperature=temperature,
        )
        async with Conductor(config) as conductor:
            result = await conductor.run(f"{gene} {variant}", tumor_type=tumor_type)

        # Build response
        return _build_response(result)

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


def _build_response(result) -> Dict[str, Any]:
    """Build the standard response dict from a Result object."""
    llm = result.llm  # May be None if LLM was not enabled

    return {
        "variant": {
            "gene": result.identifiers.gene,
            "variant": result.identifiers.variant,
            "tumor_type": result.clinical.tumor_type,
        },
        "insight": {
            "summary": result.get_summary(),  # Always 1-line summary
            "llm_narrative": llm.llm_summary if llm else None,  # LLM insight when available
            "rationale": llm.rationale if llm else None,
        },
        "identifiers": {
            "cosmic_id": result.identifiers.cosmic_id,
            "ncbi_gene_id": result.identifiers.ncbi_gene_id,
            "dbsnp_id": result.identifiers.dbsnp_id,
            "clinvar_id": result.identifiers.clinvar_id,
        },
        "hgvs": {
            "genomic": result.identifiers.hgvs_genomic,
            "protein": result.identifiers.hgvs_protein,
            "transcript": result.identifiers.hgvs_transcript,
        },
        "clinvar": {
            "clinical_significance": result.clinical.clinvar_clinical_significance,
            "accession": result.clinical.clinvar_accession,
        },
        "annotations": {
            "snpeff_effect": result.functional.snpeff_effect,
            "polyphen2_prediction": result.functional.polyphen2_prediction,
            "cadd_score": result.functional.cadd_score,
            "gnomad_exome_af": result.functional.gnomad_exome_af,
            "alphamissense_score": result.functional.alphamissense_score,
            "alphamissense_prediction": result.functional.alphamissense_prediction,
        },
        "transcript": {
            "id": result.identifiers.transcript_id,
            "consequence": result.identifiers.transcript_consequence,
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
                    for t in result.clinical.clinical_trials
                ],
                "fda_approvals": [
                    {
                        "drug_name": a.drug_name,
                        "brand_name": a.brand_name,
                        "indication": a.indication,
                    }
                    for a in result.clinical.fda_approvals
                ],
            },
        },
        "recommended_therapies": [
            {
                "drug_name": t.drug_name,
                "evidence_level": t.evidence_level,
                "approval_status": t.approval_status,
                "clinical_context": t.clinical_context,
            }
            for t in (llm.recommended_therapies if llm else result.evidence.get_recommended_therapies())
        ],
        "result_data": result.model_dump(mode="json"),
    }

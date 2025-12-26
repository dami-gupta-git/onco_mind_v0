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
    """Build the standard response dict from a Result object.

    Uses the flat Evidence structure where each source is a simple list:
    - evidence.fda_approvals, .civic_assertions, .civic_evidence, etc.
    """
    llm = result.llm  # May be None if LLM was not enabled
    evidence = result.evidence  # Evidence object with flat list structure

    return {
        "variant": {
            "gene": result.identifiers.gene,
            "variant": result.identifiers.variant,
            "tumor_type": result.context.tumor_type,
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
            "clinical_significance": evidence.clinvar_significance,
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
        # Per-source evidence (flat lists - frontend decides how to display)
        "fda_approvals": [
            {
                "drug_name": a.drug_name,
                "brand_name": a.brand_name,
                "generic_name": a.generic_name,
                "indication": a.indication,
                "approval_date": a.approval_date,
                "companion_diagnostic": a.companion_diagnostic,
                "black_box_warning": a.black_box_warning,
                "dosing_for_variant": a.dosing_for_variant,
            }
            for a in evidence.fda_approvals
        ],
        "civic_assertions": [
            {
                "id": a.assertion_id,
                "therapies": a.therapies,
                "disease": a.disease,
                "significance": a.significance,
                "amp_level": a.amp_level,
                "amp_tier": a.amp_tier,
                "description": a.description,
            }
            for a in evidence.civic_assertions
        ],
        "civic_evidence": [
            {
                "evidence_type": e.evidence_type,
                "evidence_level": e.evidence_level,
                "clinical_significance": e.clinical_significance,
                "disease": e.disease,
                "drugs": e.drugs,
                "description": e.description,
                "pmid": e.pmid,
                "source_url": e.source_url,
                "trust_rating": e.trust_rating or e.rating,  # Use trust_rating if available, else rating
            }
            for e in evidence.civic_evidence
        ],
        "vicc_evidence": [
            {
                "source": v.source,
                "drugs": v.drugs,
                "disease": v.disease,
                "response_type": v.response_type,
                "evidence_level": v.evidence_level,
                "molecular_profile": v.molecular_profile,
                "molecular_profile_score": v.molecular_profile_score,
            }
            for v in evidence.vicc_evidence
        ],
        "cgi_biomarkers": [
            {
                "drug": b.drug,
                "association": b.association,
                "tumor_type": b.tumor_type,
                "evidence_level": b.evidence_level,
                "fda_approved": b.fda_approved,
            }
            for b in evidence.cgi_biomarkers
        ],
        "clinvar_entries": [
            {
                "clinical_significance": c.clinical_significance,
                "conditions": c.conditions,
                "review_status": c.review_status,
            }
            for c in evidence.clinvar_entries
        ],
        "cosmic_entries": [
            {
                "mutation_id": c.mutation_id,
                "primary_site": c.primary_site,
                "sample_count": c.sample_count,
            }
            for c in evidence.cosmic_entries
        ],
        "clinical_trials": [
            {
                "nct_id": t.nct_id,
                "title": t.title,
                "phase": t.phase,
                "status": t.status,
                "drugs": t.interventions,
                "conditions": t.conditions,
                "url": t.url,
                "variant_specific": t.variant_specific,
            }
            for t in evidence.clinical_trials
        ],
        "pubmed_articles": [
            {
                "pmid": a.pmid,
                "title": a.title,
                "year": a.year,
                "journal": a.journal,
                "signal_type": a.signal_type,
            }
            for a in evidence.pubmed_articles
        ],
        "literature_knowledge": evidence.literature_knowledge.summary if evidence.literature_knowledge else None,
        "preclinical_biomarkers": [
            {
                "drug": b.drug,
                "association": b.association,
                "evidence_level": b.evidence_level,
            }
            for b in evidence.preclinical_biomarkers
        ],
        "early_phase_biomarkers": [
            {
                "drug": b.drug,
                "association": b.association,
                "evidence_level": b.evidence_level,
            }
            for b in evidence.early_phase_biomarkers
        ],
        "cbioportal_evidence": {
            "gene": evidence.cbioportal_evidence.gene,
            "variant": evidence.cbioportal_evidence.variant,
            "tumor_type": evidence.cbioportal_evidence.tumor_type,
            "study_id": evidence.cbioportal_evidence.study_id,
            "total_samples": evidence.cbioportal_evidence.total_samples,
            "samples_with_gene_mutation": evidence.cbioportal_evidence.samples_with_gene_mutation,
            "samples_with_exact_variant": evidence.cbioportal_evidence.samples_with_exact_variant,
            "gene_prevalence_pct": evidence.cbioportal_evidence.gene_prevalence_pct,
            "variant_prevalence_pct": evidence.cbioportal_evidence.variant_prevalence_pct,
            "co_occurring": [c.model_dump() for c in evidence.cbioportal_evidence.co_occurring],
            "mutually_exclusive": [m.model_dump() for m in evidence.cbioportal_evidence.mutually_exclusive],
        } if evidence.cbioportal_evidence else None,
        "recommended_therapies": [
            {
                "drug_name": t.drug_name,
                "evidence_level": t.evidence_level,
                "approval_status": t.approval_status,
                "clinical_context": t.clinical_context,
            }
            for t in (llm.recommended_therapies if llm else evidence.get_recommended_therapies())
        ],
        "result_data": result.model_dump(mode="json"),
    }

"""
Backend logic for OncoMind Streamlit app.

Integrates with the oncomind package for:
- Variant insight generation (single and batch)
- Evidence gathering from multiple databases

ARCHITECTURE:
    All paths use EvidenceBuilder → EvidencePanel
    With optional LLMService → LLMInsight for narrative generation

Modes:
    - Lite mode: EvidenceBuilder only (fast, no LLM)
    - Default mode: EvidenceBuilder + LLMService narrative
    - Full mode: EvidenceBuilder + Literature + LLMService narrative
"""

from typing import Any, Dict, List, Optional, Callable

from oncomind.evidence import EvidenceBuilder, EvidenceBuilderConfig
from oncomind.llm.service import LLMService
from oncomind.models.evidence import EvidenceForLLM


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
        # Step 1: Build evidence panel using EvidenceBuilder
        config = EvidenceBuilderConfig(enable_literature=enable_literature)
        async with EvidenceBuilder(config) as builder:
            panel = await builder.build_evidence_panel(f"{gene} {variant}", tumor_type=tumor_type)

        # Step 2: Generate LLM narrative if enabled
        llm_insight = None
        if enable_llm:
            # Convert EvidencePanel to EvidenceForLLM for the LLM service
            evidence = _panel_to_evidence_for_llm(panel)

            # Add literature if available
            if panel.literature.pubmed_articles:
                evidence.pubmed_articles = panel.literature.pubmed_articles
                evidence.literature_knowledge = panel.literature.literature_knowledge

            # Generate LLM insight
            llm_service = LLMService(model=model, temperature=temperature)
            llm_insight = await llm_service.get_llm_insight(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                evidence=evidence,
            )

        # Build response
        return _build_response(panel, llm_insight)

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


def _panel_to_evidence_for_llm(panel) -> EvidenceForLLM:
    """Convert an EvidencePanel to EvidenceForLLM for LLM consumption."""
    return EvidenceForLLM(
        variant_id=panel.identifiers.variant_id,
        gene=panel.identifiers.gene,
        variant=panel.identifiers.variant,
        cosmic_id=panel.identifiers.cosmic_id,
        ncbi_gene_id=panel.identifiers.ncbi_gene_id,
        dbsnp_id=panel.identifiers.dbsnp_id,
        clinvar_id=panel.identifiers.clinvar_id,
        hgvs_genomic=panel.identifiers.hgvs_genomic,
        hgvs_protein=panel.identifiers.hgvs_protein,
        hgvs_transcript=panel.identifiers.hgvs_transcript,
        transcript_id=panel.identifiers.transcript_id,
        transcript_consequence=panel.identifiers.transcript_consequence,
        civic=panel.kb.civic,
        civic_assertions=panel.kb.civic_assertions,
        clinvar=panel.kb.clinvar,
        cosmic=panel.kb.cosmic,
        cgi_biomarkers=panel.kb.cgi_biomarkers,
        vicc=panel.kb.vicc,
        fda_approvals=panel.clinical.fda_approvals,
        clinical_trials=panel.clinical.clinical_trials,
        alphamissense_score=panel.functional.alphamissense_score,
        alphamissense_prediction=panel.functional.alphamissense_prediction,
        cadd_score=panel.functional.cadd_score,
        polyphen2_prediction=panel.functional.polyphen2_prediction,
        sift_prediction=panel.functional.sift_prediction,
        gnomad_exome_af=panel.functional.gnomad_exome_af,
        clinvar_clinical_significance=panel.clinical.clinvar_clinical_significance,
    )


def _build_response(panel, llm_insight=None) -> Dict[str, Any]:
    """Build the standard response dict from panel and optional LLM insight."""
    return {
        "variant": {
            "gene": panel.identifiers.gene,
            "variant": panel.identifiers.variant,
            "tumor_type": panel.clinical.tumor_type,
        },
        "insight": {
            "summary": panel.get_summary(),  # Always 1-line summary
            "llm_narrative": llm_insight.llm_summary if llm_insight else None,  # LLM insight when available
            "rationale": llm_insight.rationale if llm_insight else None,
            "evidence_strength": llm_insight.evidence_strength if llm_insight else panel.meta.evidence_strength,
        },
        "identifiers": {
            "cosmic_id": panel.identifiers.cosmic_id,
            "ncbi_gene_id": panel.identifiers.ncbi_gene_id,
            "dbsnp_id": panel.identifiers.dbsnp_id,
            "clinvar_id": panel.identifiers.clinvar_id,
        },
        "hgvs": {
            "genomic": panel.identifiers.hgvs_genomic,
            "protein": panel.identifiers.hgvs_protein,
            "transcript": panel.identifiers.hgvs_transcript,
        },
        "clinvar": {
            "clinical_significance": panel.clinical.clinvar_clinical_significance,
            "accession": panel.clinical.clinvar_accession,
        },
        "annotations": {
            "snpeff_effect": panel.functional.snpeff_effect,
            "polyphen2_prediction": panel.functional.polyphen2_prediction,
            "cadd_score": panel.functional.cadd_score,
            "gnomad_exome_af": panel.functional.gnomad_exome_af,
            "alphamissense_score": panel.functional.alphamissense_score,
            "alphamissense_prediction": panel.functional.alphamissense_prediction,
        },
        "transcript": {
            "id": panel.identifiers.transcript_id,
            "consequence": panel.identifiers.transcript_consequence,
        },
        "recommended_therapies": (
            [
                {
                    "drug_name": t.drug_name,
                    "evidence_level": t.evidence_level,
                    "approval_status": t.approval_status,
                    "clinical_context": t.clinical_context,
                }
                for t in llm_insight.recommended_therapies
            ]
            if llm_insight
            else _extract_therapies(panel)
        ),
        "evidence_panel": panel.model_dump(mode="json"),
    }


def _extract_therapies(panel) -> List[Dict[str, Any]]:
    """Extract therapy recommendations from an EvidencePanel."""
    therapies = []

    # From FDA approvals
    for approval in panel.clinical.fda_approvals:
        drug_name = approval.brand_name or approval.generic_name or approval.drug_name
        if drug_name:
            therapies.append({
                "drug_name": drug_name,
                "evidence_level": "A",
                "approval_status": "FDA Approved",
                "clinical_context": approval.indication if approval.indication else "",
            })

    # From CGI biomarkers
    for biomarker in panel.kb.cgi_biomarkers:
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

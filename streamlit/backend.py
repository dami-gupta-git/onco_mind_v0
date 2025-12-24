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

import asyncio
from typing import Any, Dict, List, Optional, Callable

from oncomind import get_insight, get_insights, InsightConfig
from oncomind.evidence import EvidenceBuilder, EvidenceBuilderConfig
from oncomind.llm.service import LLMService
from oncomind.models.evidence import EvidenceForLLM


async def get_variant_annotation(
    gene: str,
    variant: str,
    tumor_type: Optional[str] = None,
    enable_llm: bool = False,
    enable_literature: bool = True,
    model: str = "gpt-4o-mini",
    temperature: float = 0.1
) -> Dict[str, Any]:
    """
    Generate insight for a single variant.

    This is the primary function for variant insight in the UI.

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
        result = {
            "variant": {
                "gene": panel.identifiers.gene,
                "variant": panel.identifiers.variant,
                "tumor_type": panel.clinical.tumor_type,
            },
            "insight": {
                "summary": llm_insight.summary if llm_insight else _generate_summary(panel),
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

        return result

    except Exception as e:
        return {"error": f"Insight generation failed: {str(e)}"}


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


def _generate_summary(panel) -> str:
    """Generate a summary string from an EvidencePanel."""
    parts = []

    gene = panel.identifiers.gene
    variant = panel.identifiers.variant
    parts.append(f"{gene} {variant}")

    if panel.clinical.gene_role:
        parts.append(f"({panel.clinical.gene_role})")

    drugs = panel.clinical.get_approved_drugs()
    if drugs:
        parts.append(f"has FDA-approved therapies: {', '.join(drugs[:3])}")
    else:
        parts.append("has no FDA-approved targeted therapies")

    if panel.clinical.clinical_trials:
        variant_specific = sum(1 for t in panel.clinical.clinical_trials if t.variant_specific)
        if variant_specific:
            parts.append(f"with {variant_specific} variant-specific clinical trials")
        else:
            parts.append(f"with {len(panel.clinical.clinical_trials)} relevant clinical trials")

    return " ".join(parts)


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
                "clinical_context": approval.indication[:200] if approval.indication else "",
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


# Backward compatibility alias
async def get_single_variant_insight(
    gene: str,
    variant: str,
    tumor_type: Optional[str] = None,
    enable_llm: bool = False,
    enable_literature: bool = True,
    model: str = "gpt-4o-mini",
    temperature: float = 0.1
) -> Dict[str, Any]:
    """Alias for get_variant_annotation (backward compatibility)."""
    return await get_variant_annotation(
        gene, variant, tumor_type, enable_llm, enable_literature, model, temperature
    )


# Backward compatibility alias - routes to get_variant_annotation with enable_llm=True
async def get_variant_insight(
    gene: str,
    variant: str,
    tumor_type: Optional[str] = None,
    model: str = "gpt-4o-mini",
    temperature: float = 0.1
) -> Dict[str, Any]:
    """
    Generate insight with LLM narrative.

    This is an alias for get_variant_annotation with enable_llm=True.
    Maintained for backward compatibility.
    """
    return await get_variant_annotation(
        gene=gene,
        variant=variant,
        tumor_type=tumor_type,
        enable_llm=True,
        enable_literature=True,
        model=model,
        temperature=temperature,
    )


async def batch_get_variant_annotations(
    variants: List[Dict[str, str]],
    model: str = "gpt-4o-mini",
    temperature: float = 0.1,
    enable_llm: bool = False,
    progress_callback: Optional[Callable[[int, int], None]] = None
) -> List[Dict[str, Any]]:
    """
    Generate annotations for multiple variants using the public API.

    Args:
        variants: List of dicts with 'gene', 'variant', and optional 'tumor_type'
        model: LLM model to use (if enable_llm=True)
        temperature: LLM temperature (0.0-1.0)
        enable_llm: Enable LLM enhancement for literature analysis
        progress_callback: Optional callback(current, total) for progress updates

    Returns:
        List of annotation results
    """
    try:
        config = InsightConfig(
            enable_llm=enable_llm,
            llm_model=model,
            llm_temperature=temperature,
        )

        # Convert to variant strings
        variant_strs = []
        tumor_types = []
        for v in variants:
            variant_strs.append(f"{v['gene']} {v['variant']}")
            tumor_types.append(v.get('tumor_type'))

        # Process variants
        results = []
        panels = await get_insights(
            variant_strs,
            config=config,
            progress_callback=progress_callback,
        )

        for i, panel in enumerate(panels):
            # Apply tumor type if provided
            if tumor_types[i]:
                panel.clinical.tumor_type = tumor_types[i]

            results.append({
                "variant": {
                    "gene": panel.identifiers.gene,
                    "variant": panel.identifiers.variant,
                    "tumor_type": panel.clinical.tumor_type,
                },
                "insight": {
                    "summary": _generate_summary(panel),
                    "evidence_strength": panel.meta.evidence_strength,
                },
                "identifiers": {
                    "cosmic_id": panel.identifiers.cosmic_id,
                    "ncbi_gene_id": panel.identifiers.ncbi_gene_id,
                    "dbsnp_id": panel.identifiers.dbsnp_id,
                    "clinvar_id": panel.identifiers.clinvar_id,
                },
                "recommended_therapies": _extract_therapies(panel),
                "evidence_panel": panel.model_dump(mode="json"),
            })

        return results

    except Exception as e:
        return [{"error": f"Batch insight generation failed: {str(e)}"}]


async def batch_get_variant_insights(
    variants: List[Dict[str, str]],
    model: str = "gpt-4o-mini",
    temperature: float = 0.1,
    progress_callback: Optional[Callable[[int, int], None]] = None,
    enable_llm: bool = False,
    enable_literature: bool = True
) -> List[Dict[str, Any]]:
    """
    Generate annotations for multiple variants.

    Uses the public API with configurable LLM and literature options.

    Args:
        variants: List of dicts with 'gene', 'variant', and optional 'tumor_type'
        model: LLM model to use
        temperature: LLM temperature (0.0-1.0)
        progress_callback: Optional callback(current, total) for progress updates
        enable_llm: Enable LLM enhancement
        enable_literature: Enable literature search (disable for fast mode)

    Returns:
        List of annotation results
    """
    results = []
    total = len(variants)

    for i, v in enumerate(variants):
        if progress_callback:
            progress_callback(i + 1, total)

        try:
            result = await get_variant_annotation(
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


# Future feature placeholders

async def fetch_esmfold_structure(gene: str) -> str:
    """
    Fetch protein structure from ESMFold API.

    TODO: Implement ESMFold integration
    - Call ESMFold API with gene sequence
    - Return PDB format string for visualization
    - Cache results to avoid repeated API calls

    Returns:
        PDB format string
    """
    raise NotImplementedError("ESMFold integration coming soon")


async def predict_splice_impact(
    gene: str,
    variant: str,
    hgvs_genomic: str
) -> Dict[str, Any]:
    """
    Run SpliceAI predictions for variant.

    TODO: Implement SpliceAI integration
    - Load SpliceAI model
    - Run predictions for genomic position
    - Return acceptor/donor gain/loss scores

    Returns:
        Dict with splice scores and positions
    """
    raise NotImplementedError("SpliceAI integration coming soon")


async def run_agent_workflow(
    gene: str,
    variant: str,
    tumor_type: str
) -> Dict[str, Any]:
    """
    Execute multi-agent analysis workflow.

    TODO: Implement LangGraph/CrewAI agentic workflow
    - Evidence gathering agent (MyVariant.info)
    - Structure prediction agent (ESMFold)
    - Splice analysis agent (SpliceAI)
    - Literature search agent (PubMed)
    - Synthesis agent (LLM summary)

    Returns:
        Comprehensive analysis report
    """
    raise NotImplementedError("Agentic workflow coming soon")

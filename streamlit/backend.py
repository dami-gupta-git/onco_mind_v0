"""
Backend logic for OncoMind Streamlit app.

Integrates with the oncomind package for:
- Variant insight generation (single and batch)
- Evidence gathering from multiple databases

Supports two modes:
- Fast insight mode: Uses new public API (no LLM, returns EvidencePanel)
- Full insight mode: Uses InsightEngine (with LLM narrative)
"""

import asyncio
from typing import Any, Dict, List, Optional, Callable

from oncomind import get_insight, get_insights, AnnotationConfig


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
    Generate annotation for a single variant using the new public API.

    This is the recommended function for variant annotation.
    Returns EvidencePanel data.

    Args:
        gene: Gene symbol (e.g., BRAF)
        variant: Variant notation (e.g., V600E)
        tumor_type: Optional tumor type (e.g., Melanoma)
        enable_llm: Enable LLM enhancement for literature analysis
        enable_literature: Enable literature search (disable for fast mode)
        model: LLM model to use (if enable_llm=True)
        temperature: LLM temperature (0.0-1.0)

    Returns:
        Dict containing annotation results with identifiers, evidence, etc.
    """
    try:
        config = AnnotationConfig(
            enable_llm=enable_llm,
            enable_literature=enable_literature,
            llm_model=model,
            llm_temperature=temperature,
        )

        # Build variant string
        variant_str = f"{gene} {variant}"

        # Get evidence panel
        panel = await get_insight(variant_str, tumor_type=tumor_type, config=config)

        # Convert to dict for JSON serialization
        return {
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
            "recommended_therapies": _extract_therapies(panel),
            "evidence_panel": panel.model_dump(mode="json"),
        }

    except Exception as e:
        return {"error": f"Annotation generation failed: {str(e)}"}


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


async def get_variant_insight(
    gene: str,
    variant: str,
    tumor_type: Optional[str] = None,
    model: str = "gpt-4o-mini",
    temperature: float = 0.1
) -> Dict[str, Any]:
    """
    Generate annotation insight for a single variant with LLM narrative.

    This uses the legacy InsightEngine for full LLM narrative generation.
    For faster annotation without LLM, use get_variant_annotation instead.

    Args:
        gene: Gene symbol (e.g., BRAF)
        variant: Variant notation (e.g., V600E)
        tumor_type: Optional tumor type (e.g., Melanoma)
        model: LLM model to use
        temperature: LLM temperature (0.0-1.0)

    Returns:
        Dict containing annotation results with identifiers, evidence, etc.
    """
    try:
        # Import legacy engine
        from oncomind.engine import InsightEngine
        from oncomind.models.variant import VariantInput

        # Create insight engine
        engine = InsightEngine(llm_model=model, llm_temperature=temperature)

        # Create variant input
        variant_input = VariantInput(
            gene=gene,
            variant=variant,
            tumor_type=tumor_type
        )

        # Generate insight
        async with engine:
            insight = await engine.get_insight(variant_input)

        # Convert to dict for JSON serialization
        return {
            "variant": {
                "gene": insight.gene,
                "variant": insight.variant,
                "tumor_type": insight.tumor_type,
            },
            "insight": {
                "summary": insight.summary,
                "rationale": insight.rationale,
                "evidence_strength": insight.evidence_strength,
            },
            "identifiers": {
                "cosmic_id": insight.cosmic_id,
                "ncbi_gene_id": insight.ncbi_gene_id,
                "dbsnp_id": insight.dbsnp_id,
                "clinvar_id": insight.clinvar_id,
            },
            "hgvs": {
                "genomic": insight.hgvs_genomic,
                "protein": insight.hgvs_protein,
                "transcript": insight.hgvs_transcript,
            },
            "clinvar": {
                "clinical_significance": insight.clinvar_clinical_significance,
                "accession": insight.clinvar_accession,
            },
            "annotations": {
                "snpeff_effect": insight.snpeff_effect,
                "polyphen2_prediction": insight.polyphen2_prediction,
                "cadd_score": insight.cadd_score,
                "gnomad_exome_af": insight.gnomad_exome_af,
                "alphamissense_score": insight.alphamissense_score,
                "alphamissense_prediction": insight.alphamissense_prediction,
            },
            "transcript": {
                "id": insight.transcript_id,
                "consequence": insight.transcript_consequence,
            },
            "recommended_therapies": [
                {
                    "drug_name": therapy.drug_name,
                    "evidence_level": therapy.evidence_level,
                    "approval_status": therapy.approval_status,
                    "clinical_context": therapy.clinical_context,
                }
                for therapy in insight.recommended_therapies
            ],
        }

    except Exception as e:
        return {"error": f"Annotation generation failed: {str(e)}"}


async def batch_get_variant_annotations(
    variants: List[Dict[str, str]],
    model: str = "gpt-4o-mini",
    temperature: float = 0.1,
    enable_llm: bool = False,
    progress_callback: Optional[Callable[[int, int], None]] = None
) -> List[Dict[str, Any]]:
    """
    Generate annotations for multiple variants using the new public API.

    This is the recommended function for batch annotation.

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
        config = AnnotationConfig(
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
        return [{"error": f"Batch annotation generation failed: {str(e)}"}]


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

    Uses the new public API with configurable LLM and literature options.

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

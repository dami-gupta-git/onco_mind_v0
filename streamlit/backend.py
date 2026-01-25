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

from datetime import datetime
from zoneinfo import ZoneInfo
from typing import Any, Dict, List, Optional, Callable

import streamlit as st

from oncomind.insight_builder import Conductor, ConductorConfig
from oncomind.config.debug import get_logger
from oncomind.config.constants import (
    LLM_DEFAULT_MODEL,
    LLM_DEFAULT_TEMPERATURE,
    is_biomarker_selection_drug,
)
from oncomind.models.evidence.base import tumor_types_match, extract_variant_position
from oncomind.utils import dedupe_civic_evidence, dedupe_vicc_evidence

logger = get_logger(__name__)


def _variant_matches_indication(queried_variant: str, indication) -> bool:
    """Check if queried variant matches the FDA biomarker indication.

    Returns True if:
    - Specificity is 'gene' level (any mutation in gene is approved)
    - Queried variant exactly matches one of the specified_variants
    - Queried variant is at the same codon (position) as the indication

    Args:
        queried_variant: The variant being queried (e.g., "L858R")
        indication: FDABiomarkerEvidence object

    Returns:
        True if the variant matches this indication
    """
    specificity = indication.specificity.value if indication.specificity else None

    # Gene-level approval covers all variants
    if specificity == "gene":
        return True

    # Pathway-level (e.g., HRR-deficient) - include for now
    if specificity == "pathway":
        return True

    # Check for exact variant match
    specified = [v.upper() for v in (indication.specified_variants or [])]
    if queried_variant in specified:
        return True

    # Check for codon-level match (same position)
    if specificity in ("variant", "codon"):
        queried_pos = extract_variant_position(queried_variant)
        if queried_pos and indication.codon:
            if queried_pos == indication.codon:
                return True
        # Also check positions of specified variants
        for spec_var in specified:
            spec_pos = extract_variant_position(spec_var)
            if queried_pos and spec_pos and queried_pos == spec_pos:
                return True

    return False


async def get_variant_insight(
    gene: str,
    variant: str,
    tumor_type: Optional[str] = None,
    enable_llm: bool = False,
    enable_literature: bool = False,
    literature_source: str = "pubmed",
    model: str = LLM_DEFAULT_MODEL,
    temperature: float = LLM_DEFAULT_TEMPERATURE,
    enable_timing: bool = False,
) -> Dict[str, Any]:
    """
    Generate insight for a single variant.

    Args:
        gene: Gene symbol (e.g., BRAF)
        variant: Variant notation (e.g., V600E)
        tumor_type: Optional tumor type (e.g., Melanoma)
        enable_llm: Enable LLM narrative generation
        enable_literature: Enable literature search
        literature_source: "none", "pubmed", or "semantic_scholar"
        model: LLM model to use (if enable_llm=True)
        temperature: LLM temperature (0.0-1.0)
        enable_timing: Print timing breakdown to console

    Returns:
        Dict containing insight results with identifiers, evidence, etc.
    """
    try:
        # Get current time in EST for human-readable logs
        est = ZoneInfo("America/New_York")
        timestamp = datetime.now(est).strftime("%b %d %Y %I:%M:%S %p EST")
        # Get client IP from Streamlit headers (X-Forwarded-For for proxied requests, or Host for local)
        client_ip = "local"
        try:
            headers = st.context.headers
            # Check for proxy headers first (deployed environment)
            forwarded_for = headers.get("X-Forwarded-For") or headers.get("X-Real-Ip")
            if forwarded_for:
                # Take first IP if multiple (original client)
                client_ip = forwarded_for.split(",")[0].strip()
            else:
                # Local development - extract host without port
                host = headers.get("Host", "")
                if host and not host.startswith("localhost"):
                    client_ip = host.split(":")[0]
        except Exception:
            pass
        print(f"******** [{timestamp}] [{client_ip}] Retrieving data for {gene} {variant} (tumor: {tumor_type or 'not specified'})...", flush=True)

        logger.debug(f"get_variant_insight: {gene} {variant} (tumor={tumor_type})")
        logger.debug(f"  enable_llm={enable_llm}, enable_literature={enable_literature}, model={model}")

        # Configure and run the Conductor
        config = ConductorConfig(
            enable_literature=enable_literature,
            literature_source=literature_source,
            enable_llm=enable_llm,
            llm_model=model,
            llm_temperature=temperature,
            enable_timing=enable_timing,
        )
        async with Conductor(config) as conductor:
            result = await conductor.run(f"{gene} {variant}", tumor_type=tumor_type)

        logger.debug(f"  Result: {len(result.evidence.fda_biomarker_evidence)} FDA, "
                    f"{len(result.evidence.civic_assertions)} CIViC")

        # Build response
        return _build_response(result)

    except Exception as e:
        logger.error(f"Insight generation failed for {gene} {variant}: {e}")
        return {"error": f"Insight generation failed: {str(e)}"}


async def batch_get_variant_insights(
    variants: List[Dict[str, str]],
    enable_llm: bool = False,
    enable_literature: bool = False,
    literature_source: str = "pubmed",
    model: str = LLM_DEFAULT_MODEL,
    temperature: float = LLM_DEFAULT_TEMPERATURE,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> List[Dict[str, Any]]:
    """
    Generate insights for multiple variants.

    Args:
        variants: List of dicts with 'gene', 'variant', and optional 'tumor_type'
        enable_llm: Enable LLM narrative generation
        enable_literature: Enable literature search
        literature_source: "none", "pubmed", or "semantic_scholar"
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
                literature_source=literature_source,
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


def _build_response(result) -> Dict[str, Any]:
    """Build the standard response dict from a Result object.

    Uses the flat Evidence structure where each source is a simple list:
    - evidence.fda_biomarker_evidence, .civic_assertions, .civic_evidence, etc.
    """
    llm = result.llm  # May be None if LLM was not enabled
    evidence = result.evidence  # Evidence object with flat list structure
    queried_gene = result.identifiers.gene.upper() if result.identifiers.gene else ""
    queried_variant = result.identifiers.variant.upper() if result.identifiers.variant else ""
    queried_tumor = result.context.tumor_type

    # Get therapeutic evidence - always use evidence.get_therapeutic_evidence() as primary source
    # LLM therapeutic_evidence is typically empty since the LLM doesn't populate it
    therapeutic_list = evidence.get_therapeutic_evidence()

    # Filter FDA biomarker evidence for display (already deduplicated in evidence_aggregator)
    # Use include_level="same_codon" to show related drugs at same codon (even if not matching)
    filtered_fda_evidence = evidence.get_filtered_fda_evidence(
        queried_gene,
        queried_variant,
        queried_tumor,
        include_level="same_codon",
    )

    return {
        "variant": {
            "gene": result.identifiers.gene,
            "variant": result.identifiers.variant,
            "tumor_type": result.context.tumor_type,
        },
        "insight": {
            "summary": result.get_summary(),  # Always 1-line summary
            "llm_narrative": llm.llm_summary if llm else None,  # Plain text summary (no formatting)
            "rationale": llm.rationale if llm else None,
            # Raw component fields for UI formatting
            "functional_summary": llm.functional_summary if llm else None,
            "biological_context": llm.biological_context if llm else None,
            "therapeutic_summary": llm.therapeutic_summary if llm else None,
            "therapeutic_landscape": llm.therapeutic_landscape if llm else None,
            # Research-focused fields
            "evidence_quality": llm.evidence_quality if llm else None,
            "knowledge_gaps": llm.knowledge_gaps if llm else [],
            "well_characterized": llm.well_characterized if llm else [],
            "conflicting_evidence": llm.conflicting_evidence if llm else [],
            "research_implications": llm.research_implications if llm else None,
            "references": llm.references if llm else [],
            "evidence_tags": llm.evidence_tags if llm else [],
            "research_hypotheses": llm.research_hypotheses if llm else [],
        },
        # Cross-source analysis (disabled)
        # "cross_source_analysis": result.cross_source_analysis,
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
        # FDA biomarker evidence (pre-computed above for consistency with gap analysis)
        # Include fda_label_url property since model_dump() doesn't include @property fields
        "fda_biomarker_evidence": [{**e.model_dump(), "fda_label_url": e.fda_label_url} for e in filtered_fda_evidence],
        "civic_assertions": [
            {
                "id": a.assertion_id,
                "aid": a.aid,  # Formatted ID (e.g., "AID20")
                "civic_url": a.civic_url,  # Direct link to CIViC
                "therapies": a.therapies,
                "disease": a.disease,
                "significance": a.significance,
                "amp_level": a.amp_level,
                "amp_tier": a.amp_tier,
                "description": a.description,
                # Match specificity tracking
                "locus_match": a.locus_match,
                "matched_profile": a.matched_profile,
                "tumor_match": a.tumor_match,
            }
            for a in evidence.civic_assertions
        ],
        "civic_evidence": dedupe_civic_evidence(evidence.civic_evidence),
        # Use get_vicc_unique() to exclude CIViC/CGI sources (avoid double-counting)
        # Then deduplicate by (drugs, disease, response_type) to avoid repeated entries
        "vicc_evidence": dedupe_vicc_evidence(evidence.get_vicc_unique()),
        "cgi_biomarkers": [
            {
                "drug": b.drug,
                "association": b.association,
                "tumor_type": b.tumor_type,
                "evidence_level": b.evidence_level,
                "fda_approved": b.fda_approved,
                "fda_url": b.fda_url,
                # Match specificity tracking
                "locus_match": b.locus_match,
                "matched_alteration": b.matched_alteration,
                "tumor_match": b.tumor_match,
            }
            for b in evidence.cgi_biomarkers
        ],
        "clinvar_entries": [
            {
                "variation_id": c.variation_id,
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
                "variant_specific": t.locus_match == "variant",
                "matched_biomarker": t.matched_biomarker,
                "match_scope": t.match_scope,
                "tumor_match": t.tumor_match,
                "locus_match": t.locus_match,
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
                "url": a.url,
                "locus_match": a.locus_match,
                "cancer_specificity": a.cancer_specificity,
            }
            for a in evidence.pubmed_articles
        ],
        "literature_knowledge": evidence.literature_knowledge.summary if evidence.literature_knowledge else None,
        "preclinical_biomarkers": [
            {
                "drug": b.drug,
                "association": b.association,
                "evidence_level": b.evidence_level,
                "tumor_type": b.tumor_type,
                "locus_match": b.locus_match,
            }
            for b in evidence.preclinical_biomarkers
        ],
        "early_phase_biomarkers": [
            {
                "drug": b.drug,
                "association": b.association,
                "evidence_level": b.evidence_level,
                "tumor_type": b.tumor_type,
                "locus_match": b.locus_match,
            }
            for b in evidence.early_phase_biomarkers
        ],
        "cbioportal_evidence": {
            "gene": evidence.cbioportal_evidence.gene,
            "variant": evidence.cbioportal_evidence.variant,
            "tumor_type": evidence.cbioportal_evidence.tumor_type,
            "study_id": evidence.cbioportal_evidence.study_id,
            "study_name": evidence.cbioportal_evidence.study_name,
            "total_samples": evidence.cbioportal_evidence.total_samples,
            "samples_with_gene_mutation": evidence.cbioportal_evidence.samples_with_gene_mutation,
            "samples_with_exact_variant": evidence.cbioportal_evidence.samples_with_exact_variant,
            "gene_prevalence_pct": evidence.cbioportal_evidence.gene_prevalence_pct,
            "variant_prevalence_pct": evidence.cbioportal_evidence.variant_prevalence_pct,
            "co_occurring": [c.model_dump() for c in evidence.cbioportal_evidence.co_occurring],
            "mutually_exclusive": [m.model_dump() for m in evidence.cbioportal_evidence.mutually_exclusive],
        } if evidence.cbioportal_evidence else None,
        "depmap_evidence": {
            "gene": evidence.depmap_evidence.gene,
            "variant": evidence.depmap_evidence.variant,
            "gene_dependency": {
                "gene": evidence.depmap_evidence.gene_dependency.gene,
                "mean_dependency_score": evidence.depmap_evidence.gene_dependency.mean_dependency_score,
                "n_dependent_lines": evidence.depmap_evidence.gene_dependency.n_dependent_lines,
                "n_total_lines": evidence.depmap_evidence.gene_dependency.n_total_lines,
                "dependency_pct": evidence.depmap_evidence.gene_dependency.dependency_pct,
            } if evidence.depmap_evidence.gene_dependency else None,
            "drug_sensitivities": [
                {
                    "drug_name": ds.drug_name,
                    "mean_log2fc": ds.mean_log2fc,
                    "n_cell_lines": ds.n_cell_lines,
                    "sensitive_lines": ds.sensitive_lines,
                }
                for ds in evidence.depmap_evidence.drug_sensitivities
            ],
            "cell_line_models": [
                {
                    "name": cl.name,
                    "depmap_id": cl.depmap_id,
                    "primary_disease": cl.primary_disease,
                    "subtype": cl.subtype,
                    "has_mutation": cl.has_mutation,
                    "mutation_details": cl.mutation_details,
                }
                for cl in evidence.depmap_evidence.cell_line_models
            ],
            "is_essential": evidence.depmap_evidence.is_essential(),
            "n_cell_lines_screened": evidence.depmap_evidence.n_cell_lines_screened,
        } if evidence.depmap_evidence and evidence.depmap_evidence.has_data() else None,
        "hotspots_evidence": evidence.hotspots_evidence.to_dict() if evidence.hotspots_evidence and (
            evidence.hotspots_evidence.has_data() or evidence.hotspots_evidence.is_adjacent_to_hotspot()
        ) else None,
        "recommended_therapies": [
            {
                "drug_name": t.drug_name,
                "evidence_level": t.evidence_level,
                "approval_status": t.approval_status,
                "clinical_context": t.clinical_context,
                # New research-focused fields
                "response_type": t.response_type,
                "mechanism": t.mechanism,
                "tumor_types_tested": t.tumor_types_tested,
                "source": t.source,
                "source_url": t.source_url,
                "confidence": t.confidence,
                # Match specificity tracking
                "locus_match": t.locus_match,
                # Cancer type specificity
                "cancer_specificity": t.cancer_specificity,
            }
            for t in therapeutic_list
        ],
        # Structured evidence gaps (from gap_detector, not LLM)
        "evidence_gaps": _build_evidence_gaps(evidence),
        "result_data": result.model_dump(mode="json"),
    }


def _build_evidence_gaps(evidence) -> dict:
    """Build structured evidence gaps dict from Evidence model.

    This is the deterministic gap analysis, not LLM-generated.
    """
    # Compute gaps if not already computed
    gaps = evidence.evidence_gaps
    if gaps is None:
        gaps = evidence.compute_evidence_gaps()

    return {
        "overall_quality": gaps.overall_evidence_quality,
        "research_priority": gaps.research_priority,
        "well_characterized": gaps.well_characterized,
        "well_characterized_detailed": [
            {
                "category": wc.category.value if wc.category else None,
                "aspect": wc.aspect,
                "basis": wc.basis,
                "matches_on": wc.matches_on,
                "tumor_match": wc.tumor_match,
            }
            for wc in gaps.well_characterized_detailed
        ],
        "poorly_characterized": gaps.poorly_characterized,
        "gaps": [
            {
                "category": g.category.value,
                "severity": g.severity.value,
                "description": g.description,
                "suggested_studies": g.suggested_studies,
                "addressable_with": g.addressable_with,
            }
            for g in gaps.gaps
        ],
        "has_critical_gaps": gaps.has_critical_gaps(),
    }

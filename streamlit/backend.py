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
from typing import Any, Dict, List, Optional, Callable

import streamlit as st

from oncomind.insight_builder import Conductor, ConductorConfig
from oncomind.config.debug import get_logger
from oncomind.config.constants import (
    LLM_DEFAULT_MODEL,
    LLM_DEFAULT_TEMPERATURE,
    is_biomarker_selection_drug,
    FDA_KIT_EXCLUSION_PATTERNS,
)

logger = get_logger(__name__)


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
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # Get client IP from Streamlit headers (X-Forwarded-For for proxied requests, or direct)
        client_ip = "unknown"
        try:
            headers = st.context.headers
            # Debug: print all available headers
            print(f"DEBUG headers: {dict(headers)}")
            client_ip = headers.get("X-Forwarded-For", headers.get("X-Real-Ip", headers.get("Host", "unknown")))
            if client_ip and "," in client_ip:
                client_ip = client_ip.split(",")[0].strip()
        except Exception as e:
            print(f"DEBUG: Exception getting headers: {e}")
        print(f"[{timestamp}] [{client_ip}] Retrieving data for {gene} {variant} (tumor: {tumor_type or 'not specified'})...")
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

        logger.debug(f"  Result: {len(result.evidence.fda_approvals)} FDA, "
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


# === Private helper functions ===


def _sort_fda_labels(labels: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Sort FDA labels by match level (variant > codon > gene).

    Args:
        labels: List of FDA label dicts with biomarker_match info

    Returns:
        Sorted list with variant-level matches first
    """
    def sort_key(label: Dict[str, Any]) -> tuple:
        # Get match level from biomarker_match
        biomarker_match = label.get("biomarker_match") or {}
        match_level = biomarker_match.get("match_level", "gene")

        # Priority: variant (0) > codon (1) > gene (2) > unknown (3)
        level_priority = {
            "variant": 0,
            "codon": 1,
            "gene": 2,
        }
        priority = level_priority.get(match_level, 3)

        # Secondary sort by drug name for stability
        drug_name = (label.get("drug") or "").lower()

        return (priority, drug_name)

    return sorted(labels, key=sort_key)


def _is_kit_false_positive(indication_text: str | None) -> bool:
    """Check if an FDA label is a KIT false positive (diagnostic kit, not KIT gene).

    Args:
        indication_text: The indications_and_usage text from FDA label

    Returns:
        True if this appears to be a diagnostic/preparation kit, not a KIT oncogene drug
    """
    if not indication_text:
        return False

    text_lower = indication_text.lower()

    # Check for KIT exclusion patterns (diagnostic kit, test kit, etc.)
    has_exclusion = any(pattern in text_lower for pattern in FDA_KIT_EXCLUSION_PATTERNS)

    # Check for oncology context
    oncology_terms = ["cancer", "tumor", "malignant", "neoplasm", "carcinoma", "leukemia",
                      "lymphoma", "melanoma", "sarcoma", "gist", "mastocytosis"]
    has_oncology = any(term in text_lower for term in oncology_terms)

    # If it has exclusion patterns and no oncology context, it's a false positive
    return has_exclusion and not has_oncology


def _dedupe_civic_evidence(civic_evidence_list) -> List[Dict[str, Any]]:
    """Deduplicate CIViC evidence items by evidence_id.

    Defensive measure to ensure no duplicate EIDs appear in the UI,
    even if the API returns them.
    """
    seen_ids = set()
    deduped = []
    for e in civic_evidence_list:
        if e.evidence_id in seen_ids:
            continue
        seen_ids.add(e.evidence_id)
        deduped.append({
            "evidence_id": e.evidence_id,
            "eid": e.eid,  # Formatted ID (e.g., "EID5586")
            "civic_url": e.civic_url,  # Direct link to CIViC
            "evidence_type": e.evidence_type,
            "evidence_level": e.evidence_level,
            "clinical_significance": e.clinical_significance,
            "disease": e.disease,
            "drugs": e.drugs,
            "description": e.description,
            "pmid": e.pmid,
            "source_url": e.source_url,
            "trust_rating": e.trust_rating or e.rating,  # Use trust_rating if available, else rating
            "evidence_direction": e.evidence_direction,
            # Match specificity tracking
            "locus_match": e.locus_match,
            "matched_profile": e.matched_profile,
            "tumor_match": e.tumor_match,
        })
    return deduped


def _build_response(result) -> Dict[str, Any]:
    """Build the standard response dict from a Result object.

    Uses the flat Evidence structure where each source is a simple list:
    - evidence.fda_approvals, .civic_assertions, .civic_evidence, etc.
    """
    llm = result.llm  # May be None if LLM was not enabled
    evidence = result.evidence  # Evidence object with flat list structure
    queried_gene = result.identifiers.gene.upper() if result.identifiers.gene else ""

    # Get therapeutic evidence - always use evidence.get_therapeutic_evidence() as primary source
    # LLM therapeutic_evidence is typically empty since the LLM doesn't populate it
    therapeutic_list = evidence.get_therapeutic_evidence()

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
        # Cross-source analysis (separate LLM call)
        "cross_source_analysis": result.cross_source_analysis,
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
                "companion_diagnostic": a.companion_diagnostic,
                "black_box_warning": a.black_box_warning,
                "dosing_for_variant": a.dosing_for_variant,
                # Match specificity tracking
                "locus_match": a.locus_match,
            }
            for a in evidence.fda_approvals
        ],
        "fda_labels": _sort_fda_labels([
            {
                "drug": l.drug,
                "gene": l.gene,
                "brand_name": l.brand_name,
                "generic_name": l.generic_name,
                "manufacturer": l.manufacturer,
                "indications_and_usage": l.indications_and_usage,
                "clinical_studies": {
                    "trial_name": l.clinical_studies.trial_name,
                    "nct_id": l.clinical_studies.nct_id,
                    "patients_n": l.clinical_studies.patients_n,
                    "pfs_months_treatment": l.clinical_studies.pfs_months_treatment,
                    "pfs_months_control": l.clinical_studies.pfs_months_control,
                    "hazard_ratio": l.clinical_studies.hazard_ratio,
                    "hazard_ratio_ci": l.clinical_studies.hazard_ratio_ci,
                    "orr_treatment": l.clinical_studies.orr_treatment,
                    "orr_control": l.clinical_studies.orr_control,
                    "biomarker_breakdown": l.clinical_studies.biomarker_breakdown,
                } if l.clinical_studies else None,
                "mechanism_of_action": {
                    "targets": l.mechanism_of_action.targets,
                    "mechanism": l.mechanism_of_action.mechanism,
                    "preclinical": l.mechanism_of_action.preclinical,
                } if l.mechanism_of_action else None,
                "adverse_reactions": {
                    "common_toxicities": l.adverse_reactions.common_toxicities,
                    "serious_rate": l.adverse_reactions.serious_rate,
                    "discontinuation_rate": l.adverse_reactions.discontinuation_rate,
                } if l.adverse_reactions else None,
                "effective_time": l.effective_time,
                "approved_indications": l.approved_indications or [],
                "last_label_update": l.last_label_update,
                "update_reason": l.update_reason,
                "clinical_studies_text": l.clinical_studies_text,
                "mechanism_of_action_text": l.mechanism_of_action_text,
                "adverse_reactions_text": l.adverse_reactions_text,
                # Match specificity tracking
                "locus_match": l.locus_match,
                "biomarker_match": {
                    "matched": l.biomarker_match.matched,
                    "match_level": l.biomarker_match.match_level,
                    "tumor_matched": l.biomarker_match.tumor_matched,
                    "tumor_match_type": l.biomarker_match.tumor_match_type,
                    "combination_partners": l.biomarker_match.combination_partners,
                } if l.biomarker_match else None,
            }
            for l in evidence.fda_labels
            # Filter: gene must match AND tumor must match (or be pan-cancer)
            # Also filter out biomarker selection drugs (e.g., datopotamab targets TROP2, not EGFR)
            # And filter out KIT false positives (diagnostic kits, not KIT oncogene)
            if (l.gene and l.gene.upper() == queried_gene
                and l.biomarker_match and l.biomarker_match.tumor_matched
                and not is_biomarker_selection_drug(l.drug or l.generic_name or l.brand_name or "", queried_gene)
                and not (queried_gene == "KIT" and _is_kit_false_positive(l.indications_and_usage)))
        ]),
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
        "civic_evidence": _dedupe_civic_evidence(evidence.civic_evidence),
        # Use get_vicc_unique() to exclude CIViC/CGI sources (avoid double-counting)
        "vicc_evidence": [
            {
                "source": v.source,
                "drugs": v.drugs,
                "disease": v.disease,
                "response_type": v.response_type,
                "evidence_level": v.evidence_level,
                "molecular_profile": v.molecular_profile,
                "molecular_profile_score": v.molecular_profile_score,
                "publication_url": v.publication_url[0] if isinstance(v.publication_url, list) and v.publication_url else v.publication_url,
                # Match specificity tracking
                "locus_match": v.locus_match,
                "matched_profile": v.matched_profile,
                "tumor_match": v.tumor_match,
            }
            for v in evidence.get_vicc_unique()
        ],
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

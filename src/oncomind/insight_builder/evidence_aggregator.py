"""Evidence builder for aggregating variant evidence into Evidence.

This module provides the core evidence aggregation logic as a reusable,
LLM-independent component.

ARCHITECTURE:
    ParsedVariant → build_evidence() → Evidence

    The builder:
    1. Fetches evidence from all API clients in parallel
    2. Handles API failures gracefully (logs warnings, continues)
    3. Assembles results into strongly-typed Evidence
    4. Does NOT call LLM services (that's a separate layer)

Key Design:
- Async context manager for HTTP session lifecycle
- Parallel fetching via asyncio.gather
- Graceful exception handling per-source
- No LLM dependency - pure evidence aggregation
"""

import asyncio
import os
import time
from dataclasses import dataclass, field
from typing import Any

from dotenv import load_dotenv

load_dotenv()  # Load .env file before accessing environment variables

from oncomind.api.myvariant import MyVariantClient
from oncomind.api.cgi import CGIClient
from oncomind.api.oncotree import OncoTreeClient
from oncomind.api.vicc import VICCClient
from oncomind.api.civic import CIViCClient
from oncomind.api.cbioportal import CBioPortalClient
from oncomind.api.depmap import DepMapClient
from oncomind.api.hotspots import HotspotsClient
from oncomind.api.clinicaltrials import ClinicalTrialsClient, ClinicalTrialsRateLimitError
from oncomind.api.pubmed import PubMedClient, PubMedRateLimitError
from oncomind.api.semantic_scholar import SemanticScholarClient, SemanticScholarRateLimitError
from oncomind.api.fda_drugs import get_fda_labels_for_biomarker
from oncomind.api.fda_label_parser import FDALabelParser

from oncomind.models.evidence import (
    Evidence,
    VariantIdentifiers,
    FunctionalScores,
    VariantContext,
    CBioPortalEvidence,
    CoMutationEntry,
    CGIBiomarkerEvidence,
    CIViCAssertionEvidence,
    CIViCEvidence,
    ClinicalTrialEvidence,
    FDABiomarkerEvidence,
    HotspotsEvidence,
    PubMedEvidence,
    VICCEvidence,
)
from oncomind.models.evidence.depmap import DepMapEvidence, CellLineModel
from oncomind.models.evidence.base import EvidenceLevel, is_pan_cancer_term

from oncomind.normalization import ParsedVariant
from oncomind.models.gene_context import get_gene_context, is_variant_not_actionable
from oncomind.config.debug import get_logger
from oncomind.config.constants import (
    MAX_CIVIC_ASSERTIONS,
    MAX_CLINICAL_TRIALS,
    MAX_LITERATURE_RESULTS,
    MAX_VICC_RESULTS,
)

from oncomind.utils import normalize_drug_name

logger = get_logger(__name__)


# =============================================================================
# FDA EVIDENCE DEDUPLICATION
# =============================================================================


def _dedupe_fda_biomarker_evidence(
    evidence_list: list[FDABiomarkerEvidence],
) -> list[FDABiomarkerEvidence]:
    """Deduplicate FDA biomarker evidence by drug + tumor type + combination partners.

    Multiple entries for the same drug with the same tumor types (e.g., different
    salt formulations like ERLOTINIB vs ERLOTINIB HYDROCHLORIDE) are deduplicated,
    keeping the one with the most specific match (exact > codon > gene).

    Entries for the same drug but DIFFERENT tumor types are kept separate because
    each has its own tumor_match flag that tracks whether it matches the queried tumor.

    Different combination regimens (e.g., encorafenib + cetuximab vs
    encorafenib + binimetinib) are considered separate and are NOT deduplicated.

    Args:
        evidence_list: List of FDABiomarkerEvidence objects

    Returns:
        Deduplicated list with one entry per unique (drug, tumor_types, combination)
    """
    # Group by (normalized_drug_name + tumor_types + sorted combination partners)
    best_by_key: dict[tuple, FDABiomarkerEvidence] = {}

    # Priority: exact > codon > gene (lower is better)
    match_priority = {
        "exact": 0,
        "codon": 1,
        "gene": 2,
    }

    for ev in evidence_list:
        # Build drug key: normalized drug name + sorted combination partners
        all_drugs = [normalize_drug_name(ev.drug_name or "")]
        if ev.combination_partners:
            all_drugs.extend([normalize_drug_name(p) for p in ev.combination_partners])
        drug_key = " + ".join(sorted(all_drugs))

        # Include tumor types in the key to keep separate entries per tumor
        tumor_key = tuple(sorted(t.lower() for t in ev.tumor_types)) if ev.tumor_types else ()

        dedup_key = (drug_key, tumor_key)

        match_type = getattr(ev, "variant_match_result", None) or "gene"
        current_priority = match_priority.get(match_type, 3)

        if dedup_key not in best_by_key:
            best_by_key[dedup_key] = ev
        else:
            existing = best_by_key[dedup_key]
            existing_match = getattr(existing, "variant_match_result", None) or "gene"
            existing_priority = match_priority.get(existing_match, 3)

            # Keep the more specific match
            if current_priority < existing_priority:
                best_by_key[dedup_key] = ev

    return list(best_by_key.values())


# =============================================================================
# CONFIGURATION
# =============================================================================

# Literature source options
LITERATURE_SOURCE_NONE = "none"
LITERATURE_SOURCE_PUBMED = "pubmed"
LITERATURE_SOURCE_SEMANTIC_SCHOLAR = "semantic_scholar"


@dataclass
class EvidenceAggregatorConfig:
    """Configuration for the insight builder.

    Controls which data sources to query and resource limits.
    """

    # Source toggles
    enable_vicc: bool = True
    enable_civic_assertions: bool = True
    enable_clinical_trials: bool = True
    enable_literature: bool = True

    # Literature source: "none", "pubmed", or "semantic_scholar"
    literature_source: str = LITERATURE_SOURCE_PUBMED

    # Semantic Scholar: filter to recent papers (last N years, 0 = no filter)
    # Default is 5 years to get recent, relevant literature
    semantic_scholar_recent_years: int = 5

    # Result limits (from constants.py)
    max_vicc_results: int = MAX_VICC_RESULTS
    max_civic_assertions: int = MAX_CIVIC_ASSERTIONS
    max_clinical_trials: int = MAX_CLINICAL_TRIALS
    max_literature_results: int = MAX_LITERATURE_RESULTS

    # Concurrency control for rate-limited APIs
    literature_concurrency: int = 1

    # API keys (from environment by default)
    semantic_scholar_api_key: str | None = field(
        default_factory=lambda: os.environ.get("SEMANTIC_SCHOLAR_API_KEY")
    )


# =============================================================================
# RESULT TRACKING
# =============================================================================

@dataclass
class FetchResults:
    """Container for tracking fetch results and source status."""
    sources_queried: list[str] = field(default_factory=list)
    sources_with_data: list[str] = field(default_factory=list)
    sources_failed: list[str] = field(default_factory=list)

    def handle_result(self, result: Any, source_name: str) -> Any:
        """Handle an API result with error tracking.

        Returns:
            The result if successful, empty list/None if failed
        """
        if isinstance(result, Exception):
            logger.warning(f"{source_name} fetch failed: {str(result)}")
            self.sources_failed.append(source_name)
            return [] if source_name not in ("cBioPortal", "DepMap", "MyVariant") else None
        elif result:
            self.sources_with_data.append(source_name)
            return result
        return [] if source_name not in ("cBioPortal", "DepMap", "MyVariant") else None


# =============================================================================
# EVIDENCE AGGREGATOR
# =============================================================================

class EvidenceAggregator:
    """Builder for aggregating variant evidence from multiple sources.

    Use as an async context manager to ensure proper HTTP session lifecycle:

        async with EvidenceAggregator() as aggregator:
            evidence = await aggregator.build_evidence(parsed_variant, tumor_type)
    """

    def __init__(self, config: EvidenceAggregatorConfig | None = None):
        self.config = config or EvidenceAggregatorConfig()
        self._init_clients()
        self._literature_semaphore = asyncio.Semaphore(self.config.literature_concurrency)

    def _init_clients(self) -> None:
        """Initialize all API clients based on configuration."""
        # Core clients (always enabled)
        self.myvariant_client = MyVariantClient()
        self.cgi_client = CGIClient()
        self.oncotree_client = OncoTreeClient()
        self.cbioportal_client = CBioPortalClient()
        self.depmap_client = DepMapClient()
        self.hotspots_client = HotspotsClient()  # File-based, no async needed

        # Optional clients based on config
        self.vicc_client = VICCClient() if self.config.enable_vicc else None
        self.civic_client = CIViCClient() if self.config.enable_civic_assertions else None
        self.clinical_trials_client = (
            ClinicalTrialsClient() if self.config.enable_clinical_trials else None
        )

        # Literature clients
        if self.config.enable_literature:
            self.semantic_scholar_client = SemanticScholarClient(
                api_key=self.config.semantic_scholar_api_key
            )
            self.pubmed_client = PubMedClient()
        else:
            self.semantic_scholar_client = None
            self.pubmed_client = None

    def _get_all_clients(self) -> list[Any]:
        """Get list of all active clients for context management."""
        # Note: cgi_client is synchronous, not an async context manager
        clients = [
            self.myvariant_client,
            self.oncotree_client,
            self.cbioportal_client,
            self.depmap_client,
        ]
        optional = [
            self.vicc_client,
            self.civic_client,
            self.clinical_trials_client,
            self.semantic_scholar_client,
            self.pubmed_client,
        ]
        return clients + [c for c in optional if c is not None]

    async def __aenter__(self):
        """Initialize HTTP client sessions."""
        for client in self._get_all_clients():
            await client.__aenter__()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Close HTTP client sessions."""
        for client in self._get_all_clients():
            await client.__aexit__(exc_type, exc_val, exc_tb)

    # -------------------------------------------------------------------------
    # Tumor Type Resolution
    # -------------------------------------------------------------------------

    async def _resolve_tumor_type(self, tumor_type: str | None) -> str | None:
        """Resolve tumor type using OncoTree."""
        if not tumor_type:
            return None

        try:
            resolved = await self.oncotree_client.resolve_tumor_type(tumor_type)
            if resolved != tumor_type:
                logger.debug(f"Resolved tumor type: {tumor_type} → {resolved}")
            return resolved
        except Exception as e:
            logger.warning(f"OncoTree resolution failed: {str(e)}")
            return tumor_type

    # -------------------------------------------------------------------------
    # Literature Fetching
    # -------------------------------------------------------------------------

    async def _fetch_literature(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
    ) -> list[PubMedEvidence]:
        """Fetch literature based on configured source.

        Sources:
        - pubmed: Fast PubMed E-utilities with relevance ranking
        - semantic_scholar: Richer metadata (citations, TLDR) but slower
        - none: Skip literature fetching
        """
        source = self.config.literature_source

        if source == LITERATURE_SOURCE_NONE:
            return []

        async with self._literature_semaphore:
            if source == LITERATURE_SOURCE_SEMANTIC_SCHOLAR:
                return await self._fetch_semantic_scholar_literature(gene, variant, tumor_type)
            else:  # Default to PubMed
                return await self._fetch_pubmed_literature(gene, variant, tumor_type)

    async def _fetch_pubmed_literature(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
    ) -> list[PubMedEvidence]:
        """Fetch literature from PubMed with relevance ranking."""
        if not self.pubmed_client:
            return []

        try:
            return await self.pubmed_client.search_pubmed_evidence(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                max_results=self.config.max_literature_results,
            )
        except PubMedRateLimitError:
            return []

    async def _fetch_semantic_scholar_literature(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
    ) -> list[PubMedEvidence]:
        """Fetch literature from Semantic Scholar with citation quality filter.

        Uses Semantic Scholar's minCitationCount filter to get higher-quality papers.
        Also filters to recent papers (configurable, default 5 years).
        """
        if not self.semantic_scholar_client:
            return []

        # Compute year range for recency filter
        year_range = None
        if self.config.semantic_scholar_recent_years > 0:
            from datetime import datetime
            current_year = datetime.now().year
            start_year = current_year - self.config.semantic_scholar_recent_years
            year_range = (start_year, current_year)

        try:
            # Semantic Scholar API filters by min_citations=5 for quality
            return await self.semantic_scholar_client.search_pubmed_evidence(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                max_results=self.config.max_literature_results,
                year_range=year_range,
            )

        except SemanticScholarRateLimitError:
            # Fall back to PubMed if rate limited
            return await self._fetch_pubmed_literature(gene, variant, tumor_type)

    # -------------------------------------------------------------------------
    # Result Processing Helpers
    # -------------------------------------------------------------------------

    def _process_cgi_biomarkers(
        self, all_biomarkers: list[CGIBiomarkerEvidence]
    ) -> tuple[list[CGIBiomarkerEvidence], list[CGIBiomarkerEvidence], list[CGIBiomarkerEvidence]]:
        """Split CGI biomarkers into FDA-approved, preclinical, and early-phase."""
        fda_approved = [b for b in all_biomarkers if b.fda_approved]
        preclinical = [
            b for b in all_biomarkers
            if not b.fda_approved and b.evidence_level in ("Pre-clinical", "Cell line")
        ]
        early_phase = [
            b for b in all_biomarkers
            if not b.fda_approved and b.evidence_level not in ("Pre-clinical", "Cell line", None)
        ]
        return fda_approved, preclinical, early_phase

    def _process_cbioportal_result(
        self, result: Any, tracker: FetchResults
    ) -> CBioPortalEvidence | None:
        """Convert cBioPortal result to CBioPortalEvidence."""
        if isinstance(result, Exception):
            tracker.sources_failed.append("cBioPortal")
            return None

        if not result:
            return None

        try:
            evidence = CBioPortalEvidence(
                gene=result.gene,
                variant=result.variant,
                tumor_type=result.tumor_type,
                study_id=result.study_id,
                study_name=result.study_name,
                total_samples=result.total_samples,
                samples_with_gene_mutation=result.samples_with_gene_mutation,
                samples_with_exact_variant=result.samples_with_exact_variant,
                gene_prevalence_pct=result.gene_prevalence_pct,
                variant_prevalence_pct=result.variant_prevalence_pct,
                co_occurring=[CoMutationEntry(**c) for c in result.co_occurring],
                mutually_exclusive=[CoMutationEntry(**m) for m in result.mutually_exclusive],
            )
            if evidence.has_data():
                tracker.sources_with_data.append("cBioPortal")
            return evidence
        except Exception:
            tracker.sources_failed.append("cBioPortal")
            return None

    def _process_cell_lines(
        self, result: Any, tracker: FetchResults
    ) -> list[CellLineModel]:
        """Convert cBioPortal cell line result to CellLineModel list."""
        if isinstance(result, Exception) or not result:
            return []

        models = []
        for cl in result:
            models.append(CellLineModel(
                name=cl.get("name", ""),
                ccle_name=cl.get("sample_id"),
                primary_disease=cl.get("tissue"),
                has_mutation=True,
                mutation_details=cl.get("protein_change"),
            ))

        if models:
            tracker.sources_with_data.append("CCLE")
        return models

    def _process_depmap_result(
        self,
        result: Any,
        cell_line_models: list[CellLineModel],
        gene: str,
        variant: str,
        tracker: FetchResults,
    ) -> DepMapEvidence | None:
        """Process DepMap result, merging with cBioPortal cell lines if available."""
        if isinstance(result, Exception):
            tracker.sources_failed.append("DepMap")
            # Even if DepMap failed, create evidence from cBioPortal cell lines
            if cell_line_models:
                return DepMapEvidence(gene=gene, variant=variant, cell_line_models=cell_line_models)
            return None

        if not result:
            return None

        # Merge cell lines from both sources (DepMap local CSV + cBioPortal CCLE API)
        # DepMap local data is often more complete, so we union both sets
        merged_cell_lines = self._merge_cell_lines(result.cell_line_models, cell_line_models)

        result = DepMapEvidence(
            gene=result.gene,
            variant=result.variant,
            gene_dependency=result.gene_dependency,
            co_dependencies=result.co_dependencies,
            drug_sensitivities=result.drug_sensitivities,
            cell_line_models=merged_cell_lines,
            data_version=result.data_version,
            n_cell_lines_screened=result.n_cell_lines_screened,
        )

        if result.has_data():
            tracker.sources_with_data.append("DepMap")
        return result

    def _merge_cell_lines(
        self,
        depmap_lines: list[CellLineModel],
        cbioportal_lines: list[CellLineModel],
    ) -> list[CellLineModel]:
        """Merge cell lines from DepMap and cBioPortal, deduplicating by name."""
        # Use normalized name as key for deduplication
        seen: dict[str, CellLineModel] = {}

        # Add DepMap lines first (primary source with more complete data)
        for cl in depmap_lines or []:
            key = cl.name.upper().replace("-", "").replace("_", "")
            if key not in seen:
                seen[key] = cl

        # Add cBioPortal lines (may have additional metadata)
        for cl in cbioportal_lines or []:
            key = cl.name.upper().replace("-", "").replace("_", "")
            if key not in seen:
                seen[key] = cl

        return list(seen.values())

    def _extract_myvariant_data(
        self, evidence: Any, variant: ParsedVariant, normalized_variant: str
    ) -> tuple[FunctionalScores, dict[str, Any], str | None, list, list, list]:
        """Extract structured data from MyVariant evidence."""
        identifiers_data = {
            "variant_id": f"{variant.gene}:{normalized_variant}",
            "gene": variant.gene,
            "variant": variant.variant,
            "variant_normalized": normalized_variant,
            "variant_type": variant.variant_type,
        }

        if not evidence:
            return FunctionalScores(), identifiers_data, None, [], [], []

        functional_scores = FunctionalScores(
            alphamissense_score=getattr(evidence, 'alphamissense_score', None),
            alphamissense_prediction=getattr(evidence, 'alphamissense_prediction', None),
            cadd_score=getattr(evidence, 'cadd_score', None),
            cadd_raw=getattr(evidence, 'cadd_raw', None),
            polyphen2_prediction=getattr(evidence, 'polyphen2_prediction', None),
            polyphen2_score=getattr(evidence, 'polyphen2_score', None),
            sift_prediction=getattr(evidence, 'sift_prediction', None),
            sift_score=getattr(evidence, 'sift_score', None),
            snpeff_effect=getattr(evidence, 'snpeff_effect', None),
            snpeff_impact=getattr(evidence, 'snpeff_impact', None),
            gnomad_exome_af=getattr(evidence, 'gnomad_exome_af', None),
            gnomad_genome_af=getattr(evidence, 'gnomad_genome_af', None),
        )

        identifiers_data.update({
            "cosmic_id": getattr(evidence, 'cosmic_id', None),
            "ncbi_gene_id": getattr(evidence, 'ncbi_gene_id', None),
            "dbsnp_id": getattr(evidence, 'dbsnp_rsid', None),
            "clinvar_id": getattr(evidence, 'clinvar_variation_id', None),
            "hgvs_genomic": getattr(evidence, 'hgvs_genomic', None),
            "hgvs_protein": getattr(evidence, 'hgvs_protein', None),
            "hgvs_transcript": getattr(evidence, 'hgvs_coding', None),
        })

        clinvar_significance = getattr(evidence, 'clinvar_clinical_significance', None)
        clinvar_entries = getattr(evidence, 'clinvar', []) or []
        cosmic_entries = getattr(evidence, 'cosmic', []) or []
        civic_entries = getattr(evidence, 'civic', []) or []

        return functional_scores, identifiers_data, clinvar_significance, clinvar_entries, cosmic_entries, civic_entries

    def _get_gene_context_data(self, gene: str) -> tuple[str | None, str | None, str | None]:
        """Get gene role, class, and pathway information."""
        gene_context = get_gene_context(gene)
        if not gene_context or not gene_context.role:
            return None, None, None

        gene_role = gene_context.role.value
        gene_class = gene_context.role.value

        # Check for pathway-actionable TSGs
        from oncomind.models.gene_context import get_pathway_actionable_info
        pathway_info = get_pathway_actionable_info(gene)
        pathway = pathway_info.get("pathway") if pathway_info else None

        return gene_role, gene_class, pathway

    # -------------------------------------------------------------------------
    # Main Build Method
    # -------------------------------------------------------------------------

    async def build_evidence(
        self,
        variant: ParsedVariant | str,
        tumor_type: str | None = None,
        enable_timing: bool = False,
    ) -> Evidence:
        """Build an Evidence for a single variant.

        Args:
            variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
            tumor_type: Optional tumor type for context
            enable_timing: If True, print timing breakdown to console

        Returns:
            Evidence with all aggregated evidence
        """
        total_start = time.time()
        source_timings: dict[str, float] = {}

        # Parse string input if needed
        if isinstance(variant, str):
            from oncomind.normalization import parse_variant_input
            variant = parse_variant_input(variant, tumor_type)

        gene = variant.gene
        normalized_variant = variant.variant_normalized or variant.variant
        tumor = tumor_type or variant.tumor_type

        # Initialize tracking
        tracker = FetchResults()
        tracker.sources_queried = ["MyVariant", "FDA", "CGI"]

        # Resolve tumor type
        t0 = time.time()
        resolved_tumor = await self._resolve_tumor_type(tumor)
        source_timings["OncoTree"] = time.time() - t0

        # Parallel fetch from all sources
        t0 = time.time()
        results, fetch_timings = await self._fetch_all_sources(gene, normalized_variant, resolved_tumor, tracker)
        source_timings.update(fetch_timings)
        parallel_fetch_time = time.time() - t0

        # Process results
        t0 = time.time()
        evidence = self._assemble_evidence(
            results, variant, normalized_variant, tumor, resolved_tumor, tracker
        )
        source_timings["assemble"] = time.time() - t0

        # Note: fda_labels has been removed - using fda_biomarker_evidence exclusively
        # which is populated from FDALabelParser with proper negation detection

        total_time = time.time() - total_start

        # Print timing breakdown if enabled
        if enable_timing:
            print(f"\n⏱️  Evidence aggregation timing for {gene} {normalized_variant}:")
            print(f"   Total: {total_time:.2f}s (parallel fetch: {parallel_fetch_time:.2f}s)")
            print("   Per-source breakdown:")
            # Sort by time descending to show slowest first
            sorted_timings = sorted(source_timings.items(), key=lambda x: x[1], reverse=True)
            for source, t in sorted_timings:
                if t >= 0.01:  # Only show sources that took >= 10ms
                    print(f"     {source:20s} {t:.2f}s")

        return evidence

    async def _fetch_all_sources(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
        tracker: FetchResults,
    ) -> tuple[tuple, dict[str, float]]:
        """Fetch data from all sources in parallel.

        Returns:
            Tuple of (results, timings) where timings maps source name to seconds
        """
        timings: dict[str, float] = {}

        async def timed_fetch(name: str, coro):
            """Wrap a coroutine with timing."""
            t0 = time.time()
            try:
                result = await coro
                timings[name] = time.time() - t0
                return result
            except Exception as e:
                timings[name] = time.time() - t0
                raise e

        async def fetch_vicc():
            if self.vicc_client:
                tracker.sources_queried.append("VICC")
                return await self.vicc_client.fetch_vicc_evidence(
                    gene=gene, variant=variant, tumor_type=tumor_type,
                    max_results=self.config.max_vicc_results,
                )
            return []

        async def fetch_civic_assertions():
            if self.civic_client:
                tracker.sources_queried.append("CIViC")
                return await self.civic_client.fetch_assertion_evidence(
                    gene=gene, variant=variant, tumor_type=tumor_type,
                    max_results=self.config.max_civic_assertions,
                )
            return []

        async def fetch_civic_evidence():
            if self.civic_client:
                return await self.civic_client.fetch_evidence_items(
                    gene=gene, variant=variant, tumor_type=tumor_type,
                    max_results=self.config.max_civic_assertions,
                )
            return []

        async def fetch_trials():
            if not self.clinical_trials_client:
                return []
            tracker.sources_queried.append("ClinicalTrials.gov")
            try:
                return await self.clinical_trials_client.search_trial_evidence(
                    gene=gene, variant=variant, tumor_type=tumor_type,
                    recruiting_only=True, max_results=self.config.max_clinical_trials,
                )
            except ClinicalTrialsRateLimitError:
                return []

        async def fetch_literature():
            if self.config.enable_literature:
                tracker.sources_queried.append("Literature")
                return await self._fetch_literature(gene, variant, tumor_type)
            return []

        async def fetch_cbioportal():
            if self.cbioportal_client:
                tracker.sources_queried.append("cBioPortal")
                return await self.cbioportal_client.fetch_co_mutation_data(
                    gene=gene, variant=variant, tumor_type=tumor_type,
                )
            return None

        async def fetch_cell_lines():
            if self.cbioportal_client:
                try:
                    return await self.cbioportal_client.fetch_cell_lines_with_mutation(
                        gene=gene, variant=variant,
                    )
                except Exception as e:
                    logger.warning(f"Cell line fetch failed for {gene} {variant}: {e}")
                    return []
            return []

        async def fetch_depmap():
            if self.depmap_client:
                tracker.sources_queried.append("DepMap")
                return await self.depmap_client.fetch_depmap_evidence(
                    gene=gene, variant=variant,
                )
            return None

        async def fetch_fda_biomarker_parsed():
            """Fetch and parse FDA labels using FDALabelParser.

            This uses the new parser to extract structured biomarker-indication
            associations with proper negation handling and variant specificity.
            """
            try:
                # Fetch raw labels from OpenFDA
                raw_labels = await get_fda_labels_for_biomarker(
                    gene=gene,
                    variant=variant,
                    tumor_type=tumor_type,
                    limit=500,
                )

                if not raw_labels:
                    return []

                # Parse labels into BiomarkerIndication objects
                parser = FDALabelParser()
                all_indications: list[FDABiomarkerEvidence] = []
                for label_data in raw_labels:
                    try:
                        indications = parser.parse_label(label_data)
                        # Extract spl_set_id for linking (OpenFDA field name)
                        set_id = label_data.get("openfda", {}).get("spl_set_id", [None])[0]
                        spl_version = label_data.get("spl_version")

                        # Convert to FDABiomarkerEvidence models
                        for ind in indications:
                            evidence = FDABiomarkerEvidence.from_biomarker_indication(
                                indication=ind,
                                set_id=set_id,
                                spl_version=spl_version,
                            )
                            all_indications.append(evidence)
                    except Exception as parse_err:
                        logger.warning(f"Failed to parse FDA label: {parse_err}")
                        continue

                # Match each indication against the query variant
                matched_indications: list[FDABiomarkerEvidence] = []
                for evidence in all_indications:
                    match_result = evidence.matches_variant(gene, variant)

                    # Skip indications for different genes
                    if match_result["match_type"] is None and "Different gene" in match_result.get("reason", ""):
                        continue

                    # Populate match result fields on the evidence
                    evidence.variant_match_result = match_result["match_type"]
                    evidence.variant_match_reason = match_result["reason"]

                    # Only include:
                    # 1. True matches (exact, codon, gene level approvals)
                    # 2. Exclusions (drugs contraindicated for this variant - important safety info)
                    # 3. Same codon different variant (e.g., G12C drug when querying G12D)
                    # Skip: different_variant, different_codon
                    if match_result["matches"] or match_result["match_type"] in ("excluded", "same_codon_different_variant"):
                        matched_indications.append(evidence)

                # Sort: matches first, then exclusions
                matched_indications.sort(key=lambda x: (
                    not (x.variant_match_result in ("exact", "codon", "gene")),  # Matches first
                    x.variant_match_result != "exact",  # Exact matches before partial
                    x.variant_match_result == "excluded",  # Exclusions last
                ))

                return matched_indications
            except Exception as e:
                logger.warning(f"FDA biomarker parsed fetch failed for {gene}: {e}")
                return []

        results = await asyncio.gather(
            timed_fetch("MyVariant", self.myvariant_client.fetch_evidence(gene=gene, variant=variant)),
            timed_fetch("CGI", asyncio.to_thread(self.cgi_client.fetch_biomarker_evidence, gene, variant, tumor_type)),
            timed_fetch("VICC", fetch_vicc()),
            timed_fetch("CIViC_assertions", fetch_civic_assertions()),
            timed_fetch("CIViC_evidence", fetch_civic_evidence()),
            timed_fetch("ClinicalTrials", fetch_trials()),
            timed_fetch("Literature", fetch_literature()),
            timed_fetch("cBioPortal", fetch_cbioportal()),
            timed_fetch("CellLines", fetch_cell_lines()),
            timed_fetch("DepMap", fetch_depmap()),
            timed_fetch("FDA_biomarker_parsed", fetch_fda_biomarker_parsed()),
            return_exceptions=True,
        )

        return results, timings

    def _assemble_evidence(
        self,
        results: tuple,
        variant: ParsedVariant,
        normalized_variant: str,
        tumor: str | None,
        resolved_tumor: str | None,
        tracker: FetchResults,
    ) -> Evidence:
        """Assemble Evidence from fetch results."""
        (
            myvariant_result, cgi_result, vicc_result, civic_assertions_result,
            civic_evidence_result, trials_result, literature_result, cbioportal_result,
            cell_lines_result, depmap_result, fda_biomarker_parsed_result,
        ) = results

        gene = variant.gene

        # Process standard results
        myvariant_evidence = tracker.handle_result(myvariant_result, "MyVariant")
        all_cgi: list[CGIBiomarkerEvidence] = tracker.handle_result(cgi_result, "CGI") or []
        vicc_evidence: list[VICCEvidence] = tracker.handle_result(vicc_result, "VICC") or []
        civic_assertions: list[CIViCAssertionEvidence] = tracker.handle_result(civic_assertions_result, "CIViC") or []
        civic_evidence: list[CIViCEvidence] = tracker.handle_result(civic_evidence_result, "CIViC Evidence") or []
        clinical_trials: list[ClinicalTrialEvidence] = tracker.handle_result(trials_result, "ClinicalTrials.gov") or []
        pubmed_articles: list[PubMedEvidence] = tracker.handle_result(literature_result, "Literature") or []

        # FDA biomarker evidence - structured indications from FDALabelParser
        fda_biomarker_evidence_raw: list[FDABiomarkerEvidence] = tracker.handle_result(
            fda_biomarker_parsed_result, "FDA Biomarker Parsed"
        ) or []

        # Deduplicate FDA evidence (e.g., ERLOTINIB vs ERLOTINIB HYDROCHLORIDE)
        # keeping the most specific match for each drug+gene+tumor+combination
        fda_biomarker_evidence = _dedupe_fda_biomarker_evidence(fda_biomarker_evidence_raw)
        if fda_biomarker_evidence:
            logger.info(
                f"FDA biomarker parser found {len(fda_biomarker_evidence_raw)} indications, "
                f"{len(fda_biomarker_evidence)} after deduplication for {gene}"
            )

        # Process CGI biomarkers (split by evidence level)
        cgi_biomarkers, preclinical_biomarkers, early_phase_biomarkers = self._process_cgi_biomarkers(all_cgi)

        # Process complex results
        cbioportal_evidence = self._process_cbioportal_result(cbioportal_result, tracker)
        cell_line_models = self._process_cell_lines(cell_lines_result, tracker)
        depmap_evidence = self._process_depmap_result(
            depmap_result, cell_line_models, gene, normalized_variant, tracker
        )

        # Fetch hotspots data (synchronous - reads from local file)
        hotspots_evidence = self.hotspots_client.fetch_hotspot_evidence(gene, normalized_variant)
        if hotspots_evidence and hotspots_evidence.has_data():
            tracker.sources_with_data.append("CancerHotspots")

        # Extract MyVariant data (CIViC evidence now comes from CIViC client directly)
        (
            functional_scores, identifiers_data, clinvar_significance,
            clinvar_entries, cosmic_entries, _civic_entries_unused
        ) = self._extract_myvariant_data(myvariant_evidence, variant, normalized_variant)

        # Check if variant is not actionable (benign polymorphism)
        # If so, skip FDA trial matching to avoid false positives
        gnomad_af = functional_scores.gnomad_exome_af
        not_actionable, not_actionable_reason = is_variant_not_actionable(
            gene=gene,
            variant=normalized_variant,
            clinvar_significance=clinvar_significance,
            gnomad_af=gnomad_af,
        )

        if not_actionable:
            logger.debug(f"Skipping FDA matching for {gene} {normalized_variant}: {not_actionable_reason}")

        # Note: fda_approvals has been removed - using fda_biomarker_evidence exclusively

        # Get gene context
        gene_role, gene_class, pathway = self._get_gene_context_data(gene)

        # Build Evidence
        return Evidence(
            identifiers=VariantIdentifiers(**identifiers_data),
            functional=functional_scores,
            context=VariantContext(
                tumor_type=tumor,
                tumor_type_resolved=resolved_tumor,
                gene_role=gene_role,
                gene_class=gene_class,
                mutation_class=None,
                pathway=pathway,
            ),
            fda_biomarker_evidence=fda_biomarker_evidence,
            civic_assertions=civic_assertions,
            civic_evidence=civic_evidence,  # Now from CIViC GraphQL API
            vicc_evidence=vicc_evidence,
            cgi_biomarkers=cgi_biomarkers,
            clinvar_entries=clinvar_entries,
            clinvar_significance=clinvar_significance,
            cosmic_entries=cosmic_entries,
            clinical_trials=clinical_trials,
            pubmed_articles=pubmed_articles,
            literature_knowledge=None,  # Set by LLM layer later
            preclinical_biomarkers=preclinical_biomarkers,
            early_phase_biomarkers=early_phase_biomarkers,
            cbioportal_evidence=cbioportal_evidence,
            depmap_evidence=depmap_evidence,
            hotspots_evidence=hotspots_evidence,
            literature_searched=self.config.enable_literature,
        )

    async def build_evidences(
        self,
        variants: list[ParsedVariant | str],
        tumor_type: str | None = None,
    ) -> list[Evidence]:
        """Build Evidence for multiple variants in parallel."""
        tasks = [self.build_evidence(v, tumor_type) for v in variants]
        results = await asyncio.gather(*tasks, return_exceptions=True)

        evidences = []
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                logger.error(f"Failed to process variant {i}: {str(result)}")
            else:
                evidences.append(result)
        return evidences


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

async def build_evidence(
    variant: ParsedVariant | str,
    tumor_type: str | None = None,
    config: EvidenceAggregatorConfig | None = None,
) -> Evidence:
    """Convenience function to build Evidence for a single variant.

    This is the recommended entry point for single-variant annotation.

    Args:
        variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
        tumor_type: Optional tumor type for context
        config: Optional configuration for the builder

    Returns:
        Evidence with all aggregated evidence

    Example:
        >>> evidence = await build_evidence("BRAF V600E", tumor_type="Melanoma")
        >>> print(evidence.identifiers.gene)
        BRAF
    """
    async with EvidenceAggregator(config) as builder:
        return await builder.build_evidence(variant, tumor_type)


async def build_evidences(
    variants: list[ParsedVariant | str],
    tumor_type: str | None = None,
    config: EvidenceAggregatorConfig | None = None,
) -> list[Evidence]:
    """Convenience function to build Evidence for multiple variants."""
    async with EvidenceAggregator(config) as builder:
        return await builder.build_evidences(variants, tumor_type)


__all__ = [
    "EvidenceAggregator",
    "EvidenceAggregatorConfig",
    "build_evidence",
    "build_evidences",
]

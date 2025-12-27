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
from dataclasses import dataclass, field
from typing import Any

from dotenv import load_dotenv

load_dotenv()  # Load .env file before accessing environment variables

from oncomind.api.myvariant import MyVariantClient
from oncomind.api.fda import FDAClient
from oncomind.api.cgi import CGIClient
from oncomind.api.oncotree import OncoTreeClient
from oncomind.api.vicc import VICCClient
from oncomind.api.civic import CIViCClient
from oncomind.api.cbioportal import CBioPortalClient
from oncomind.api.depmap import DepMapClient
from oncomind.api.clinicaltrials import ClinicalTrialsClient, ClinicalTrialsRateLimitError
from oncomind.api.pubmed import PubMedClient, PubMedRateLimitError
from oncomind.api.semantic_scholar import SemanticScholarClient, SemanticScholarRateLimitError

from oncomind.models.evidence import (
    Evidence,
    VariantIdentifiers,
    FunctionalScores,
    VariantContext,
    CBioPortalEvidence,
    CoMutationEntry,
    CGIBiomarkerEvidence,
    CIViCAssertionEvidence,
    ClinicalTrialEvidence,
    FDAApproval,
    PubMedEvidence,
    VICCEvidence,
)
from oncomind.models.evidence.depmap import DepMapEvidence, CellLineModel

from oncomind.normalization import ParsedVariant
from oncomind.utils import normalize_variant
from oncomind.models.gene_context import get_gene_context, GeneRole


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

    # Result limits
    max_vicc_results: int = 50
    max_civic_assertions: int = 20
    max_clinical_trials: int = 10
    max_literature_results: int = 6

    # Concurrency control for rate-limited APIs
    # Semantic Scholar: 1 RPS without key, ~10 RPS with key
    # PubMed: 3 RPS without key, 10 RPS with key
    literature_concurrency: int = 1  # Limit concurrent literature requests

    # API keys (from environment by default)
    semantic_scholar_api_key: str | None = field(
        default_factory=lambda: os.environ.get("SEMANTIC_SCHOLAR_API_KEY")
    )


class EvidenceAggregator:
    """Builder for aggregating variant evidence from multiple sources.

    Use as an async context manager to ensure proper HTTP session lifecycle:

        async with EvidenceAggregator() as aggregator:
            evidence = await aggregator.build_evidence(parsed_variant, tumor_type)

    Or for batch processing:

        async with EvidenceAggregator() as aggregator:
            evidences = await aggregator.build_evidences(variants, tumor_type)
    """

    def __init__(self, config: EvidenceAggregatorConfig | None = None):
        self.config = config or EvidenceAggregatorConfig()

        # Initialize API clients
        self.myvariant_client = MyVariantClient()
        self.fda_client = FDAClient()
        self.cgi_client = CGIClient()
        self.oncotree_client = OncoTreeClient()

        # Optional clients based on config
        self.vicc_client = VICCClient() if self.config.enable_vicc else None
        self.civic_client = CIViCClient() if self.config.enable_civic_assertions else None
        self.cbioportal_client = CBioPortalClient()  # Always enabled - fast API
        self.depmap_client = DepMapClient()  # Always enabled - provides preclinical context
        self.clinical_trials_client = (
            ClinicalTrialsClient() if self.config.enable_clinical_trials else None
        )

        # Literature clients with fallback
        if self.config.enable_literature:
            self.semantic_scholar_client = SemanticScholarClient(
                api_key=self.config.semantic_scholar_api_key
            )
            self.pubmed_client = PubMedClient()
        else:
            self.semantic_scholar_client = None
            self.pubmed_client = None

        # Semaphore to limit concurrent literature API requests
        self._literature_semaphore = asyncio.Semaphore(self.config.literature_concurrency)

    async def __aenter__(self):
        """Initialize HTTP client sessions."""
        await self.myvariant_client.__aenter__()
        await self.fda_client.__aenter__()
        await self.oncotree_client.__aenter__()

        if self.vicc_client:
            await self.vicc_client.__aenter__()
        if self.civic_client:
            await self.civic_client.__aenter__()
        if self.cbioportal_client:
            await self.cbioportal_client.__aenter__()
        if self.depmap_client:
            await self.depmap_client.__aenter__()
        if self.clinical_trials_client:
            await self.clinical_trials_client.__aenter__()
        if self.semantic_scholar_client:
            await self.semantic_scholar_client.__aenter__()
        if self.pubmed_client:
            await self.pubmed_client.__aenter__()

        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Close HTTP client sessions."""
        await self.myvariant_client.__aexit__(exc_type, exc_val, exc_tb)
        await self.fda_client.__aexit__(exc_type, exc_val, exc_tb)
        await self.oncotree_client.__aexit__(exc_type, exc_val, exc_tb)

        if self.vicc_client:
            await self.vicc_client.__aexit__(exc_type, exc_val, exc_tb)
        if self.civic_client:
            await self.civic_client.__aexit__(exc_type, exc_val, exc_tb)
        if self.cbioportal_client:
            await self.cbioportal_client.__aexit__(exc_type, exc_val, exc_tb)
        if self.depmap_client:
            await self.depmap_client.__aexit__(exc_type, exc_val, exc_tb)
        if self.clinical_trials_client:
            await self.clinical_trials_client.__aexit__(exc_type, exc_val, exc_tb)
        if self.semantic_scholar_client:
            await self.semantic_scholar_client.__aexit__(exc_type, exc_val, exc_tb)
        if self.pubmed_client:
            await self.pubmed_client.__aexit__(exc_type, exc_val, exc_tb)

    async def _resolve_tumor_type(self, tumor_type: str | None) -> str | None:
        """Resolve tumor type using OncoTree."""
        if not tumor_type:
            return None

        try:
            resolved = await self.oncotree_client.resolve_tumor_type(tumor_type)
            if resolved != tumor_type:
                print(f"  Resolved tumor type: {tumor_type} → {resolved}")
            return resolved
        except Exception as e:
            print(f"  Warning: OncoTree resolution failed: {str(e)}")
            return tumor_type

    def _handle_result(
        self,
        result: Any,
        source_name: str,
        sources_with_data: list[str],
        sources_failed: list[str],
    ) -> Any:
        """Handle an API result with error logging.

        Args:
            result: The result from asyncio.gather (could be Exception)
            source_name: Name of the source for logging
            sources_with_data: List to append to if data found
            sources_failed: List to append to if failed

        Returns:
            The result if successful, empty list if failed
        """
        if isinstance(result, Exception):
            print(f"  Warning: {source_name} failed: {str(result)}")
            sources_failed.append(source_name)
            return []
        elif result:
            sources_with_data.append(source_name)
            return result
        return []

    async def _fetch_literature(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
    ) -> list[PubMedEvidence]:
        """Fetch literature from Semantic Scholar with PubMed fallback.

        Uses semaphore to limit concurrent requests and avoid rate limits.
        Returns PubMedEvidence directly using the new *_evidence() methods.

        Returns:
            List of PubMedEvidence objects
        """
        if not self.semantic_scholar_client:
            return []

        # Use semaphore to limit concurrent literature API requests
        async with self._literature_semaphore:
            try:
                # Use new search_pubmed_evidence method
                return await self.semantic_scholar_client.search_pubmed_evidence(
                    gene=gene,
                    variant=variant,
                    tumor_type=tumor_type,
                    max_results=self.config.max_literature_results,
                )

            except SemanticScholarRateLimitError:
                # Fall back to PubMed (still within semaphore)
                return await self._fetch_pubmed_literature(gene, variant, tumor_type)

    async def _fetch_pubmed_literature(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
    ) -> list[PubMedEvidence]:
        """Fetch literature from PubMed (fallback from Semantic Scholar).

        Returns PubMedEvidence directly using the new *_evidence() method.
        Retry is handled by tenacity in PubMedClient.
        """
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
            # After tenacity retries exhausted, return empty
            return []

    async def build_evidence(
        self,
        variant: ParsedVariant | str,
        tumor_type: str | None = None,
    ) -> Evidence:
        """Build an Evidence for a single variant.

        Args:
            variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
            tumor_type: Optional tumor type for context

        Returns:
            Evidence with all aggregated evidence
        """
        # Handle string input
        if isinstance(variant, str):
            from oncomind.normalization import parse_variant_input
            variant = parse_variant_input(variant, tumor_type)

        gene = variant.gene
        normalized_variant = variant.variant_normalized or variant.variant
        tumor = tumor_type or variant.tumor_type

        # Track sources and failures
        sources_queried = ["MyVariant", "FDA", "CGI"]
        sources_with_data: list[str] = []
        sources_failed: list[str] = []
        processing_notes: list[str] = []

        # Resolve tumor type
        resolved_tumor = await self._resolve_tumor_type(tumor)

        # Build async fetch functions (using *_evidence methods)
        async def fetch_vicc():
            if self.vicc_client:
                sources_queried.append("VICC")
                return await self.vicc_client.fetch_vicc_evidence(
                    gene=gene,
                    variant=normalized_variant,
                    tumor_type=resolved_tumor,
                    max_results=self.config.max_vicc_results,
                )
            return []

        async def fetch_civic_assertions():
            if self.civic_client:
                sources_queried.append("CIViC")
                return await self.civic_client.fetch_assertion_evidence(
                    gene=gene,
                    variant=normalized_variant,
                    tumor_type=resolved_tumor,
                    max_results=self.config.max_civic_assertions,
                )
            return []

        async def fetch_clinical_trials():
            """Fetch clinical trials as evidence. Retry is handled by tenacity in client."""
            if not self.clinical_trials_client:
                return []

            sources_queried.append("ClinicalTrials.gov")

            try:
                return await self.clinical_trials_client.search_trial_evidence(
                    gene=gene,
                    variant=normalized_variant,
                    tumor_type=resolved_tumor,
                    recruiting_only=True,
                    max_results=self.config.max_clinical_trials,
                )
            except ClinicalTrialsRateLimitError:
                # After tenacity retries exhausted, return empty
                return []

        async def fetch_literature():
            if self.config.enable_literature:
                sources_queried.append("Literature")
                return await self._fetch_literature(gene, normalized_variant, resolved_tumor)
            return []

        async def fetch_cbioportal():
            if self.cbioportal_client:
                sources_queried.append("cBioPortal")
                return await self.cbioportal_client.fetch_co_mutation_data(
                    gene=gene,
                    variant=normalized_variant,
                    tumor_type=resolved_tumor,
                )
            return None

        async def fetch_cbioportal_cell_lines():
            """Fetch cell lines from cBioPortal CCLE (primary source for cell lines)."""
            if self.cbioportal_client:
                try:
                    return await self.cbioportal_client.fetch_cell_lines_with_mutation(
                        gene=gene,
                        variant=normalized_variant,
                    )
                except Exception:
                    return []
            return []

        async def fetch_depmap():
            if self.depmap_client:
                sources_queried.append("DepMap")
                return await self.depmap_client.fetch_depmap_evidence(
                    gene=gene,
                    variant=normalized_variant,
                )
            return None

        # Parallel fetch from all sources (using *_evidence methods where available)
        results = await asyncio.gather(
            self.myvariant_client.fetch_evidence(gene=gene, variant=normalized_variant),
            self.fda_client.fetch_approval_evidence(gene=gene, variant=normalized_variant),
            asyncio.to_thread(
                self.cgi_client.fetch_biomarker_evidence,
                gene,
                normalized_variant,
                resolved_tumor,
            ),
            fetch_vicc(),
            fetch_civic_assertions(),
            fetch_clinical_trials(),
            fetch_literature(),
            fetch_cbioportal(),
            fetch_cbioportal_cell_lines(),  # Primary source for cell lines
            fetch_depmap(),
            return_exceptions=True,
        )

        # Unpack results with error handling
        (
            myvariant_result,
            fda_result,
            cgi_result,
            vicc_result,
            civic_result,
            trials_result,
            literature_result,
            cbioportal_result,
            cbioportal_cell_lines_result,  # Cell lines from cBioPortal CCLE
            depmap_result,
        ) = results

        # Handle results using helper method
        myvariant_evidence = self._handle_result(
            myvariant_result, "MyVariant", sources_with_data, sources_failed
        ) or None  # Downstream code checks `if myvariant_evidence:`
        fda_approvals: list[FDAApproval] = self._handle_result(
            fda_result, "FDA", sources_with_data, sources_failed
        )
        all_cgi_biomarkers: list[CGIBiomarkerEvidence] = self._handle_result(
            cgi_result, "CGI", sources_with_data, sources_failed
        )

        # Split CGI biomarkers: FDA-approved go to kb, preclinical go to research
        cgi_biomarkers = [b for b in all_cgi_biomarkers if b.fda_approved]
        preclinical_biomarkers = [b for b in all_cgi_biomarkers if not b.fda_approved and b.evidence_level in ("Pre-clinical", "Cell line")]
        early_phase_biomarkers = [b for b in all_cgi_biomarkers if not b.fda_approved and b.evidence_level not in ("Pre-clinical", "Cell line", None)]

        vicc_evidence: list[VICCEvidence] = self._handle_result(
            vicc_result, "VICC", sources_with_data, sources_failed
        )
        civic_assertions: list[CIViCAssertionEvidence] = self._handle_result(
            civic_result, "CIViC", sources_with_data, sources_failed
        )
        clinical_trials: list[ClinicalTrialEvidence] = self._handle_result(
            trials_result, "ClinicalTrials.gov", sources_with_data, sources_failed
        )
        pubmed_articles: list[PubMedEvidence] = self._handle_result(
            literature_result, "Literature", sources_with_data, sources_failed
        )

        # Handle cBioPortal result - convert CoMutationData to CBioPortalEvidence
        cbioportal_evidence: CBioPortalEvidence | None = None
        if cbioportal_result and not isinstance(cbioportal_result, Exception):
            try:
                cbioportal_evidence = CBioPortalEvidence(
                    gene=cbioportal_result.gene,
                    variant=cbioportal_result.variant,
                    tumor_type=cbioportal_result.tumor_type,
                    study_id=cbioportal_result.study_id,
                    study_name=cbioportal_result.study_name,
                    total_samples=cbioportal_result.total_samples,
                    samples_with_gene_mutation=cbioportal_result.samples_with_gene_mutation,
                    samples_with_exact_variant=cbioportal_result.samples_with_exact_variant,
                    gene_prevalence_pct=cbioportal_result.gene_prevalence_pct,
                    variant_prevalence_pct=cbioportal_result.variant_prevalence_pct,
                    co_occurring=[CoMutationEntry(**c) for c in cbioportal_result.co_occurring],
                    mutually_exclusive=[CoMutationEntry(**m) for m in cbioportal_result.mutually_exclusive],
                )
                if cbioportal_evidence.has_data():
                    sources_with_data.append("cBioPortal")
            except Exception:
                sources_failed.append("cBioPortal")
        elif isinstance(cbioportal_result, Exception):
            sources_failed.append("cBioPortal")

        # Handle cell lines - cBioPortal CCLE is primary source, DepMap is fallback
        cell_line_models: list[CellLineModel] = []
        if cbioportal_cell_lines_result and not isinstance(cbioportal_cell_lines_result, Exception):
            # Convert cBioPortal cell line dicts to CellLineModel
            for cl in cbioportal_cell_lines_result:
                cell_line_models.append(CellLineModel(
                    name=cl.get("name", ""),
                    ccle_name=cl.get("sample_id"),
                    primary_disease=cl.get("tissue"),
                    has_mutation=True,
                    mutation_details=cl.get("protein_change"),
                ))
            if cell_line_models:
                sources_with_data.append("CCLE")

        # Handle DepMap result - already returns DepMapEvidence or None
        depmap_evidence = None
        if depmap_result and not isinstance(depmap_result, Exception):
            depmap_evidence = depmap_result
            # If we got cell lines from cBioPortal, use those (primary source)
            # Otherwise fall back to DepMap's cell line data
            if cell_line_models:
                # Replace DepMap cell lines with cBioPortal data
                depmap_evidence = DepMapEvidence(
                    gene=depmap_evidence.gene,
                    variant=depmap_evidence.variant,
                    gene_dependency=depmap_evidence.gene_dependency,
                    co_dependencies=depmap_evidence.co_dependencies,
                    drug_sensitivities=depmap_evidence.drug_sensitivities,
                    cell_line_models=cell_line_models,  # Use cBioPortal data
                    data_version=depmap_evidence.data_version,
                    n_cell_lines_screened=depmap_evidence.n_cell_lines_screened,
                )
            if depmap_evidence and depmap_evidence.has_data():
                sources_with_data.append("DepMap")
        elif isinstance(depmap_result, Exception):
            sources_failed.append("DepMap")
            # Even if DepMap failed, we might have cell lines from cBioPortal
            if cell_line_models:
                depmap_evidence = DepMapEvidence(
                    gene=gene,
                    variant=normalized_variant,
                    cell_line_models=cell_line_models,
                )

        # Extract data from MyVariant evidence
        functional_scores = FunctionalScores()
        identifiers_data: dict[str, Any] = {
            "variant_id": f"{gene}:{normalized_variant}",
            "gene": gene,
            "variant": variant.variant,
            "variant_normalized": normalized_variant,
            "variant_type": variant.variant_type,
        }

        clinvar_significance = None
        clinvar_entries = []
        cosmic_entries = []
        civic_entries = []

        if myvariant_evidence:
            # Extract functional scores
            functional_scores = FunctionalScores(
                alphamissense_score=getattr(myvariant_evidence, 'alphamissense_score', None),
                alphamissense_prediction=getattr(myvariant_evidence, 'alphamissense_prediction', None),
                cadd_score=getattr(myvariant_evidence, 'cadd_phred', None),
                cadd_raw=getattr(myvariant_evidence, 'cadd_raw', None),
                polyphen2_prediction=getattr(myvariant_evidence, 'polyphen2_prediction', None),
                polyphen2_score=getattr(myvariant_evidence, 'polyphen2_score', None),
                sift_prediction=getattr(myvariant_evidence, 'sift_prediction', None),
                sift_score=getattr(myvariant_evidence, 'sift_score', None),
                snpeff_effect=getattr(myvariant_evidence, 'snpeff_effect', None),
                snpeff_impact=getattr(myvariant_evidence, 'snpeff_impact', None),
                gnomad_exome_af=getattr(myvariant_evidence, 'gnomad_exome_af', None),
                gnomad_genome_af=getattr(myvariant_evidence, 'gnomad_genome_af', None),
            )

            # Extract identifiers
            identifiers_data.update({
                "cosmic_id": getattr(myvariant_evidence, 'cosmic_id', None),
                "ncbi_gene_id": getattr(myvariant_evidence, 'ncbi_gene_id', None),
                "dbsnp_id": getattr(myvariant_evidence, 'dbsnp_rsid', None),
                "clinvar_id": getattr(myvariant_evidence, 'clinvar_variation_id', None),
                "hgvs_genomic": getattr(myvariant_evidence, 'hgvs_genomic', None),
                "hgvs_protein": getattr(myvariant_evidence, 'hgvs_protein', None),
                "hgvs_transcript": getattr(myvariant_evidence, 'hgvs_coding', None),
            })

            # Extract ClinVar significance
            clinvar_significance = getattr(myvariant_evidence, 'clinvar_clinical_significance', None)

            # Extract evidence lists if available
            clinvar_entries = getattr(myvariant_evidence, 'clinvar', []) or []
            cosmic_entries = getattr(myvariant_evidence, 'cosmic', []) or []
            civic_entries = getattr(myvariant_evidence, 'civic', []) or []

        # Get gene context
        gene_context = get_gene_context(gene)
        gene_role = None
        gene_class = None
        pathway = None
        mutation_class = None

        if gene_context:
            if gene_context.role:
                gene_role = gene_context.role.value
            # GeneContext doesn't have functional_class/pathway - derive from role
            if gene_context.role:
                gene_class = gene_context.role.value  # Use role as the class
            # Check for pathway-actionable TSGs
            from oncomind.models.gene_context import get_pathway_actionable_info
            pathway_info = get_pathway_actionable_info(gene_context.gene)
            if pathway_info:
                pathway = pathway_info.get("pathway")

        # Build the Evidence with flat list structure
        evidence = Evidence(
            identifiers=VariantIdentifiers(**identifiers_data),
            functional=functional_scores,
            context=VariantContext(
                tumor_type=tumor,
                tumor_type_resolved=resolved_tumor,
                gene_role=gene_role,
                gene_class=gene_class,
                mutation_class=mutation_class,
                pathway=pathway,
            ),
            # Evidence lists (flat structure - frontend decides how to display)
            fda_approvals=fda_approvals,
            civic_assertions=civic_assertions,
            civic_evidence=civic_entries,
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
            # Track what was searched (for accurate gap detection)
            literature_searched=self.config.enable_literature,
        )

        return evidence

    async def build_evidences(
        self,
        variants: list[ParsedVariant | str],
        tumor_type: str | None = None,
    ) -> list[Evidence]:
        """Build Evidence for multiple variants in parallel.

        Args:
            variants: List of ParsedVariant objects or variant strings
            tumor_type: Optional tumor type (applied to all variants)

        Returns:
            List of Evidence objects
        """
        tasks = [
            self.build_evidence(v, tumor_type)
            for v in variants
        ]

        results = await asyncio.gather(*tasks, return_exceptions=True)

        # Filter out exceptions
        panels = []
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                print(f"  Warning: Failed to process variant {i}: {str(result)}")
            else:
                panels.append(result)

        return panels


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
        >>> print(evidence.get_approved_drugs())
        ['Dabrafenib', 'Vemurafenib']
    """
    async with EvidenceAggregator(config) as builder:
        return await builder.build_evidence(variant, tumor_type)


async def build_evidences(
    variants: list[ParsedVariant | str],
    tumor_type: str | None = None,
    config: EvidenceAggregatorConfig | None = None,
) -> list[Evidence]:
    """Convenience function to build Evidence for multiple variants.

    Args:
        variants: List of ParsedVariant objects or variant strings
        tumor_type: Optional tumor type (applied to all variants)
        config: Optional configuration for the builder

    Returns:
        List of Evidence objects
    """
    async with EvidenceAggregator(config) as builder:
        return await builder.build_evidences(variants, tumor_type)


__all__ = [
    "EvidenceAggregator",
    "EvidenceAggregatorConfig",
    "build_evidence",
    "build_evidences",
]

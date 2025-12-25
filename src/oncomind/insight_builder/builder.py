"""Evidence builder for aggregating variant evidence into Insight.

This module provides the core evidence aggregation logic as a reusable,
LLM-independent component.

ARCHITECTURE:
    ParsedVariant → build_insight() → Insight

    The builder:
    1. Fetches evidence from all API clients in parallel
    2. Handles API failures gracefully (logs warnings, continues)
    3. Assembles results into strongly-typed Insight
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
from oncomind.api.clinicaltrials import ClinicalTrialsClient, ClinicalTrialsRateLimitError
from oncomind.api.pubmed import PubMedClient, PubMedRateLimitError
from oncomind.api.semantic_scholar import SemanticScholarClient, SemanticScholarRateLimitError

from oncomind.models.insight import (
    Insight,
    VariantIdentifiers,
    KnowledgebaseEvidence,
    FunctionalScores,
    ClinicalContext,
    LiteratureEvidence,
)
from oncomind.models.insight.cgi import CGIBiomarkerEvidence
from oncomind.models.insight.civic import CIViCAssertionEvidence
from oncomind.models.insight.clinical_trials import ClinicalTrialEvidence
from oncomind.models.insight.fda import FDAApproval
from oncomind.models.insight.pubmed import PubMedEvidence

from oncomind.normalization import ParsedVariant
from oncomind.utils import normalize_variant
from oncomind.models.gene_context import get_gene_context, GeneRole


@dataclass
class InsightBuilderConfig:
    """Configuration for the insight builder.

    Controls which data sources to query and resource limits.
    """

    # Source toggles
    enable_vicc: bool = True
    enable_civic_assertions: bool = True
    enable_clinical_trials: bool = True
    enable_literature: bool = True

    # Result limits
    max_vicc_results: int = 15
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


class InsightBuilder:
    """Builder for aggregating variant evidence from multiple sources.

    Use as an async context manager to ensure proper HTTP session lifecycle:

        async with InsightBuilder() as builder:
            insight = await builder.build_insight(parsed_variant, tumor_type)

    Or for batch processing:

        async with InsightBuilder() as builder:
            insights = await builder.build_insights(variants, tumor_type)
    """

    def __init__(self, config: InsightBuilderConfig | None = None):
        self.config = config or InsightBuilderConfig()

        # Initialize API clients
        self.myvariant_client = MyVariantClient()
        self.fda_client = FDAClient()
        self.cgi_client = CGIClient()
        self.oncotree_client = OncoTreeClient()

        # Optional clients based on config
        self.vicc_client = VICCClient() if self.config.enable_vicc else None
        self.civic_client = CIViCClient() if self.config.enable_civic_assertions else None
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

    async def build_insight(
        self,
        variant: ParsedVariant | str,
        tumor_type: str | None = None,
    ) -> Insight:
        """Build an Insight for a single variant.

        Args:
            variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
            tumor_type: Optional tumor type for context

        Returns:
            Insight with all aggregated evidence
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
        ) = results

        # Handle MyVariant result
        myvariant_evidence = None
        if isinstance(myvariant_result, Exception):
            print(f"  Warning: MyVariant API failed: {str(myvariant_result)}")
            sources_failed.append("MyVariant")
        else:
            myvariant_evidence = myvariant_result
            if myvariant_evidence:
                sources_with_data.append("MyVariant")

        # Handle FDA result (now returns FDAApproval directly)
        fda_approvals: list[FDAApproval] = []
        if isinstance(fda_result, Exception):
            print(f"  Warning: FDA API failed: {str(fda_result)}")
            sources_failed.append("FDA")
        elif fda_result:
            fda_approvals = fda_result
            if fda_approvals:
                sources_with_data.append("FDA")

        # Handle CGI result (now returns CGIBiomarkerEvidence directly)
        cgi_biomarkers: list[CGIBiomarkerEvidence] = []
        if isinstance(cgi_result, Exception):
            print(f"  Warning: CGI biomarkers failed: {str(cgi_result)}")
            sources_failed.append("CGI")
        elif cgi_result:
            cgi_biomarkers = cgi_result
            if cgi_biomarkers:
                sources_with_data.append("CGI")

        # Handle VICC result (now returns VICCEvidence directly)
        vicc_evidence = []
        if isinstance(vicc_result, Exception):
            print(f"  Warning: VICC MetaKB API failed: {str(vicc_result)}")
            sources_failed.append("VICC")
        elif vicc_result:
            vicc_evidence = vicc_result
            if vicc_evidence:
                sources_with_data.append("VICC")

        # Handle CIViC result (now returns CIViCAssertionEvidence directly)
        civic_assertions: list[CIViCAssertionEvidence] = []
        if isinstance(civic_result, Exception):
            print(f"  Warning: CIViC Assertions API failed: {str(civic_result)}")
            sources_failed.append("CIViC")
        elif civic_result:
            civic_assertions = civic_result
            if civic_assertions:
                sources_with_data.append("CIViC")

        # Handle clinical trials result (now returns ClinicalTrialEvidence directly)
        clinical_trials: list[ClinicalTrialEvidence] = []
        if isinstance(trials_result, Exception):
            print(f"  Warning: ClinicalTrials.gov API failed: {str(trials_result)}")
            sources_failed.append("ClinicalTrials.gov")
        elif trials_result:
            clinical_trials = trials_result
            if clinical_trials:
                sources_with_data.append("ClinicalTrials.gov")

        # Handle literature result (now returns PubMedEvidence directly)
        pubmed_articles: list[PubMedEvidence] = []
        if isinstance(literature_result, Exception):
            print(f"  Warning: Literature search failed: {str(literature_result)}")
            sources_failed.append("Literature")
        elif literature_result:
            pubmed_articles = literature_result
            if pubmed_articles:
                sources_with_data.append("Literature")

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

        # Build the Insight
        panel = Insight(
            identifiers=VariantIdentifiers(**identifiers_data),
            kb=KnowledgebaseEvidence(
                civic=civic_entries,
                civic_assertions=civic_assertions,
                clinvar=clinvar_entries,
                cosmic=cosmic_entries,
                cgi_biomarkers=cgi_biomarkers,
                vicc=vicc_evidence,
            ),
            functional=functional_scores,
            clinical=ClinicalContext(
                tumor_type=tumor,
                tumor_type_resolved=resolved_tumor,
                fda_approvals=fda_approvals,
                clinical_trials=clinical_trials,
                gene_role=gene_role,
                gene_class=gene_class,
                mutation_class=mutation_class,
                pathway=pathway,
                clinvar_clinical_significance=clinvar_significance,
            ),
            literature=LiteratureEvidence(
                pubmed_articles=pubmed_articles,
                literature_knowledge=None,  # Set by LLM layer later
                key_pmids=[a.pmid for a in pubmed_articles[:5] if a.pmid],
            ),
        )

        return panel

    async def build_insights(
        self,
        variants: list[ParsedVariant | str],
        tumor_type: str | None = None,
    ) -> list[Insight]:
        """Build Insights for multiple variants in parallel.

        Args:
            variants: List of ParsedVariant objects or variant strings
            tumor_type: Optional tumor type (applied to all variants)

        Returns:
            List of Insight objects
        """
        tasks = [
            self.build_insight(v, tumor_type)
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


async def build_insight(
    variant: ParsedVariant | str,
    tumor_type: str | None = None,
    config: InsightBuilderConfig | None = None,
) -> Insight:
    """Convenience function to build an Insight for a single variant.

    This is the recommended entry point for single-variant annotation.

    Args:
        variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
        tumor_type: Optional tumor type for context
        config: Optional configuration for the builder

    Returns:
        Insight with all aggregated evidence

    Example:
        >>> panel = await build_insight("BRAF V600E", tumor_type="Melanoma")
        >>> print(panel.identifiers.gene)
        BRAF
        >>> print(panel.clinical.get_approved_drugs())
        ['Dabrafenib', 'Vemurafenib']
    """
    async with InsightBuilder(config) as builder:
        return await builder.build_insight(variant, tumor_type)


async def build_insights(
    variants: list[ParsedVariant | str],
    tumor_type: str | None = None,
    config: InsightBuilderConfig | None = None,
) -> list[Insight]:
    """Convenience function to build Insights for multiple variants.

    Args:
        variants: List of ParsedVariant objects or variant strings
        tumor_type: Optional tumor type (applied to all variants)
        config: Optional configuration for the builder

    Returns:
        List of Insight objects
    """
    async with InsightBuilder(config) as builder:
        return await builder.build_insights(variants, tumor_type)


__all__ = [
    "InsightBuilder",
    "InsightBuilderConfig",
    "build_insight",
    "build_insights",
]

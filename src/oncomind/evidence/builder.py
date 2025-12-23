"""Evidence builder for aggregating variant evidence into EvidencePanel.

This module provides the core evidence aggregation logic, extracting it from
the InsightEngine into a reusable, LLM-independent component.

ARCHITECTURE:
    ParsedVariant → build_evidence_panel() → EvidencePanel

    The builder:
    1. Fetches evidence from all API clients in parallel
    2. Handles API failures gracefully (logs warnings, continues)
    3. Assembles results into strongly-typed EvidencePanel
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

from oncomind.api.myvariant import MyVariantClient
from oncomind.api.fda import FDAClient
from oncomind.api.cgi import CGIClient
from oncomind.api.oncotree import OncoTreeClient
from oncomind.api.vicc import VICCClient
from oncomind.api.civic import CIViCClient
from oncomind.api.clinicaltrials import ClinicalTrialsClient
from oncomind.api.pubmed import PubMedClient, PubMedRateLimitError
from oncomind.api.semantic_scholar import SemanticScholarClient, SemanticScholarRateLimitError

from oncomind.models.evidence.evidence_panel import (
    EvidencePanel,
    VariantIdentifiers,
    KnowledgebaseEvidence,
    FunctionalScores,
    ClinicalContext,
    LiteratureEvidence,
    EvidenceMeta,
)
from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
from oncomind.models.evidence.civic import CIViCAssertionEvidence
from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
from oncomind.models.evidence.fda import FDAApproval
from oncomind.models.evidence.pubmed import PubMedEvidence

from oncomind.normalization import ParsedVariant
from oncomind.utils import normalize_variant
from oncomind.models.gene_context import get_gene_context, GeneRole


@dataclass
class EvidenceBuilderConfig:
    """Configuration for the evidence builder.

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


class EvidenceBuilder:
    """Builder for aggregating variant evidence from multiple sources.

    Use as an async context manager to ensure proper HTTP session lifecycle:

        async with EvidenceBuilder() as builder:
            panel = await builder.build_evidence_panel(parsed_variant, tumor_type)

    Or for batch processing:

        async with EvidenceBuilder() as builder:
            panels = await builder.build_evidence_panels(variants, tumor_type)
    """

    def __init__(self, config: EvidenceBuilderConfig | None = None):
        self.config = config or EvidenceBuilderConfig()

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
    ) -> tuple[list, str | None]:
        """Fetch literature from Semantic Scholar with PubMed fallback.

        Uses semaphore to limit concurrent requests and avoid rate limits.

        Returns:
            Tuple of (papers/articles, source_name)
        """
        if not self.semantic_scholar_client:
            return [], None

        # Use semaphore to limit concurrent literature API requests
        async with self._literature_semaphore:
            try:
                # Search both resistance AND general variant literature
                resistance_papers, variant_papers = await asyncio.gather(
                    self.semantic_scholar_client.search_resistance_literature(
                        gene=gene,
                        variant=variant,
                        tumor_type=tumor_type,
                        max_results=self.config.max_literature_results // 2,
                    ),
                    self.semantic_scholar_client.search_variant_literature(
                        gene=gene,
                        variant=variant,
                        tumor_type=tumor_type,
                        max_results=self.config.max_literature_results // 2,
                    ),
                )

                # Merge and deduplicate by paper_id
                seen_ids = set()
                merged_papers = []
                for paper in resistance_papers + variant_papers:
                    if paper.paper_id not in seen_ids:
                        seen_ids.add(paper.paper_id)
                        merged_papers.append(paper)

                return merged_papers[:self.config.max_literature_results], "semantic_scholar"

            except SemanticScholarRateLimitError:
                # Fall back to PubMed (still within semaphore)
                return await self._fetch_pubmed_literature(gene, variant, tumor_type)

    async def _fetch_pubmed_literature(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None,
    ) -> tuple[list, str | None]:
        """Fetch literature from PubMed with fast retry.

        Uses short backoff: 0.5s, 1s between retries (max 2 retries).
        """
        if not self.pubmed_client:
            return [], None

        import random

        max_retries = 2  # Reduced from 4 for speed
        for attempt in range(max_retries):
            try:
                resistance_articles, variant_articles = await asyncio.gather(
                    self.pubmed_client.search_resistance_literature(
                        gene=gene,
                        variant=variant,
                        tumor_type=tumor_type,
                        max_results=self.config.max_literature_results // 2,
                    ),
                    self.pubmed_client.search_variant_literature(
                        gene=gene,
                        variant=variant,
                        tumor_type=tumor_type,
                        max_results=self.config.max_literature_results // 2,
                    ),
                )

                # Merge and deduplicate by pmid
                seen_pmids = set()
                merged_articles = []
                for article in resistance_articles + variant_articles:
                    if article.pmid not in seen_pmids:
                        seen_pmids.add(article.pmid)
                        merged_articles.append(article)

                return merged_articles[:self.config.max_literature_results], "pubmed"

            except PubMedRateLimitError:
                if attempt < max_retries - 1:
                    # Short backoff: 0.5s, 1s + small jitter
                    wait_time = 0.5 * (attempt + 1) + random.uniform(0, 0.3)
                    await asyncio.sleep(wait_time)

        # Silently return empty - meta.sources_failed will indicate the issue
        return [], None

    async def build_evidence_panel(
        self,
        variant: ParsedVariant | str,
        tumor_type: str | None = None,
    ) -> EvidencePanel:
        """Build an EvidencePanel for a single variant.

        Args:
            variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
            tumor_type: Optional tumor type for context

        Returns:
            EvidencePanel with all aggregated evidence
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

        # Build async fetch functions
        async def fetch_vicc():
            if self.vicc_client:
                sources_queried.append("VICC")
                return await self.vicc_client.fetch_associations(
                    gene=gene,
                    variant=normalized_variant,
                    tumor_type=resolved_tumor,
                    max_results=self.config.max_vicc_results,
                )
            return []

        async def fetch_civic_assertions():
            if self.civic_client:
                sources_queried.append("CIViC")
                return await self.civic_client.fetch_assertions(
                    gene=gene,
                    variant=normalized_variant,
                    tumor_type=resolved_tumor,
                    max_results=self.config.max_civic_assertions,
                )
            return []

        async def fetch_clinical_trials():
            if self.clinical_trials_client:
                sources_queried.append("ClinicalTrials.gov")
                return await self.clinical_trials_client.search_trials(
                    gene=gene,
                    variant=normalized_variant,
                    tumor_type=resolved_tumor,
                    recruiting_only=True,
                    max_results=self.config.max_clinical_trials,
                )
            return []

        async def fetch_literature():
            if self.config.enable_literature:
                sources_queried.append("Literature")
                return await self._fetch_literature(gene, normalized_variant, resolved_tumor)
            return [], None

        # Parallel fetch from all sources
        results = await asyncio.gather(
            self.myvariant_client.fetch_evidence(gene=gene, variant=normalized_variant),
            self.fda_client.fetch_drug_approvals(gene=gene, variant=normalized_variant),
            asyncio.to_thread(
                self.cgi_client.fetch_biomarkers,
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

        # Handle FDA result
        fda_approvals: list[FDAApproval] = []
        if isinstance(fda_result, Exception):
            print(f"  Warning: FDA API failed: {str(fda_result)}")
            sources_failed.append("FDA")
        elif fda_result:
            for approval_record in fda_result:
                parsed = self.fda_client.parse_approval_data(
                    approval_record, gene, normalized_variant
                )
                if parsed:
                    fda_approvals.append(FDAApproval(**parsed))
            if fda_approvals:
                sources_with_data.append("FDA")

        # Handle CGI result
        cgi_biomarkers: list[CGIBiomarkerEvidence] = []
        if isinstance(cgi_result, Exception):
            print(f"  Warning: CGI biomarkers failed: {str(cgi_result)}")
            sources_failed.append("CGI")
        elif cgi_result:
            for biomarker in cgi_result:
                cgi_biomarkers.append(CGIBiomarkerEvidence(
                    gene=biomarker.gene,
                    alteration=biomarker.alteration,
                    drug=biomarker.drug,
                    drug_status=biomarker.drug_status,
                    association=biomarker.association,
                    evidence_level=biomarker.evidence_level,
                    source=biomarker.source,
                    tumor_type=biomarker.tumor_type,
                    fda_approved=biomarker.is_fda_approved(),
                ))
            if cgi_biomarkers:
                sources_with_data.append("CGI")

        # Handle VICC result
        vicc_evidence = []
        if isinstance(vicc_result, Exception):
            print(f"  Warning: VICC MetaKB API failed: {str(vicc_result)}")
            sources_failed.append("VICC")
        elif vicc_result:
            from oncomind.models.evidence.vicc import VICCEvidence
            for assoc in vicc_result:
                vicc_evidence.append(VICCEvidence(
                    description=assoc.description,
                    gene=assoc.gene,
                    variant=assoc.variant,
                    disease=assoc.disease,
                    drugs=assoc.drugs,
                    evidence_level=assoc.evidence_level,
                    response_type=assoc.response_type,
                    source=assoc.source,
                    publication_url=assoc.publication_url,
                    oncogenic=assoc.oncogenic,
                    is_sensitivity=assoc.is_sensitivity(),
                    is_resistance=assoc.is_resistance(),
                    oncokb_level=assoc.get_oncokb_level(),
                ))
            if vicc_evidence:
                sources_with_data.append("VICC")

        # Handle CIViC result
        civic_assertions: list[CIViCAssertionEvidence] = []
        if isinstance(civic_result, Exception):
            print(f"  Warning: CIViC Assertions API failed: {str(civic_result)}")
            sources_failed.append("CIViC")
        elif civic_result:
            for assertion in civic_result:
                civic_assertions.append(CIViCAssertionEvidence(
                    assertion_id=assertion.assertion_id,
                    name=assertion.name,
                    amp_level=assertion.amp_level,
                    amp_tier=assertion.get_amp_tier(),
                    amp_level_letter=assertion.get_amp_level(),
                    assertion_type=assertion.assertion_type,
                    significance=assertion.significance,
                    status=assertion.status,
                    molecular_profile=assertion.molecular_profile,
                    disease=assertion.disease,
                    therapies=assertion.therapies,
                    fda_companion_test=assertion.fda_companion_test,
                    nccn_guideline=assertion.nccn_guideline,
                    description=assertion.description,
                    is_sensitivity=assertion.is_sensitivity(),
                    is_resistance=assertion.is_resistance(),
                ))
            if civic_assertions:
                sources_with_data.append("CIViC")

        # Handle clinical trials result
        clinical_trials: list[ClinicalTrialEvidence] = []
        if isinstance(trials_result, Exception):
            print(f"  Warning: ClinicalTrials.gov API failed: {str(trials_result)}")
            sources_failed.append("ClinicalTrials.gov")
        elif trials_result:
            for trial in trials_result:
                variant_specific = trial.mentions_variant(normalized_variant, gene=gene)
                clinical_trials.append(ClinicalTrialEvidence(
                    nct_id=trial.nct_id,
                    title=trial.title,
                    status=trial.status,
                    phase=trial.phase,
                    conditions=trial.conditions,
                    interventions=trial.interventions,
                    sponsor=trial.sponsor,
                    url=trial.url,
                    variant_specific=variant_specific,
                ))
            if clinical_trials:
                sources_with_data.append("ClinicalTrials.gov")

        # Handle literature result
        pubmed_articles: list[PubMedEvidence] = []
        literature_source = None
        if isinstance(literature_result, Exception):
            print(f"  Warning: Literature search failed: {str(literature_result)}")
            sources_failed.append("Literature")
        elif literature_result:
            literature_items, literature_source = literature_result
            if literature_items:
                # Convert to PubMedEvidence (basic conversion without LLM scoring)
                for item in literature_items:
                    if literature_source == "semantic_scholar":
                        url = (
                            f"https://pubmed.ncbi.nlm.nih.gov/{item.pmid}/"
                            if item.pmid
                            else f"https://www.semanticscholar.org/paper/{item.paper_id}"
                        )
                        pubmed_articles.append(PubMedEvidence(
                            pmid=item.pmid or item.paper_id,
                            title=item.title,
                            abstract=item.abstract or "",
                            authors=[],
                            journal=item.venue or "",
                            year=str(item.year) if item.year else None,
                            doi=None,
                            url=url,
                            signal_type=None,  # Set by LLM layer later
                            drugs_mentioned=[],
                            citation_count=item.citation_count,
                            influential_citation_count=item.influential_citation_count,
                            tldr=item.tldr,
                            is_open_access=item.is_open_access,
                            open_access_pdf_url=item.open_access_pdf_url,
                            semantic_scholar_id=item.paper_id,
                        ))
                    else:
                        pubmed_articles.append(PubMedEvidence(
                            pmid=item.pmid,
                            title=item.title,
                            abstract=item.abstract,
                            authors=item.authors,
                            journal=item.journal,
                            year=item.year,
                            doi=item.doi,
                            url=item.url,
                            signal_type=None,
                            drugs_mentioned=[],
                            citation_count=None,
                            influential_citation_count=None,
                            tldr=None,
                            is_open_access=None,
                            open_access_pdf_url=None,
                            semantic_scholar_id=None,
                        ))
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

        # Build the EvidencePanel
        panel = EvidencePanel(
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
            meta=EvidenceMeta(
                sources_queried=sources_queried,
                sources_with_data=sources_with_data,
                sources_failed=sources_failed,
                processing_notes=processing_notes,
            ),
        )

        return panel

    async def build_evidence_panels(
        self,
        variants: list[ParsedVariant | str],
        tumor_type: str | None = None,
    ) -> list[EvidencePanel]:
        """Build EvidencePanels for multiple variants in parallel.

        Args:
            variants: List of ParsedVariant objects or variant strings
            tumor_type: Optional tumor type (applied to all variants)

        Returns:
            List of EvidencePanel objects
        """
        tasks = [
            self.build_evidence_panel(v, tumor_type)
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


async def build_evidence_panel(
    variant: ParsedVariant | str,
    tumor_type: str | None = None,
    config: EvidenceBuilderConfig | None = None,
) -> EvidencePanel:
    """Convenience function to build an EvidencePanel for a single variant.

    This is the recommended entry point for single-variant annotation.

    Args:
        variant: ParsedVariant object or variant string (e.g., "BRAF V600E")
        tumor_type: Optional tumor type for context
        config: Optional configuration for the builder

    Returns:
        EvidencePanel with all aggregated evidence

    Example:
        >>> panel = await build_evidence_panel("BRAF V600E", tumor_type="Melanoma")
        >>> print(panel.identifiers.gene)
        BRAF
        >>> print(panel.clinical.get_approved_drugs())
        ['Dabrafenib', 'Vemurafenib']
    """
    async with EvidenceBuilder(config) as builder:
        return await builder.build_evidence_panel(variant, tumor_type)


async def build_evidence_panels(
    variants: list[ParsedVariant | str],
    tumor_type: str | None = None,
    config: EvidenceBuilderConfig | None = None,
) -> list[EvidencePanel]:
    """Convenience function to build EvidencePanels for multiple variants.

    Args:
        variants: List of ParsedVariant objects or variant strings
        tumor_type: Optional tumor type (applied to all variants)
        config: Optional configuration for the builder

    Returns:
        List of EvidencePanel objects
    """
    async with EvidenceBuilder(config) as builder:
        return await builder.build_evidence_panels(variants, tumor_type)


__all__ = [
    "EvidenceBuilder",
    "EvidenceBuilderConfig",
    "build_evidence_panel",
    "build_evidence_panels",
]

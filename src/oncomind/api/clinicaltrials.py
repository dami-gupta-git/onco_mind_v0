"""ClinicalTrials.gov API client for fetching clinical trial data.

ARCHITECTURE:
    Gene + Variant + Tumor Type → ClinicalTrials.gov API v2 → Active trials

Fetches recruiting clinical trials for cancer biomarkers to support
Tier II classification for variants with active investigational therapies.

Key Design:
- Async HTTP with connection pooling (httpx.AsyncClient)
- Searches by gene+variant keyword and cancer type
- Filters for recruiting/active trials only
- Returns structured trial information
- Rate limiting: ~50 requests/minute, exponential backoff on 429/403
"""

from typing import Any
from dataclasses import dataclass

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential_jitter,
)

from oncomind.models.evidence.base import tumor_types_match


class ClinicalTrialsRateLimitError(Exception):
    """Raised when rate limited by ClinicalTrials.gov API.

    Attributes:
        retry_after: Seconds to wait before retrying (from Retry-After header)
    """

    def __init__(self, message: str, retry_after: float | None = None):
        super().__init__(message)
        self.retry_after = retry_after


class ClinicalTrialsError(Exception):
    """Exception raised for ClinicalTrials.gov API errors."""

    pass


@dataclass
class ClinicalTrial:
    """A clinical trial from ClinicalTrials.gov."""

    nct_id: str
    title: str
    status: str  # RECRUITING, ACTIVE_NOT_RECRUITING, etc.
    phase: str | None
    conditions: list[str]
    interventions: list[str]
    brief_summary: str | None
    eligibility_criteria: str | None
    sponsor: str | None
    url: str

    def is_recruiting(self) -> bool:
        """Check if trial is actively recruiting."""
        return self.status in ["RECRUITING", "ENROLLING_BY_INVITATION"]

    def is_active(self) -> bool:
        """Check if trial is active (recruiting or not)."""
        return self.status in [
            "RECRUITING",
            "ENROLLING_BY_INVITATION",
            "ACTIVE_NOT_RECRUITING",
            "NOT_YET_RECRUITING",
        ]

    def mentions_variant(
        self, variant: str | None, gene: str | None = None
    ) -> tuple[str, str | None]:
        """Check if trial mentions the gene and/or variant with codon-level matching.

        Uses parse_biomarker_specificity() to extract what variant the trial describes,
        then is_variant_covered() to compare against the query. This enables proper
        codon-level matching (e.g., G12C trial matches G12D query at codon level).

        Args:
            variant: Variant notation (e.g., "G12D", "V600E"). Can be None to check gene only.
            gene: Gene symbol (e.g., "NRAS", "BRAF"). If provided, ensures
                  the variant is mentioned in context of this gene to avoid
                  false positives (e.g., KRAS G12D trial matching NRAS G12D query).

        Returns:
            Tuple of (match_type, matched_biomarker):
            - ("variant", "KRAS G12C") - exact variant match
            - ("codon", "KRAS G12") - same codon, different variant (e.g., G12C trial, G12D query)
            - ("ambiguous", "KRAS G12") - gene found and ambiguous variant in BROAD_VARIANTS
            - ("gene", "KRAS") - only gene found in text
            - ("none", None) - no match found
        """
        from oncomind.config.constants import BROAD_VARIANTS
        from oncomind.insight_builder.fda_processor import (
            parse_biomarker_specificity,
            is_variant_covered,
        )

        import re

        gene_upper = gene.upper() if gene else None

        # Build gene pattern for word boundary matching
        gene_pattern = rf"\b{re.escape(gene_upper)}\b" if gene_upper else None

        # Priority check: gene must appear in title OR eligibility criteria
        # If gene only appears in summary (background info), it's not a real match
        # e.g., a trial for "IDH1-Mutated AML" may mention "IDH1 or IDH2" in background
        # but is only actually studying IDH1
        title_upper = self.title.upper() if self.title else ""
        eligibility_upper = (
            self.eligibility_criteria.upper() if self.eligibility_criteria else ""
        )
        primary_text = f"{title_upper} {eligibility_upper}"

        gene_in_primary = bool(gene_pattern and re.search(gene_pattern, primary_text))

        if not gene_in_primary:
            return ("none", None)

        # Combine all text for deeper analysis (variant matching, etc.)
        search_texts = [self.title]
        if self.eligibility_criteria:
            search_texts.append(self.eligibility_criteria)
        if self.brief_summary:
            search_texts.append(self.brief_summary)

        full_text = " ".join(search_texts)

        # Use parse_biomarker_specificity to extract what the trial describes
        biomarker_spec = parse_biomarker_specificity(full_text, gene)

        if biomarker_spec and variant:
            # Use is_variant_covered for proper codon-level matching
            covered, match_level = is_variant_covered(variant, biomarker_spec)

            if match_level == "variant":
                # Extract the matched variant from the spec
                matched_var = biomarker_spec.get("specified_variant")
                if not matched_var and biomarker_spec.get("specified_variants"):
                    # For variant_list, find the matching one
                    matched_var = next(
                        (
                            v
                            for v in biomarker_spec["specified_variants"]
                            if v.upper() == variant.upper()
                        ),
                        biomarker_spec["specified_variants"][0],
                    )
                biomarker = (
                    f"{gene} {matched_var}" if matched_var else f"{gene} {variant}"
                )
                return ("variant", biomarker)

            elif match_level == "codon":
                # Codon-level match (e.g., trial has G12C, query is G12D)
                codon = biomarker_spec.get("codon")
                biomarker = f"{gene} {codon}" if codon else gene
                return ("codon", biomarker)

            elif match_level == "gene":
                return ("gene", gene)

        # Fallback: check for ambiguous variants (BROAD_VARIANTS like G12, V600)
        # Use primary_text (title + eligibility) for consistency
        ambig_variants = [
            (g, v)
            for (g, v) in BROAD_VARIANTS
            if g.upper() in primary_text and v.upper() in primary_text
        ]
        matched_ambig = next(
            ((g, v) for (g, v) in ambig_variants if g.upper() == gene_upper), None
        )
        if matched_ambig:
            g, v = matched_ambig
            return ("ambiguous", f"{g} {v}")

        # Gene-level only (we already verified gene_in_primary above)
        return ("gene", gene)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "nct_id": self.nct_id,
            "title": self.title,
            "status": self.status,
            "phase": self.phase,
            "conditions": self.conditions,
            "interventions": self.interventions,
            "brief_summary": self.brief_summary[:500] if self.brief_summary else None,
            "sponsor": self.sponsor,
            "url": self.url,
        }


class ClinicalTrialsClient:
    """Client for ClinicalTrials.gov API v2.

    Uses the modern v2 API (launched March 2024) to search for
    clinical trials by gene, variant, and cancer type.

    API Documentation: https://clinicaltrials.gov/data-api/api

    Rate Limiting:
    - ~50 requests/minute per IP
    - Exponential backoff with jitter on 429/403 errors
    - Request only necessary fields to reduce payload
    """

    BASE_URL = "https://clinicaltrials.gov/api/v2/studies"
    DEFAULT_TIMEOUT = 30.0
    DEFAULT_PAGE_SIZE = 100  # Increased to reduce total requests
    USER_AGENT = "OncoMind/0.1.0 (contact: oncomind-research@example.com)"

    # Only request fields we actually use (reduces payload size)
    FIELDS = [
        "protocolSection.identificationModule.nctId",
        "protocolSection.identificationModule.briefTitle",
        "protocolSection.identificationModule.officialTitle",
        "protocolSection.statusModule.overallStatus",
        "protocolSection.designModule.phases",
        "protocolSection.conditionsModule.conditions",
        "protocolSection.armsInterventionsModule.interventions",
        "protocolSection.descriptionModule.briefSummary",
        "protocolSection.eligibilityModule.eligibilityCriteria",
        "protocolSection.sponsorCollaboratorsModule.leadSponsor",
    ]

    def __init__(self, timeout: float = DEFAULT_TIMEOUT):
        """Initialize the ClinicalTrials client.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout
        self._client: httpx.AsyncClient | None = None
        self._headers = {
            "User-Agent": self.USER_AGENT,
            "Accept": "application/json",
        }

    async def __aenter__(self) -> "ClinicalTrialsClient":
        """Async context manager entry."""
        self._client = httpx.AsyncClient(timeout=self.timeout, headers=self._headers)
        return self

    async def __aexit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Async context manager exit."""
        if self._client:
            await self._client.aclose()
            self._client = None

    def _get_client(self) -> httpx.AsyncClient:
        """Get or create the HTTP client."""
        if self._client is None:
            self._client = httpx.AsyncClient(
                timeout=self.timeout, headers=self._headers
            )
        return self._client

    def _build_search_query(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
    ) -> str:
        """Build search query string.

        Args:
            gene: Gene symbol (e.g., "KRAS")
            variant: Optional variant notation (e.g., "G12D")
            tumor_type: Optional tumor type (e.g., "Pancreatic Cancer")

        Returns:
            Search query string
        """
        # Build query parts
        parts = [gene.upper()]

        if variant:
            # Add variant in multiple formats for better matching
            variant_clean = variant.upper().replace("P.", "")
            parts.append(variant_clean)
            # Also try gene+variant concatenated (e.g., "KRASG12D")
            parts.append(f"{gene.upper()}{variant_clean}")

        # Join with OR for broader matching
        query = " OR ".join(parts)

        return query

    def _parse_study(self, study: dict[str, Any]) -> ClinicalTrial | None:
        """Parse a study from the API response.

        Args:
            study: Raw study data from API

        Returns:
            ClinicalTrial object or None if parsing fails
        """
        try:
            protocol = study.get("protocolSection", {})

            # Identification
            id_module = protocol.get("identificationModule", {})
            nct_id = id_module.get("nctId", "")
            title = id_module.get("briefTitle", "") or id_module.get(
                "officialTitle", ""
            )

            # Status
            status_module = protocol.get("statusModule", {})
            status = status_module.get("overallStatus", "UNKNOWN")

            # Design - phase
            design_module = protocol.get("designModule", {})
            phases = design_module.get("phases", [])
            phase = phases[0] if phases else None

            # Conditions
            conditions_module = protocol.get("conditionsModule", {})
            conditions = conditions_module.get("conditions", [])

            # Interventions
            arms_module = protocol.get("armsInterventionsModule", {})
            interventions_list = arms_module.get("interventions", [])
            interventions = [
                i.get("name", "") for i in interventions_list if i.get("name")
            ]

            # Description
            desc_module = protocol.get("descriptionModule", {})
            brief_summary = desc_module.get("briefSummary", "")

            # Eligibility
            eligibility_module = protocol.get("eligibilityModule", {})
            eligibility_criteria = eligibility_module.get("eligibilityCriteria", "")

            # Sponsor
            sponsor_module = protocol.get("sponsorCollaboratorsModule", {})
            lead_sponsor = sponsor_module.get("leadSponsor", {})
            sponsor = lead_sponsor.get("name", "")

            # Build URL
            url = f"https://clinicaltrials.gov/study/{nct_id}"

            return ClinicalTrial(
                nct_id=nct_id,
                title=title,
                status=status,
                phase=phase,
                conditions=conditions,
                interventions=interventions,
                brief_summary=brief_summary,
                eligibility_criteria=eligibility_criteria,
                sponsor=sponsor,
                url=url,
            )
        except Exception:
            return None

    def _parse_retry_after(self, response: httpx.Response) -> float | None:
        """Parse Retry-After header from response.

        Args:
            response: HTTP response object

        Returns:
            Seconds to wait, or None if header not present/parseable
        """
        retry_after = response.headers.get("Retry-After")
        if not retry_after:
            return None

        try:
            # Retry-After can be seconds (integer) or HTTP-date
            # Try parsing as seconds first (most common)
            return float(retry_after)
        except ValueError:
            # Could be HTTP-date format, but for simplicity return None
            # and let tenacity handle the backoff
            return None

    @retry(
        retry=retry_if_exception_type(ClinicalTrialsRateLimitError),
        stop=stop_after_attempt(3),
        wait=wait_exponential_jitter(initial=2, max=15, jitter=1),
        reraise=True,
    )
    async def _make_request(
        self,
        params: dict[str, Any],
    ) -> dict[str, Any]:
        """Make API request with retry on rate limit.

        Args:
            params: Query parameters for the API request

        Returns:
            Parsed JSON response

        Raises:
            ClinicalTrialsRateLimitError: If rate limited after retries
            ClinicalTrialsError: For other API errors
        """
        client = self._get_client()

        try:
            response = await client.get(self.BASE_URL, params=params)

            # Check for rate limit errors (429 or 403)
            if response.status_code in (429, 403):
                retry_after = self._parse_retry_after(response)
                raise ClinicalTrialsRateLimitError(
                    f"Rate limited by ClinicalTrials.gov (status {response.status_code})",
                    retry_after=retry_after,
                )

            response.raise_for_status()
            return response.json()

        except httpx.HTTPStatusError as e:
            if e.response.status_code in (429, 403):
                retry_after = self._parse_retry_after(e.response)
                raise ClinicalTrialsRateLimitError(
                    f"Rate limited by ClinicalTrials.gov (status {e.response.status_code})",
                    retry_after=retry_after,
                )
            raise ClinicalTrialsError(f"HTTP error: {e}")

        except httpx.HTTPError as e:
            raise ClinicalTrialsError(f"Network error: {e}")

    async def search_trials(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        recruiting_only: bool = True,
        max_results: int = 10,
    ) -> list[ClinicalTrial]:
        """Search for clinical trials by gene/variant/tumor type.

        Args:
            gene: Gene symbol (e.g., "KRAS")
            variant: Optional variant notation (e.g., "G12D")
            tumor_type: Optional tumor type filter
            recruiting_only: If True, only return recruiting trials
            max_results: Maximum number of results

        Returns:
            List of ClinicalTrial objects
        """
        # Build query
        query = self._build_search_query(gene, variant, tumor_type)

        # Filter to last 10 years
        from datetime import datetime

        current_year = datetime.now().year
        min_year = current_year - 10

        # Build parameters - request only needed fields to reduce payload
        params: dict[str, Any] = {
            "query.term": query,
            "pageSize": min(max_results * 2, self.DEFAULT_PAGE_SIZE),
            "countTotal": "true",
            "format": "json",
            "fields": "|".join(self.FIELDS),  # Only fetch fields we need
            "filter.advanced": f"AREA[StartDate]RANGE[01/01/{min_year}, MAX]",  # Started in last 10 years
        }

        # Add condition filter if tumor type specified
        if tumor_type:
            params["query.cond"] = tumor_type

        # Filter by status if recruiting only
        if recruiting_only:
            params["filter.overallStatus"] = (
                "RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING"
            )
        else:
            params["filter.overallStatus"] = (
                "RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING|ACTIVE_NOT_RECRUITING"
            )

        # Make request - let errors propagate for caller to handle retry
        data = await self._make_request(params)

        # Parse studies
        studies = data.get("studies", [])
        trials = []

        for study in studies:
            trial = self._parse_study(study)
            if trial is None:
                continue

            # Filter by variant mention if variant specified
            # Pass gene to avoid false positives (e.g., KRAS G12D matching NRAS G12D query)
            if variant:
                match_type, _ = trial.mentions_variant(variant, gene=gene)
                if match_type == "none":
                    # Still include if it mentions the gene prominently
                    if gene.upper() not in trial.title.upper():
                        continue

            trials.append(trial)

            if len(trials) >= max_results:
                break

        return trials

    async def search_trials_by_disease(
        self,
        tumor_type: str,
        recruiting_only: bool = True,
        max_results: int = 10,
    ) -> list[ClinicalTrial]:
        """Search for clinical trials by disease/tumor type only.

        This finds disease-based trials (e.g., "Chemotherapy for NSCLC")
        that may not mention specific genes or variants.

        Args:
            tumor_type: Tumor type / disease (e.g., "NSCLC", "Melanoma")
            recruiting_only: If True, only return recruiting trials
            max_results: Maximum number of results

        Returns:
            List of ClinicalTrial objects
        """
        # Filter to last 10 years
        from datetime import datetime

        current_year = datetime.now().year
        min_year = current_year - 10

        # Build parameters - search by condition only
        params: dict[str, Any] = {
            "query.cond": tumor_type,
            "pageSize": min(max_results * 2, self.DEFAULT_PAGE_SIZE),
            "countTotal": "true",
            "format": "json",
            "fields": "|".join(self.FIELDS),
            "filter.advanced": f"AREA[StartDate]RANGE[01/01/{min_year}, MAX]",  # Started in last 10 years
        }

        # Filter by status
        if recruiting_only:
            params["filter.overallStatus"] = (
                "RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING"
            )
        else:
            params["filter.overallStatus"] = (
                "RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING|ACTIVE_NOT_RECRUITING"
            )

        # Make request
        data = await self._make_request(params)

        # Parse studies
        studies = data.get("studies", [])
        trials = []

        for study in studies:
            trial = self._parse_study(study)
            if trial is None:
                continue

            trials.append(trial)

            if len(trials) >= max_results:
                break

        return trials

    async def search_trial_evidence(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None = None,
        recruiting_only: bool = True,
        max_results: int = 10,
    ) -> list["ClinicalTrialEvidence"]:
        """Search for clinical trials and convert to evidence model.

        This method searches for trials and converts them to the
        ClinicalTrialEvidence model for use in Insight.

        Args:
            gene: Gene symbol (e.g., "KRAS")
            variant: Variant notation (e.g., "G12D")
            tumor_type: Optional tumor type filter
            recruiting_only: If True, only return recruiting trials
            max_results: Maximum number of results

        Returns:
            List of ClinicalTrialEvidence objects
        """
        from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
        from oncomind.models.evidence.base import EvidenceLevel

        trials = await self.search_trials(
            gene=gene,
            variant=variant,
            tumor_type=tumor_type,
            recruiting_only=recruiting_only,
            max_results=max_results,
        )
        evidence_list = []

        for trial in trials:
            # mentions_variant returns: (match_type, matched_biomarker)
            # match_type: "variant", "codon", "ambiguous", "gene", or "none"
            match_type, matched_biomarker = trial.mentions_variant(variant, gene=gene)

            # Determine level and scope based on match type
            if match_type == "variant":
                level = "variant"
                scope = "specific"
            elif match_type == "codon":
                level = "codon"
                scope = "specific"
            elif match_type == "ambiguous":
                level = "variant"
                scope = "ambiguous"
            elif match_type == "gene":
                level = "gene"
                scope = "unspecified"
            else:  # "none"
                level = "gene"
                scope = "unspecified"

            # Build locus_variant_match based on match type
            locus_variant_match = EvidenceLevel(
                level=level,
                scope=scope,
                origin="trial",
            )

            # Build cancer_type_match based on whether trial targets specific tumor type
            # We still fetch trials broadly, but track whether each matches the queried tumor
            cancer_type_match = None
            if tumor_type:
                cancer_matches = any(
                    tumor_types_match(condition, tumor_type)
                    for condition in (trial.conditions or [])
                )
                cancer_type_match = EvidenceLevel(
                    level="cancer_specific" if cancer_matches else None,
                    scope="specific" if cancer_matches else "unspecified",
                    origin="trial",
                )

            evidence_list.append(
                ClinicalTrialEvidence(
                    nct_id=trial.nct_id,
                    title=trial.title,
                    status=trial.status,
                    phase=trial.phase,
                    conditions=trial.conditions,
                    interventions=trial.interventions,
                    sponsor=trial.sponsor,
                    url=trial.url,
                    locus_variant_match=locus_variant_match,
                    cancer_type_match=cancer_type_match,
                    matched_biomarker=matched_biomarker,
                )
            )

        return evidence_list

    async def search_trial_evidence_by_disease(
        self,
        tumor_type: str,
        gene: str | None = None,
        variant: str | None = None,
        recruiting_only: bool = True,
        max_results: int = 10,
    ) -> list["ClinicalTrialEvidence"]:
        """Search for clinical trials by disease and convert to evidence model.

        This searches by disease/tumor type only. If gene/variant are provided,
        they are used to determine the locus_variant_match (whether trial mentions
        the biomarker), but the search itself is disease-based.

        Args:
            tumor_type: Tumor type / disease (e.g., "NSCLC", "Melanoma")
            gene: Optional gene symbol to check for mentions
            variant: Optional variant notation to check for mentions
            recruiting_only: If True, only return recruiting trials
            max_results: Maximum number of results

        Returns:
            List of ClinicalTrialEvidence objects
        """
        from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
        from oncomind.models.evidence.base import EvidenceLevel

        trials = await self.search_trials_by_disease(
            tumor_type=tumor_type,
            recruiting_only=recruiting_only,
            max_results=max_results,
        )
        evidence_list = []

        for trial in trials:
            # Determine locus_variant_match using mentions_variant
            # mentions_variant returns: (match_type, matched_biomarker)
            # match_type: "variant", "codon", "ambiguous", "gene", or "none"
            locus_variant_match = None
            matched_biomarker = None
            if gene:
                match_type, matched_biomarker = trial.mentions_variant(
                    variant, gene=gene
                )
                if match_type == "variant":
                    locus_variant_match = EvidenceLevel(
                        level="variant",
                        scope="specific",
                        origin="trial",
                    )
                elif match_type == "codon":
                    locus_variant_match = EvidenceLevel(
                        level="codon",
                        scope="specific",
                        origin="trial",
                    )
                elif match_type == "ambiguous":
                    locus_variant_match = EvidenceLevel(
                        level="variant",
                        scope="ambiguous",
                        origin="trial",
                    )
                elif match_type == "gene":
                    locus_variant_match = EvidenceLevel(
                        level="gene",
                        scope="unspecified",
                        origin="trial",
                    )
                # else: match_type == "none", locus_variant_match stays None

            # cancer_type_match - we searched by disease so it should match
            # Args: source_disease (condition from trial), queried_tumor (user's query)
            cancer_matches = any(
                tumor_types_match(condition, tumor_type)
                for condition in (trial.conditions or [])
            )

            cancer_type_match = EvidenceLevel(
                level="cancer_specific" if cancer_matches else None,
                scope="specific" if cancer_matches else "unspecified",
                origin="trial",
            )

            evidence_list.append(
                ClinicalTrialEvidence(
                    nct_id=trial.nct_id,
                    title=trial.title,
                    status=trial.status,
                    phase=trial.phase,
                    conditions=trial.conditions,
                    interventions=trial.interventions,
                    sponsor=trial.sponsor,
                    url=trial.url,
                    locus_variant_match=locus_variant_match,
                    cancer_type_match=cancer_type_match,
                    matched_biomarker=matched_biomarker,
                )
            )

        return evidence_list

    async def search_all_trial_evidence(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None = None,
        recruiting_only: bool = True,
        max_results: int = 10,
    ) -> list["ClinicalTrialEvidence"]:
        """Search for clinical trials using both biomarker and disease searches, then merge.

        This performs two searches:
        1. Gene/variant-based search (finds biomarker-specific trials)
        2. Disease-based search if tumor_type provided (finds disease-based trials)

        Results are merged and deduplicated by NCT ID.

        Args:
            gene: Gene symbol (e.g., "KRAS")
            variant: Variant notation (e.g., "G12D")
            tumor_type: Optional tumor type for disease-based search
            recruiting_only: If True, only return recruiting trials
            max_results: Maximum number of results (after merging)

        Returns:
            List of ClinicalTrialEvidence objects, deduplicated
        """
        import asyncio

        # Run both searches in parallel if tumor_type provided
        if tumor_type:
            biomarker_task = self.search_trial_evidence(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                recruiting_only=recruiting_only,
                max_results=max_results,
            )
            disease_task = self.search_trial_evidence_by_disease(
                tumor_type=tumor_type,
                gene=gene,
                variant=variant,
                recruiting_only=recruiting_only,
                max_results=max_results,
            )
            biomarker_results, disease_results = await asyncio.gather(
                biomarker_task, disease_task
            )
        else:
            # No tumor_type, only do biomarker search
            biomarker_results = await self.search_trial_evidence(
                gene=gene,
                variant=variant,
                tumor_type=tumor_type,
                recruiting_only=recruiting_only,
                max_results=max_results,
            )
            disease_results = []

        # Merge and deduplicate by NCT ID
        # Prefer biomarker results (they have more specific locus_variant_match)
        seen_nct_ids: set[str] = set()
        merged: list["ClinicalTrialEvidence"] = []

        # Add biomarker results first (higher priority)
        for trial in biomarker_results:
            if trial.nct_id not in seen_nct_ids:
                seen_nct_ids.add(trial.nct_id)
                merged.append(trial)

        # Add disease results that weren't in biomarker results
        for trial in disease_results:
            if trial.nct_id not in seen_nct_ids:
                seen_nct_ids.add(trial.nct_id)
                merged.append(trial)

        # Limit to max_results
        return merged[:max_results]

    async def close(self) -> None:
        """Close the HTTP client."""
        if self._client:
            await self._client.aclose()
            self._client = None

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
        return self.status in ['RECRUITING', 'ENROLLING_BY_INVITATION']

    def is_active(self) -> bool:
        """Check if trial is active (recruiting or not)."""
        return self.status in [
            'RECRUITING',
            'ENROLLING_BY_INVITATION',
            'ACTIVE_NOT_RECRUITING',
            'NOT_YET_RECRUITING'
        ]

    def mentions_variant(self, variant: str | None, gene: str | None = None) -> tuple[str, str | None]:
        """Check if trial mentions the gene and/or variant.

        Args:
            variant: Variant notation (e.g., "G12D", "V600E"). Can be None to check gene only.
            gene: Gene symbol (e.g., "NRAS", "BRAF"). If provided, ensures
                  the variant is mentioned in context of this gene to avoid
                  false positives (e.g., KRAS G12D trial matching NRAS G12D query).

        Returns:
            Tuple of (match_type, matched_biomarker):
            - ("specific", "KRAS G12D") - both gene and variant found in text
            - ("ambiguous", "KRAS G12") - gene found and ambiguous variant in BROAD_VARIANTS
            - ("gene", "KRAS") - only gene found in text
            - ("none", None) - no match found
        """
        from oncomind.config.constants import BROAD_VARIANTS

        gene_upper = gene.upper() if gene else None
        variant_upper = variant.upper() if variant else None

        # Combine all text to search
        search_texts = [self.title.upper()]
        if self.eligibility_criteria:
            search_texts.append(self.eligibility_criteria.upper())
        if self.brief_summary:
            search_texts.append(self.brief_summary.upper())

        full_text = " ".join(search_texts)

        gene_found = gene_upper and gene_upper in full_text
        variant_found = variant_upper and variant_upper in full_text

        # Find ambiguous variants where both gene and variant are in the text
        ambig_variants = [
            (g, v) for (g, v) in BROAD_VARIANTS
            if g.upper() in full_text and v.upper() in full_text
        ]
        # Find the ambiguous variant that matches the queried gene
        matched_ambig = next(
            ((g, v) for (g, v) in ambig_variants if g.upper() == gene_upper),
            None
        )

        if variant_found:
            biomarker = f"{gene} {variant}" if gene else variant
            return ('specific', biomarker)
        elif matched_ambig:
            g, v = matched_ambig
            return ('ambiguous', f"{g} {v}")
        elif gene_found:
            return ('gene', gene)
        else:
            return ('none', None)

       

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            'nct_id': self.nct_id,
            'title': self.title,
            'status': self.status,
            'phase': self.phase,
            'conditions': self.conditions,
            'interventions': self.interventions,
            'brief_summary': self.brief_summary[:500] if self.brief_summary else None,
            'sponsor': self.sponsor,
            'url': self.url,
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
            self._client = httpx.AsyncClient(timeout=self.timeout, headers=self._headers)
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
            variant_clean = variant.upper().replace('P.', '')
            parts.append(variant_clean)
            # Also try gene+variant concatenated (e.g., "KRASG12D")
            parts.append(f"{gene.upper()}{variant_clean}")

        # Join with OR for broader matching
        query = ' OR '.join(parts)

        return query

    def _parse_study(self, study: dict[str, Any]) -> ClinicalTrial | None:
        """Parse a study from the API response.

        Args:
            study: Raw study data from API

        Returns:
            ClinicalTrial object or None if parsing fails
        """
        try:
            protocol = study.get('protocolSection', {})

            # Identification
            id_module = protocol.get('identificationModule', {})
            nct_id = id_module.get('nctId', '')
            title = id_module.get('briefTitle', '') or id_module.get('officialTitle', '')

            # Status
            status_module = protocol.get('statusModule', {})
            status = status_module.get('overallStatus', 'UNKNOWN')

            # Design - phase
            design_module = protocol.get('designModule', {})
            phases = design_module.get('phases', [])
            phase = phases[0] if phases else None

            # Conditions
            conditions_module = protocol.get('conditionsModule', {})
            conditions = conditions_module.get('conditions', [])

            # Interventions
            arms_module = protocol.get('armsInterventionsModule', {})
            interventions_list = arms_module.get('interventions', [])
            interventions = [
                i.get('name', '')
                for i in interventions_list
                if i.get('name')
            ]

            # Description
            desc_module = protocol.get('descriptionModule', {})
            brief_summary = desc_module.get('briefSummary', '')

            # Eligibility
            eligibility_module = protocol.get('eligibilityModule', {})
            eligibility_criteria = eligibility_module.get('eligibilityCriteria', '')

            # Sponsor
            sponsor_module = protocol.get('sponsorCollaboratorsModule', {})
            lead_sponsor = sponsor_module.get('leadSponsor', {})
            sponsor = lead_sponsor.get('name', '')

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

        # Build parameters - request only needed fields to reduce payload
        params: dict[str, Any] = {
            'query.term': query,
            'pageSize': min(max_results * 2, self.DEFAULT_PAGE_SIZE),
            'countTotal': 'true',
            'format': 'json',
            'fields': '|'.join(self.FIELDS),  # Only fetch fields we need
        }

        # Add condition filter if tumor type specified
        if tumor_type:
            params['query.cond'] = tumor_type

        # Filter by status if recruiting only
        if recruiting_only:
            params['filter.overallStatus'] = 'RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING'
        else:
            params['filter.overallStatus'] = 'RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING|ACTIVE_NOT_RECRUITING'

        # Make request - let errors propagate for caller to handle retry
        data = await self._make_request(params)

        # Parse studies
        studies = data.get('studies', [])
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
        # Build parameters - search by condition only
        params: dict[str, Any] = {
            'query.cond': tumor_type,
            'pageSize': min(max_results * 2, self.DEFAULT_PAGE_SIZE),
            'countTotal': 'true',
            'format': 'json',
            'fields': '|'.join(self.FIELDS),
        }

        # Filter by status
        if recruiting_only:
            params['filter.overallStatus'] = 'RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING'
        else:
            params['filter.overallStatus'] = 'RECRUITING|ENROLLING_BY_INVITATION|NOT_YET_RECRUITING|ACTIVE_NOT_RECRUITING'

        # Make request
        data = await self._make_request(params)

        # Parse studies
        studies = data.get('studies', [])
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
            match_type, matched_biomarker = trial.mentions_variant(variant, gene=gene)

            # Determine level and scope based on match type
            if match_type == "specific":
                level = "variant"
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

            # Build variant_level based on match type
            variant_level = EvidenceLevel(
                level=level,
                scope=scope,
                origin="trial",
            )

            # Build cancer_type_level based on whether trial targets specific tumor type
            # Only set if tumor_type was queried; otherwise leave as None (unknown)
            cancer_type_level = None
            if tumor_type:
                tumor_type_lower = tumor_type.lower()
                cancer_matches = False
                if trial.conditions:
                    for condition in trial.conditions:
                        if tumor_type_lower in condition.lower() or condition.lower() in tumor_type_lower:
                            cancer_matches = True
                            break

                cancer_type_level = EvidenceLevel(
                    level="cancer_specific" if cancer_matches else "pan_cancer",
                    scope="specific" if cancer_matches else "unspecified",
                    origin="trial",
                )

            evidence_list.append(ClinicalTrialEvidence(
                nct_id=trial.nct_id,
                title=trial.title,
                status=trial.status,
                phase=trial.phase,
                conditions=trial.conditions,
                interventions=trial.interventions,
                sponsor=trial.sponsor,
                url=trial.url,
                variant_level=variant_level,
                cancer_type_level=cancer_type_level,
                matched_biomarker=matched_biomarker,
            ))

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
        they are used to determine the variant_level (whether trial mentions
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
            # Determine variant_level using mentions_variant
            # mentions_variant returns: (match_type, matched_biomarker)
            variant_level = None
            matched_biomarker = None
            if gene:
                match_type, matched_biomarker = trial.mentions_variant(variant, gene=gene)
                if match_type == "specific":
                    variant_level = EvidenceLevel(
                        level="variant",
                        scope="specific",
                        origin="trial",
                    )
                elif match_type == "ambiguous":
                    variant_level = EvidenceLevel(
                        level="variant",
                        scope="ambiguous",
                        origin="trial",
                    )
                elif match_type == "gene":
                    variant_level = EvidenceLevel(
                        level="gene",
                        scope="unspecified",
                        origin="trial",
                    )
                # else: match_type == "none", variant_level stays None

            # cancer_type_level - we searched by disease so it should match
            tumor_type_lower = tumor_type.lower()
            cancer_matches = False
            if trial.conditions:
                for condition in trial.conditions:
                    if tumor_type_lower in condition.lower() or condition.lower() in tumor_type_lower:
                        cancer_matches = True
                        break

            cancer_type_level = EvidenceLevel(
                level="cancer_specific" if cancer_matches else "pan_cancer",
                scope="specific" if cancer_matches else "unspecified",
                origin="trial",
            )

            evidence_list.append(ClinicalTrialEvidence(
                nct_id=trial.nct_id,
                title=trial.title,
                status=trial.status,
                phase=trial.phase,
                conditions=trial.conditions,
                interventions=trial.interventions,
                sponsor=trial.sponsor,
                url=trial.url,
                variant_level=variant_level,
                cancer_type_level=cancer_type_level,
                matched_biomarker=matched_biomarker,
            ))

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
        # Prefer biomarker results (they have more specific variant_level)
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

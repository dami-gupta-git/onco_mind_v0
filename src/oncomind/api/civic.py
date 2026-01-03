"""CIViC (Clinical Interpretation of Variants in Cancer) GraphQL API client.

ARCHITECTURE:
    Gene + Variant → CIViC GraphQL API → AMP/ASCO/CAP Assertions with NCCN guidelines

CIViC provides curated clinical interpretations with:
- AMP/ASCO/CAP tier classifications (Tier I/II Level A/B/C/D)
- FDA companion test status
- NCCN guideline references
- Assertion types: PREDICTIVE, PROGNOSTIC, DIAGNOSTIC, ONCOGENIC

Key Design:
- GraphQL API for flexible querying
- Assertions are curated summaries with AMP tier assignments
- Complements CGI and VICC by providing guideline-backed tiers
- Free and open source (no license required unlike OncoKB)
"""

from typing import Any

import httpx

from oncomind.config.constants import TUMOR_TYPE_MAPPINGS
from oncomind.config.debug import get_logger

logger = get_logger(__name__)


class CIViCError(Exception):
    """Exception raised for CIViC API errors."""
    pass


class CIViCAssertion:
    """A curated assertion from CIViC database.

    Assertions represent the clinical significance of a molecular profile
    in a specific disease context, with AMP/ASCO/CAP tier classification.
    """

    def __init__(
        self,
        assertion_id: int,
        name: str,
        amp_level: str | None,
        assertion_type: str,  # PREDICTIVE, PROGNOSTIC, DIAGNOSTIC, ONCOGENIC
        assertion_direction: str,  # SUPPORTS, DOES_NOT_SUPPORT
        significance: str,  # SENSITIVITYRESPONSE, RESISTANCE, ONCOGENIC, etc.
        status: str,  # ACCEPTED, SUBMITTED, REJECTED
        molecular_profile: str,
        disease: str,
        therapies: list[str],
        fda_companion_test: bool | None,
        nccn_guideline: str | None,
        description: str | None = None,
    ):
        self.assertion_id = assertion_id
        self.name = name
        self.amp_level = amp_level
        self.assertion_type = assertion_type
        self.assertion_direction = assertion_direction
        self.significance = significance
        self.status = status
        self.molecular_profile = molecular_profile
        self.disease = disease
        self.therapies = therapies
        self.fda_companion_test = fda_companion_test
        self.nccn_guideline = nccn_guideline
        self.description = description

    def get_amp_tier(self) -> str | None:
        """Extract AMP tier from amp_level (e.g., 'TIER_I_LEVEL_A' -> 'Tier I')."""
        if not self.amp_level:
            return None
        if "TIER_I" in self.amp_level:
            return "Tier I"
        elif "TIER_II" in self.amp_level:
            return "Tier II"
        elif "TIER_III" in self.amp_level:
            return "Tier III"
        elif "TIER_IV" in self.amp_level:
            return "Tier IV"
        return None

    def get_amp_level(self) -> str | None:
        """Extract AMP level from amp_level (e.g., 'TIER_I_LEVEL_A' -> 'A')."""
        if not self.amp_level:
            return None
        if "LEVEL_A" in self.amp_level:
            return "A"
        elif "LEVEL_B" in self.amp_level:
            return "B"
        elif "LEVEL_C" in self.amp_level:
            return "C"
        elif "LEVEL_D" in self.amp_level:
            return "D"
        return None

    def is_sensitivity(self) -> bool:
        """Check if this represents a sensitivity/response assertion."""
        if not self.significance:
            return False
        sig_upper = self.significance.upper()
        return any(term in sig_upper for term in ["SENSITIV", "RESPONSE"])

    def is_resistance(self) -> bool:
        """Check if this represents a resistance assertion."""
        if not self.significance:
            return False
        return "RESIST" in self.significance.upper()

    def is_accepted(self) -> bool:
        """Check if assertion has been accepted (vs submitted/rejected)."""
        return self.status == "ACCEPTED"

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "assertion_id": self.assertion_id,
            "name": self.name,
            "amp_level": self.amp_level,
            "amp_tier": self.get_amp_tier(),
            "amp_level_letter": self.get_amp_level(),
            "assertion_type": self.assertion_type,
            "assertion_direction": self.assertion_direction,
            "significance": self.significance,
            "status": self.status,
            "molecular_profile": self.molecular_profile,
            "disease": self.disease,
            "therapies": self.therapies,
            "fda_companion_test": self.fda_companion_test,
            "nccn_guideline": self.nccn_guideline,
            "description": self.description,
            "is_sensitivity": self.is_sensitivity(),
            "is_resistance": self.is_resistance(),
        }


class CIViCClient:
    """Client for CIViC GraphQL API.

    CIViC (Clinical Interpretation of Variants in Cancer) provides
    curated clinical interpretations with AMP/ASCO/CAP tier assignments.

    API Documentation: https://griffithlab.github.io/civic-v2/
    GraphQL Endpoint: https://civicdb.org/api/graphql
    """

    GRAPHQL_URL = "https://civicdb.org/api/graphql"
    DEFAULT_TIMEOUT = 30.0

    # GraphQL query for assertions
    ASSERTIONS_QUERY = """
    query GetAssertions($molecularProfileName: String, $first: Int) {
        assertions(molecularProfileName: $molecularProfileName, first: $first) {
            nodes {
                id
                name
                ampLevel
                assertionType
                assertionDirection
                significance
                status
                description
                therapies {
                    name
                }
                disease {
                    name
                }
                molecularProfile {
                    name
                }
                fdaCompanionTest
                nccnGuideline {
                    name
                }
            }
        }
    }
    """

    def __init__(self, timeout: float = DEFAULT_TIMEOUT):
        """Initialize the CIViC client.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout
        self._client: httpx.AsyncClient | None = None

    async def __aenter__(self):
        """Initialize HTTP client session."""
        self._client = httpx.AsyncClient(timeout=self.timeout)
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Close HTTP client session."""
        if self._client:
            await self._client.aclose()
            self._client = None

    def _get_client(self) -> httpx.AsyncClient:
        """Get the HTTP client, creating one if needed."""
        if self._client is None:
            self._client = httpx.AsyncClient(timeout=self.timeout)
        return self._client

    def _determine_match_level(
        self,
        molecular_profile: str,
        gene: str,
        variant: str | None,
    ) -> str:
        """Determine the match specificity level for a molecular profile.

        Args:
            molecular_profile: The molecular profile string from CIViC (e.g., "EGFR L858R")
            gene: The queried gene symbol
            variant: The queried variant (e.g., "L858R")

        Returns:
            Match level: 'variant' (exact), 'codon' (same position), or 'gene' (gene-only)
        """
        import re

        profile_upper = molecular_profile.upper() if molecular_profile else ""
        gene_upper = gene.upper()

        # Check if gene is even in the profile
        if gene_upper not in profile_upper:
            return "gene"  # Shouldn't happen but fallback

        if not variant:
            return "gene"

        # Clean variant for comparison
        clean_variant = variant.replace("p.", "").upper()

        # Exact variant match
        if clean_variant in profile_upper:
            return "variant"

        # Check for codon-level match (same position, different amino acid change)
        # Extract position from variant (e.g., L858R -> 858, V600E -> 600)
        pos_match = re.search(r'[A-Z](\d+)', clean_variant)
        if pos_match:
            position = pos_match.group(1)
            # Check if profile contains same position with any amino acid
            codon_pattern = rf'[A-Z]{position}[A-Z]?'
            if re.search(codon_pattern, profile_upper):
                # Make sure it's actually a different variant at same position
                return "codon"

        return "gene"

    def _is_pan_cancer(self, disease: str | None) -> bool:
        """Check if disease is a generic pan-cancer term.

        Pan-cancer diseases like "Cancer" or "Solid Tumor" apply broadly
        and should not be considered a specific tumor type match.

        Args:
            disease: Disease string from CIViC

        Returns:
            True if disease is a generic pan-cancer term
        """
        if not disease:
            return False

        disease_lower = disease.lower().strip()

        # Generic pan-cancer terms
        pan_cancer_terms = {
            "cancer",
            "solid tumor",
            "solid tumors",
            "solid tumour",
            "solid tumours",
            "advanced solid tumor",
            "advanced solid tumors",
            "malignant neoplasm",
            "malignant neoplasms",
            "neoplasm",
            "tumor",
            "tumour",
        }

        return disease_lower in pan_cancer_terms

    def _tumor_matches(self, civic_disease: str, tumor_type: str | None) -> bool:
        """Check if CIViC disease matches user tumor type.

        Args:
            civic_disease: Disease string from CIViC
            tumor_type: User-provided tumor type

        Returns:
            True if tumor types match. Returns False for pan-cancer diseases
            when a specific tumor type is queried (they should be included
            but marked as disease_match=False).
        """
        if not tumor_type:
            return True  # No filter

        # Pan-cancer diseases don't count as a specific tumor match
        if self._is_pan_cancer(civic_disease):
            return False

        civic_lower = civic_disease.lower() if civic_disease else ""
        tumor_lower = tumor_type.lower()

        # Direct substring match (but not if civic_disease is too generic)
        # e.g., "breast cancer" in "her2-negative breast cancer" = True
        # but "cancer" in "breast cancer" would be caught by _is_pan_cancer above
        if tumor_lower in civic_lower or civic_lower in tumor_lower:
            return True

        # Check tumor type mappings
        for abbrev, full_names in TUMOR_TYPE_MAPPINGS.items():
            # Check if tumor matches this mapping (either as abbrev or substring match)
            tumor_matches_mapping = (
                tumor_lower == abbrev or
                any(tumor_lower in name for name in full_names) or
                any(name in tumor_lower for name in full_names)
            )
            if tumor_matches_mapping:
                # Check if civic disease matches any full name
                if any(name in civic_lower for name in full_names):
                    return True

        return False

    def _parse_assertion(self, node: dict[str, Any]) -> CIViCAssertion | None:
        """Parse a GraphQL assertion node into an assertion object.

        Args:
            node: Raw assertion node from GraphQL response

        Returns:
            CIViCAssertion or None if parsing fails
        """
        try:
            # Extract therapies
            therapies = []
            for therapy in node.get("therapies", []):
                if therapy.get("name"):
                    therapies.append(therapy["name"])

            # Extract disease
            disease = node.get("disease", {}).get("name", "")

            # Extract molecular profile
            molecular_profile = node.get("molecularProfile", {}).get("name", "")

            # Extract NCCN guideline
            nccn = node.get("nccnGuideline")
            nccn_guideline = nccn.get("name") if nccn else None

            return CIViCAssertion(
                assertion_id=node.get("id"),
                name=node.get("name", ""),
                amp_level=node.get("ampLevel"),
                assertion_type=node.get("assertionType", ""),
                assertion_direction=node.get("assertionDirection", ""),
                significance=node.get("significance", ""),
                status=node.get("status", ""),
                molecular_profile=molecular_profile,
                disease=disease,
                therapies=therapies,
                fda_companion_test=node.get("fdaCompanionTest"),
                nccn_guideline=nccn_guideline,
                description=node.get("description"),
            )

        except Exception:
            return None

    async def fetch_assertions(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        max_results: int = 50,
    ) -> list[CIViCAssertion]:
        """Fetch CIViC assertions for a gene/variant combination.

        Args:
            gene: Gene symbol (e.g., "EGFR")
            variant: Optional variant notation (e.g., "L858R")
            tumor_type: Optional tumor type to filter results
            max_results: Maximum number of results to return

        Returns:
            List of CIViCAssertion objects
        """
        client = self._get_client()

        # Build search term for molecular profile
        search_term = gene.upper()
        if variant:
            clean_variant = variant.replace("p.", "").upper()
            search_term = f"{gene.upper()} {clean_variant}"

        # Execute GraphQL query
        variables = {
            "molecularProfileName": search_term,
            "first": max_results * 2,  # Fetch extra to account for filtering
        }

        try:
            response = await client.post(
                self.GRAPHQL_URL,
                json={
                    "query": self.ASSERTIONS_QUERY,
                    "variables": variables,
                },
            )
            response.raise_for_status()
            data = response.json()
        except httpx.HTTPError as e:
            logger.error(f"CIViC API request failed: {e}")
            raise CIViCError(f"CIViC API request failed: {e}")
        except Exception as e:
            logger.error(f"Failed to parse CIViC response: {e}")
            raise CIViCError(f"Failed to parse CIViC response: {e}")

        # Check for GraphQL errors
        if "errors" in data:
            raise CIViCError(f"GraphQL errors: {data['errors']}")

        # Parse assertions
        assertions = []
        nodes = data.get("data", {}).get("assertions", {}).get("nodes", [])

        for node in nodes:
            assertion = self._parse_assertion(node)
            if assertion is None:
                continue

            # Filter by tumor type if specified
            # Include pan-cancer assertions (e.g., "Cancer") but filter out non-matching specific cancers
            if tumor_type and assertion.disease:
                is_pan_cancer = self._is_pan_cancer(assertion.disease)
                tumor_matches = self._tumor_matches(assertion.disease, tumor_type)
                if not is_pan_cancer and not tumor_matches:
                    continue

            # Check if molecular profile contains our variant
            if variant:
                clean_variant = variant.replace("p.", "").upper()
                if clean_variant not in assertion.molecular_profile.upper():
                    continue

            assertions.append(assertion)

            if len(assertions) >= max_results:
                break

        return assertions

    async def fetch_predictive_assertions(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        max_results: int = 25,
    ) -> list[CIViCAssertion]:
        """Fetch only PREDICTIVE assertions (therapy response).

        Args:
            gene: Gene symbol
            variant: Optional variant notation
            tumor_type: Optional tumor type filter
            max_results: Maximum results

        Returns:
            List of predictive assertions
        """
        all_assertions = await self.fetch_assertions(gene, variant, tumor_type, max_results * 2)
        return [a for a in all_assertions if a.assertion_type == "PREDICTIVE"][:max_results]

    async def fetch_tier_i_assertions(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
    ) -> list[CIViCAssertion]:
        """Fetch only Tier I (Level A or B) assertions.

        These represent the strongest clinical evidence for actionability.

        Args:
            gene: Gene symbol
            variant: Optional variant notation
            tumor_type: Optional tumor type filter

        Returns:
            List of Tier I assertions
        """
        all_assertions = await self.fetch_assertions(gene, variant, tumor_type, max_results=50)
        return [a for a in all_assertions if a.get_amp_tier() == "Tier I"]

    # GraphQL query for evidence items
    EVIDENCE_ITEMS_QUERY = """
    query GetEvidenceItems($molecularProfileName: String, $first: Int) {
        molecularProfiles(name: $molecularProfileName, first: 20) {
            nodes {
                id
                name
                evidenceItems(first: $first) {
                    nodes {
                        id
                        evidenceType
                        evidenceLevel
                        evidenceDirection
                        significance
                        status
                        description
                        evidenceRating
                        therapies {
                            name
                        }
                        disease {
                            name
                        }
                        source {
                            id
                            sourceType
                            citation
                            pmid: citationId
                        }
                    }
                }
            }
        }
    }
    """

    async def fetch_evidence_items(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        max_results: int = 50,
    ) -> list["CIViCEvidence"]:
        """Fetch CIViC evidence items (EIDs) directly from GraphQL API.

        Evidence items are the individual pieces of evidence (linked to publications)
        that support clinical interpretations. Each has a unique EID.

        Args:
            gene: Gene symbol (e.g., "EGFR")
            variant: Optional variant notation (e.g., "L858R")
            tumor_type: Optional tumor type to filter results
            max_results: Maximum number of results to return

        Returns:
            List of CIViCEvidence objects
        """
        from oncomind.models.evidence.civic import CIViCEvidence

        client = self._get_client()

        # Build search terms - query specific variant and codon-level
        # Note: We intentionally do NOT search for "{gene} MUTATION" because:
        # 1. Gene-level evidence adds noise to the UI
        # 2. The LLM only uses civic_assertions (not civic_evidence)
        # 3. Gap detector can rely on assertions + VICC for gene-level signals
        search_terms = []
        gene_upper = gene.upper()

        if variant:
            clean_variant = variant.replace("p.", "").upper()
            search_terms.append(f"{gene_upper} {clean_variant}")

            # Also try codon-level (e.g., V600E -> V600)
            import re
            codon_match = re.match(r'^([A-Z])(\d+)[A-Z]*$', clean_variant)
            if codon_match:
                codon_variant = f"{codon_match.group(1)}{codon_match.group(2)}"
                if codon_variant != clean_variant:
                    search_terms.append(f"{gene_upper} {codon_variant}")
        else:
            # Gene-only query (no variant specified)
            search_terms.append(gene_upper)

        evidence_list = []
        seen_evidence_ids = set()

        for search_term in search_terms:
            variables = {
                "molecularProfileName": search_term,
                "first": max_results * 2,
            }

            try:
                response = await client.post(
                    self.GRAPHQL_URL,
                    json={
                        "query": self.EVIDENCE_ITEMS_QUERY,
                        "variables": variables,
                    },
                )
                response.raise_for_status()
                data = response.json()
            except httpx.HTTPError as e:
                logger.warning(f"CIViC evidence query failed for {search_term}: {e}")
                continue
            except Exception as e:
                logger.warning(f"Failed to parse CIViC evidence response for {search_term}: {e}")
                continue

            # Check for GraphQL errors
            if "errors" in data:
                logger.warning(f"GraphQL errors for {search_term}: {data['errors']}")
                continue

            # Parse evidence items from molecular profiles
            profiles = data.get("data", {}).get("molecularProfiles", {}).get("nodes", [])

            for profile in profiles:
                profile_name = profile.get("name", "")
                evidence_nodes = profile.get("evidenceItems", {}).get("nodes", [])

                for node in evidence_nodes:
                    evidence_id = node.get("id")

                    # Skip duplicates
                    if evidence_id in seen_evidence_ids:
                        continue
                    seen_evidence_ids.add(evidence_id)

                    # Extract disease
                    disease = node.get("disease", {}).get("name") if node.get("disease") else None

                    # Filter by tumor type if specified
                    # Include pan-cancer evidence (e.g., "Cancer") but filter out non-matching specific cancers
                    if tumor_type and disease:
                        is_pan_cancer = self._is_pan_cancer(disease)
                        tumor_matches = self._tumor_matches(disease, tumor_type)
                        if not is_pan_cancer and not tumor_matches:
                            continue

                    # Extract therapies/drugs
                    drugs = []
                    for therapy in node.get("therapies", []):
                        if therapy.get("name"):
                            drugs.append(therapy["name"])

                    # Extract source info
                    source_data = node.get("source", {})
                    source = source_data.get("citation") if source_data else None
                    pmid = source_data.get("pmid") if source_data else None

                    # Determine match level
                    match_level = self._determine_match_level(profile_name, gene, variant)

                    # Determine disease match
                    disease_match = self._tumor_matches(disease, tumor_type) if tumor_type and disease else True

                    evidence_list.append(CIViCEvidence(
                        evidence_id=evidence_id,
                        evidence_type=node.get("evidenceType"),
                        evidence_level=node.get("evidenceLevel"),
                        evidence_direction=node.get("evidenceDirection"),
                        clinical_significance=node.get("significance"),
                        disease=disease,
                        drugs=drugs,
                        description=node.get("description"),
                        source=source,
                        rating=node.get("evidenceRating"),
                        pmid=str(pmid) if pmid else None,
                        source_url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None,
                        trust_rating=node.get("evidenceRating"),
                        match_level=match_level,
                        matched_profile=profile_name,
                        disease_match=disease_match,
                    ))

                    if len(evidence_list) >= max_results:
                        return evidence_list

        return evidence_list

    async def fetch_assertion_evidence(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        max_results: int = 50,
    ) -> list["CIViCAssertionEvidence"]:
        """Fetch CIViC assertions and convert to evidence model.

        This method fetches assertions and converts them to the
        CIViCAssertionEvidence model for use in Insight.

        Args:
            gene: Gene symbol (e.g., "EGFR")
            variant: Optional variant notation (e.g., "L858R")
            tumor_type: Optional tumor type to filter results
            max_results: Maximum number of results to return

        Returns:
            List of CIViCAssertionEvidence objects
        """
        from oncomind.models.evidence.civic import CIViCAssertionEvidence

        assertions = await self.fetch_assertions(gene, variant, tumor_type, max_results)
        evidence_list = []

        for assertion in assertions:
            # Determine match specificity
            match_level = self._determine_match_level(
                assertion.molecular_profile,
                gene,
                variant,
            )
            # Check if disease matches (tumor type filtering already applied, but track it)
            disease_match = self._tumor_matches(assertion.disease, tumor_type) if tumor_type else True

            evidence_list.append(CIViCAssertionEvidence(
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
                match_level=match_level,
                matched_profile=assertion.molecular_profile,
                disease_match=disease_match,
            ))

        return evidence_list

    async def close(self) -> None:
        """Close the HTTP client."""
        if self._client:
            await self._client.aclose()
            self._client = None

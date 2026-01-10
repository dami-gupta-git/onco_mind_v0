"""VICC MetaKB client for harmonized cancer variant interpretations.

ARCHITECTURE:
    Gene + Variant → VICC MetaKB API → Harmonized evidence from CIViC, CGI, JAX, OncoKB, PMKB

The VICC (Variant Interpretation for Cancer Consortium) Meta-Knowledgebase aggregates
and harmonizes clinical interpretations from multiple cancer variant knowledgebases:
- CIViC (Clinical Interpretations of Variants in Cancer)
- Cancer Genome Interpreter (CGI)
- JAX-CKB (Jackson Laboratory Clinical Knowledgebase)
- OncoKB
- PMKB (Precision Medicine Knowledgebase)
- MolecularMatch

Key Design:
- Lucene query syntax for flexible variant matching
- Evidence levels: A (validated), B (clinical), C (case study), D (preclinical)
- Response types: Responsive/Sensitivity, Resistant, 1A/1B/2A/2B/etc (OncoKB-style)
- Sources attributed to original KB for provenance tracking
"""
from typing import Any

import httpx

from oncomind.models.evidence.base import is_pan_cancer_term, tumor_types_match, determine_locus_match


class VICCError(Exception):
    """Exception raised for VICC MetaKB-related errors."""

    pass


class VICCAssociation:
    """A clinical interpretation association from VICC MetaKB.

    Represents harmonized evidence linking a variant to a drug/disease with
    clinical significance and evidence level.
    """

    def __init__(
        self,
        description: str,
        gene: str,
        variant: str | None,
        disease: str,
        drugs: list[str],
        evidence_level: str | None,  # A, B, C, D
        response_type: str | None,  # Responsive, Resistant, Sensitivity, 1A, etc.
        source: str,  # civic, cgi, jax, oncokb, pmkb
        publication_url: str | list[str] | None = None,  # Can be string or list
        oncogenic: str | None = None,
    ):
        self.description = description
        self.gene = gene
        self.variant = variant
        self.disease = disease
        self.drugs = drugs
        self.evidence_level = evidence_level
        self.response_type = response_type
        self.source = source
        self.publication_url = publication_url
        self.oncogenic = oncogenic

    def is_sensitivity(self) -> bool:
        """Check if this represents a sensitivity/response association."""
        if not self.response_type:
            return False
        rt_upper = self.response_type.upper()
        return any(term in rt_upper for term in ["SENSITIV", "RESPONSE", "RESPONSIVE"])

    def is_resistance(self) -> bool:
        """Check if this represents a resistance association."""
        if not self.response_type:
            return False
        return "RESIST" in self.response_type.upper()

    def get_oncokb_level(self) -> str | None:
        """Extract OncoKB-style level if present (1A, 1B, 2A, 2B, 3A, 3B, 4, R1, R2)."""
        if not self.response_type:
            return None
        # OncoKB uses levels like 1A, 1B, 2A, 2B, 3A, 3B, 4, R1, R2
        import re
        match = re.match(r'^([1234][AB]?|R[12])$', self.response_type.upper())
        if match:
            return match.group(1)
        return None

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "description": self.description,
            "gene": self.gene,
            "variant": self.variant,
            "disease": self.disease,
            "drugs": self.drugs,
            "evidence_level": self.evidence_level,
            "response_type": self.response_type,
            "source": self.source,
            "publication_url": self.publication_url,
            "oncogenic": self.oncogenic,
            "is_sensitivity": self.is_sensitivity(),
            "is_resistance": self.is_resistance(),
            "oncokb_level": self.get_oncokb_level(),
        }


class VICCClient:
    """Client for VICC MetaKB API.

    The VICC Meta-Knowledgebase provides harmonized clinical interpretations
    from multiple cancer variant databases, enabling comprehensive evidence
    lookup across the ecosystem.

    API Documentation: https://search.cancervariants.org/api/v1/ui/
    """

    BASE_URL = "https://search.cancervariants.org/api/v1"
    DEFAULT_TIMEOUT = 30.0
    DEFAULT_SIZE = 50  # Number of results per query

    def __init__(self, timeout: float = DEFAULT_TIMEOUT):
        """Initialize the VICC client.

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

    def _get_kit_exon(self, variant: str) -> int | None:
        """Get the KIT exon number from a variant string.

        KIT exon boundaries (codon ranges):
        - Exon 9: codons 503-510 (extracellular domain)
        - Exon 11: codons 550-591 (juxtamembrane domain) - most common in GIST
        - Exon 13: codons 627-664 (ATP binding)
        - Exon 17: codons 788-828 (activation loop)

        Args:
            variant: Variant string (e.g., "W557_K558del", "V560D")

        Returns:
            Exon number or None if not determinable
        """
        import re
        # Extract first position from variant (e.g., W557_K558del → 557, V560D → 560)
        pos_match = re.search(r'[A-Z](\d+)', variant.upper())
        if pos_match:
            position = int(pos_match.group(1))
            if 503 <= position <= 510:
                return 9
            elif 550 <= position <= 591:
                return 11
            elif 627 <= position <= 664:
                return 13
            elif 788 <= position <= 828:
                return 17
        return None

    def _build_query(
        self, gene: str, variant: str | None = None, tumor_type: str | None = None
    ) -> str:
        """Build a Lucene query string for the VICC API.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional variant notation (e.g., "V600E")
            tumor_type: Optional tumor type (e.g., "breast cancer")

        Returns:
            Lucene query string
        """
        # Gene symbol search in feature names
        query_parts = [gene.upper()]

        if variant:
            # Remove p. prefix if present
            clean_variant = variant.replace("p.", "").upper()
            query_parts.append(clean_variant)

        if tumor_type:
            # Normalize tumor type for query
            # Use quotes for multi-word tumor types to ensure phrase matching
            tumor_normalized = tumor_type.strip().lower()
            # VICC uses disease field - add tumor type as quoted phrase
            query_parts.append(f'"{tumor_normalized}"')

        return " AND ".join(query_parts)

    def _build_exon_query(
        self, gene: str, exon: int, tumor_type: str | None = None
    ) -> str:
        """Build a Lucene query for exon-level search.

        Args:
            gene: Gene symbol (e.g., "KIT")
            exon: Exon number (e.g., 11)
            tumor_type: Optional tumor type (e.g., "GIST")

        Returns:
            Lucene query string for exon-level search
        """
        query = f'{gene.upper()} AND "exon {exon}"'
        if tumor_type:
            tumor_normalized = tumor_type.strip().lower()
            query += f' AND "{tumor_normalized}"'
        return query

    def _determine_locus_match(
        self,
        assoc_variant: str | None,
        queried_gene: str,
        queried_variant: str | None,
    ) -> str:
        """Determine the match specificity level for a VICC association.

        Delegates to the core determine_locus_match() function.

        Args:
            assoc_variant: The variant from the VICC association (e.g., "V600E")
            queried_gene: The gene being queried (unused, kept for API compatibility)
            queried_variant: The variant being queried (e.g., "V600E")

        Returns:
            Match level: 'variant' (exact), 'codon' (same position), or 'gene' (gene-only)
        """
        return determine_locus_match(assoc_variant, queried_variant)

    def _parse_association(self, hit: dict[str, Any]) -> VICCAssociation | None:
        """Parse a VICC API hit into an association object.

        Args:
            hit: Raw hit from VICC API

        Returns:
            VICCAssociation or None if parsing fails
        """
        try:
            assoc = hit.get("association", {})
            features = hit.get("features", [])

            # Extract gene and variant from features
            gene = ""
            variant = None
            for feature in features:
                if feature.get("geneSymbol"):
                    gene = feature["geneSymbol"]
                if feature.get("name"):
                    # VICC returns variant in 'name' field (e.g., "V600E", not "BRAF V600E")
                    name = feature["name"]
                    if gene and gene in name:
                        # Handle cases like "BRAF V600E" -> "V600E"
                        variant = name.replace(gene, "").strip()
                    elif name and not variant:
                        # Use name directly as variant (e.g., "V600E")
                        variant = name.strip()

            # Extract disease
            disease = hit.get("diseases", "") or ""

            # Extract drugs
            drugs = []
            drug_str = hit.get("drugs", "")
            if drug_str:
                # VICC concatenates drugs with spaces/commas
                drugs = [d.strip() for d in drug_str.replace(",", " ").split() if d.strip()]

            # Extract evidence info
            evidence_level = hit.get("evidence_label")  # A, B, C, D
            response_type = assoc.get("response_type")

            # Extract source from evidence
            source = "vicc"
            evidence_list = assoc.get("evidence", [])
            if evidence_list:
                for ev in evidence_list:
                    source_name = ev.get("evidenceType", {}).get("sourceName")
                    if source_name:
                        source = source_name.lower()
                        break

            # Extract publication
            publication_url = assoc.get("publication_url")

            # Extract oncogenic info
            oncogenic = assoc.get("oncogenic")

            # Extract description
            description = assoc.get("description", "")

            return VICCAssociation(
                description=description,
                gene=gene,
                variant=variant,
                disease=disease,
                drugs=drugs,
                evidence_level=evidence_level,
                response_type=response_type,
                source=source,
                publication_url=publication_url,
                oncogenic=oncogenic,
            )

        except Exception:
            return None

    def _is_compound_mutation_resistance(self, assoc: "VICCAssociation", variant: str | None) -> bool:
        """Check if resistance is due to a compound/secondary mutation, not the queried variant.

        Args:
            assoc: The VICC association to check
            variant: The variant being queried (e.g., "V560D")

        Returns:
            True if this is resistance from a secondary mutation, not the queried variant
        """
        if not variant or not assoc.is_resistance():
            return False

        desc_lower = assoc.description.lower() if assoc.description else ""

        # Check for indicators of secondary/compound mutations
        secondary_indicators = [
            "secondary mutation",
            "acquired mutation",
            "harboring " + variant.lower() + " and ",
            variant.lower() + " and " + assoc.gene.lower() if assoc.gene else "",
            "developed resistance",
            "resistance developed",
        ]

        for indicator in secondary_indicators:
            if indicator and indicator in desc_lower:
                return True

        return False

    async def _fetch_with_query(
        self,
        query: str,
        variant: str | None,
        tumor_type: str | None,
        max_results: int,
    ) -> list[VICCAssociation]:
        """Internal method to fetch associations with a specific query.

        Args:
            query: Lucene query string
            variant: Optional variant for filtering compound mutations
            tumor_type: Optional tumor type filter
            max_results: Maximum results

        Returns:
            List of VICCAssociation objects
        """
        client = self._get_client()
        url = f"{self.BASE_URL}/associations"
        params = {"q": query, "size": max_results}

        try:
            response = await client.get(url, params=params)
            response.raise_for_status()
            data = response.json()
        except httpx.HTTPError as e:
            raise VICCError(f"VICC API request failed: {e}")
        except Exception as e:
            raise VICCError(f"Failed to parse VICC response: {e}")

        associations = []
        hits = data.get("hits", {}).get("hits", [])

        for hit in hits:
            assoc = self._parse_association(hit)
            if assoc is None:
                continue

            # Filter out MolecularMatch data - it has unreliable response_type values
            # (returns AMP tier codes like "2C" instead of clinical response types)
            if assoc.source and assoc.source.lower() == "molecularmatch":
                continue

            # Filter out CIViC data from VICC - we get CIViC directly with evidence_direction
            # VICC doesn't pass through evidence_direction, so CIViC via VICC loses the
            # SUPPORTS/DOES_NOT_SUPPORT distinction (e.g., MK-2206 shows "Sensitivity"
            # when paper actually says "did not show sensitivity")
            if assoc.source and assoc.source.lower() == "civic":
                continue

            if tumor_type and not is_pan_cancer_term(assoc.disease) and not tumor_types_match(assoc.disease, tumor_type):
                continue

            if self._is_compound_mutation_resistance(assoc, variant):
                continue

            associations.append(assoc)

        return associations

    async def fetch_associations(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        max_results: int = 50,
    ) -> list[VICCAssociation]:
        """Fetch clinical interpretation associations from VICC MetaKB.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional variant notation (e.g., "V600E")
            tumor_type: Optional tumor type to filter results
            max_results: Maximum number of results to return

        Returns:
            List of VICCAssociation objects
        """
        # Build primary query with exact variant and tumor type
        query = self._build_query(gene, variant, tumor_type)
        associations = await self._fetch_with_query(query, variant, tumor_type, max_results)

        # For KIT variants, also try exon-level search if we didn't find much
        # This helps find evidence for deletions like W557_K558del which may be
        # catalogued as "KIT exon 11 deletion" in databases
        if gene.upper() == "KIT" and variant and len(associations) < 5:
            exon = self._get_kit_exon(variant)
            if exon:
                exon_query = self._build_exon_query(gene, exon, tumor_type)
                exon_associations = await self._fetch_with_query(
                    exon_query, variant, tumor_type, max_results
                )
                # Add unique exon-level associations
                seen_descriptions = {a.description for a in associations}
                for assoc in exon_associations:
                    if assoc.description not in seen_descriptions:
                        associations.append(assoc)
                        seen_descriptions.add(assoc.description)

        return associations

    async def fetch_sensitivity_associations(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        max_results: int = 25,
    ) -> list[VICCAssociation]:
        """Fetch only sensitivity/response associations.

        Args:
            gene: Gene symbol
            variant: Optional variant notation
            tumor_type: Optional tumor type filter
            max_results: Maximum results

        Returns:
            List of sensitivity associations
        """
        all_assocs = await self.fetch_associations(gene, variant, tumor_type, max_results * 2)
        return [a for a in all_assocs if a.is_sensitivity()][:max_results]

    async def fetch_resistance_associations(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        max_results: int = 25,
    ) -> list[VICCAssociation]:
        """Fetch only resistance associations.

        Args:
            gene: Gene symbol
            variant: Optional variant notation
            tumor_type: Optional tumor type filter
            max_results: Maximum results

        Returns:
            List of resistance associations
        """
        all_assocs = await self.fetch_associations(gene, variant, tumor_type, max_results * 2)
        return [a for a in all_assocs if a.is_resistance()][:max_results]

    async def fetch_vicc_evidence(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
        max_results: int = 50,
    ) -> list["VICCEvidence"]:
        """Fetch VICC associations and convert to evidence model.

        This method fetches associations and converts them to the
        VICCEvidence model for use in Insight.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional variant notation (e.g., "V600E")
            tumor_type: Optional tumor type to filter results
            max_results: Maximum number of results to return

        Returns:
            List of VICCEvidence objects
        """
        from oncomind.models.evidence.vicc import VICCEvidence

        associations = await self.fetch_associations(gene, variant, tumor_type, max_results)
        evidence_list = []

        for assoc in associations:
            # Determine match specificity
            locus_match = self._determine_locus_match(assoc.variant, gene, variant)
            # Build matched profile string
            matched_profile = f"{assoc.gene} {assoc.variant}" if assoc.gene and assoc.variant else assoc.gene

            # Determine tumor match using centralized function
            # Args: source_disease (from KB), queried_tumor (user's query)
            tumor_matches = tumor_types_match(assoc.disease, tumor_type) if tumor_type else False

            # Build EvidenceLevel objects for consistency with other models
            from oncomind.models.evidence.base import EvidenceLevel
            locus_variant_match = EvidenceLevel(
                level=locus_match,
                scope="specific" if locus_match == "variant" else "unspecified",
                origin="kb",
            )
            cancer_type_match = None
            if tumor_type:
                cancer_type_match = EvidenceLevel(
                    level="cancer_specific" if tumor_matches else (assoc.disease or "pan_cancer"),
                    scope="specific" if tumor_matches else "unspecified",
                    origin="kb",
                )

            evidence_list.append(VICCEvidence(
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
                matched_profile=matched_profile,
                locus_variant_match=locus_variant_match,
                cancer_type_match=cancer_type_match,
            ))

        return evidence_list

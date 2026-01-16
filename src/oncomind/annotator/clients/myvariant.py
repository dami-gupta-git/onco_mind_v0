"""Annotator MyVariant.info API client."""

from typing import Any

import httpx

from oncomind.annotator.models import MyVariantAnnotation
from oncomind.config.debug import get_logger
from oncomind.utils import to_hgvs_protein_three_letter

logger = get_logger(__name__)


class AnnotatorMyVariantClient:
    """MyVariant.info client for the annotator workflow."""

    BASE_URL = "https://myvariant.info/v1"
    DEFAULT_TIMEOUT = 30.0

    def __init__(self, timeout: float = DEFAULT_TIMEOUT) -> None:
        """Initialize the client.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout
        self._client: httpx.AsyncClient | None = None

    async def __aenter__(self) -> "AnnotatorMyVariantClient":
        """Async context manager entry."""
        self._client = httpx.AsyncClient(timeout=self.timeout)
        return self

    async def __aexit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Async context manager exit."""
        if self._client:
            await self._client.aclose()
            self._client = None

    def _get_client(self) -> httpx.AsyncClient:
        """Get or create the HTTP client."""
        if self._client is None:
            self._client = httpx.AsyncClient(timeout=self.timeout)
        return self._client

    async def _query(self, query: str, fields: list[str] | None = None) -> dict[str, Any]:
        """Execute a query against MyVariant API.

        Args:
            query: Query string (e.g., "BRAF:V600E" or "chr7:140453136")
            fields: Specific fields to retrieve

        Returns:
            API response as dictionary
        """
        client = self._get_client()
        params: dict[str, str] = {"q": query}

        if fields:
            params["fields"] = ",".join(fields)

        response = await client.get(f"{self.BASE_URL}/query", params=params)
        response.raise_for_status()
        return response.json()

    async def fetch_annotation(self, gene: str, variant: str) -> MyVariantAnnotation:
        """Fetch annotation and return only selected fields.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Variant notation (e.g., "p.V600E")

        Returns:
            MyVariantAnnotation with validated fields
        """
        annot = await self.build_annotation(gene, variant)

        # Extract selected fields from the first hit
        hits = annot.get("hits", [])
        if not hits:
            return MyVariantAnnotation()

        first_hit = hits[0]
        return MyVariantAnnotation.from_hit(first_hit)

    async def build_annotation(self, gene: str, variant: str) -> dict[str, Any]:
        """Fetch evidence for a variant from multiple sources.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Variant notation (e.g., "V600E")

        Returns:
            Aggregated evidence from all sources
        """
        # Request specific fields from ClinVar, COSMIC, and other annotation sources
        fields = [
            "clinvar",
            "cosmic",
            "dbsnp",
            "cadd",
            "entrezgene",  # NCBI Gene ID
            "cosmic.cosmic_id",  # COSMIC mutation ID
            "clinvar.variant_id",  # ClinVar variation ID
            "clinvar.rcv",  # ClinVar RCV records (contains clinical_significance and accession)
            "dbsnp.rsid",  # dbSNP rs number
            "hgvs",  # HGVS notations (genomic, protein, transcript)
            "snpeff",  # SnpEff effect prediction
            "dbnsfp.polyphen2.hdiv.pred",  # PolyPhen2 prediction
            "dbnsfp.polyphen2.hdiv.score",  # PolyPhen2 score
            "dbnsfp.sift.pred",  # SIFT prediction
            "dbnsfp.sift.score",  # SIFT score
            "dbnsfp.cadd.phred",  # CADD phred score
            "dbnsfp.alphamissense",  # AlphaMissense pathogenicity prediction
            "gnomad_exome.af.af",  # gnomAD exome allele frequency
            "vcf.alt",  # VCF alternative allele
            "vcf.ref",  # VCF reference allele
        ]

        # Try multiple query strategies to find the variant
        # Strategy 1: Gene with protein notation (e.g., "BRAF p.V600E")
        query = f"{gene} {variant}"
        result = await self._query(query, fields=fields)

        # Strategy 2: If no hits, try simple gene:variant (e.g., "BRAF:p.V600E")
        if result.get("total", 0) == 0:
            query = f"{gene}:{variant}"
            result = await self._query(query, fields=fields)

        # Strategy 3: If no hits, try three-letter amino acid codes (e.g., "EGFR p.Leu858Arg")
        if result.get("total", 0) == 0:
            three_letter = to_hgvs_protein_three_letter(variant)
            if three_letter:
                query = f"{gene} {three_letter}"
                result = await self._query(query, fields=fields)

        return result

    async def close(self) -> None:
        """Close the HTTP client."""
        if self._client:
            await self._client.aclose()
            self._client = None

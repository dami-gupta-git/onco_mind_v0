"""MyVariant.info API client for fetching variant evidence.

ARCHITECTURE:
    Gene + Variant → MyVariant.info API → Evidence (ClinVar/COSMIC/AlphaMissense/CADD/gnomAD)

Aggregates variant information from multiple databases for LLM assessment.
CIViC evidence is fetched separately via CIViCClient (GraphQL API).

Key Design:
- Async HTTP with connection pooling (httpx.AsyncClient)
- Retry with exponential backoff (tenacity)
- Structured parsing to typed Evidence models
- Context manager for session cleanup
"""
import re
from typing import Any

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from oncomind.models.myvariant import MyVariantHit, MyVariantResponse

from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.cosmic import COSMICEvidence
from oncomind.models.evidence.myvariant_evidence import MyVariantEvidence
from oncomind.config.debug import get_logger

logger = get_logger(__name__)


class MyVariantAPIError(Exception):
    """Exception raised for MyVariant API errors."""

    pass


# Amino acid single-letter to three-letter code mapping
AA_1_TO_3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter',  # Stop codon
}


def _convert_to_hgvs_p_three_letter(variant: str) -> str | None:
    """Convert single-letter amino acid variant to three-letter HGVS protein notation.

    Args:
        variant: Variant notation like "G12D", "V600E", "R248*"

    Returns:
        Three-letter HGVS notation like "p.Gly12Asp", or None if cannot convert
    """
    # Strip p. prefix if present
    v = variant.upper()
    if v.startswith("P."):
        v = v[2:]

    # Match pattern: RefAA + Position + AltAA (e.g., G12D, V600E, R248*)
    match = re.match(r'^([A-Z])(\d+)([A-Z*])$', v)
    if not match:
        return None

    ref_aa, pos, alt_aa = match.groups()
    ref_three = AA_1_TO_3.get(ref_aa)
    alt_three = AA_1_TO_3.get(alt_aa)

    if not ref_three or not alt_three:
        return None

    return f"p.{ref_three}{pos}{alt_three}"


class MyVariantClient:
    """Client for MyVariant.info API.

    MyVariant.info aggregates variant annotations from multiple sources
    including ClinVar, COSMIC, AlphaMissense, CADD, gnomAD and more.

    Note: CIViC evidence is fetched separately via CIViCClient (GraphQL API)
    to get the most up-to-date clinical interpretations.
    """

    BASE_URL = "https://myvariant.info/v1"
    DEFAULT_TIMEOUT = 30.0

    def __init__(
        self,
        timeout: float = DEFAULT_TIMEOUT,
        max_retries: int = 3,
    ) -> None:
        """Initialize the MyVariant client.

        Args:
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
        """
        self.timeout = timeout
        self.max_retries = max_retries
        self._client: httpx.AsyncClient | None = None

    async def __aenter__(self) -> "MyVariantClient":
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

    @retry(
        retry=retry_if_exception_type((httpx.HTTPError, httpx.TimeoutException)),
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=2, max=10),
    )
    async def _query(self, query: str, fields: list[str] | None = None) -> dict[str, Any]:
        """Execute a query against MyVariant API.

        Args:
            query: Query string (e.g., "BRAF:V600E" or "chr7:140453136")
            fields: Specific fields to retrieve

        Returns:
            API response as dictionary

        Raises:
            MyVariantAPIError: If the API request fails
        """
        client = self._get_client()
        params: dict[str, str] = {"q": query}

        if fields:
            params["fields"] = ",".join(fields)

        response = await client.get(f"{self.BASE_URL}/query", params=params)
        response.raise_for_status()
        data = response.json()

        if "error" in data:
            raise MyVariantAPIError(f"API error: {data['error']}")

        return data

    async def get_variant(self, variant_id: str) -> dict[str, Any]:
        """Get variant by ID.

        Args:
            variant_id: Variant identifier (HGVS, dbSNP, etc.)

        Returns:
            Variant data
        """
        client = self._get_client()
        response = await client.get(f"{self.BASE_URL}/variant/{variant_id}")
        response.raise_for_status()
        return response.json()

    def _parse_clinvar_evidence(
        self, clinvar_data: dict[str, Any] | list[Any]
    ) -> list[ClinVarEvidence]:
        """Parse ClinVar data into evidence objects.

        Args:
            clinvar_data: Raw ClinVar data from API

        Returns:
            List of ClinVar evidence objects
        """
        evidence_list: list[ClinVarEvidence] = []

        # Handle both single dict and list of dicts
        items = clinvar_data if isinstance(clinvar_data, list) else [clinvar_data]

        for item in items:
            if not isinstance(item, dict):
                continue

            variant_id = item.get("variant_id")

            # MyVariant ClinVar format: rcv array contains the clinical interpretations
            rcv_data = item.get("rcv", [])
            if rcv_data:
                rcv_list = rcv_data if isinstance(rcv_data, list) else [rcv_data]
                # Deduplicate by clinical significance to avoid showing 6x "Pathogenic"
                seen_significances: set[str] = set()
                for rcv in rcv_list:
                    if not isinstance(rcv, dict):
                        continue
                    clin_sig = rcv.get("clinical_significance")
                    sig_key = str(clin_sig).lower() if clin_sig else ""
                    if sig_key in seen_significances:
                        continue  # Skip duplicates
                    seen_significances.add(sig_key)

                    conditions = []
                    # RCV entries may have conditions
                    if "conditions" in rcv:
                        cond_data = rcv["conditions"]
                        if isinstance(cond_data, list):
                            for cond in cond_data:
                                if isinstance(cond, dict):
                                    conditions.append(cond.get("name", ""))
                                else:
                                    conditions.append(str(cond))
                        elif isinstance(cond_data, dict):
                            conditions.append(cond_data.get("name", ""))

                    evidence_list.append(
                        ClinVarEvidence(
                            clinical_significance=str(clin_sig) if clin_sig else None,
                            review_status=rcv.get("review_status"),
                            conditions=conditions,
                            last_evaluated=rcv.get("last_evaluated"),
                            variation_id=str(variant_id) if variant_id else None,
                        )
                    )
            else:
                # Fallback: old format with clinical_significance at top level
                clin_sig = item.get("clinical_significance")
                if isinstance(clin_sig, list):
                    clin_sig = ", ".join(str(s) for s in clin_sig)

                conditions = []
                if "conditions" in item:
                    cond_data = item["conditions"]
                    if isinstance(cond_data, list):
                        for cond in cond_data:
                            if isinstance(cond, dict):
                                conditions.append(cond.get("name", ""))
                            else:
                                conditions.append(str(cond))
                    elif isinstance(cond_data, dict):
                        conditions.append(cond_data.get("name", ""))

                evidence_list.append(
                    ClinVarEvidence(
                        clinical_significance=str(clin_sig) if clin_sig else None,
                        review_status=item.get("review_status"),
                        conditions=conditions,
                        last_evaluated=item.get("last_evaluated"),
                        variation_id=str(variant_id) if variant_id else None,
                    )
                )

        return evidence_list

    def _parse_cosmic_evidence(
        self, cosmic_data: dict[str, Any] | list[Any]
    ) -> list[COSMICEvidence]:
        """Parse COSMIC data into evidence objects.

        Args:
            cosmic_data: Raw COSMIC data from API

        Returns:
            List of COSMIC evidence objects
        """
        evidence_list: list[COSMICEvidence] = []

        # Handle both single dict and list of dicts
        items = cosmic_data if isinstance(cosmic_data, list) else [cosmic_data]

        for item in items:
            if not isinstance(item, dict):
                continue

            evidence_list.append(
                COSMICEvidence(
                    mutation_id=item.get("mutation_id"),
                    primary_site=item.get("primary_site"),
                    site_subtype=item.get("site_subtype"),
                    primary_histology=item.get("primary_histology"),
                    histology_subtype=item.get("histology_subtype"),
                    sample_count=item.get("sample_count"),
                    mutation_somatic_status=item.get("mutation_somatic_status"),
                )
            )

        return evidence_list

    def _extract_from_hit(
        self, hit: MyVariantHit, gene: str, variant: str
    ) -> MyVariantEvidence:
        """Extract Evidence fields from a parsed MyVariantHit using Pydantic models.

        This method uses Pydantic's automatic parsing instead of manual nested
        dictionary navigation, making it cleaner and more maintainable.

        Args:
            hit: Parsed MyVariant API hit
            gene: Gene symbol
            variant: Variant notation

        Returns:
            Evidence object with all extracted fields
        """
        # Extract database identifiers using Pydantic models
        cosmic_id = None
        if hit.cosmic:
            cosmic_data = hit.cosmic if isinstance(hit.cosmic, list) else [hit.cosmic]
            if cosmic_data and cosmic_data[0].cosmic_id:
                cosmic_id = cosmic_data[0].cosmic_id

        ncbi_gene_id = None
        if hit.entrezgene:
            ncbi_gene_id = str(hit.entrezgene)
        elif hit.dbsnp and hit.dbsnp.gene:
            # gene can be a single DbSNPGene or a list of them
            gene_data = hit.dbsnp.gene
            if isinstance(gene_data, list):
                # Take the first gene if it's a list
                if gene_data and gene_data[0].geneid:
                    ncbi_gene_id = str(gene_data[0].geneid)
            elif gene_data.geneid:
                ncbi_gene_id = str(gene_data.geneid)

        dbsnp_id = None
        if hit.dbsnp and hit.dbsnp.rsid:
            rsid = hit.dbsnp.rsid
            dbsnp_id = f"rs{rsid}" if not rsid.startswith("rs") else rsid

        # Extract ClinVar data
        clinvar_id = None
        clinvar_clinical_significance = None
        clinvar_accession = None
        if hit.clinvar:
            clinvar_list = hit.clinvar if isinstance(hit.clinvar, list) else [hit.clinvar]
            if clinvar_list:
                first_clinvar = clinvar_list[0]
                if first_clinvar.variant_id:
                    clinvar_id = str(first_clinvar.variant_id)
                # Extract from rcv (can be single object or list)
                if first_clinvar.rcv:
                    rcv_list = first_clinvar.rcv if isinstance(first_clinvar.rcv, list) else [first_clinvar.rcv]
                    if rcv_list:
                        first_rcv = rcv_list[0]
                        if first_rcv.clinical_significance:
                            clinvar_clinical_significance = first_rcv.clinical_significance
                        if first_rcv.accession:
                            clinvar_accession = first_rcv.accession

        # Extract HGVS notations
        hgvs_genomic = None
        hgvs_protein = None
        hgvs_transcript = None

        # Use variant id as genomic HGVS if it looks like HGVS
        if hit.id and (hit.id.startswith("chr") or hit.id.startswith("NC_")):
            hgvs_genomic = hit.id

        if hit.hgvs:
            hgvs_list = [hit.hgvs] if isinstance(hit.hgvs, str) else hit.hgvs
            for hgvs in hgvs_list:
                if hgvs.startswith("chr") or hgvs.startswith("NC_"):
                    hgvs_genomic = hgvs
                elif ":p." in hgvs and not hgvs_protein:
                    hgvs_protein = hgvs
                elif ":c." in hgvs and not hgvs_transcript:
                    hgvs_transcript = hgvs

        # Extract functional annotations using Pydantic models
        snpeff_effect = None
        transcript_id = None
        transcript_consequence = None
        if hit.snpeff and hit.snpeff.ann:
            ann = hit.snpeff.ann
            ann_data = ann if isinstance(ann, list) else [ann]
            if ann_data:
                first_ann = ann_data[0]
                snpeff_effect = first_ann.effect
                transcript_id = first_ann.feature_id
                transcript_consequence = first_ann.effect

        polyphen2_prediction = None
        polyphen2_score = None
        if hit.dbnsfp and hit.dbnsfp.polyphen2 and hit.dbnsfp.polyphen2.hdiv:
            hdiv = hit.dbnsfp.polyphen2.hdiv
            # Handle both string and list[str] from API
            if hdiv.pred is not None:
                if isinstance(hdiv.pred, list):
                    polyphen2_prediction = hdiv.pred[0] if hdiv.pred else None
                else:
                    polyphen2_prediction = hdiv.pred
            # Extract score
            if hdiv.score is not None:
                try:
                    if isinstance(hdiv.score, list):
                        polyphen2_score = float(hdiv.score[0]) if hdiv.score else None
                    else:
                        polyphen2_score = float(hdiv.score)
                except (ValueError, TypeError):
                    pass

        # Extract SIFT prediction
        sift_prediction = None
        sift_score = None
        if hit.dbnsfp and hit.dbnsfp.sift:
            sift = hit.dbnsfp.sift
            # Handle both string and list[str] from API
            if sift.pred is not None:
                if isinstance(sift.pred, list):
                    sift_prediction = sift.pred[0] if sift.pred else None
                else:
                    sift_prediction = sift.pred
            # Extract score (lower = more deleterious)
            if sift.score is not None:
                try:
                    if isinstance(sift.score, list):
                        sift_score = float(sift.score[0]) if sift.score else None
                    else:
                        sift_score = float(sift.score)
                except (ValueError, TypeError):
                    pass

        cadd_score = None
        # Try dbnsfp first, then top-level cadd
        if hit.dbnsfp and hit.dbnsfp.cadd and hit.dbnsfp.cadd.phred:
            try:
                cadd_score = float(hit.dbnsfp.cadd.phred)
            except (ValueError, TypeError):
                pass
        if cadd_score is None and hit.cadd and hit.cadd.phred:
            try:
                cadd_score = float(hit.cadd.phred)
            except (ValueError, TypeError):
                pass

        gnomad_exome_af = None
        if hit.gnomad_exome and hit.gnomad_exome.af and hit.gnomad_exome.af.af:
            try:
                gnomad_exome_af = float(hit.gnomad_exome.af.af)
            except (ValueError, TypeError):
                pass

        # Extract AlphaMissense prediction
        alphamissense_score = None
        alphamissense_prediction = None
        if hit.dbnsfp and hit.dbnsfp.alphamissense:
            am = hit.dbnsfp.alphamissense
            # Handle score (can be float or list[float])
            if am.score is not None:
                try:
                    if isinstance(am.score, list):
                        alphamissense_score = float(am.score[0]) if am.score else None
                    else:
                        alphamissense_score = float(am.score)
                except (ValueError, TypeError):
                    pass
            # Handle prediction (can be str or list[str])
            if am.pred is not None:
                if isinstance(am.pred, list):
                    alphamissense_prediction = am.pred[0] if am.pred else None
                else:
                    alphamissense_prediction = am.pred

        # Parse evidence using existing parsers (ClinVar, COSMIC)
        # Note: CIViC evidence is fetched separately via CIViCClient (GraphQL API)
        clinvar_evidence = []
        if hit.clinvar:
            # Convert back to dict for existing parser
            clinvar_data = hit.clinvar
            if isinstance(clinvar_data, list):
                clinvar_evidence = self._parse_clinvar_evidence([c.model_dump() for c in clinvar_data])
            else:
                clinvar_evidence = self._parse_clinvar_evidence(clinvar_data.model_dump())

        cosmic_evidence = []
        if hit.cosmic:
            # Convert back to dict for existing parser
            cosmic_data = hit.cosmic
            if isinstance(cosmic_data, list):
                cosmic_evidence = self._parse_cosmic_evidence([c.model_dump() for c in cosmic_data])
            else:
                cosmic_evidence = self._parse_cosmic_evidence(cosmic_data.model_dump())

        return MyVariantEvidence(
            variant_id=hit.id,
            gene=gene,
            variant=variant,
            cosmic_id=cosmic_id,
            ncbi_gene_id=ncbi_gene_id,
            dbsnp_id=dbsnp_id,
            clinvar_id=clinvar_id,
            clinvar_clinical_significance=clinvar_clinical_significance,
            clinvar_accession=clinvar_accession,
            hgvs_genomic=hgvs_genomic,
            hgvs_transcript=hgvs_transcript,
            hgvs_protein=hgvs_protein,
            snpeff_effect=snpeff_effect,
            polyphen2_prediction=polyphen2_prediction,
            polyphen2_score=polyphen2_score,
            sift_prediction=sift_prediction,
            sift_score=sift_score,
            cadd_score=cadd_score,
            gnomad_exome_af=gnomad_exome_af,
            alphamissense_score=alphamissense_score,
            alphamissense_prediction=alphamissense_prediction,
            transcript_id=transcript_id,
            transcript_consequence=transcript_consequence,
            civic=[],  # CIViC evidence is fetched separately via CIViCClient
            clinvar=clinvar_evidence,
            cosmic=cosmic_evidence,
            raw_data=hit.model_dump(by_alias=True),
        )

    async def _fetch_clinvar_fallback(self, gene: str, variant: str) -> dict[str, Any] | None:
        """
        Fallback to fetch ClinVar data directly from NCBI E-utilities when MyVariant doesn't have it.

        Uses NCBI's E-utilities API to search ClinVar for the variant.

        Args:
            gene: Gene symbol (e.g., "EGFR")
            variant: Variant notation (e.g., "L858R")

        Returns:
            ClinVar data dict with variant_id, clinical_significance, and accession if found
        """
        import re

        try:
            client = self._get_client()
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

            # Build search terms - try multiple formats for better matching
            # Convert X/Ter notation: K3326X -> Lys3326Ter, V600E -> Val600Glu
            aa_map = {
                'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
                'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
                'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
                'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
                'X': 'Ter', '*': 'Ter',
            }

            search_terms = []

            # Parse variant like V600E, K3326X, etc.
            match = re.match(r'^([A-Z])(\d+)([A-Z*])$', variant.upper())
            if match:
                ref_aa, pos, alt_aa = match.groups()
                ref_long = aa_map.get(ref_aa, ref_aa)
                alt_long = aa_map.get(alt_aa, alt_aa)
                # Try three-letter format FIRST - more reliable for ClinVar
                search_terms.append(f"{ref_long}{pos}{alt_long}")
                # Then try position-only for more lenient matching
                search_terms.append(f"{ref_long}{pos}")

            # Finally try original notation
            search_terms.append(variant)

            id_list = []
            for search_variant in search_terms:
                search_term = f"{gene}[gene] AND {search_variant}"
                search_params = {
                    "db": "clinvar",
                    "term": search_term,
                    "retmode": "json",
                    "retmax": 5  # Get multiple results to find best match
                }

                search_response = await client.get(search_url, params=search_params)
                if search_response.status_code == 200:
                    search_data = search_response.json()
                    id_list = search_data.get("esearchresult", {}).get("idlist", [])
                    if id_list:
                        break  # Found results, stop searching

            if not id_list:
                return None

            # Fetch summaries for all candidates to find best match
            summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            summary_params = {
                "db": "clinvar",
                "id": ",".join(id_list),  # Fetch all at once
                "retmode": "json"
            }

            summary_response = await client.get(summary_url, params=summary_params)
            if summary_response.status_code != 200:
                return None

            summary_data = summary_response.json()

            # Find the best matching result by checking the title contains our variant
            # Build expected patterns to match in title
            expected_patterns = []
            if match:
                # For K3326X, look for "Lys3326Ter" or "p.Lys3326Ter" in title
                expected_patterns.append(f"{ref_long}{pos}{alt_long}".lower())
                expected_patterns.append(f"p.{ref_long}{pos}{alt_long}".lower())
            expected_patterns.append(variant.lower())

            best_result = None
            best_variant_id = None

            for vid in id_list:
                result = summary_data.get("result", {}).get(vid, {})
                title = result.get("title", "").lower()

                # Check if any expected pattern is in the title
                is_match = any(pat in title for pat in expected_patterns)
                if is_match:
                    best_result = result
                    best_variant_id = vid
                    break

            # Fall back to first result if no exact match found
            if not best_result:
                best_variant_id = id_list[0]
                best_result = summary_data.get("result", {}).get(best_variant_id, {})

            # Extract relevant fields - NCBI ClinVar API has multiple classification fields
            # Try clinical_significance first (legacy), then germline_classification, then somatic_classification
            clinical_significance = best_result.get("clinical_significance", {}).get("description")
            if not clinical_significance:
                clinical_significance = best_result.get("germline_classification", {}).get("description")
            if not clinical_significance:
                clinical_significance = best_result.get("somatic_classification", {}).get("description")

            # Also get review status
            review_status = None
            if best_result.get("germline_classification", {}).get("review_status"):
                review_status = best_result["germline_classification"]["review_status"]
            elif best_result.get("somatic_classification", {}).get("review_status"):
                review_status = best_result["somatic_classification"]["review_status"]

            accession = best_result.get("accession")

            if clinical_significance or accession:
                return {
                    "variant_id": best_variant_id,
                    "clinical_significance": clinical_significance,
                    "accession": accession,
                    "review_status": review_status,
                }

            return None

        except Exception:
            return None

    async def fetch_evidence(self, gene: str, variant: str) -> MyVariantEvidence:
        """Fetch evidence for a variant from multiple sources.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Variant notation (e.g., "V600E")

        Returns:
            Aggregated evidence from all sources

        Raises:
            MyVariantAPIError: If the API request fails
        """
        # Request specific fields from ClinVar, COSMIC, and other annotation sources
        # Note: CIViC evidence is fetched separately via CIViCClient (GraphQL API)
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

        try:
            # Try multiple query strategies to find the variant
            # Strategy 1: Gene with protein notation (e.g., "BRAF p.V600E")
            # This works best with MyVariant API
            protein_notation = f"p.{variant}" if not variant.startswith("p.") else variant
            query = f"{gene} {protein_notation}"
            result = await self._query(query, fields=fields)

            # Strategy 2: If no hits, try simple gene:variant (e.g., "BRAF:V600E")
            if result.get("total", 0) == 0:
                query = f"{gene}:{variant}"
                result = await self._query(query, fields=fields)

            # Strategy 3: If still no hits, try searching by gene name and variant without prefix
            if result.get("total", 0) == 0:
                query = f"{gene} {variant}"
                result = await self._query(query, fields=fields)

            # Strategy 4: Use dbnsfp/snpeff field query with three-letter amino acid codes
            # This is more reliable for hotspot variants like KRAS G12D
            # Query format: dbnsfp.genename:KRAS AND snpeff.ann.hgvs_p:p.Gly12Asp
            if result.get("total", 0) == 0:
                hgvs_p_three = _convert_to_hgvs_p_three_letter(variant)
                if hgvs_p_three:
                    query = f"dbnsfp.genename:{gene} AND snpeff.ann.hgvs_p:{hgvs_p_three}"
                    result = await self._query(query, fields=fields)

            # Strategy 5: Use VEP to get genomic coordinates, then re-query MyVariant
            # VEP converts protein notation to HGVS genomic, which MyVariant indexes better
            vep_annotation = None
            if result.get("total", 0) == 0:
                try:
                    from oncomind.api.vep import VEPClient
                    vep_client = VEPClient()
                    vep_annotation = await vep_client.annotate_variant(gene, variant)

                    if vep_annotation and vep_annotation.myvariant_query:
                        # Re-query MyVariant with genomic notation
                        result = await self._query(vep_annotation.myvariant_query, fields=fields)
                except Exception:
                    # VEP failed - continue with existing result
                    pass

            # Parse response using Pydantic
            parsed_response = MyVariantResponse(**result)

            if not parsed_response.hits:
                # No data found in MyVariant - try ClinVar fallback
                # Note: CIViC evidence is fetched separately via CIViCClient
                clinvar_fallback = await self._fetch_clinvar_fallback(gene, variant)

                # Extract ClinVar data from fallback
                clinvar_id = None
                clinvar_significance = None
                clinvar_accession = None
                if clinvar_fallback:
                    clinvar_id = clinvar_fallback.get("variant_id")
                    clinvar_significance = clinvar_fallback.get("clinical_significance")
                    clinvar_accession = clinvar_fallback.get("accession")

                # Use VEP predictions if available (from Strategy 5)
                vep_hgvs_genomic = None
                vep_hgvs_transcript = None
                vep_polyphen = None
                vep_cadd = None
                vep_alphamissense_score = None
                vep_alphamissense_pred = None
                if vep_annotation:
                    vep_hgvs_genomic = vep_annotation.hgvs_genomic
                    vep_hgvs_transcript = vep_annotation.hgvs_transcript
                    vep_polyphen = vep_annotation.polyphen_prediction
                    vep_cadd = vep_annotation.cadd_phred
                    vep_alphamissense_score = vep_annotation.alphamissense_score
                    vep_alphamissense_pred = vep_annotation.alphamissense_prediction

                # Return evidence with fallback data (ClinVar and VEP predictions)
                # CIViC evidence should be fetched separately via CIViCClient
                return MyVariantEvidence(
                    variant_id=f"{gene}:{variant}",
                    gene=gene,
                    variant=variant,
                    cosmic_id=None,
                    ncbi_gene_id=None,
                    dbsnp_id=None,
                    clinvar_id=clinvar_id,
                    clinvar_clinical_significance=clinvar_significance,
                    clinvar_accession=clinvar_accession,
                    hgvs_genomic=vep_hgvs_genomic,
                    hgvs_protein=None,
                    hgvs_transcript=vep_hgvs_transcript,
                    snpeff_effect=None,
                    polyphen2_prediction=vep_polyphen,
                    cadd_score=vep_cadd,
                    gnomad_exome_af=None,
                    alphamissense_score=vep_alphamissense_score,
                    alphamissense_prediction=vep_alphamissense_pred,
                    transcript_id=None,
                    transcript_consequence=None,
                    civic=[],  # CIViC evidence is fetched separately via CIViCClient
                    clinvar=[],
                    cosmic=[],
                    raw_data=result,
                )

            # Use the first hit (most relevant) and extract using Pydantic models
            first_hit = parsed_response.hits[0]
            evidence = self._extract_from_hit(first_hit, gene, variant)

            # If MyVariant returned no ClinVar data, try ClinVar fallback
            # Note: CIViC evidence is fetched separately via CIViCClient
            if not evidence.clinvar_id and not evidence.clinvar_clinical_significance:
                clinvar_fallback = await self._fetch_clinvar_fallback(gene, variant)
                if clinvar_fallback:
                    # Update evidence with ClinVar fallback data
                    evidence.clinvar_id = clinvar_fallback.get("variant_id")
                    evidence.clinvar_clinical_significance = clinvar_fallback.get("clinical_significance")
                    evidence.clinvar_accession = clinvar_fallback.get("accession")

            return evidence

        except MyVariantAPIError:
            raise
        except Exception as e:
            logger.error(f"Failed to parse evidence for {gene} {variant}: {str(e)}")
            raise MyVariantAPIError(f"Failed to parse evidence: {str(e)}")

    async def close(self) -> None:
        """Close the HTTP client."""
        if self._client:
            await self._client.aclose()
            self._client = None

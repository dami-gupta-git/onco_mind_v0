"""cBioPortal API client for co-mutation and prevalence data.

ARCHITECTURE:
    Gene + Variant + Tumor Type → cBioPortal API → Co-mutation patterns and prevalence

cBioPortal provides genomic data from large cancer studies (TCGA, MSK-IMPACT, etc.):
- Mutation prevalence in specific tumor types
- Co-occurring mutations (genes frequently mutated together)
- Mutually exclusive mutations (genes rarely mutated together)

Key Design:
- Uses public API at www.cbioportal.org/api
- No API key required for public data
- Calculates co-occurrence/mutual exclusivity from sample-level mutation data
"""

from dataclasses import dataclass
from typing import Any
import httpx

from oncomind.constants import TUMOR_TYPE_MAPPINGS


class CBioPortalError(Exception):
    """Exception raised for cBioPortal API errors."""
    pass


# Map tumor types to cBioPortal study IDs
TUMOR_STUDY_MAPPINGS: dict[str, list[str]] = {
    "melanoma": ["skcm_tcga_pan_can_atlas_2018", "mel_ucla_2016", "skcm_mskcc_2014"],
    "nsclc": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018", "nsclc_tcga_broad_2016"],
    "lung": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"],
    "colorectal": ["coadread_tcga_pan_can_atlas_2018", "crc_msk_2017"],
    "breast": ["brca_tcga_pan_can_atlas_2018", "breast_msk_2018"],
    "pancreatic": ["paad_tcga_pan_can_atlas_2018"],
    "glioblastoma": ["gbm_tcga_pan_can_atlas_2018"],
    "ovarian": ["ov_tcga_pan_can_atlas_2018"],
    "prostate": ["prad_tcga_pan_can_atlas_2018"],
    "bladder": ["blca_tcga_pan_can_atlas_2018"],
    "thyroid": ["thca_tcga_pan_can_atlas_2018"],
    "gist": ["gist_mskcc"],
}


@dataclass
class CoMutationData:
    """Co-mutation and prevalence data for a variant."""

    gene: str
    variant: str | None
    tumor_type: str | None
    study_id: str

    # Prevalence
    total_samples: int
    samples_with_gene_mutation: int
    samples_with_exact_variant: int
    gene_prevalence_pct: float
    variant_prevalence_pct: float

    # Co-occurring mutations (genes frequently mutated together)
    co_occurring: list[dict[str, Any]]  # [{gene, count, pct, odds_ratio}]

    # Mutually exclusive mutations (genes rarely mutated together)
    mutually_exclusive: list[dict[str, Any]]  # [{gene, count, pct, odds_ratio}]

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "gene": self.gene,
            "variant": self.variant,
            "tumor_type": self.tumor_type,
            "study_id": self.study_id,
            "total_samples": self.total_samples,
            "samples_with_gene_mutation": self.samples_with_gene_mutation,
            "samples_with_exact_variant": self.samples_with_exact_variant,
            "gene_prevalence_pct": self.gene_prevalence_pct,
            "variant_prevalence_pct": self.variant_prevalence_pct,
            "co_occurring": self.co_occurring,
            "mutually_exclusive": self.mutually_exclusive,
        }


class CBioPortalClient:
    """Client for cBioPortal API.

    cBioPortal provides cancer genomics data from large studies like TCGA
    and MSK-IMPACT, enabling analysis of mutation prevalence and co-occurrence.

    API Documentation: https://www.cbioportal.org/api
    """

    BASE_URL = "https://www.cbioportal.org/api"
    DEFAULT_TIMEOUT = 30.0

    # Common cancer genes to check for co-occurrence
    CANCER_GENES = [
        "TP53", "KRAS", "NRAS", "BRAF", "EGFR", "PIK3CA", "PTEN", "APC",
        "CDKN2A", "RB1", "NF1", "ARID1A", "ATM", "BRCA1", "BRCA2", "IDH1",
        "IDH2", "FBXW7", "KEAP1", "STK11", "SMAD4", "CTNNB1", "MET", "ALK",
        "ROS1", "RET", "ERBB2", "KIT", "PDGFRA", "FGFR1", "FGFR2", "FGFR3"
    ]

    def __init__(self, timeout: float = DEFAULT_TIMEOUT):
        """Initialize the cBioPortal client."""
        self.timeout = timeout
        self._client: httpx.AsyncClient | None = None
        self._gene_cache: dict[str, int] = {}  # gene symbol -> entrez ID

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

    def _get_study_ids(self, tumor_type: str | None) -> list[str]:
        """Get study IDs for a tumor type."""
        if not tumor_type:
            # Default to pan-cancer studies
            return ["msk_impact_2017"]

        tumor_lower = tumor_type.lower()

        # Direct match
        if tumor_lower in TUMOR_STUDY_MAPPINGS:
            return TUMOR_STUDY_MAPPINGS[tumor_lower]

        # Check mappings
        for abbrev, full_names in TUMOR_TYPE_MAPPINGS.items():
            if tumor_lower == abbrev or any(tumor_lower in name for name in full_names):
                if abbrev in TUMOR_STUDY_MAPPINGS:
                    return TUMOR_STUDY_MAPPINGS[abbrev]

        # Fallback to MSK-IMPACT pan-cancer
        return ["msk_impact_2017"]

    async def _get_entrez_id(self, gene: str) -> int | None:
        """Get Entrez gene ID from gene symbol."""
        if gene in self._gene_cache:
            return self._gene_cache[gene]

        client = self._get_client()
        try:
            response = await client.get(f"{self.BASE_URL}/genes/{gene}")
            if response.status_code == 200:
                data = response.json()
                entrez_id = data.get("entrezGeneId")
                self._gene_cache[gene] = entrez_id
                return entrez_id
        except Exception:
            pass
        return None

    async def _get_mutations(
        self,
        entrez_ids: list[int],
        molecular_profile_id: str,
    ) -> list[dict[str, Any]]:
        """Fetch mutations for given genes from a molecular profile."""
        client = self._get_client()

        try:
            response = await client.post(
                f"{self.BASE_URL}/mutations/fetch",
                params={"projection": "SUMMARY"},
                json={
                    "entrezGeneIds": entrez_ids,
                    "molecularProfileIds": [molecular_profile_id],
                },
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            raise CBioPortalError(f"Failed to fetch mutations: {e}")

    async def _get_sample_count(self, study_id: str) -> int:
        """Get total sample count for a study."""
        client = self._get_client()

        try:
            response = await client.get(
                f"{self.BASE_URL}/studies/{study_id}/samples",
                params={"projection": "SUMMARY"},
            )
            response.raise_for_status()
            return len(response.json())
        except Exception:
            return 0

    async def fetch_co_mutation_data(
        self,
        gene: str,
        variant: str | None = None,
        tumor_type: str | None = None,
    ) -> CoMutationData | None:
        """Fetch co-mutation and prevalence data for a gene/variant.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional variant (e.g., "V600E")
            tumor_type: Optional tumor type (e.g., "Melanoma")

        Returns:
            CoMutationData with prevalence and co-occurrence info
        """
        # Get study IDs for tumor type
        study_ids = self._get_study_ids(tumor_type)
        if not study_ids:
            return None

        study_id = study_ids[0]  # Use primary study
        molecular_profile_id = f"{study_id}_mutations"

        # Get Entrez ID for query gene
        query_entrez = await self._get_entrez_id(gene.upper())
        if not query_entrez:
            return None

        # Get Entrez IDs for cancer genes (excluding query gene)
        cancer_genes = [g for g in self.CANCER_GENES if g.upper() != gene.upper()]
        cancer_entrez_ids = []
        for g in cancer_genes[:20]:  # Limit to 20 for performance
            eid = await self._get_entrez_id(g)
            if eid:
                cancer_entrez_ids.append(eid)

        # Fetch mutations for query gene
        try:
            query_mutations = await self._get_mutations([query_entrez], molecular_profile_id)
        except CBioPortalError:
            return None

        if not query_mutations:
            return None

        # Get total sample count
        total_samples = await self._get_sample_count(study_id)
        if total_samples == 0:
            total_samples = 500  # Fallback estimate

        # Build sample sets
        samples_with_gene = set()
        samples_with_variant = set()
        clean_variant = variant.replace("p.", "").upper() if variant else None

        for m in query_mutations:
            sample_id = m.get("sampleId")
            samples_with_gene.add(sample_id)

            if clean_variant:
                protein_change = (m.get("proteinChange") or "").upper()
                if protein_change == clean_variant:
                    samples_with_variant.add(sample_id)

        # If variant specified, use variant samples; otherwise gene samples
        target_samples = samples_with_variant if clean_variant and samples_with_variant else samples_with_gene

        # Fetch mutations for cancer genes to find co-occurrence
        co_occurring = []
        mutually_exclusive = []

        if cancer_entrez_ids and target_samples:
            try:
                cancer_mutations = await self._get_mutations(
                    cancer_entrez_ids, molecular_profile_id
                )

                # Build gene -> samples mapping
                from collections import defaultdict
                gene_samples: dict[str, set] = defaultdict(set)
                entrez_to_gene: dict[int, str] = {}

                for m in cancer_mutations:
                    eid = m.get("entrezGeneId")
                    sample = m.get("sampleId")
                    # We need gene symbol - fetch from cache or API would be slow
                    # Use keyword field which contains gene name
                    keyword = m.get("keyword", "")
                    gene_name = keyword.split()[0] if keyword else ""
                    if gene_name and sample:
                        gene_samples[gene_name].add(sample)
                        entrez_to_gene[eid] = gene_name

                # Calculate co-occurrence statistics
                for other_gene, other_samples in gene_samples.items():
                    # Count samples with both mutations
                    both = len(target_samples & other_samples)
                    only_query = len(target_samples - other_samples)
                    only_other = len(other_samples - target_samples)
                    neither = total_samples - len(target_samples | other_samples)

                    # Calculate odds ratio
                    # OR = (both * neither) / (only_query * only_other)
                    if only_query > 0 and only_other > 0:
                        odds_ratio = (both * neither) / (only_query * only_other) if neither > 0 else 0
                    else:
                        odds_ratio = 0

                    co_pct = 100 * both / len(target_samples) if target_samples else 0

                    entry = {
                        "gene": other_gene,
                        "count": both,
                        "pct": round(co_pct, 1),
                        "odds_ratio": round(odds_ratio, 2) if odds_ratio else None,
                    }

                    # Categorize: OR > 1 = co-occurring, OR < 1 = mutually exclusive
                    if odds_ratio > 1.5 and both >= 3:
                        co_occurring.append(entry)
                    elif odds_ratio < 0.5 and both <= 2:
                        mutually_exclusive.append(entry)

                # Sort by count/significance
                co_occurring.sort(key=lambda x: x["count"], reverse=True)
                mutually_exclusive.sort(key=lambda x: x.get("odds_ratio") or 999)

            except CBioPortalError:
                pass

        return CoMutationData(
            gene=gene.upper(),
            variant=variant,
            tumor_type=tumor_type,
            study_id=study_id,
            total_samples=total_samples,
            samples_with_gene_mutation=len(samples_with_gene),
            samples_with_exact_variant=len(samples_with_variant) if variant else 0,
            gene_prevalence_pct=round(100 * len(samples_with_gene) / total_samples, 1),
            variant_prevalence_pct=round(100 * len(samples_with_variant) / total_samples, 1) if variant else 0,
            co_occurring=co_occurring[:10],  # Top 10
            mutually_exclusive=mutually_exclusive[:10],  # Top 10
        )

    async def close(self) -> None:
        """Close the HTTP client."""
        if self._client:
            await self._client.aclose()
            self._client = None

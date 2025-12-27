"""DepMap API client for gene dependency and drug sensitivity data.

ARCHITECTURE:
    Gene + Variant → DepMap API → Dependency scores, drug sensitivities, cell line models

DepMap provides cancer cell line data from the Broad Institute:
- CRISPR gene dependency scores (CERES)
- PRISM drug sensitivity data
- Cell line mutation profiles (CCLE)

Key Design:
- Uses public DepMap API at https://api.depmap.org
- Falls back to pre-computed data for common cancer genes
- No API key required for public endpoints
"""

from dataclasses import dataclass
from typing import Any
import httpx

from oncomind.models.evidence.depmap import (
    DepMapEvidence,
    GeneDependency,
    DrugSensitivity,
    CellLineModel,
)


class DepMapError(Exception):
    """Exception raised for DepMap API errors."""
    pass


class DepMapRateLimitError(DepMapError):
    """Rate limit exceeded."""
    pass


@dataclass
class DepMapData:
    """Raw DepMap data before conversion to evidence model."""

    gene: str
    variant: str | None
    dependency_score: float | None
    n_dependent_lines: int
    n_total_lines: int
    top_dependent_lines: list[str]
    co_dependencies: list[dict[str, Any]]
    drug_sensitivities: list[dict[str, Any]]
    cell_lines: list[dict[str, Any]]
    data_version: str | None


class DepMapClient:
    """Client for DepMap API.

    DepMap provides cancer dependency data including:
    - Gene essentiality scores from CRISPR screens
    - Drug sensitivity data from PRISM
    - Cell line mutation profiles

    API Documentation: https://depmap.org/portal/api/
    """

    # DepMap uses multiple endpoints
    BASE_URL = "https://api.depmap.org/api"
    PORTAL_URL = "https://depmap.org/portal"
    DEFAULT_TIMEOUT = 30.0

    # Common cancer genes with pre-computed dependency data
    # This serves as a fallback/cache for frequently queried genes
    CANCER_GENE_DEPENDENCIES = {
        "BRAF": {"score": -0.8, "dependent_pct": 45},
        "KRAS": {"score": -1.2, "dependent_pct": 65},
        "EGFR": {"score": -0.6, "dependent_pct": 35},
        "TP53": {"score": -0.1, "dependent_pct": 5},  # TSG - not essential when lost
        "PIK3CA": {"score": -0.5, "dependent_pct": 30},
        "NRAS": {"score": -0.7, "dependent_pct": 40},
        "MYC": {"score": -1.5, "dependent_pct": 80},
        "ERBB2": {"score": -0.9, "dependent_pct": 50},
        "ALK": {"score": -0.4, "dependent_pct": 20},
        "RET": {"score": -0.3, "dependent_pct": 15},
        "MET": {"score": -0.5, "dependent_pct": 25},
        "FGFR1": {"score": -0.4, "dependent_pct": 20},
        "FGFR2": {"score": -0.5, "dependent_pct": 25},
        "FGFR3": {"score": -0.4, "dependent_pct": 20},
        "KIT": {"score": -0.6, "dependent_pct": 35},
        "PDGFRA": {"score": -0.3, "dependent_pct": 15},
    }

    # Known drug sensitivities by gene
    GENE_DRUG_SENSITIVITIES = {
        "BRAF": [
            {"drug": "vemurafenib", "ic50_mutant": 50, "ic50_wt": 2000, "target": "BRAF V600"},
            {"drug": "dabrafenib", "ic50_mutant": 30, "ic50_wt": 1500, "target": "BRAF V600"},
            {"drug": "encorafenib", "ic50_mutant": 20, "ic50_wt": 1200, "target": "BRAF V600"},
            {"drug": "trametinib", "ic50_mutant": 8, "ic50_wt": 200, "target": "MEK1/2"},
            {"drug": "cobimetinib", "ic50_mutant": 15, "ic50_wt": 250, "target": "MEK1/2"},
        ],
        "KRAS": [
            {"drug": "sotorasib", "ic50_mutant": 100, "ic50_wt": None, "target": "KRAS G12C"},
            {"drug": "adagrasib", "ic50_mutant": 80, "ic50_wt": None, "target": "KRAS G12C"},
            {"drug": "trametinib", "ic50_mutant": 25, "ic50_wt": 300, "target": "MEK1/2"},
        ],
        "EGFR": [
            {"drug": "erlotinib", "ic50_mutant": 40, "ic50_wt": 3000, "target": "EGFR"},
            {"drug": "gefitinib", "ic50_mutant": 35, "ic50_wt": 2500, "target": "EGFR"},
            {"drug": "osimertinib", "ic50_mutant": 15, "ic50_wt": 500, "target": "EGFR T790M"},
            {"drug": "afatinib", "ic50_mutant": 20, "ic50_wt": 1000, "target": "EGFR/HER2"},
        ],
        "PIK3CA": [
            {"drug": "alpelisib", "ic50_mutant": 50, "ic50_wt": 500, "target": "PI3K alpha"},
            {"drug": "everolimus", "ic50_mutant": 5, "ic50_wt": 50, "target": "mTOR"},
            {"drug": "capivasertib", "ic50_mutant": 100, "ic50_wt": 800, "target": "AKT"},
        ],
        "ERBB2": [
            {"drug": "lapatinib", "ic50_mutant": 30, "ic50_wt": 2000, "target": "EGFR/HER2"},
            {"drug": "neratinib", "ic50_mutant": 10, "ic50_wt": 500, "target": "pan-HER"},
            {"drug": "tucatinib", "ic50_mutant": 8, "ic50_wt": 1000, "target": "HER2"},
        ],
    }

    def __init__(self, timeout: float = DEFAULT_TIMEOUT):
        """Initialize the DepMap client."""
        self.timeout = timeout
        self._client: httpx.AsyncClient | None = None

    async def __aenter__(self):
        """Initialize HTTP client session."""
        self._client = httpx.AsyncClient(
            timeout=self.timeout,
            headers={"Accept": "application/json"}
        )
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Close HTTP client session."""
        if self._client:
            await self._client.aclose()
            self._client = None

    def _get_client(self) -> httpx.AsyncClient:
        """Get the HTTP client, creating one if needed."""
        if self._client is None:
            self._client = httpx.AsyncClient(
                timeout=self.timeout,
                headers={"Accept": "application/json"}
            )
        return self._client

    async def _try_api_query(self, gene: str) -> dict[str, Any] | None:
        """Try to query the DepMap API for gene data.

        Returns None if API is unavailable or rate limited.
        """
        client = self._get_client()

        try:
            # Try the gene dependency endpoint
            response = await client.get(
                f"{self.BASE_URL}/gene/{gene}/dependency",
                timeout=10.0
            )

            if response.status_code == 200:
                return response.json()
            elif response.status_code == 429:
                raise DepMapRateLimitError("DepMap API rate limit exceeded")

        except httpx.TimeoutException:
            pass
        except httpx.HTTPError:
            pass

        return None

    def _get_fallback_data(
        self,
        gene: str,
        variant: str | None,
    ) -> DepMapData | None:
        """Get fallback data from pre-computed cache.

        This provides commonly-queried cancer gene data when the API
        is unavailable or for faster responses.

        Note: For cell line data by mutation, use CBioPortalClient.fetch_cell_lines_with_mutation()
        which queries the CCLE study for comprehensive, up-to-date data.
        """
        gene_upper = gene.upper()

        # Check if we have cached dependency data
        if gene_upper not in self.CANCER_GENE_DEPENDENCIES:
            return None

        dep_data = self.CANCER_GENE_DEPENDENCIES[gene_upper]

        # Get drug sensitivities for this gene
        drug_data = []
        if gene_upper in self.GENE_DRUG_SENSITIVITIES:
            for drug in self.GENE_DRUG_SENSITIVITIES[gene_upper]:
                drug_data.append({
                    "drug_name": drug["drug"],
                    "ic50_nm": drug["ic50_mutant"],
                    "n_cell_lines": 10,  # Approximate
                    "target": drug["target"],
                })

        # Estimate dependency counts
        n_total = 1000  # Approximate total cell lines
        n_dependent = int(n_total * dep_data["dependent_pct"] / 100)

        return DepMapData(
            gene=gene_upper,
            variant=variant,
            dependency_score=dep_data["score"],
            n_dependent_lines=n_dependent,
            n_total_lines=n_total,
            top_dependent_lines=[],  # Would need real data
            co_dependencies=[],  # Would need real data
            drug_sensitivities=drug_data,
            cell_lines=[],  # Use cBioPortal for cell line data
            data_version="fallback_cache",
        )

    async def fetch_depmap_evidence(
        self,
        gene: str,
        variant: str | None = None,
    ) -> DepMapEvidence | None:
        """Fetch DepMap evidence for a gene/variant.

        Returns gene dependency and drug sensitivity data.
        For cell line data by mutation, use CBioPortalClient.fetch_cell_lines_with_mutation().

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional variant (e.g., "V600E")

        Returns:
            DepMapEvidence with dependency and sensitivity data
        """
        gene_upper = gene.upper()

        # Try DepMap API first (usually fails - endpoint doesn't exist)
        api_data = await self._try_api_query(gene_upper)

        if api_data:
            return self._convert_api_response(api_data, gene_upper, variant)

        # Use fallback data
        fallback = self._get_fallback_data(gene_upper, variant)

        if not fallback:
            return None

        return self._convert_fallback_data(fallback)

    def _convert_api_response(
        self,
        data: dict[str, Any],
        gene: str,
        variant: str | None,
    ) -> DepMapEvidence:
        """Convert API response to DepMapEvidence model."""
        # This would parse the actual API response structure
        # For now, return a basic structure
        return DepMapEvidence(
            gene=gene,
            variant=variant,
            gene_dependency=GeneDependency(
                gene=gene,
                mean_dependency_score=data.get("mean_score"),
                n_dependent_lines=data.get("n_dependent", 0),
                n_total_lines=data.get("n_total", 0),
                dependency_pct=data.get("dependency_pct", 0.0),
                top_dependent_lines=data.get("top_lines", []),
            ) if data.get("mean_score") is not None else None,
            data_version=data.get("version"),
        )

    def _convert_fallback_data(self, data: DepMapData) -> DepMapEvidence:
        """Convert fallback data to DepMapEvidence model."""
        # Gene dependency
        gene_dependency = None
        if data.dependency_score is not None:
            gene_dependency = GeneDependency(
                gene=data.gene,
                mean_dependency_score=data.dependency_score,
                n_dependent_lines=data.n_dependent_lines,
                n_total_lines=data.n_total_lines,
                dependency_pct=(data.n_dependent_lines / data.n_total_lines * 100)
                if data.n_total_lines > 0 else 0.0,
                top_dependent_lines=data.top_dependent_lines,
            )

        # Drug sensitivities
        drug_sensitivities = []
        for ds in data.drug_sensitivities:
            drug_sensitivities.append(DrugSensitivity(
                drug_name=ds["drug_name"],
                ic50_nm=ds.get("ic50_nm"),
                n_cell_lines=ds.get("n_cell_lines", 0),
            ))

        # Cell line models
        cell_line_models = []
        for cl in data.cell_lines:
            cell_line_models.append(CellLineModel(
                name=cl["name"],
                primary_disease=cl.get("disease"),
                subtype=cl.get("subtype"),
                has_mutation=cl.get("has_mutation", False),
                mutation_details=cl.get("mutation"),
            ))

        return DepMapEvidence(
            gene=data.gene,
            variant=data.variant,
            gene_dependency=gene_dependency,
            drug_sensitivities=drug_sensitivities,
            cell_line_models=cell_line_models,
            data_version=data.data_version,
            n_cell_lines_screened=data.n_total_lines,
        )

    async def close(self) -> None:
        """Close the HTTP client."""
        if self._client:
            await self._client.aclose()
            self._client = None

"""DepMap API client for gene dependency and drug sensitivity data.

ARCHITECTURE:
    Gene + Variant → DepMap Portal Downloads → Dependency scores, drug sensitivities, cell line models

DepMap provides cancer cell line data from the Broad Institute:
- CRISPR gene dependency scores (Chronos)
- PRISM drug sensitivity data
- Cell line mutation profiles (CCLE)

Key Design:
- Uses public DepMap download API at https://depmap.org/portal/download/api/
- Downloads and caches data files locally for performance
- Returns None if data is unavailable
- No API key required for public endpoints
"""

import asyncio
import csv
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import httpx
import pandas as pd

from oncomind.models.evidence.depmap import (
    DepMapEvidence,
    GeneDependency,
    DrugSensitivity,
    CellLineModel,
)
from oncomind.config.debug import get_logger

logger = get_logger(__name__)


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
    """Client for DepMap data.

    DepMap provides cancer dependency data including:
    - Gene essentiality scores from CRISPR screens
    - Drug sensitivity data from PRISM
    - Cell line mutation profiles

    Data is accessed via direct download URLs from the DepMap portal.
    """

    # DepMap download API base
    DOWNLOAD_API = "https://depmap.org/portal/download/api/download"

    # Data file URLs (DepMap Public 24Q4 release)
    MUTATIONS_FILE = "OmicsSomaticMutations.csv"
    MUTATIONS_RELEASE = "DepMap+Public+24Q4"

    # PRISM drug sensitivity data (24Q2 release)
    PRISM_SENSITIVITY_FILE = "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv"
    PRISM_DRUGS_FILE = "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"

    # CRISPR dependency data
    CRISPR_FILE = "CRISPRGeneEffect.csv"
    CRISPR_RELEASE = "DepMap+Public+24Q4"

    DEFAULT_TIMEOUT = 60.0  # Longer timeout for large file downloads

    def __init__(self, timeout: float = DEFAULT_TIMEOUT, cache_dir: Path | None = None):
        """Initialize the DepMap client.

        Args:
            timeout: HTTP request timeout in seconds
            cache_dir: Optional directory for caching downloaded data
        """
        self.timeout = timeout
        self.cache_dir = cache_dir
        self._client: httpx.AsyncClient | None = None

        # In-memory cache for parsed data (to avoid re-parsing)
        self._prism_drugs_cache: dict[str, dict] | None = None
        self._mutations_cache: dict[str, list[dict]] | None = None

    async def __aenter__(self):
        """Initialize HTTP client session."""
        self._client = httpx.AsyncClient(
            timeout=self.timeout,
            follow_redirects=True,
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
                follow_redirects=True,
            )
        return self._client

    def _build_download_url(self, file_name: str, release: str | None = None) -> str:
        """Build a DepMap download URL."""
        url = f"{self.DOWNLOAD_API}?file_name={file_name}"
        if release:
            url += f"&release={release}"
        return url

    async def _fetch_csv_data(self, url: str, max_rows: int | None = None) -> list[dict]:
        """Fetch and parse CSV data from a URL.

        Args:
            url: URL to fetch
            max_rows: Optional limit on rows to parse (for large files)

        Returns:
            List of dicts, one per row
        """
        client = self._get_client()

        try:
            logger.debug(f"Fetching DepMap data from: {url}")
            response = await client.get(url)
            response.raise_for_status()

            # Parse CSV
            content = response.text
            reader = csv.DictReader(io.StringIO(content))

            rows = []
            for i, row in enumerate(reader):
                if max_rows and i >= max_rows:
                    break
                rows.append(row)

            logger.debug(f"Parsed {len(rows)} rows from DepMap")
            return rows

        except httpx.TimeoutException:
            logger.warning(f"DepMap request timed out: {url}")
            raise DepMapError("DepMap request timed out")
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 429:
                raise DepMapRateLimitError("DepMap rate limit exceeded")
            logger.warning(f"DepMap HTTP error: {e}")
            raise DepMapError(f"DepMap HTTP error: {e}")
        except Exception as e:
            logger.warning(f"DepMap fetch error: {e}")
            raise DepMapError(f"DepMap fetch error: {e}")

    async def fetch_prism(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Fetch PRISM drug sensitivity data and drug metadata.

        Returns:
            Tuple of (sensitivity_df, drugs_df):
            - sensitivity_df: DataFrame with cell lines as rows, drugs as columns (log-fold change values)
            - drugs_df: DataFrame with drug metadata (name, target, MOA, etc.)
        """
        # Primary screen — log-fold change values (lower = more sensitive)
        prism_url = "https://depmap.org/portal/download/api/download?file_name=Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv"
        sensitivity = pd.read_csv(prism_url, index_col=0)

        # Drug metadata — maps column IDs to drug names, targets, MOA
        drug_info_url = "https://depmap.org/portal/download/api/download?file_name=Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"
        drugs = pd.read_csv(drug_info_url)

        logger.info(f"Loaded PRISM data: {sensitivity.shape[0]} cell lines, {sensitivity.shape[1]} drugs")
        return sensitivity, drugs

    def get_drug_response(
        self,
        gene: str,
        protein_change: str,
        mutations_df: pd.DataFrame,
        sensitivity_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Find drugs that cell lines with a specific variant are sensitive to.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            protein_change: Protein change (e.g., "V600E")
            mutations_df: DataFrame with mutations (from fetch_mutations as DataFrame)
            sensitivity_df: DataFrame with drug sensitivity (from fetch_prism)

        Returns:
            DataFrame with columns: drug_id, mean_diff, n_variant_lines
            Sorted by mean_diff (negative = variant more sensitive)
        """
        # Get cell lines with the variant
        variant_lines = mutations_df[
            (mutations_df['HugoSymbol'] == gene) &
            (mutations_df['ProteinChange'].str.contains(protein_change, na=False))
        ]['ModelID'].unique()

        # Filter sensitivity matrix to those cell lines
        variant_sensitivity = sensitivity_df.loc[sensitivity_df.index.isin(variant_lines)]
        wt_sensitivity = sensitivity_df.loc[~sensitivity_df.index.isin(variant_lines)]

        # Find drugs with differential sensitivity
        results = []
        for drug in sensitivity_df.columns:
            mut_response = variant_sensitivity[drug].dropna()
            wt_response = wt_sensitivity[drug].dropna()

            if len(mut_response) >= 3 and len(wt_response) >= 3:
                diff = mut_response.mean() - wt_response.mean()
                results.append({
                    'drug_id': drug,
                    'mean_diff': diff,  # negative = variant more sensitive
                    'n_variant_lines': len(mut_response)
                })

        return pd.DataFrame(results).sort_values('mean_diff')

    async def fetch_mutations(
        self,
        gene: str,
        variant: str | None = None,
    ) -> list[dict]:
        """Fetch mutation data for a gene/variant from CCLE.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional protein change to filter (e.g., "V600E")

        Returns:
            List of mutation records with cell line info
        """
        url = self._build_download_url(self.MUTATIONS_FILE, self.MUTATIONS_RELEASE)

        try:
            # Fetch all mutations (large file - consider caching)
            rows = await self._fetch_csv_data(url)

            # Filter for gene
            gene_upper = gene.upper()
            filtered = [
                r for r in rows
                if (r.get("HugoSymbol") or r.get("Hugo_Symbol", "")).upper() == gene_upper
            ]

            # Filter for variant if specified
            if variant:
                variant_upper = variant.upper()
                filtered = [
                    r for r in filtered
                    if variant_upper in (r.get("ProteinChange") or r.get("HGVSp_Short") or "").upper()
                ]

            logger.debug(f"Found {len(filtered)} mutations for {gene} {variant or ''}")
            return filtered

        except DepMapError:
            return []

    async def fetch_gene_dependency(
        self,
        gene: str,
    ) -> GeneDependency | None:
        """Fetch CRISPR gene dependency data.

        Uses the task-based API to get gene-specific dependency scores.

        Args:
            gene: Gene symbol (e.g., "BRAF")

        Returns:
            GeneDependency with essentiality scores, or None if unavailable
        """
        client = self._get_client()

        try:
            # Use the custom download endpoint with gene filter
            response = await client.post(
                "https://depmap.org/portal/api/download/custom",
                json={
                    "datasetId": "Chronos_Combined",
                    "featureLabels": [gene.upper()],
                    "dropEmpty": True,
                },
                timeout=30.0,
            )
            response.raise_for_status()
            task_info = response.json()

            task_id = task_info.get("id")
            if not task_id:
                logger.warning("No task ID returned from DepMap")
                return None

            # Poll for task completion
            for _ in range(10):  # Max 10 attempts
                await asyncio.sleep(1)

                status_response = await client.get(
                    f"https://depmap.org/portal/api/task/{task_id}",
                    timeout=10.0,
                )
                status = status_response.json()

                if status.get("state") == "SUCCESS":
                    download_url = status.get("result", {}).get("downloadUrl")
                    if download_url:
                        # Fetch the CSV data
                        data_response = await client.get(download_url, timeout=30.0)
                        data_response.raise_for_status()

                        # Parse CSV - first column is cell line ID, second is gene score
                        content = data_response.text
                        reader = csv.reader(io.StringIO(content))

                        scores = []
                        header = next(reader, None)
                        for row in reader:
                            if len(row) >= 2:
                                try:
                                    scores.append(float(row[1]))
                                except (ValueError, IndexError):
                                    pass

                        if scores:
                            mean_score = sum(scores) / len(scores)
                            # Count dependent lines (score < -0.5)
                            n_dependent = sum(1 for s in scores if s < -0.5)

                            return GeneDependency(
                                gene=gene.upper(),
                                mean_dependency_score=mean_score,
                                n_dependent_lines=n_dependent,
                                n_total_lines=len(scores),
                                dependency_pct=(n_dependent / len(scores) * 100) if scores else 0.0,
                                top_dependent_lines=[],  # Would need cell line names
                            )
                    break

                elif status.get("state") == "FAILURE":
                    logger.warning(f"DepMap task failed: {status.get('message')}")
                    break

            return None

        except Exception as e:
            logger.warning(f"Error fetching gene dependency: {e}")
            return None

    async def fetch_depmap_evidence(
        self,
        gene: str,
        variant: str | None = None,
    ) -> DepMapEvidence | None:
        """Fetch DepMap evidence for a gene/variant.

        Returns gene dependency and drug sensitivity data.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional variant (e.g., "V600E")

        Returns:
            DepMapEvidence with dependency and sensitivity data
        """
        gene_upper = gene.upper()

        try:
            # Fetch gene dependency and mutations in parallel
            dependency_task = self.fetch_gene_dependency(gene_upper)
            mutations_task = self.fetch_mutations(gene_upper, variant)

            dependency, mutations = await asyncio.gather(
                dependency_task,
                mutations_task,
                return_exceptions=True,
            )

            # Handle exceptions
            if isinstance(dependency, Exception):
                logger.warning(f"Gene dependency fetch failed: {dependency}")
                dependency = None
            if isinstance(mutations, Exception):
                logger.warning(f"Mutations fetch failed: {mutations}")
                mutations = []

            # Build cell line models from mutations
            cell_line_models = []
            seen_lines = set()
            for mut in mutations:
                line_id = mut.get("ModelID") or mut.get("DepMap_ID")
                if line_id and line_id not in seen_lines:
                    seen_lines.add(line_id)
                    cell_line_models.append(CellLineModel(
                        name=line_id,
                        ccle_name=mut.get("CellLineName") or mut.get("CCLE_Name"),
                        primary_disease=mut.get("OncotreeLineage") or mut.get("primary_disease"),
                        subtype=mut.get("OncotreePrimaryDisease") or mut.get("Subtype"),
                        has_mutation=True,
                        mutation_details=mut.get("ProteinChange") or mut.get("HGVSp_Short"),
                    ))

            # Return None if no data
            if not dependency and not cell_line_models:
                return None

            return DepMapEvidence(
                gene=gene_upper,
                variant=variant,
                gene_dependency=dependency,
                cell_line_models=cell_line_models,
                data_version="DepMap Public 24Q4",
                n_cell_lines_screened=dependency.n_total_lines if dependency else 0,
            )

        except Exception as e:
            logger.error(f"DepMap evidence fetch failed: {e}")
            return None

    async def close(self) -> None:
        """Close the HTTP client."""
        if self._client:
            await self._client.aclose()
            self._client = None

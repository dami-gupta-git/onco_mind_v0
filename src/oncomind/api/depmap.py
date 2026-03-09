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
from oncomind.config.constants import (
    DEPMAP_DOWNLOAD_API,
    DEPMAP_CUSTOM_DOWNLOAD_API,
    DEPMAP_TASK_STATUS_API,
    DEPMAP_PRISM_SENSITIVITY_URL,
    DEPMAP_PRISM_DRUG_INFO_URL,
    DEPMAP_MUTATIONS_FILE,
    DEPMAP_CRISPR_FILE,
    DEPMAP_RELEASE,
    DEPMAP_DATA_VERSION,
    DEPMAP_DEFAULT_TIMEOUT,
    DEPMAP_TASK_POLL_TIMEOUT,
    DEPMAP_DOWNLOAD_TIMEOUT,
    DEPMAP_MAX_TASK_POLL_ATTEMPTS,
    DEPMAP_TASK_POLL_INTERVAL,
    DEPMAP_SENSITIVITY_THRESHOLD,
    DEPMAP_MIN_CELL_LINES,
    DEPMAP_TOP_DRUGS_LIMIT,
    DEPMAP_SENSITIVE_LINES_DISPLAY,
    DEPMAP_DEPENDENCY_THRESHOLD,
)

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

    def __init__(
        self, timeout: float = DEPMAP_DEFAULT_TIMEOUT, cache_dir: Path | None = None
    ):
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
        url = f"{DEPMAP_DOWNLOAD_API}?file_name={file_name}"
        if release:
            url += f"&release={release}"
        return url

    async def _fetch_csv_data(
        self, url: str, max_rows: int | None = None
    ) -> list[dict]:
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
        sensitivity = pd.read_csv(DEPMAP_PRISM_SENSITIVITY_URL, index_col=0)

        # Drug metadata — maps column IDs to drug names, targets, MOA
        drugs = pd.read_csv(DEPMAP_PRISM_DRUG_INFO_URL)

        logger.info(
            f"Loaded PRISM data: {sensitivity.shape[0]} cell lines, {sensitivity.shape[1]} drugs"
        )
        return sensitivity, drugs

    def fetch_drug_sensitivities(
        self,
        variant_cell_lines: list[str],
        sensitivity_threshold: float = DEPMAP_SENSITIVITY_THRESHOLD,
        min_cell_lines: int = DEPMAP_MIN_CELL_LINES,
        top_n: int = DEPMAP_TOP_DRUGS_LIMIT,
    ) -> list[DrugSensitivity]:
        """Find drugs that cell lines with a specific mutation are sensitive to.

        Uses locally cached PRISM data files.

        Args:
            variant_cell_lines: List of DepMap IDs for cell lines with the variant
            sensitivity_threshold: Log2FC threshold for sensitivity (default -1.7)
            min_cell_lines: Minimum cell lines required for a drug to be included
            top_n: Maximum number of drugs to return

        Returns:
            List of DrugSensitivity objects, sorted by mean_log2fc (most sensitive first)
        """
        if not variant_cell_lines:
            return []

        data_dir = Path(__file__).parent.parent.parent.parent / "data" / "depmap"
        sensitivity_file = (
            data_dir / "primary-screen-replicate-collapsed-logfold-change.csv"
        )
        treatment_file = (
            data_dir / "primary-screen-replicate-collapsed-treatment-info.csv"
        )

        if not sensitivity_file.exists() or not treatment_file.exists():
            logger.warning("PRISM data files not found")
            return []

        try:
            # Load treatment info to get drug names
            treatment_df = pd.read_csv(treatment_file)
            # Create mapping from column_name to drug name
            # The column names in sensitivity file are like "BRD-A00077618-236-07-6::2.5::HTS"
            drug_name_map = {}
            for _, row in treatment_df.iterrows():
                col_name = row.get("column_name")
                drug_name = row.get("name")
                if col_name and drug_name:
                    drug_name_map[col_name] = drug_name

            # Load sensitivity data (cell lines as rows, drugs as columns)
            sensitivity_df = pd.read_csv(sensitivity_file, index_col=0)

            # Get cell lines that are in our variant set AND in the sensitivity data
            variant_lines_set = set(variant_cell_lines)
            available_variant_lines = [
                line for line in sensitivity_df.index if line in variant_lines_set
            ]

            if len(available_variant_lines) < min_cell_lines:
                logger.debug(
                    f"Only {len(available_variant_lines)} variant lines in PRISM data"
                )
                return []

            # Calculate mean response for variant cell lines for each drug
            variant_df = sensitivity_df.loc[available_variant_lines]
            sensitivities = []

            for drug_col in variant_df.columns:
                responses = variant_df[drug_col].dropna()
                if len(responses) < min_cell_lines:
                    continue

                mean_log2fc = responses.mean()

                # Only include if mean response is below threshold (sensitive)
                if mean_log2fc <= sensitivity_threshold:
                    # Get drug name from mapping, or extract from column name
                    drug_name = drug_name_map.get(drug_col)
                    if not drug_name:
                        # Try to match by broad_id prefix
                        broad_id = (
                            drug_col.split("::")[0] if "::" in drug_col else drug_col
                        )
                        matches = treatment_df[treatment_df["broad_id"] == broad_id]
                        if not matches.empty:
                            drug_name = matches.iloc[0].get("name")

                    if not drug_name:
                        drug_name = drug_col.split("::")[0]  # Use broad_id as fallback

                    # Get list of sensitive cell lines
                    sensitive_lines = [
                        line
                        for line in available_variant_lines
                        if pd.notna(variant_df.loc[line, drug_col])
                        and variant_df.loc[line, drug_col] <= sensitivity_threshold
                    ]

                    sensitivities.append(
                        DrugSensitivity(
                            drug_name=drug_name,
                            mean_log2fc=round(mean_log2fc, 3),
                            n_cell_lines=len(responses),
                            sensitive_lines=sensitive_lines[
                                :DEPMAP_SENSITIVE_LINES_DISPLAY
                            ],
                        )
                    )

            # Sort by mean_log2fc (most negative = most sensitive)
            sensitivities.sort(
                key=lambda x: x.mean_log2fc if x.mean_log2fc is not None else 0
            )

            logger.debug(f"Found {len(sensitivities)} sensitive drugs for variant")
            return sensitivities[:top_n]

        except Exception as e:
            logger.warning(f"Error calculating drug sensitivities: {e}")
            return []

    def fetch_mutations(
        self,
        gene: str,
        variant: str | None = None,
    ) -> list[dict]:
        """Fetch mutation data for a gene/variant from local CCLE file.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Optional protein change to filter (e.g., "V600E")

        Returns:
            List of mutation records with cell line info
        """
        # Path to local data files
        data_dir = Path(__file__).parent.parent.parent.parent / "data" / "depmap"
        mutations_file = data_dir / DEPMAP_MUTATIONS_FILE
        sample_info_file = data_dir / "sample_info.csv"

        if not mutations_file.exists():
            logger.warning(f"Mutations file not found: {mutations_file}")
            return []

        try:
            gene_upper = gene.upper()

            # Read mutations file, filtering by gene during read for efficiency
            # Use chunked reading for the large file
            chunks = pd.read_csv(
                mutations_file,
                usecols=[
                    "ModelID",
                    "HugoSymbol",
                    "ProteinChange",
                    "VariantType",
                    "DNAChange",
                ],
                chunksize=100000,
            )

            # Filter for gene in each chunk
            gene_mutations = []
            for chunk in chunks:
                chunk_filtered = chunk[chunk["HugoSymbol"].str.upper() == gene_upper]
                if not chunk_filtered.empty:
                    gene_mutations.append(chunk_filtered)

            if not gene_mutations:
                logger.debug(f"No mutations found for {gene}")
                return []

            mutations_df = pd.concat(gene_mutations, ignore_index=True)

            # Filter for variant if specified
            # Use precise matching to avoid substring false positives:
            # - S1982R should NOT match S1982RfsTer22 (frameshift)
            # - V600E should match p.V600E but not V600Efs
            if variant:
                import re

                variant_upper = variant.upper()
                # Escape regex special chars, then require word boundary or end after variant
                # This ensures S1982R matches "p.S1982R" but not "p.S1982RfsTer22"
                pattern = re.escape(variant_upper) + r"(?:[^A-Z]|$)"
                mutations_df = mutations_df[
                    mutations_df["ProteinChange"]
                    .fillna("")
                    .str.upper()
                    .str.contains(pattern, regex=True)
                ]

            if mutations_df.empty:
                logger.debug(f"No mutations found for {gene} {variant or ''}")
                return []

            # Join with sample_info for cell line metadata
            if sample_info_file.exists():
                sample_info = pd.read_csv(
                    sample_info_file,
                    usecols=[
                        "DepMap_ID",
                        "cell_line_name",
                        "CCLE_Name",
                        "primary_disease",
                        "Subtype",
                    ],
                )
                mutations_df = mutations_df.merge(
                    sample_info,
                    left_on="ModelID",
                    right_on="DepMap_ID",
                    how="left",
                )

            # Convert to list of dicts
            result = mutations_df.to_dict("records")
            logger.debug(f"Found {len(result)} mutations for {gene} {variant or ''}")
            return result

        except Exception as e:
            logger.warning(f"Error reading mutations file: {e}")
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
                DEPMAP_CUSTOM_DOWNLOAD_API,
                json={
                    "datasetId": "Chronos_Combined",
                    "featureLabels": [gene.upper()],
                    "dropEmpty": True,
                },
                timeout=DEPMAP_DOWNLOAD_TIMEOUT,
            )
            response.raise_for_status()
            task_info = response.json()

            task_id = task_info.get("id")
            if not task_id:
                logger.warning("No task ID returned from DepMap")
                return None

            # Poll for task completion
            for _ in range(DEPMAP_MAX_TASK_POLL_ATTEMPTS):
                await asyncio.sleep(DEPMAP_TASK_POLL_INTERVAL)

                status_response = await client.get(
                    f"{DEPMAP_TASK_STATUS_API}/{task_id}",
                    timeout=DEPMAP_TASK_POLL_TIMEOUT,
                )
                status = status_response.json()

                if status.get("state") == "SUCCESS":
                    download_url = status.get("result", {}).get("downloadUrl")
                    if download_url:
                        # Fetch the CSV data
                        data_response = await client.get(
                            download_url, timeout=DEPMAP_DOWNLOAD_TIMEOUT
                        )
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
                            # Count dependent lines (score below threshold)
                            n_dependent = sum(
                                1 for s in scores if s < DEPMAP_DEPENDENCY_THRESHOLD
                            )

                            return GeneDependency(
                                gene=gene.upper(),
                                mean_dependency_score=mean_score,
                                n_dependent_lines=n_dependent,
                                n_total_lines=len(scores),
                                dependency_pct=(
                                    (n_dependent / len(scores) * 100) if scores else 0.0
                                ),
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
            # Fetch mutations from local file (synchronous)
            mutations = self.fetch_mutations(gene_upper, variant)

            # Fetch gene dependency asynchronously
            dependency = await self.fetch_gene_dependency(gene_upper)
            if isinstance(dependency, Exception):
                logger.warning(f"Gene dependency fetch failed: {dependency}")
                dependency = None

            # Build cell line models from mutations
            cell_line_models = []
            variant_cell_line_ids = []
            seen_lines = set()

            # Helper to convert pandas NaN to None
            def clean_val(val) -> str | None:
                if val is None or (isinstance(val, float) and pd.isna(val)):
                    return None
                return str(val) if val else None

            for mut in mutations:
                line_id = mut.get("ModelID") or mut.get("DepMap_ID")
                if line_id and line_id not in seen_lines:
                    seen_lines.add(line_id)
                    variant_cell_line_ids.append(line_id)

                    cell_line_name = clean_val(mut.get("cell_line_name"))
                    cell_line_models.append(
                        CellLineModel(
                            name=cell_line_name or line_id,
                            depmap_id=line_id,
                            ccle_name=clean_val(mut.get("CCLE_Name")),
                            primary_disease=clean_val(mut.get("primary_disease")),
                            subtype=clean_val(mut.get("Subtype")),
                            has_mutation=True,
                            mutation_details=clean_val(mut.get("ProteinChange")),
                        )
                    )

            # Fetch drug sensitivities for variant cell lines
            drug_sensitivities = self.fetch_drug_sensitivities(variant_cell_line_ids)

            # Return None if no data
            if not dependency and not cell_line_models:
                return None

            return DepMapEvidence(
                gene=gene_upper,
                variant=variant,
                gene_dependency=dependency,
                drug_sensitivities=drug_sensitivities,
                cell_line_models=cell_line_models,
                data_version=DEPMAP_DATA_VERSION,
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

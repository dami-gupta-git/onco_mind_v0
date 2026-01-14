"""OncoKB API client for cancer gene list.

OncoKB is a precision oncology knowledge base that contains information about
the effects and treatment implications of specific cancer gene alterations.

This module fetches the curated cancer gene list from OncoKB's public API.
The list is used to determine if a variant is in a known cancer gene (Tier III-B).

Source: https://www.oncokb.org/cancerGenes
API: https://www.oncokb.org/api/v1/utils/cancerGeneList
"""

import logging
import aiohttp
from functools import lru_cache
from typing import Set

logger = logging.getLogger(__name__)

ONCOKB_CANCER_GENE_LIST_URL = "https://www.oncokb.org/api/v1/utils/cancerGeneList"

# Cache for the cancer gene set (module-level to persist across calls)
_cancer_gene_cache: Set[str] | None = None


async def fetch_cancer_gene_list() -> Set[str]:
    """Fetch the list of known cancer genes from OncoKB API.

    Returns:
        Set of gene symbols (Hugo symbols) that are known cancer genes.
        Returns empty set if API call fails.
    """
    global _cancer_gene_cache

    # Return cached result if available
    if _cancer_gene_cache is not None:
        return _cancer_gene_cache

    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(
                ONCOKB_CANCER_GENE_LIST_URL,
                timeout=aiohttp.ClientTimeout(total=10)
            ) as response:
                if response.status != 200:
                    logger.warning(f"OncoKB API returned status {response.status}")
                    return set()

                data = await response.json()

                # Extract hugoSymbol from each gene entry
                genes = set()
                for gene_entry in data:
                    if isinstance(gene_entry, dict) and 'hugoSymbol' in gene_entry:
                        genes.add(gene_entry['hugoSymbol'].upper())

                logger.info(f"Fetched {len(genes)} cancer genes from OncoKB")
                _cancer_gene_cache = genes
                return genes

    except aiohttp.ClientError as e:
        logger.warning(f"Failed to fetch OncoKB cancer gene list: {e}")
        return set()
    except Exception as e:
        logger.error(f"Unexpected error fetching OncoKB cancer gene list: {e}")
        return set()


def is_known_cancer_gene_sync(gene: str) -> bool:
    """Synchronous check if a gene is in the OncoKB cancer gene list.

    Uses cached data if available. For async contexts, use fetch_cancer_gene_list().

    Args:
        gene: Gene symbol (e.g., 'EGFR', 'BRAF')

    Returns:
        True if gene is in the OncoKB cancer gene list, False otherwise.
        Returns False if cache is not populated (need to call fetch_cancer_gene_list first).
    """
    if _cancer_gene_cache is None:
        # Cache not populated - return False (conservative)
        # In practice, the engine should populate the cache during initialization
        return False

    return gene.upper() in _cancer_gene_cache


def get_cached_cancer_genes() -> Set[str]:
    """Get the cached cancer gene set.

    Returns:
        Set of gene symbols, or empty set if cache not populated.
    """
    return _cancer_gene_cache or set()



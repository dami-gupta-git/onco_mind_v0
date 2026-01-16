"""Annotator for orchestrating variant annotation.

This module provides the Annotator class that coordinates
annotation from various sources.
"""

import time
from pathlib import Path
from typing import Any

from oncomind.annotator.config import AnnotatorConfig
from oncomind.annotator.result import AnnotationResult
from oncomind.annotator.vcf_record import VCFRecord
from oncomind.config.debug import get_logger

from oncomind.annotator.clients import (
    AnnotatorMyVariantClient,
    AnnotatorFDAClient,
    AnnotatorCGIClient,
    AnnotatorOncoTreeClient,
    AnnotatorVICCClient,
    AnnotatorCIViCClient,
    AnnotatorCBioPortalClient,
    AnnotatorDepMapClient,
    AnnotatorHotspotsClient,
    AnnotatorClinicalTrialsClient,
    AnnotatorPubMedClient,
    AnnotatorSemanticScholarClient,
)

logger = get_logger(__name__)


class Annotator:
    """Annotates variants from VCF files.

    Use as an async context manager:

        async with Annotator(config) as annotator:
            result = await annotator.run(vcf_path)
    """

    def __init__(self, config: AnnotatorConfig | None = None):
        self.config = config or AnnotatorConfig()
        self._init_clients()

    def _init_clients(self) -> None:
        """Initialize all API clients based on configuration."""
        # Core clients (always enabled)
        self.myvariant_client = AnnotatorMyVariantClient()
        # self.fda_client = AnnotatorFDAClient()
        # self.cgi_client = AnnotatorCGIClient()
        # self.oncotree_client = AnnotatorOncoTreeClient()
        # self.cbioportal_client = AnnotatorCBioPortalClient()
        # self.depmap_client = AnnotatorDepMapClient()
        # self.hotspots_client = AnnotatorHotspotsClient()
        #
        # # Optional clients based on config
        # self.vicc_client = AnnotatorVICCClient()
        # self.civic_client = AnnotatorCIViCClient()
        # self.clinical_trials_client = AnnotatorClinicalTrialsClient()
        #
        # # Literature clients
        # self.semantic_scholar_client = AnnotatorSemanticScholarClient()
        # self.pubmed_client = AnnotatorPubMedClient()

    def _get_all_clients(self) -> list[Any]:
        """Get list of all active clients for context management."""
        clients = [
            self.myvariant_client,
        ]
        return [c for c in clients if c is not None]

    async def __aenter__(self) -> "Annotator":
        """Initialize HTTP client sessions."""
        for client in self._get_all_clients():
            await client.__aenter__()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb) -> None:
        """Close HTTP client sessions."""
        for client in self._get_all_clients():
            await client.__aexit__(exc_type, exc_val, exc_tb)

    async def run(self, vcf_path: Path) -> AnnotationResult:
        """Annotate variants from a VCF file.

        Args:
            vcf_path: Path to the VCF file

        Returns:
            AnnotationResult with annotated variants
        """
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

        logger.debug(f"Parsing VCF file: {vcf_path}")

        vcf_records = self._parse_vcf(vcf_path)
        parsed_variants = [record.to_parsed_variant() for record in vcf_records]

        logger.info(f"Parsed {len(parsed_variants)} variants from {vcf_path}")

        # Annotate each variant
        annotated_variants = []
        for variant in parsed_variants:
            annotation = await self.build_annotations(variant)
            annotated_variants.append(annotation)

        logger.info(f"Annotated {len(annotated_variants)} variants")

        return AnnotationResult(variants=annotated_variants)

    async def build_annotations(self, variant, tumor_type: str | None = None) -> dict[str, Any]:
        """Build annotations for a parsed variant.

        Args:
            variant: AnnParsedVariant to annotate
            tumor_type: Optional tumor type override

        Returns:
            Annotated variant data
        """
        total_start = time.time()
        source_timings: dict[str, float] = {}

        gene = variant.gene
        protein = variant.protein or variant.variant
        tumor = tumor_type or variant.tumor_type

        # Fetch evidence from MyVariant
        start = time.time()
        myvariant_annotation = await self.myvariant_client.fetch_annotation(gene, protein)
        source_timings["myvariant"] = time.time() - start

        total_time = time.time() - total_start
        logger.debug(f"build_annotations completed in {total_time:.2f}s")

        return {
            "gene": gene,
            "protein": protein,
            "tumor_type": tumor,
            "myvariant": myvariant_annotation.model_dump(),
            "timings": source_timings,
        }

    def _parse_vcf(self, vcf_path: Path) -> list[VCFRecord]:
        """Parse VCF file into VCFRecord objects."""
        records = []

        with open(vcf_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) < 8:
                    logger.warning(f"Skipping malformed line: {line[:50]}...")
                    continue

                chrom, pos, variant_id, ref, alt, qual, filt, info_str = fields[:8]

                info = self._parse_info(info_str)

                record = VCFRecord(
                    chrom=chrom,
                    pos=int(pos),
                    id=variant_id if variant_id != "." else None,
                    ref=ref,
                    alt=alt,
                    info=info,
                )
                records.append(record)

        return records

    def _parse_info(self, info_str: str) -> dict[str, str]:
        """Parse VCF INFO field into a dictionary."""
        if not info_str or info_str == ".":
            return {}

        info = {}
        for item in info_str.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                info[key] = value
            else:
                info[item] = "true"
        return info


__all__ = [
    "Annotator",
]

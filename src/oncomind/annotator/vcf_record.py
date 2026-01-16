"""VCF record data model."""

from dataclasses import dataclass, field

from oncomind.annotator.ann_parsed_variant import AnnParsedVariant


@dataclass
class VCFRecord:
    """A raw VCF record before conversion to ParsedVariant."""

    chrom: str
    pos: int
    id: str | None
    ref: str
    alt: str
    info: dict[str, str] = field(default_factory=dict)

    def to_parsed_variant(self) -> AnnParsedVariant:
        """Convert VCF record to AnnParsedVariant."""
        # Basic variant representation
        variant_id = f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"

        # Check INFO for gene annotation
        gene = "UNKNOWN"
        if self.info:
            gene = self.info.get("GENE", self.info.get("Gene", self.info.get("gene", "UNKNOWN")))

        # Check INFO for protein change (HGVSP or PROTEIN)
        protein = None
        if self.info:
            protein = self.info.get("HGVSP", self.info.get("PROTEIN", self.info.get("protein")))

        # Check INFO for tumor type
        tumor_type = None
        if self.info:
            tumor_type = self.info.get("TUMOR_TYPE", self.info.get("tumor_type"))

        return AnnParsedVariant(
            gene=gene,
            variant=variant_id,
            protein=protein,
            variant_normalized=None,
            variant_type="snv" if len(self.ref) == 1 and len(self.alt) == 1 else "indel",
            tumor_type=tumor_type,
            raw_input=variant_id,
            parse_confidence=0.5 if gene == "UNKNOWN" else 0.8,
            parse_warnings=["VCF parsing is basic - use VEP for full annotation"],
        )


__all__ = [
    "VCFRecord",
]

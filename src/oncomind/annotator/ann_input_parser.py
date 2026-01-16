"""Input parsing for the annotator workflow."""

from oncomind.annotator.ann_parsed_variant import AnnParsedVariant
from oncomind.utils.variant_normalization import VariantNormalizer
from oncomind.config.debug import get_logger

logger = get_logger(__name__)


class AnnInputParser:
    """Parser for converting variant inputs to AnnParsedVariant."""

    def __init__(
        self,
        gene: str,
        variant: str,
        tumor_type: str | None = None,
    ):
        """Initialize the parser with variant information.

        Args:
            gene: Gene symbol (e.g., "BRAF")
            variant: Variant notation (e.g., "V600E")
            tumor_type: Optional tumor type
        """
        self.gene = gene.upper()
        self.variant = variant
        self.tumor_type = tumor_type

    def to_parsed_variant(self) -> AnnParsedVariant:
        """Convert input to AnnParsedVariant."""
        # Get protein notation using VariantNormalizer
        protein_result = VariantNormalizer.normalize_protein_change(self.variant)
        protein = protein_result.get("hgvs_protein")

        # Get variant type
        variant_type = VariantNormalizer.classify_variant_type(self.variant)

        return AnnParsedVariant(
            gene=self.gene,
            variant=self.variant,
            protein=protein,
            variant_normalized=protein_result.get("short_form"),
            variant_type=variant_type,
            tumor_type=self.tumor_type,
            raw_input=f"{self.gene} {self.variant}",
            parse_confidence=1.0,
            parse_warnings=[],
        )


__all__ = [
    "AnnInputParser",
]

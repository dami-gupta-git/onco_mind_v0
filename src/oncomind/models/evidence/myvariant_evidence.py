"""Evidence data from MyVariant.info API.

This is a simple data class that holds the parsed response from
the MyVariant.info API. It's used internally by MyVariantClient
and EvidenceAggregator.
"""

from typing import Any

from pydantic import BaseModel, Field

from oncomind.models.evidence.civic import CIViCEvidence
from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.cosmic import COSMICEvidence


class MyVariantEvidence(BaseModel):
    """Parsed evidence from MyVariant.info API response.

    This is a simple data container for the MyVariant API response.
    Used by EvidenceAggregator to extract evidence for the Insight model.
    """

    # Core identifiers
    variant_id: str
    gene: str
    variant: str

    # Database identifiers
    cosmic_id: str | None = None
    ncbi_gene_id: str | None = None
    dbsnp_id: str | None = None
    clinvar_id: str | None = None
    clinvar_clinical_significance: str | None = None
    clinvar_accession: str | None = None

    # HGVS notations
    hgvs_genomic: str | None = None
    hgvs_protein: str | None = None
    hgvs_transcript: str | None = None

    # Functional annotations
    snpeff_effect: str | None = None
    polyphen2_prediction: str | None = None
    polyphen2_score: float | None = None
    sift_prediction: str | None = None
    sift_score: float | None = None
    cadd_score: float | None = None
    gnomad_exome_af: float | None = None
    alphamissense_score: float | None = None
    alphamissense_prediction: str | None = None
    transcript_id: str | None = None
    transcript_consequence: str | None = None

    # Evidence lists from databases
    civic: list[CIViCEvidence] = Field(default_factory=list)
    clinvar: list[ClinVarEvidence] = Field(default_factory=list)
    cosmic: list[COSMICEvidence] = Field(default_factory=list)

    # Raw data for debugging
    raw_data: dict[str, Any] = Field(default_factory=dict)

    def has_evidence(self) -> bool:
        """Check if any evidence was found."""
        return bool(self.civic or self.clinvar or self.cosmic)

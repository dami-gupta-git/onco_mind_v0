from pydantic import Field

from oncomind.models.evidence.base import EvidenceItemBase


class COSMICEvidence(EvidenceItemBase):
    """Evidence from COSMIC (Catalogue of Somatic Mutations in Cancer)."""

    mutation_id: str | None = None
    primary_site: str | None = None
    site_subtype: str | None = None
    primary_histology: str | None = None
    histology_subtype: str | None = None
    sample_count: int | None = None
    mutation_somatic_status: str | None = None

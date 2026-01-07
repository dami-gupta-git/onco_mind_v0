from pydantic import Field

from oncomind.models.evidence.base import EvidenceItemBase, EvidenceLevel


class CGIBiomarkerEvidence(EvidenceItemBase):
    """Evidence from Cancer Genome Interpreter biomarkers database."""

    gene: str | None = None
    alteration: str | None = None
    drug: str | None = None
    drug_status: str | None = None
    association: str | None = None
    evidence_level: str | None = None
    source: str | None = None
    tumor_type: str | None = None
    fda_approved: bool = False
    # Additional match tracking fields (not in base class)
    matched_alteration: str | None = Field(
        default=None,
        description="The alteration that was actually matched"
    )

    def to_fda_approval(self) -> "FDAApproval":
        """Convert this CGI biomarker to an FDAApproval object.

        Only meaningful for CGI entries where fda_approved=True
        (evidence_level contains "FDA guidelines" or "NCCN guidelines").

        Returns:
            FDAApproval object with data from this CGI entry
        """
        # Import here to avoid circular imports
        from oncomind.models.evidence.fda import FDAApproval

        # Build indication text from CGI fields
        # Format: "[Gene Alteration] [Association] in [Tumor Type]"
        indication_parts = []
        if self.gene and self.alteration:
            indication_parts.append(f"{self.gene} {self.alteration}")
        if self.association:
            indication_parts.append(f"({self.association})")
        if self.tumor_type:
            indication_parts.append(f"in {self.tumor_type}")

        indication = " ".join(indication_parts) if indication_parts else None

        return FDAApproval(
            drug_name=self.drug,
            generic_name=self.drug,
            indication=indication,
            gene=self.gene,
            variant_in_indications=True,  # CGI explicitly maps gene/variant to drug
            # Preserve locus and cancer match from CGI
            locus_variant_match=self.locus_variant_match,
            cancer_type_match=self.cancer_type_match,
        )


# Forward reference for type hints
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from oncomind.models.evidence.fda import FDAApproval


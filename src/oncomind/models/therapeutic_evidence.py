"""Therapeutic evidence model for research context.

Expands on RecommendedTherapy to include preclinical data, mechanisms,
and research-relevant context beyond FDA approval status.
"""

from pydantic import BaseModel, Field


class TherapeuticEvidence(BaseModel):
    """Evidence for therapeutic relevance at any evidence level.

    This model merges RecommendedTherapy fields with additional research context.

    Original RecommendedTherapy fields (preserved):
    - drug_name
    - evidence_level
    - approval_status
    - clinical_context

    New research fields:
    - response_type, mechanism, source, etc.
    """

    # === Original RecommendedTherapy fields (PRESERVED) ===
    drug_name: str = Field(..., description="Name of the therapeutic agent")
    evidence_level: str | None = Field(
        None,
        description="Level of evidence: FDA-approved | Phase 3 | Phase 2 | Phase 1 | Preclinical | In vitro | Computational | Case report"
    )
    approval_status: str | None = Field(
        None,
        description="FDA approval status: Approved in indication | Approved different histology | Investigational | Off-label | Not approved"
    )
    clinical_context: str | None = Field(
        None,
        description="Clinical context (e.g., first-line, second-line, resistance setting, maintenance)"
    )

    # === New research-focused fields ===

    # Response type
    response_type: str | None = Field(
        None,
        description="Sensitivity | Resistance | Mixed | Unknown"
    )

    # Mechanism (critical for research)
    mechanism: str | None = Field(
        None,
        description="Mechanism of action or resistance (e.g., 'Constitutive kinase activation', 'EGFR feedback reactivation')"
    )

    # Tumor context
    tumor_types_tested: list[str] = Field(
        default_factory=list,
        description="Tumor types where this drug-variant relationship was tested"
    )

    # Preclinical context
    cell_lines_tested: list[str] = Field(
        default_factory=list,
        description="Cell lines where response was observed (e.g., 'A375', 'SK-MEL-28')"
    )

    # Quantitative data (when available)
    ic50_nm: float | None = Field(
        None,
        description="IC50 value in nanomolar (from cell line studies)"
    )
    response_rate_pct: float | None = Field(
        None,
        description="Clinical response rate percentage (from trials)"
    )

    # Source attribution
    source: str | None = Field(
        None,
        description="Evidence source: CIViC | CGI | VICC | FDA | GDSC | CTRP | Literature"
    )
    pmids: list[str] = Field(
        default_factory=list,
        description="PubMed IDs supporting this evidence"
    )
    source_url: str | None = Field(
        None,
        description="URL to source evidence (e.g., CIViC assertion page)"
    )

    # Confidence assessment
    confidence: str = Field(
        "low",
        description="Confidence level: high | moderate | low"
    )

    # Match specificity tracking
    match_level: str | None = Field(
        None,
        description="Level of match: 'variant' (exact), 'codon' (same position), 'gene' (gene-only)"
    )

    # Cancer type specificity tracking
    cancer_specificity: str | None = Field(
        None,
        description="Cancer type specificity: 'cancer_specific' (matches queried tumor), "
                    "'pan_cancer' (tumor-agnostic), or specific cancer name (e.g., 'ovarian cancer') "
                    "when evidence is for a different cancer than queried"
    )

    # === Helper methods ===

    def is_sensitivity(self) -> bool:
        """Check if this represents drug sensitivity."""
        if self.response_type:
            return self.response_type.lower() in ("sensitivity", "response", "responsive")
        return False

    def is_resistance(self) -> bool:
        """Check if this represents drug resistance."""
        if self.response_type:
            return "resist" in self.response_type.lower()
        return False

    def is_fda_approved(self) -> bool:
        """Check if this is FDA-approved in any indication."""
        if self.approval_status:
            return "approved" in self.approval_status.lower()
        if self.evidence_level:
            return self.evidence_level.lower() == "fda-approved"
        return False

    def is_clinical_evidence(self) -> bool:
        """Check if evidence is from clinical studies (vs preclinical)."""
        if self.evidence_level:
            level = self.evidence_level.lower()
            return any(term in level for term in ["fda", "phase", "clinical", "case"])
        return False

    def is_preclinical_evidence(self) -> bool:
        """Check if evidence is from preclinical studies."""
        if self.evidence_level:
            level = self.evidence_level.lower()
            return any(term in level for term in ["preclinical", "in vitro", "cell line", "computational"])
        return False

    def get_evidence_tier(self) -> int:
        """Get numeric tier for sorting (1 = highest evidence).

        Returns:
            1: FDA-approved
            2: Phase 3
            3: Phase 2
            4: Phase 1 / Case reports
            5: Preclinical
            6: In vitro / Computational
            7: Unknown
        """
        if not self.evidence_level:
            return 7

        level = self.evidence_level.lower()

        if "fda" in level:
            return 1
        elif "phase 3" in level or "phase3" in level:
            return 2
        elif "phase 2" in level or "phase2" in level:
            return 3
        elif "phase 1" in level or "phase1" in level or "case" in level:
            return 4
        elif "preclinical" in level:
            return 5
        elif "in vitro" in level or "computational" in level:
            return 6
        else:
            return 7

    def to_display_string(self) -> str:
        """Format for display."""
        parts = [self.drug_name]

        if self.evidence_level:
            parts.append(f"[{self.evidence_level}]")

        if self.response_type:
            parts.append(f"({self.response_type})")

        if self.mechanism:
            parts.append(f"- {self.mechanism}")

        return " ".join(parts)


# Backwards compatibility alias
RecommendedTherapy = TherapeuticEvidence

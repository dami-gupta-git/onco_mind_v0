"""Experimental AMP/ASCO/CAP tier computation.

DISCLAIMER: This is experimental code. The tier assignments produced by
this module are NOT authoritative and should NOT be used for clinical
decision-making. Authoritative tiering requires expert human review.

AMP/ASCO/CAP 2017 Guidelines Summary:
- Tier I: Variants with strong clinical significance (FDA-approved therapies or professional guidelines)
  - Level A: FDA-approved therapies
  - Level B: Well-powered studies
- Tier II: Variants with potential clinical significance
  - Level C: Case studies, small trials
  - Level D: Preclinical evidence
- Tier III: Variants of unknown clinical significance
  - Level A: Not actionable or conflicting evidence
  - Level B: VUS in known cancer genes
- Tier IV: Benign or likely benign variants

References:
- Li MM, et al. J Mol Diagn. 2017;19(1):4-23. doi:10.1016/j.jmoldx.2016.10.002
"""

from dataclasses import dataclass
from typing import Any

from oncomind.models.evidence import Evidence


@dataclass
class TierResult:
    """Result of experimental tier computation.

    Attributes:
        tier: Computed tier (e.g., "Tier I", "Tier II")
        level: Evidence level (e.g., "A", "B", "C", "D")
        tier_level: Combined tier-level (e.g., "Tier I-A")
        confidence: Confidence in the tier assignment (0-1)
        rationale: Explanation for the tier assignment
        supporting_evidence: Key evidence that supports the tier
        caveats: Important caveats about the tier assignment
    """

    tier: str
    level: str | None = None
    tier_level: str | None = None
    confidence: float = 0.0
    rationale: str = ""
    supporting_evidence: list[str] | None = None
    caveats: list[str] | None = None

    def __post_init__(self):
        if self.tier and self.level:
            self.tier_level = f"{self.tier}-{self.level}"
        if self.supporting_evidence is None:
            self.supporting_evidence = []
        if self.caveats is None:
            self.caveats = [
                "This tier is computed experimentally and is NOT authoritative.",
                "Clinical decisions should NOT be based solely on this tier.",
                "Consult with a clinical geneticist or molecular tumor board.",
            ]


def compute_experimental_tier(
    panel: Evidence,
    tumor_type: str | None = None,
) -> TierResult:
    """Compute an experimental AMP/ASCO/CAP tier for a variant.

    IMPORTANT: This is experimental code. The tier assignments are NOT
    authoritative and should NOT be used for clinical decision-making.

    This function uses heuristics based on:
    - FDA approvals in the Insight
    - CIViC assertion tiers
    - CGI biomarker FDA approval status
    - VICC evidence levels
    - ClinVar pathogenicity

    Args:
        panel: Evidence with aggregated evidence
        tumor_type: Optional tumor type for context-specific tiering

    Returns:
        TierResult with experimental tier, confidence, and rationale
    """
    tumor = tumor_type or panel.context.tumor_type
    supporting_evidence: list[str] = []
    rationale_parts: list[str] = []

    # Check for Tier I-A: FDA-approved targeted therapy
    if panel.fda_approvals:
        drugs = panel.get_approved_drugs()
        if drugs:
            supporting_evidence.append(f"FDA-approved drugs: {', '.join(drugs[:3])}")
            rationale_parts.append(
                "FDA-approved targeted therapy exists for this gene/variant."
            )
            return TierResult(
                tier="Tier I",
                level="A",
                confidence=0.8,
                rationale=" ".join(rationale_parts),
                supporting_evidence=supporting_evidence,
            )

    # Check for Tier I: CIViC assertions with Tier I classification
    tier_i_assertions = [
        a for a in panel.civic_assertions
        if a.amp_tier == "Tier I"
    ]
    if tier_i_assertions:
        assertion = tier_i_assertions[0]
        supporting_evidence.append(
            f"CIViC Tier I assertion: {assertion.name}"
        )
        if assertion.nccn_guideline:
            supporting_evidence.append(f"NCCN guideline: {assertion.nccn_guideline}")
            rationale_parts.append(
                "CIViC curated assertion indicates Tier I with NCCN guideline support."
            )
            return TierResult(
                tier="Tier I",
                level="B",
                confidence=0.7,
                rationale=" ".join(rationale_parts),
                supporting_evidence=supporting_evidence,
            )
        else:
            rationale_parts.append("CIViC curated assertion indicates Tier I.")
            return TierResult(
                tier="Tier I",
                level="B",
                confidence=0.6,
                rationale=" ".join(rationale_parts),
                supporting_evidence=supporting_evidence,
            )

    # Check for Tier II: CGI FDA-approved biomarkers
    fda_biomarkers = [b for b in panel.cgi_biomarkers if b.fda_approved]
    if fda_biomarkers:
        biomarker = fda_biomarkers[0]
        supporting_evidence.append(
            f"CGI FDA-approved biomarker: {biomarker.drug} ({biomarker.association})"
        )
        if biomarker.association and "RESIST" in biomarker.association.upper():
            rationale_parts.append(
                "CGI indicates FDA-approved drug with resistance association."
            )
            return TierResult(
                tier="Tier II",
                level="C",
                confidence=0.5,
                rationale=" ".join(rationale_parts),
                supporting_evidence=supporting_evidence,
            )
        else:
            rationale_parts.append("CGI indicates FDA-approved biomarker status.")
            return TierResult(
                tier="Tier I",
                level="A",
                confidence=0.6,
                rationale=" ".join(rationale_parts),
                supporting_evidence=supporting_evidence,
            )

    # Check for Tier II: CIViC Tier II assertions
    tier_ii_assertions = [
        a for a in panel.civic_assertions
        if a.amp_tier == "Tier II"
    ]
    if tier_ii_assertions:
        assertion = tier_ii_assertions[0]
        supporting_evidence.append(f"CIViC Tier II assertion: {assertion.name}")
        rationale_parts.append("CIViC curated assertion indicates Tier II.")
        return TierResult(
            tier="Tier II",
            level="C",
            confidence=0.5,
            rationale=" ".join(rationale_parts),
            supporting_evidence=supporting_evidence,
        )

    # Check for Tier II: VICC Level A/B evidence
    vicc_high = [v for v in panel.vicc_evidence if v.evidence_level in ("A", "B")]
    if vicc_high:
        ev = vicc_high[0]
        supporting_evidence.append(
            f"VICC Level {ev.evidence_level} evidence: {ev.description[:50]}..."
        )
        rationale_parts.append(f"VICC MetaKB has Level {ev.evidence_level} evidence.")
        return TierResult(
            tier="Tier II",
            level="C" if ev.evidence_level == "B" else "B",
            confidence=0.4,
            rationale=" ".join(rationale_parts),
            supporting_evidence=supporting_evidence,
        )

    # Check for clinical trials as Tier II-D (investigational)
    if panel.clinical_trials:
        variant_specific = [t for t in panel.clinical_trials if t.variant_specific]
        if variant_specific:
            trial = variant_specific[0]
            supporting_evidence.append(
                f"Variant-specific trial: {trial.nct_id}"
            )
            rationale_parts.append(
                "Active clinical trial specifically for this variant."
            )
            return TierResult(
                tier="Tier II",
                level="D",
                confidence=0.4,
                rationale=" ".join(rationale_parts),
                supporting_evidence=supporting_evidence,
            )
        else:
            trial = panel.clinical_trials[0]
            supporting_evidence.append(f"Gene-level trial: {trial.nct_id}")
            rationale_parts.append(
                "Active clinical trial for this gene (not variant-specific)."
            )
            return TierResult(
                tier="Tier III",
                level="A",
                confidence=0.3,
                rationale=" ".join(rationale_parts),
                supporting_evidence=supporting_evidence,
            )

    # Check for Tier III-B: VUS in known cancer gene
    if panel.context.gene_role in ("oncogene", "TSG", "tumor_suppressor"):
        supporting_evidence.append(f"Gene role: {panel.context.gene_role}")
        rationale_parts.append(
            f"{panel.identifiers.gene} is a known cancer gene "
            f"({panel.context.gene_role}), but variant-specific evidence is limited."
        )
        return TierResult(
            tier="Tier III",
            level="B",
            confidence=0.3,
            rationale=" ".join(rationale_parts),
            supporting_evidence=supporting_evidence,
        )

    # Check ClinVar for benign classification
    if panel.clinvar_significance:
        sig = panel.clinvar_significance.lower()
        if "benign" in sig and "pathogenic" not in sig:
            supporting_evidence.append(
                f"ClinVar: {panel.clinvar_significance}"
            )
            rationale_parts.append("ClinVar classifies this variant as benign.")
            return TierResult(
                tier="Tier IV",
                level=None,
                confidence=0.5,
                rationale=" ".join(rationale_parts),
                supporting_evidence=supporting_evidence,
            )

    # Default: Unknown significance
    rationale_parts.append(
        "Insufficient evidence to determine clinical significance."
    )
    return TierResult(
        tier="Tier III",
        level="A",
        confidence=0.2,
        rationale=" ".join(rationale_parts),
        supporting_evidence=supporting_evidence,
    )


__all__ = [
    "compute_experimental_tier",
    "TierResult",
]

"""Evidence - Top-level per-variant evidence aggregation model.

This is the primary output of the annotation pipeline, providing a strongly-typed
container for all structured evidence gathered about a variant from databases and APIs.

For LLM narrative synthesis, see LLMInsight.

ARCHITECTURE (Flat Structure):
    Each evidence source is a simple list field on Evidence:
    - fda_approvals: list[FDAApproval]
    - civic_assertions: list[CIViCAssertionEvidence]
    - civic_evidence: list[CIViCEvidence]
    - vicc_evidence: list[VICCEvidence]
    - cgi_biomarkers: list[CGIBiomarkerEvidence]
    - clinvar_entries: list[ClinVarEvidence]
    - cosmic_entries: list[COSMICEvidence]
    - clinical_trials: list[ClinicalTrialEvidence]
    - pubmed_articles: list[PubMedEvidence]
    - preclinical_biomarkers: list[CGIBiomarkerEvidence]
    - early_phase_biomarkers: list[CGIBiomarkerEvidence]

    Plus context/metadata:
    - identifiers: VariantIdentifiers
    - functional: FunctionalScores
    - context: VariantContext
    - clinvar_significance: str | None
    - literature_knowledge: LiteratureKnowledge | None
"""

from __future__ import annotations

from typing import Any
from pydantic import BaseModel, Field

from oncomind.models.evidence.cbioportal import CBioPortalEvidence
from oncomind.models.evidence.civic import CIViCEvidence, CIViCAssertionEvidence
from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.cosmic import COSMICEvidence
from oncomind.models.evidence.fda import FDAApproval
from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
from oncomind.models.evidence.vicc import VICCEvidence
from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
from oncomind.models.evidence.pubmed import PubMedEvidence
from oncomind.models.evidence.literature_knowledge import LiteratureKnowledge
from oncomind.models.evidence.evidence_gaps import EvidenceGaps
from oncomind.models.therapeutic_evidence import TherapeuticEvidence

# Backwards compatibility alias
RecommendedTherapy = TherapeuticEvidence


class VariantIdentifiers(BaseModel):
    """Core variant identifiers and normalized notation."""

    variant_id: str = Field(..., description="Unique identifier (gene:variant)")
    gene: str = Field(..., description="Gene symbol (e.g., BRAF)")
    variant: str = Field(..., description="Variant notation (e.g., V600E)")
    variant_normalized: str | None = Field(None, description="Normalized variant notation")
    variant_type: str | None = Field(None, description="Variant classification (missense, nonsense, etc.)")

    # Database identifiers
    cosmic_id: str | None = Field(None, description="COSMIC mutation ID")
    ncbi_gene_id: str | None = Field(None, description="NCBI Gene ID")
    dbsnp_id: str | None = Field(None, description="dbSNP rsID")
    clinvar_id: str | None = Field(None, description="ClinVar Variation ID")

    # HGVS notations
    hgvs_genomic: str | None = Field(None, description="HGVS genomic notation")
    hgvs_protein: str | None = Field(None, description="HGVS protein notation")
    hgvs_transcript: str | None = Field(None, description="HGVS transcript notation")

    # Transcript info
    transcript_id: str | None = Field(None, description="Transcript ID")
    transcript_consequence: str | None = Field(None, description="Transcript consequence")


class FunctionalScores(BaseModel):
    """Functional prediction scores and computational annotations."""

    # AlphaMissense
    alphamissense_score: float | None = Field(None, description="AlphaMissense pathogenicity score (0-1)")
    alphamissense_prediction: str | None = Field(None, description="AlphaMissense prediction (P/B/A)")

    # CADD
    cadd_score: float | None = Field(None, description="CADD PHRED score")
    cadd_raw: float | None = Field(None, description="CADD raw score")

    # PolyPhen2
    polyphen2_prediction: str | None = Field(None, description="PolyPhen2 prediction")
    polyphen2_score: float | None = Field(None, description="PolyPhen2 score")

    # SIFT
    sift_prediction: str | None = Field(None, description="SIFT prediction")
    sift_score: float | None = Field(None, description="SIFT score")

    # SnpEff
    snpeff_effect: str | None = Field(None, description="SnpEff effect annotation")
    snpeff_impact: str | None = Field(None, description="SnpEff impact (HIGH/MODERATE/LOW/MODIFIER)")

    # SpliceAI (future)
    spliceai_score: float | None = Field(None, description="SpliceAI delta score")
    spliceai_prediction: str | None = Field(None, description="SpliceAI prediction")

    # Population frequencies
    gnomad_exome_af: float | None = Field(None, description="gnomAD exome allele frequency")
    gnomad_genome_af: float | None = Field(None, description="gnomAD genome allele frequency")

    def get_pathogenicity_summary(self) -> str:
        """Generate a concise pathogenicity summary from available scores."""
        parts = []

        if self.alphamissense_prediction:
            pred_map = {"P": "Pathogenic", "B": "Benign", "A": "Ambiguous"}
            pred = pred_map.get(self.alphamissense_prediction, self.alphamissense_prediction)
            score_str = f" ({self.alphamissense_score:.2f})" if self.alphamissense_score else ""
            parts.append(f"AlphaMissense: {pred}{score_str}")

        if self.cadd_score is not None:
            parts.append(f"CADD: {self.cadd_score:.1f}")

        if self.polyphen2_prediction:
            parts.append(f"PolyPhen2: {self.polyphen2_prediction}")

        if self.sift_prediction:
            parts.append(f"SIFT: {self.sift_prediction}")

        return " | ".join(parts) if parts else "No functional predictions available"


class VariantContext(BaseModel):
    """Context for variant interpretation (tumor type, gene role, mutation class)."""

    tumor_type: str | None = Field(None, description="Tumor type for context")
    tumor_type_resolved: str | None = Field(None, description="OncoTree-resolved tumor type")

    # Gene context
    gene_role: str | None = Field(None, description="Gene role (oncogene, TSG, DDR, MMR, etc.)")
    gene_class: str | None = Field(None, description="Gene functional class from context")
    mutation_class: str | None = Field(None, description="Mutation class (e.g., BRAF Class I/II/III)")
    pathway: str | None = Field(None, description="Associated pathway if relevant")


class Evidence(BaseModel):
    """Top-level evidence aggregation for a single variant.

    This is the primary output of the OncoMind annotation pipeline.
    Each evidence source is a simple list field - the frontend decides how to display them.

    Structure:
    - identifiers: Core variant info and database IDs
    - functional: Computational predictions (AlphaMissense, CADD, etc.)
    - context: Variant context (tumor type, gene role, mutation class)

    Evidence Lists (one per source):
    - fda_approvals: FDA drug approvals
    - civic_assertions: CIViC curated assertions
    - civic_evidence: CIViC evidence items
    - vicc_evidence: VICC MetaKB evidence
    - cgi_biomarkers: CGI FDA-approved biomarkers
    - clinvar_entries: ClinVar entries
    - cosmic_entries: COSMIC entries
    - clinical_trials: Clinical trials
    - pubmed_articles: PubMed articles
    - preclinical_biomarkers: CGI preclinical biomarkers
    - early_phase_biomarkers: CGI early phase biomarkers

    Example:
        >>> evidence = build_evidence("BRAF V600E", tumor_type="Melanoma")
        >>> print(evidence.identifiers.gene, evidence.identifiers.variant)
        BRAF V600E
        >>> print([a.brand_name for a in evidence.fda_approvals])
        ['Tafinlar', 'Zelboraf', 'Braftovi']
    """

    # Core info
    identifiers: VariantIdentifiers = Field(..., description="Variant identifiers and notation")
    functional: FunctionalScores = Field(
        default_factory=FunctionalScores, description="Functional predictions"
    )
    context: VariantContext = Field(
        default_factory=VariantContext, description="Variant context (tumor type, gene role)"
    )

    # === Evidence lists (one per source) ===

    # FDA
    fda_approvals: list[FDAApproval] = Field(default_factory=list, description="FDA drug approvals")

    # CIViC
    civic_assertions: list[CIViCAssertionEvidence] = Field(
        default_factory=list, description="CIViC curated assertions"
    )
    civic_evidence: list[CIViCEvidence] = Field(
        default_factory=list, description="CIViC evidence items"
    )

    # VICC
    vicc_evidence: list[VICCEvidence] = Field(default_factory=list, description="VICC MetaKB evidence")

    # CGI (FDA-approved biomarkers)
    cgi_biomarkers: list[CGIBiomarkerEvidence] = Field(
        default_factory=list, description="CGI FDA-approved biomarkers"
    )

    # ClinVar
    clinvar_entries: list[ClinVarEvidence] = Field(default_factory=list, description="ClinVar entries")
    clinvar_significance: str | None = Field(None, description="Primary ClinVar clinical significance")

    # COSMIC
    cosmic_entries: list[COSMICEvidence] = Field(default_factory=list, description="COSMIC entries")

    # Clinical Trials
    clinical_trials: list[ClinicalTrialEvidence] = Field(
        default_factory=list, description="Matching clinical trials"
    )

    # Literature
    pubmed_articles: list[PubMedEvidence] = Field(default_factory=list, description="PubMed articles")
    literature_knowledge: LiteratureKnowledge | None = Field(
        None, description="LLM-extracted structured knowledge from literature"
    )

    # Research (CGI preclinical/early phase)
    preclinical_biomarkers: list[CGIBiomarkerEvidence] = Field(
        default_factory=list, description="CGI preclinical biomarkers"
    )
    early_phase_biomarkers: list[CGIBiomarkerEvidence] = Field(
        default_factory=list, description="CGI early phase biomarkers"
    )

    # cBioPortal (co-mutation and prevalence)
    cbioportal_evidence: CBioPortalEvidence | None = Field(
        None, description="cBioPortal co-mutation and prevalence data"
    )

    # Evidence gaps (computed, not fetched)
    evidence_gaps: EvidenceGaps | None = Field(
        None, description="Detected evidence gaps for research prioritization"
    )

    # === Helper methods ===

    def compute_evidence_gaps(self) -> EvidenceGaps:
        """Compute and cache evidence gaps."""
        from oncomind.insight_builder.gap_detector import detect_evidence_gaps
        self.evidence_gaps = detect_evidence_gaps(self)
        return self.evidence_gaps

    def has_evidence(self) -> bool:
        """Check if any evidence was found."""
        return bool(
            self.fda_approvals or
            self.civic_assertions or
            self.civic_evidence or
            self.vicc_evidence or
            self.cgi_biomarkers or
            self.clinvar_entries or
            self.clinvar_significance or
            self.cosmic_entries or
            self.clinical_trials or
            self.pubmed_articles or
            self.preclinical_biomarkers or
            self.early_phase_biomarkers or
            (self.cbioportal_evidence and self.cbioportal_evidence.has_data())
        )

    def get_evidence_sources(self) -> list[str]:
        """Get list of sources that have evidence."""
        sources = []
        if self.fda_approvals:
            sources.append("FDA")
        if self.civic_assertions or self.civic_evidence:
            sources.append("CIViC")
        if self.vicc_evidence:
            sources.append("VICC")
        if self.cgi_biomarkers:
            sources.append("CGI")
        if self.clinvar_entries or self.clinvar_significance:
            sources.append("ClinVar")
        if self.cosmic_entries:
            sources.append("COSMIC")
        if self.clinical_trials:
            sources.append("Clinical Trials")
        if self.pubmed_articles:
            sources.append("Literature")
        if self.preclinical_biomarkers or self.early_phase_biomarkers:
            sources.append("Research")
        if self.cbioportal_evidence and self.cbioportal_evidence.has_data():
            sources.append("cBioPortal")
        return sources

    def get_approved_drugs(self) -> list[str]:
        """Get list of FDA-approved drug names."""
        drugs = []
        for approval in self.fda_approvals:
            name = approval.brand_name or approval.generic_name or approval.drug_name
            if name and name not in drugs:
                drugs.append(name)
        return drugs

    def get_summary(self) -> str:
        """Generate a summary of the variant evidence."""
        parts = []

        # Variant identity with mutation type
        variant_desc = f"{self.identifiers.gene} {self.identifiers.variant}"
        if self.identifiers.variant_type:
            variant_desc += f" is a {self.identifiers.variant_type} variant"
        parts.append(variant_desc)

        # HGVS notation
        if self.identifiers.hgvs_genomic:
            parts.append(f"({self.identifiers.hgvs_genomic})")
        parts.append(".")

        # FDA-approved therapies
        drugs = self.get_approved_drugs()
        if drugs:
            ellipsis = "..." if len(drugs) > 3 else ""
            parts.append(f"It is associated with FDA-approved therapies: {', '.join(drugs[:3])}{ellipsis}.")

        return " ".join(parts)

    def get_summary_stats(self) -> dict[str, Any]:
        """Get summary statistics about the evidence."""
        return {
            "gene": self.identifiers.gene,
            "variant": self.identifiers.variant,
            "tumor_type": self.context.tumor_type,
            "evidence_sources": self.get_evidence_sources(),
            "fda_approved_drugs": self.get_approved_drugs(),
            "clinical_trials_count": len(self.clinical_trials),
            "pubmed_articles_count": len(self.pubmed_articles),
            "has_evidence": self.has_evidence(),
        }

    def get_recommended_therapies(self, max_results: int = 10) -> list[TherapeuticEvidence]:
        """Backwards-compatible alias for get_therapeutic_evidence.

        DEPRECATED: Use get_therapeutic_evidence() instead.
        """
        return self.get_therapeutic_evidence(include_preclinical=False, max_results=max_results)

    def get_therapeutic_evidence(
        self,
        include_preclinical: bool = True,
        max_results: int = 20
    ) -> list[TherapeuticEvidence]:
        """Get all therapeutic evidence at all evidence levels.

        This replaces get_recommended_therapies() with a research-focused method
        that includes preclinical and investigational evidence.

        Args:
            include_preclinical: Include preclinical/in vitro evidence (default True)
            max_results: Maximum results to return

        Returns:
            List of TherapeuticEvidence sorted by evidence tier
        """
        evidence_list: list[TherapeuticEvidence] = []
        seen_drugs: set[str] = set()

        # From FDA approvals (Tier 1)
        for approval in self.fda_approvals:
            drug_key = (approval.generic_name or approval.brand_name or approval.drug_name or "").lower()
            if drug_key and drug_key not in seen_drugs:
                seen_drugs.add(drug_key)

                drug_name = approval.brand_name or approval.generic_name or approval.drug_name
                if approval.generic_name and approval.brand_name:
                    drug_name = f"{approval.generic_name} ({approval.brand_name})"

                evidence_list.append(TherapeuticEvidence(
                    drug_name=drug_name,
                    evidence_level="FDA-approved",
                    approval_status="Approved in indication" if approval.variant_in_indications else "Approved",
                    clinical_context=self._extract_line_of_therapy(approval),
                    response_type="Sensitivity",
                    mechanism=None,
                    tumor_types_tested=[self.context.tumor_type] if self.context.tumor_type else [],
                    source="FDA",
                    confidence="high",
                ))

        # From CIViC assertions (Tier 1-2)
        for assertion in self.civic_assertions:
            if assertion.therapies:
                for therapy in assertion.therapies:
                    drug_key = therapy.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        # Determine evidence level from AMP tier
                        if assertion.amp_tier == "Tier I":
                            evidence_level = "FDA-approved" if assertion.fda_companion_test else "Phase 3"
                        elif assertion.amp_tier == "Tier II":
                            evidence_level = "Phase 2"
                        else:
                            evidence_level = "Case report"

                        evidence_list.append(TherapeuticEvidence(
                            drug_name=therapy,
                            evidence_level=evidence_level,
                            approval_status=self._get_approval_status_from_tier(assertion.amp_tier),
                            clinical_context=assertion.disease,
                            response_type="Sensitivity" if assertion.significance and "SENSITIV" in assertion.significance.upper() else "Resistance" if assertion.significance and "RESIST" in assertion.significance.upper() else None,
                            mechanism=None,
                            tumor_types_tested=[assertion.disease] if assertion.disease else [],
                            source="CIViC",
                            confidence="high" if assertion.amp_tier == "Tier I" else "moderate",
                        ))

        # From CGI biomarkers - FDA approved
        for biomarker in self.cgi_biomarkers:
            if biomarker.fda_approved and biomarker.drug:
                drug_key = biomarker.drug.lower()
                if drug_key not in seen_drugs:
                    seen_drugs.add(drug_key)

                    response_type = None
                    if biomarker.association:
                        if "RESIST" in biomarker.association.upper():
                            response_type = "Resistance"
                        else:
                            response_type = "Sensitivity"

                    evidence_list.append(TherapeuticEvidence(
                        drug_name=biomarker.drug,
                        evidence_level="FDA-approved",
                        approval_status="FDA Approved (CGI)",
                        clinical_context=biomarker.tumor_type,
                        response_type=response_type,
                        mechanism=None,
                        tumor_types_tested=[biomarker.tumor_type] if biomarker.tumor_type else [],
                        source="CGI",
                        confidence="high",
                    ))

        # From VICC evidence
        for vicc in self.vicc_evidence:
            if vicc.drugs:
                for drug in vicc.drugs:
                    drug_key = drug.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        evidence_list.append(TherapeuticEvidence(
                            drug_name=drug,
                            evidence_level=self._vicc_level_to_evidence_level(vicc.evidence_level),
                            approval_status=self._get_approval_from_vicc(vicc),
                            clinical_context=vicc.disease,
                            response_type="Sensitivity" if vicc.is_sensitivity else "Resistance" if vicc.is_resistance else None,
                            mechanism=None,
                            tumor_types_tested=[vicc.disease] if vicc.disease else [],
                            source=f"VICC ({vicc.source})" if vicc.source else "VICC",
                            confidence="moderate" if vicc.evidence_level in ("A", "B") else "low",
                        ))

        # From preclinical biomarkers (if requested)
        if include_preclinical:
            for biomarker in self.preclinical_biomarkers:
                if biomarker.drug:
                    drug_key = biomarker.drug.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        response_type = None
                        if biomarker.association:
                            if "RESIST" in biomarker.association.upper():
                                response_type = "Resistance"
                            else:
                                response_type = "Sensitivity"

                        evidence_list.append(TherapeuticEvidence(
                            drug_name=biomarker.drug,
                            evidence_level="Preclinical",
                            approval_status="Investigational",
                            clinical_context=biomarker.tumor_type,
                            response_type=response_type,
                            mechanism=None,
                            tumor_types_tested=[biomarker.tumor_type] if biomarker.tumor_type else [],
                            source="CGI (preclinical)",
                            confidence="low",
                        ))

        # Sort by evidence tier
        evidence_list.sort(key=lambda x: x.get_evidence_tier())

        return evidence_list[:max_results]

    def get_therapeutic_evidence_by_level(self) -> dict[str, list[TherapeuticEvidence]]:
        """Group therapeutic evidence by evidence level.

        Returns:
            Dict with keys: 'fda_approved', 'clinical', 'preclinical'
        """
        all_evidence = self.get_therapeutic_evidence(include_preclinical=True)

        return {
            "fda_approved": [e for e in all_evidence if e.is_fda_approved()],
            "clinical": [e for e in all_evidence if e.is_clinical_evidence() and not e.is_fda_approved()],
            "preclinical": [e for e in all_evidence if e.is_preclinical_evidence()],
        }

    def get_resistance_evidence(self) -> list[TherapeuticEvidence]:
        """Get therapeutic evidence specifically for resistance."""
        all_evidence = self.get_therapeutic_evidence(include_preclinical=True)
        return [e for e in all_evidence if e.is_resistance()]

    def get_sensitivity_evidence(self) -> list[TherapeuticEvidence]:
        """Get therapeutic evidence specifically for sensitivity."""
        all_evidence = self.get_therapeutic_evidence(include_preclinical=True)
        return [e for e in all_evidence if e.is_sensitivity()]

    # Helper methods for get_therapeutic_evidence
    def _extract_line_of_therapy(self, approval) -> str | None:
        """Extract line of therapy from FDA approval."""
        if self.context.tumor_type and approval.indication:
            parsed = approval.parse_indication_for_tumor(self.context.tumor_type)
            if parsed.get('tumor_match'):
                return parsed.get('line_of_therapy')
        return None

    def _get_approval_status_from_tier(self, amp_tier: str | None) -> str:
        """Map AMP tier to approval status."""
        if amp_tier == "Tier I":
            return "Approved in indication"
        elif amp_tier == "Tier II":
            return "Approved different histology"
        else:
            return "Investigational"

    def _vicc_level_to_evidence_level(self, level: str | None) -> str:
        """Map VICC evidence level to standard format."""
        if not level:
            return "Unknown"
        if level == "A":
            return "FDA-approved"
        elif level == "B":
            return "Phase 3"
        elif level == "C":
            return "Phase 2"
        elif level == "D":
            return "Preclinical"
        return "Unknown"

    def _get_approval_from_vicc(self, vicc) -> str:
        """Get approval status from VICC evidence."""
        if vicc.evidence_level == "A":
            return "FDA Approved"
        elif vicc.evidence_level == "B":
            return "Approved different histology"
        return "Investigational"

    def to_flat_dict(self) -> dict[str, Any]:
        """Convert to a flat dictionary for DataFrame/tabular output."""
        return {
            # Identifiers
            "gene": self.identifiers.gene,
            "variant": self.identifiers.variant,
            "variant_normalized": self.identifiers.variant_normalized,
            "variant_type": self.identifiers.variant_type,
            "cosmic_id": self.identifiers.cosmic_id,
            "dbsnp_id": self.identifiers.dbsnp_id,
            "clinvar_id": self.identifiers.clinvar_id,
            "hgvs_protein": self.identifiers.hgvs_protein,

            # Clinical
            "tumor_type": self.context.tumor_type,
            "gene_role": self.context.gene_role,
            "clinvar_significance": self.clinvar_significance,
            "fda_approved_drugs": ", ".join(self.get_approved_drugs()),
            "clinical_trials_count": len(self.clinical_trials),

            # Functional
            "alphamissense_score": self.functional.alphamissense_score,
            "alphamissense_prediction": self.functional.alphamissense_prediction,
            "cadd_score": self.functional.cadd_score,
            "gnomad_af": self.functional.gnomad_exome_af,

            # Counts
            "civic_evidence_count": len(self.civic_assertions) + len(self.civic_evidence),
            "vicc_evidence_count": len(self.vicc_evidence),
            "pubmed_articles_count": len(self.pubmed_articles),
        }

    def get_evidence_summary_for_llm(self) -> str:
        """Generate a compact evidence summary for LLM prompt consumption."""
        lines = [f"Evidence for {self.identifiers.gene} {self.identifiers.variant}:\n"]
        tumor_type = self.context.tumor_type

        # FDA Approvals
        if self.fda_approvals:
            lines.append(f"FDA Approved Drugs ({len(self.fda_approvals)}):")
            for approval in self.fda_approvals[:5]:
                drug = approval.brand_name or approval.generic_name or approval.drug_name
                if tumor_type:
                    parsed = approval.parse_indication_for_tumor(tumor_type)
                    if parsed['tumor_match']:
                        line_info = parsed['line_of_therapy'].upper()
                        approval_info = parsed['approval_type'].upper()
                        lines.append(f"  - {drug} [FOR {tumor_type.upper()}]:")
                        lines.append(f"      Line of therapy: {line_info}")
                        lines.append(f"      Approval type: {approval_info}")
                        lines.append(f"      Excerpt: {parsed['indication_excerpt'][:200]}...")
                    else:
                        indication = (approval.indication or "")[:300]
                        lines.append(f"  - {drug} [DIFFERENT TUMOR TYPE]: {indication}...")
                else:
                    indication = (approval.indication or "")[:300]
                    date_str = f" (Approved: {approval.approval_date})" if approval.approval_date else ""
                    lines.append(f"  - {drug}{date_str}: {indication}...")
            lines.append("")

        # CGI Biomarkers
        if self.cgi_biomarkers:
            approved = [b for b in self.cgi_biomarkers if b.fda_approved]
            if approved:
                resistance = [b for b in approved if b.association and 'RESIST' in b.association.upper()]
                sensitivity = [b for b in approved if b.association and 'RESIST' not in b.association.upper()]

                if resistance:
                    lines.append(f"CGI FDA-APPROVED RESISTANCE MARKERS ({len(resistance)}):")
                    for b in resistance[:5]:
                        lines.append(f"  - {b.drug} [{b.association.upper()}] in {b.tumor_type or 'solid tumors'} - Evidence: {b.evidence_level}")
                    lines.append("")

                if sensitivity:
                    lines.append(f"CGI FDA-Approved Sensitivity Biomarkers ({len(sensitivity)}):")
                    for b in sensitivity[:5]:
                        lines.append(f"  - {b.drug} [{b.association}] in {b.tumor_type or 'solid tumors'} - Evidence: {b.evidence_level}")
                    lines.append("")

        # CIViC Assertions
        if self.civic_assertions:
            predictive_tier_i = [a for a in self.civic_assertions
                                  if a.amp_tier == "Tier I" and a.assertion_type == "PREDICTIVE"]
            if predictive_tier_i:
                lines.append(f"CIViC PREDICTIVE TIER I ASSERTIONS ({len(predictive_tier_i)}):")
                for a in predictive_tier_i[:5]:
                    therapies = ", ".join(a.therapies) if a.therapies else "N/A"
                    lines.append(f"  - {a.molecular_profile}: {therapies} [{a.significance}]")
                    lines.append(f"      AMP Level: {a.amp_level}, Disease: {a.disease}")
                lines.append("")

        # ClinVar
        if self.clinvar_significance:
            lines.append(f"ClinVar: {self.clinvar_significance}")
            lines.append("")

        # PubMed Literature
        if self.pubmed_articles:
            resistance_articles = [a for a in self.pubmed_articles if a.is_resistance_evidence()]
            if resistance_articles:
                lines.append(f"PUBMED RESISTANCE LITERATURE ({len(resistance_articles)} articles):")
                for article in resistance_articles[:3]:
                    drugs_str = f" [Drugs: {', '.join(article.drugs_mentioned[:3])}]" if article.drugs_mentioned else ""
                    lines.append(f"  - PMID {article.pmid}: {article.title[:100]}...{drugs_str}")
                lines.append("")

        return "\n".join(lines) if len(lines) > 1 else ""

    def get_biological_context_for_llm(self) -> str:
        """Get biological context formatted for LLM prompt.

        Combines cBioPortal data with gene context for research-focused synthesis.

        Returns:
            Formatted string with gene role, pathway, and cBioPortal data
        """
        lines = []

        # Gene context
        if self.context.gene_role:
            role_descriptions = {
                "oncogene": "an oncogene (gain-of-function mutations drive cancer)",
                "TSG": "a tumor suppressor gene (loss-of-function mutations drive cancer)",
                "tumor_suppressor": "a tumor suppressor gene (loss-of-function mutations drive cancer)",
                "tsg_pathway_actionable": "a pathway-actionable tumor suppressor (loss activates druggable pathway)",
                "ddr": "a DNA damage repair gene (loss creates synthetic lethality with PARP inhibitors)",
            }
            role_desc = role_descriptions.get(self.context.gene_role, self.context.gene_role)
            lines.append(f"GENE ROLE: {self.identifiers.gene} is {role_desc}")

            if self.context.pathway:
                lines.append(f"PATHWAY: {self.context.pathway}")
            lines.append("")

        # cBioPortal data
        if self.cbioportal_evidence and self.cbioportal_evidence.has_data():
            lines.append(self.cbioportal_evidence.to_prompt_context())
        else:
            lines.append("PREVALENCE: No cBioPortal data available")
            lines.append("")

        # Mutation class (for oncogenes like BRAF)
        if self.context.mutation_class:
            lines.append(f"MUTATION CLASS: {self.context.mutation_class}")
            lines.append("")

        return "\n".join(lines)

    def get_literature_summary_for_llm(self) -> str:
        """Get literature findings formatted for LLM prompt.

        Includes PubMed articles with signal types and extracted knowledge.

        Returns:
            Formatted string with literature findings for research synthesis
        """
        if not self.pubmed_articles and not self.literature_knowledge:
            return ""

        lines = []

        # PubMed articles
        if self.pubmed_articles:
            # Group by signal type
            resistance_articles = [a for a in self.pubmed_articles if a.is_resistance_evidence()]
            sensitivity_articles = [a for a in self.pubmed_articles if a.is_sensitivity_evidence()]
            other_articles = [a for a in self.pubmed_articles
                             if not a.is_resistance_evidence() and not a.is_sensitivity_evidence()]

            if resistance_articles:
                lines.append(f"RESISTANCE LITERATURE ({len(resistance_articles)} articles):")
                for article in resistance_articles[:5]:
                    drugs_str = f" [Drugs: {', '.join(article.drugs_mentioned[:3])}]" if article.drugs_mentioned else ""
                    summary = article.get_best_summary(200)
                    impact = article.get_impact_indicator()
                    impact_str = f" [{impact}]" if impact else ""
                    lines.append(f"  - PMID {article.pmid}: {article.title[:80]}...")
                    lines.append(f"    Signal: {article.signal_type}{drugs_str}{impact_str}")
                    if summary:
                        lines.append(f"    Summary: {summary}")
                lines.append("")

            if sensitivity_articles:
                lines.append(f"SENSITIVITY LITERATURE ({len(sensitivity_articles)} articles):")
                for article in sensitivity_articles[:5]:
                    drugs_str = f" [Drugs: {', '.join(article.drugs_mentioned[:3])}]" if article.drugs_mentioned else ""
                    summary = article.get_best_summary(200)
                    lines.append(f"  - PMID {article.pmid}: {article.title[:80]}...")
                    lines.append(f"    Signal: {article.signal_type}{drugs_str}")
                    if summary:
                        lines.append(f"    Summary: {summary}")
                lines.append("")

            if other_articles:
                lines.append(f"OTHER RELEVANT LITERATURE ({len(other_articles)} articles):")
                for article in other_articles[:3]:
                    summary = article.get_best_summary(150)
                    lines.append(f"  - PMID {article.pmid}: {article.title[:80]}...")
                    if summary:
                        lines.append(f"    Summary: {summary}")
                lines.append("")

        # Extracted structured knowledge
        if self.literature_knowledge:
            lk = self.literature_knowledge
            lines.append("EXTRACTED KNOWLEDGE (synthesized from literature):")

            if lk.mutation_type != "unknown":
                lines.append(f"  Mutation Type: {lk.mutation_type}")

            if lk.resistant_to:
                drugs = ", ".join(f"{r.drug} ({r.evidence})" for r in lk.resistant_to[:5])
                lines.append(f"  Resistant to: {drugs}")

            if lk.sensitive_to:
                drugs = ", ".join(f"{s.drug} ({s.evidence})" for s in lk.sensitive_to[:5])
                lines.append(f"  Potentially sensitive to: {drugs}")

            if lk.clinical_significance:
                lines.append(f"  Clinical Significance: {lk.clinical_significance}")

            if lk.key_findings:
                lines.append("  Key Findings:")
                for finding in lk.key_findings[:3]:
                    lines.append(f"    • {finding}")

            lines.append("")

        return "\n".join(lines)

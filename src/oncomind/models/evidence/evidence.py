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
from oncomind.models.evidence.depmap import DepMapEvidence
from oncomind.models.evidence.fda import FDAApproval, FDALabelEvidence
from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
from oncomind.models.evidence.hotspots import HotspotsEvidence
from oncomind.models.evidence.vicc import VICCEvidence
from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
from oncomind.models.evidence.pubmed import PubMedEvidence
from oncomind.models.extracted.literature_knowledge import LiteratureKnowledge
from oncomind.models.extracted.evidence_gaps import EvidenceGaps
from oncomind.models.extracted.therapeutic_data import TherapeuticData
from oncomind.models.evidence.base import tumor_types_match, is_pan_cancer_term
from oncomind.config.constants import is_biomarker_selection_drug, is_acquired_resistance_mutation

# Import FDALabelInfo with TYPE_CHECKING to avoid circular import
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from oncomind.api.fda_drugs import FDALabelInfo


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
    fda_labels: list[FDALabelEvidence] = Field(
        default_factory=list,
        description="FDA drug labels with clinical study data, mechanism, adverse reactions"
    )
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

    # DepMap/CCLE (gene dependencies and drug sensitivities)
    depmap_evidence: DepMapEvidence | None = Field(
        None, description="DepMap gene dependency and drug sensitivity data"
    )

    # Cancer Hotspots (MSK recurrent mutation sites)
    hotspots_evidence: HotspotsEvidence | None = Field(
        None, description="Cancer hotspots data from cancerhotspots.org"
    )

    # Evidence gaps (computed, not fetched)
    evidence_gaps: EvidenceGaps | None = Field(
        None, description="Detected evidence gaps for research prioritization"
    )

    # Configuration flags (track what was searched)
    literature_searched: bool = Field(
        False, description="Whether literature search was enabled/performed"
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
            (self.cbioportal_evidence and self.cbioportal_evidence.has_data()) or
            (self.depmap_evidence and self.depmap_evidence.has_data())
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
        if self.depmap_evidence and self.depmap_evidence.has_data():
            sources.append("DepMap")
        return sources

    def get_approved_drugs(self) -> list[str]:
        """Get list of FDA-approved drug names."""
        drugs = []
        for approval in self.fda_approvals:
            name = approval.brand_name or approval.generic_name or approval.drug_name
            if name and name not in drugs:
                drugs.append(name)
        return drugs

    def get_vicc_unique(self) -> list[VICCEvidence]:
        """Get VICC evidence, filtering duplicates intelligently.

        VICC MetaKB aggregates from CIViC, CGI, JAX-CKB, OncoKB, PMKB, MolecularMatch.

        Filtering logic:
        - If we have variant-level CIViC assertions/evidence, filter out CIViC-sourced VICC
        - If we have variant-level CGI biomarkers, filter out CGI-sourced VICC
        - If we have NO variant-level data from a source, keep its gene-level VICC entries

        This ensures gene-level evidence is shown for rare variants that have no
        variant-specific annotations.

        Returns:
            List of VICCEvidence items (deduplicated against direct source data)
        """
        result = []

        # Check if we have variant-level data from each source
        has_civic_data = len(self.civic_assertions) > 0 or len(self.civic_evidence) > 0
        has_cgi_data = len(self.cgi_biomarkers) > 0

        for v in self.vicc_evidence:
            source_lower = v.source.lower()

            # Filter CIViC-sourced entries only if we have direct CIViC data
            if source_lower == "civic" and has_civic_data:
                continue

            # Filter CGI-sourced entries only if we have direct CGI data
            if source_lower == "cgi" and has_cgi_data:
                continue

            result.append(v)

        return result

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

    def get_therapeutic_evidence(
        self,
        include_preclinical: bool = True,
    ) -> list[TherapeuticData]:
        """Get all therapeutic evidence at all evidence levels.

        This replaces get_recommended_therapies() with a research-focused method
        that includes preclinical and investigational evidence.

        Uses get_fda_approved_therapies() as the single source of truth for
        FDA-approved drugs, then adds non-FDA sources.

        Args:
            include_preclinical: Include preclinical/in vitro evidence (default True)

        Returns:
            List of TherapeuticData sorted by evidence tier
        """
        # Start with FDA-approved therapies (single source of truth)
        evidence_list: list[TherapeuticData] = list(self.get_fda_approved_therapies())
        seen_drugs: set[str] = {e.drug_name.lower() for e in evidence_list}

        # From CIViC assertions - non-FDA (Tier I without fda_companion_test, Tier II+)
        for assertion in self.civic_assertions:
            # Skip Tier I with fda_companion_test (already in FDA-approved)
            if assertion.amp_tier == "Tier I" and assertion.fda_companion_test:
                continue

            if assertion.therapies:
                for therapy in assertion.therapies:
                    drug_key = therapy.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        # Determine evidence level from AMP tier
                        if assertion.amp_tier == "Tier I":
                            evidence_level = "Phase 3"
                        elif assertion.amp_tier == "Tier II":
                            evidence_level = "Phase 2"
                        else:
                            evidence_level = "Case report"

                        # Determine cancer specificity
                        cancer_specificity = self._get_cancer_specificity_from_disease(assertion.disease)

                        evidence_list.append(TherapeuticData(
                            drug_name=therapy,
                            evidence_level=evidence_level,
                            approval_status=self._get_approval_status_from_tier(assertion.amp_tier),
                            clinical_context=assertion.disease,
                            response_type="Sensitivity" if assertion.significance and "SENSITIV" in assertion.significance.upper() else "Resistance" if assertion.significance and "RESIST" in assertion.significance.upper() else None,
                            mechanism=None,
                            tumor_types_tested=[assertion.disease] if assertion.disease else [],
                            source="CIViC",
                            source_url=assertion.civic_url,
                            confidence="high" if assertion.amp_tier == "Tier I" else "moderate",
                            locus_match=assertion.locus_match,
                            cancer_specificity=cancer_specificity,
                        ))

        # From CIViC evidence items (individual evidence, not curated assertions)
        for civic in self.civic_evidence:
            if civic.drugs:
                for drug in civic.drugs:
                    drug_key = drug.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        # Determine cancer specificity
                        cancer_specificity = self._get_cancer_specificity_from_disease(civic.disease)

                        # Map CIViC evidence level (A-E) to standard format
                        evidence_level = self._civic_level_to_evidence_level(civic.evidence_level)

                        # Determine response type using computed properties
                        # (these check both clinical_significance AND evidence_direction)
                        response_type = None
                        if civic.is_sensitivity:
                            response_type = "Sensitivity"
                        elif civic.is_resistance:
                            response_type = "Resistance"

                        evidence_list.append(TherapeuticData(
                            drug_name=drug,
                            evidence_level=evidence_level,
                            approval_status=self._get_approval_from_civic_level(civic.evidence_level),
                            clinical_context=civic.disease,
                            response_type=response_type,
                            mechanism=None,
                            tumor_types_tested=[civic.disease] if civic.disease else [],
                            source="CIViC",
                            source_url=civic.civic_url,
                            confidence="moderate" if civic.evidence_level in ("A", "B") else "low",
                            locus_match=civic.locus_match,
                            cancer_specificity=cancer_specificity,
                        ))

        # From VICC evidence
        for vicc in self.vicc_evidence:
            if vicc.drugs:
                for drug in vicc.drugs:
                    drug_key = drug.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        # Build source URL for VICC
                        vicc_url = vicc.publication_url[0] if isinstance(vicc.publication_url, list) and vicc.publication_url else vicc.publication_url

                        # Determine cancer specificity
                        cancer_specificity = self._get_cancer_specificity_from_disease(vicc.disease)

                        evidence_list.append(TherapeuticData(
                            drug_name=drug,
                            evidence_level=self._vicc_level_to_evidence_level(vicc.evidence_level),
                            approval_status=self._get_approval_from_vicc(vicc),
                            clinical_context=vicc.disease,
                            response_type="Sensitivity" if vicc.is_sensitivity else "Resistance" if vicc.is_resistance else None,
                            mechanism=None,
                            tumor_types_tested=[vicc.disease] if vicc.disease else [],
                            source=f"VICC ({vicc.source})" if vicc.source else "VICC",
                            source_url=vicc_url,
                            confidence="moderate" if vicc.evidence_level in ("A", "B") else "low",
                            locus_match=vicc.locus_match,
                            cancer_specificity=cancer_specificity,
                        ))

        # From preclinical biomarkers (if requested)
        if include_preclinical:
            for biomarker in self.preclinical_biomarkers:
                # Ensure drug is a non-empty string (not a list or empty)
                if biomarker.drug and isinstance(biomarker.drug, str):
                    drug_key = biomarker.drug.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        # Determine response type from association field
                        # CGI values: "Responsive", "Resistant", "No Responsive", "Increased Toxicity"
                        # Only mark as Sensitivity/Resistance if explicitly stated
                        response_type = None
                        if biomarker.association:
                            assoc_upper = biomarker.association.upper()
                            if "RESIST" in assoc_upper:
                                response_type = "Resistance"
                            elif assoc_upper == "RESPONSIVE":
                                # Only explicit "Responsive" (not "No Responsive" or toxicity)
                                response_type = "Sensitivity"

                        # Determine cancer specificity
                        cancer_specificity = self._get_cancer_specificity_from_disease(biomarker.tumor_type)

                        evidence_list.append(TherapeuticData(
                            drug_name=biomarker.drug,
                            evidence_level="Preclinical",
                            approval_status="Investigational",
                            clinical_context=biomarker.tumor_type,
                            response_type=response_type,
                            mechanism=None,
                            tumor_types_tested=[biomarker.tumor_type] if biomarker.tumor_type else [],
                            source="CGI (preclinical)",
                            confidence="low",
                            locus_match=biomarker.locus_match,
                            cancer_specificity=cancer_specificity,
                        ))

            # Also include early-phase biomarkers (clinical trials, case reports)
            for biomarker in self.early_phase_biomarkers:
                if biomarker.drug and isinstance(biomarker.drug, str):
                    drug_key = biomarker.drug.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        # Determine response type from association field
                        # CGI values: "Responsive", "Resistant", "No Responsive", "Increased Toxicity"
                        # Only mark as Sensitivity/Resistance if explicitly stated
                        response_type = None
                        if biomarker.association:
                            assoc_upper = biomarker.association.upper()
                            if "RESIST" in assoc_upper:
                                response_type = "Resistance"
                            elif assoc_upper == "RESPONSIVE":
                                # Only explicit "Responsive" (not "No Responsive" or toxicity)
                                response_type = "Sensitivity"

                        cancer_specificity = self._get_cancer_specificity_from_disease(biomarker.tumor_type)

                        evidence_list.append(TherapeuticData(
                            drug_name=biomarker.drug,
                            evidence_level=biomarker.evidence_level or "Early Phase",
                            approval_status="Investigational",
                            clinical_context=biomarker.tumor_type,
                            response_type=response_type,
                            mechanism=None,
                            tumor_types_tested=[biomarker.tumor_type] if biomarker.tumor_type else [],
                            source="CGI (early phase)",
                            confidence="low",
                            locus_match=biomarker.locus_match,
                            cancer_specificity=cancer_specificity,
                        ))

            # Include drugs from active clinical trials
            for trial in self.clinical_trials:
                if not trial.is_phase2_or_later():
                    continue  # Only include Phase 2+ trials for therapeutic evidence
                drugs = trial.get_drug_names()
                for drug in drugs:
                    drug_key = drug.lower()
                    if drug_key not in seen_drugs:
                        seen_drugs.add(drug_key)

                        # Determine evidence level from trial phase
                        if trial.phase and "PHASE3" in trial.phase.upper().replace(" ", ""):
                            evidence_level = "Phase 3"
                        else:
                            evidence_level = "Phase 2"

                        # Use the cancer_specificity property from the API client
                        # (returns "cancer_specific", specific disease name, or None)
                        cancer_specificity = trial.cancer_specificity

                        evidence_list.append(TherapeuticData(
                            drug_name=drug,
                            evidence_level=evidence_level,
                            approval_status="Clinical Trial",
                            clinical_context=", ".join(trial.conditions[:2]) if trial.conditions else None,
                            response_type=None,  # Trials don't specify response type
                            mechanism=None,
                            tumor_types_tested=trial.conditions,
                            source=f"ClinicalTrials.gov ({trial.nct_id})",
                            source_url=trial.url,
                            confidence="low",  # Trial data, not yet published results
                            locus_match=getattr(trial, 'match_scope', None),
                            cancer_specificity=cancer_specificity,
                        ))

        # Sort by evidence tier
        evidence_list.sort(key=lambda x: x.get_evidence_tier())

        return evidence_list

    def get_therapeutic_evidence_by_level(self) -> dict[str, list[TherapeuticData]]:
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

    def get_resistance_evidence(self) -> list[TherapeuticData]:
        """Get therapeutic evidence specifically for resistance."""
        all_evidence = self.get_therapeutic_evidence(include_preclinical=True)
        return [e for e in all_evidence if e.is_resistance()]

    def get_sensitivity_evidence(self) -> list[TherapeuticData]:
        """Get therapeutic evidence specifically for sensitivity."""
        all_evidence = self.get_therapeutic_evidence(include_preclinical=True)
        return [e for e in all_evidence if e.is_sensitivity()]

    def get_fda_approved_therapies(self) -> list[TherapeuticData]:
        """Get all FDA-approved therapies with no limit.

        Sources:
        - FDA approvals (direct)
        - CIViC assertions: Tier I with fda_companion_test=True
        - CGI biomarkers: fda_approved=True

        This is the single source of truth for FDA-approved drug counting.
        Used by both the Therapies tab and Gap Analysis.

        Returns:
            List of TherapeuticData with evidence_level="FDA-approved"
        """
        evidence_list: list[TherapeuticData] = []
        # Use (drug, association) as key to allow same drug with different responses
        # e.g., Imatinib Responsive + Imatinib Resistant should both appear
        seen_entries: set[tuple[str, str]] = set()

        # From FDA approvals
        for approval in self.fda_approvals:
            # For combinations (contain "+"), use drug_name as key; otherwise use generic/brand
            is_combination = approval.drug_name and "+" in approval.drug_name
            if is_combination:
                drug_key = approval.drug_name.lower()
            else:
                drug_key = (approval.generic_name or approval.brand_name or approval.drug_name or "").lower()
            assoc_key = (approval.association or "").lower()
            entry_key = (drug_key, assoc_key)
            if drug_key and entry_key not in seen_entries:
                seen_entries.add(entry_key)

                # For combinations, use the full combination name
                if is_combination:
                    drug_name = approval.drug_name
                elif approval.generic_name and approval.brand_name:
                    drug_name = f"{approval.generic_name} ({approval.brand_name})"
                else:
                    drug_name = approval.brand_name or approval.generic_name or approval.drug_name

                cancer_specificity = self._get_fda_cancer_specificity(approval)

                # FDA labels from OpenFDA don't contain sensitivity/resistance info
                # That signal comes from CGI/CIViC/VICC therapeutic evidence separately
                # The association field may be None for FDA label-derived approvals
                response_type = None
                if approval.association:
                    assoc_upper = approval.association.upper()
                    if "RESIST" in assoc_upper:
                        response_type = "Resistance"
                    elif "RESPONS" in assoc_upper or "SENSITIV" in assoc_upper:
                        response_type = "Sensitivity"

                # Generate DailyMed URL for FDA drug label
                # Use brand name if available, otherwise generic name
                dailymed_search = approval.brand_name or approval.generic_name or approval.drug_name
                dailymed_url = f"https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query={dailymed_search.replace(' ', '+')}" if dailymed_search else None

                evidence_list.append(TherapeuticData(
                    drug_name=drug_name,
                    evidence_level="FDA-approved",
                    approval_status="Approved in indication" if approval.variant_in_indications else "Approved",
                    clinical_context=self._extract_line_of_therapy(approval),
                    response_type=response_type,
                    mechanism=None,
                    tumor_types_tested=[self.context.tumor_type] if self.context.tumor_type else [],
                    source="FDA",
                    source_url=dailymed_url,
                    confidence="high",
                    locus_match=approval.locus_match,
                    cancer_specificity=cancer_specificity,
                ))

        # From CIViC assertions - only Tier I with fda_companion_test
        for assertion in self.civic_assertions:
            if assertion.amp_tier == "Tier I" and assertion.fda_companion_test and assertion.therapies:
                for therapy in assertion.therapies:
                    drug_key = therapy.lower()
                    # Determine response type for deduplication key
                    civic_response = "sensitivity" if assertion.significance and "SENSITIV" in assertion.significance.upper() else "resistance" if assertion.significance and "RESIST" in assertion.significance.upper() else ""
                    entry_key = (drug_key, civic_response)
                    if entry_key not in seen_entries:
                        seen_entries.add(entry_key)

                        cancer_specificity = self._get_cancer_specificity_from_disease(assertion.disease)

                        evidence_list.append(TherapeuticData(
                            drug_name=therapy,
                            evidence_level="FDA-approved",
                            approval_status="FDA Approved (CIViC)",
                            clinical_context=assertion.disease,
                            response_type="Sensitivity" if civic_response == "sensitivity" else "Resistance" if civic_response == "resistance" else None,
                            mechanism=None,
                            tumor_types_tested=[assertion.disease] if assertion.disease else [],
                            source="CIViC",
                            source_url=assertion.civic_url,
                            confidence="high",
                            locus_match=assertion.locus_match,
                            cancer_specificity=cancer_specificity,
                        ))

        # From CGI biomarkers - fda_approved=True
        for biomarker in self.cgi_biomarkers:
            if biomarker.fda_approved and biomarker.drug and isinstance(biomarker.drug, str):
                drug_key = biomarker.drug.lower()
                assoc_key = (biomarker.association or "").lower()
                entry_key = (drug_key, assoc_key)
                if entry_key not in seen_entries:
                    seen_entries.add(entry_key)

                    # Determine response type from association field
                    # CGI values: "Responsive", "Resistant", "No Responsive", "Increased Toxicity"
                    # Only mark as Sensitivity/Resistance if explicitly stated
                    response_type = None
                    if biomarker.association:
                        assoc_upper = biomarker.association.upper()
                        if "RESIST" in assoc_upper:
                            response_type = "Resistance"
                        elif assoc_upper == "RESPONSIVE":
                            # Only explicit "Responsive" (not "No Responsive" or toxicity)
                            response_type = "Sensitivity"

                    cancer_specificity = self._get_cancer_specificity_from_disease(biomarker.tumor_type)

                    # Generate DailyMed URL for CGI-sourced FDA approvals
                    dailymed_url = f"https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query={biomarker.drug.replace(' ', '+')}"

                    evidence_list.append(TherapeuticData(
                        drug_name=biomarker.drug,
                        evidence_level="FDA-approved",
                        approval_status="FDA Approved (CGI)",
                        clinical_context=biomarker.tumor_type,
                        response_type=response_type,
                        mechanism=None,
                        tumor_types_tested=[biomarker.tumor_type] if biomarker.tumor_type else [],
                        source="CGI",
                        source_url=dailymed_url,
                        confidence="high",
                        locus_match=biomarker.locus_match,
                        cancer_specificity=cancer_specificity,
                    ))

        return evidence_list

    # Helper methods for get_therapeutic_evidence
    def _extract_line_of_therapy(self, approval) -> str | None:
        """Extract line of therapy from FDA approval."""
        if self.context.tumor_type and approval.indication:
            parsed = approval.parse_indication_for_tumor(self.context.tumor_type)
            if parsed.get('tumor_match'):
                return parsed.get('line_of_therapy')
        return None

    def _get_response_type_from_evidence(self, drug_name: str) -> str | None:
        """Check CIViC, VICC, and CGI evidence for response type for a drug.

        This is used to determine if a drug should show Sensitivity or Resistance
        for FDA approvals. Even if FDA approved at gene level, variant-specific
        resistance evidence should override (e.g., KIT D816V is imatinib-resistant).

        Sources checked (with weights):
        - VICC (OncoKB weighted 2x, others 1x)
        - CIViC evidence items
        - CGI biomarkers

        Args:
            drug_name: Drug name (lowercase) to look up

        Returns:
            "Resistance" if evidence shows resistance
            "Sensitivity" if evidence shows sensitivity
            "Conflicting" if both resistance and sensitivity with similar scores
            None if no evidence found
        """
        drug_lower = drug_name.lower()
        resistance_score = 0
        sensitivity_score = 0

        def drug_matches(drug_list: list[str]) -> bool:
            """Check if drug_name matches any drug in the list."""
            for d in drug_list:
                d_lower = d.lower()
                if drug_lower in d_lower or d_lower in drug_lower:
                    return True
            return False

        # Check VICC evidence (includes OncoKB, JAX-CKB, etc.)
        for vicc in self.vicc_evidence:
            if not vicc.drugs:
                continue
            if not drug_matches(vicc.drugs):
                continue

            # Weight OncoKB higher (more curated)
            weight = 2 if vicc.source == "oncokb" else 1

            if vicc.is_resistance:
                resistance_score += weight
            elif vicc.is_sensitivity:
                sensitivity_score += weight

        # Check CIViC evidence items
        for civic in self.civic_evidence:
            if not civic.drugs:
                continue
            if not drug_matches(civic.drugs):
                continue

            # CIViC evidence items weighted by evidence level
            # Level A/B = 2, Level C/D/E = 1
            weight = 2 if civic.evidence_level in ("A", "B") else 1

            # Use computed properties that check both clinical_significance AND evidence_direction
            if civic.is_resistance:
                resistance_score += weight
            elif civic.is_sensitivity:
                sensitivity_score += weight

        # Check CGI biomarkers
        for biomarker in self.cgi_biomarkers:
            if not biomarker.drug:
                continue
            if not drug_matches([biomarker.drug]):
                continue

            # CGI FDA-approved entries weighted higher
            weight = 2 if biomarker.fda_approved else 1

            if biomarker.association:
                assoc_upper = biomarker.association.upper()
                if "RESIST" in assoc_upper:
                    resistance_score += weight
                elif "RESPONS" in assoc_upper or "SENSITIV" in assoc_upper:
                    sensitivity_score += weight

        # Determine response type based on scores
        # Resistance takes precedence if scores are close (within 2x)
        # because resistance is clinically more important to flag
        if resistance_score > 0 and sensitivity_score > 0:
            # Both exist - check if one clearly dominates
            if resistance_score >= sensitivity_score:
                return "Resistance"
            elif sensitivity_score > resistance_score * 2:
                # Sensitivity only if it strongly outweighs resistance
                return "Sensitivity"
            else:
                return "Resistance"  # Default to resistance when conflicting
        elif resistance_score > 0:
            return "Resistance"
        elif sensitivity_score > 0:
            return "Sensitivity"

        return None

    def _drug_has_sensitivity_for_gene(self, drug_name: str) -> bool:
        """Check if a drug has sensitivity evidence for ANY variant in this gene.

        Used to filter resistance evidence - only show resistance for drugs that
        are actually relevant to this gene/tumor. For example, lapatinib resistance
        for EGFR T790M in NSCLC should be filtered out because lapatinib has no
        sensitivity evidence for any EGFR variant in NSCLC.

        Args:
            drug_name: Drug name (lowercase) to check

        Returns:
            True if drug has sensitivity evidence for any variant in this gene
        """
        drug_lower = drug_name.lower()

        # Check FDA approvals - these are the strongest evidence
        for approval in self.fda_approvals:
            fda_drug = (approval.generic_name or approval.brand_name or approval.drug_name or "").lower()
            if drug_lower in fda_drug or fda_drug in drug_lower:
                return True

        # Check VICC sensitivity evidence (includes OncoKB, CIViC, CGI, etc.)
        for vicc in self.vicc_evidence:
            if vicc.is_sensitivity and vicc.drugs:
                for d in vicc.drugs:
                    if drug_lower in d.lower() or d.lower() in drug_lower:
                        return True

        # Check CIViC assertions for sensitivity
        for assertion in self.civic_assertions:
            if assertion.is_sensitivity and assertion.therapies:
                for therapy in assertion.therapies:
                    if drug_lower in therapy.lower() or therapy.lower() in drug_lower:
                        return True

        # Check CIViC evidence items for sensitivity
        for evidence in self.civic_evidence:
            if evidence.clinical_significance:
                sig = evidence.clinical_significance.upper()
                if "SENS" in sig or "RESPON" in sig:
                    if evidence.drugs:
                        for d in evidence.drugs:
                            if drug_lower in d.lower() or d.lower() in drug_lower:
                                return True

        # Check CGI biomarkers for sensitivity
        # Only explicit "Responsive" (not "No Responsive", "Increased Toxicity", etc.)
        for biomarker in self.cgi_biomarkers:
            if biomarker.association and biomarker.association.upper() == "RESPONSIVE":
                if biomarker.drug:
                    if drug_lower in biomarker.drug.lower() or biomarker.drug.lower() in drug_lower:
                        return True

        return False

    def _get_fda_cancer_specificity(self, approval) -> str:
        """Determine cancer_specificity for an FDA approval.

        Uses the cancer_specificity property from FDAApproval, which reads from
        cancer_type_match set during evidence aggregation. This ensures consistency
        with how ClinicalTrialEvidence handles tumor matching.

        Returns:
            - "cancer_specific": if the FDA indication matches the queried tumor type
            - "pan_cancer": if tumor-agnostic (e.g., solid tumors)
            - The actual cancer type (e.g., "breast cancer"): if approved for a different tumor
            - "unknown": if we can't determine the cancer type (should NOT be treated as pan_cancer)
        """
        # Use the property from FDAApproval (reads from cancer_type_match)
        if approval.cancer_specificity:
            return approval.cancer_specificity

        # Fallback for approvals without cancer_type_match set (legacy data)
        if self.context.tumor_type and approval.indication:
            parsed = approval.parse_indication_for_tumor(self.context.tumor_type)
            if parsed.get('tumor_match'):
                return "cancer_specific"

        # No match - extract what cancer the drug IS approved for
        indication_cancer = approval.extract_indication_cancer_type()
        if indication_cancer:
            if "pan-cancer" in indication_cancer.lower() or "solid tumor" in indication_cancer.lower():
                return "pan_cancer"
            return indication_cancer

        # Unknown - do NOT assume pan_cancer, return as unknown so it gets filtered out
        return "unknown"

    def _get_cancer_specificity_from_disease(self, disease: str | None) -> str:
        """Determine cancer_specificity from a disease/tumor type string.

        Compares the evidence disease against the queried tumor type.

        Args:
            disease: The disease/tumor type from the evidence (e.g., "Ovarian Cancer")

        Returns:
            - "cancer_specific": if the disease matches the queried tumor type
            - The disease name: if no match but we have a specific disease
            - "pan_cancer": if unknown or tumor-agnostic
        """
        if not disease:
            return None

        queried_tumor = self.context.tumor_type
        if not queried_tumor:
            return None

        # Check for tumor-agnostic terms using centralized function
        if is_pan_cancer_term(disease):
            return "pan_cancer"

        # Check if disease matches queried tumor using centralized function
        if tumor_types_match(disease, queried_tumor):
            return "cancer_specific"

        # No match - return the specific disease
        return disease

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

    def _civic_level_to_evidence_level(self, level: str | None) -> str:
        """Map CIViC evidence level (A-E) to standard format.

        CIViC evidence levels:
        - A: Validated association (FDA-approved or professional guidelines)
        - B: Clinical evidence (well-powered studies with consensus)
        - C: Case study (limited clinical evidence)
        - D: Preclinical (in vivo/in vitro models)
        - E: Inferential (indirect evidence)
        """
        if not level:
            return "Unknown"
        level_upper = level.upper()
        if level_upper == "A":
            return "FDA-approved"
        elif level_upper == "B":
            return "Phase 3"
        elif level_upper == "C":
            return "Case report"
        elif level_upper == "D":
            return "Preclinical"
        elif level_upper == "E":
            return "Preclinical"  # Inferential treated as preclinical
        return "Unknown"

    def _get_approval_from_civic_level(self, level: str | None) -> str:
        """Get approval status from CIViC evidence level."""
        if not level:
            return "Investigational"
        level_upper = level.upper()
        if level_upper == "A":
            return "FDA Approved"
        elif level_upper == "B":
            return "Clinical trials"
        return "Investigational"

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

    def _get_variant_level_response_drugs(self) -> tuple[set[str], set[str]]:
        """Get drugs with variant-level sensitivity or resistance evidence.

        Used to filter out conflicting codon/gene-level evidence.
        If a drug has variant-level resistance, we should not show codon-level sensitivity.
        If a drug has variant-level sensitivity, we should not show codon-level resistance.

        Returns:
            Tuple of (variant_level_sensitive_drugs, variant_level_resistant_drugs)
        """
        variant_sensitive: set[str] = set()
        variant_resistant: set[str] = set()

        # Check VICC evidence
        for vicc in self.vicc_evidence:
            if vicc.locus_match == "variant" and vicc.drugs:
                for drug in vicc.drugs:
                    drug_lower = drug.lower()
                    if vicc.is_sensitivity:
                        variant_sensitive.add(drug_lower)
                    if vicc.is_resistance:
                        variant_resistant.add(drug_lower)

        # Check CIViC assertions
        for assertion in self.civic_assertions:
            if assertion.locus_match == "variant" and assertion.therapies:
                for therapy in assertion.therapies:
                    drug_lower = therapy.lower()
                    if assertion.is_sensitivity:
                        variant_sensitive.add(drug_lower)
                    if assertion.is_resistance:
                        variant_resistant.add(drug_lower)

        # Check CIViC evidence items
        for evidence in self.civic_evidence:
            if evidence.locus_match == "variant" and evidence.drugs and evidence.clinical_significance:
                sig_upper = evidence.clinical_significance.upper()
                for drug in evidence.drugs:
                    drug_lower = drug.lower()
                    if "SENS" in sig_upper or "RESPON" in sig_upper:
                        variant_sensitive.add(drug_lower)
                    if "RESIST" in sig_upper:
                        variant_resistant.add(drug_lower)

        # Check CGI biomarkers
        for biomarker in self.cgi_biomarkers:
            if biomarker.locus_match == "variant" and biomarker.drug and biomarker.association:
                drug_lower = biomarker.drug.lower()
                assoc_upper = biomarker.association.upper()
                if "SENS" in assoc_upper or "RESPON" in assoc_upper:
                    variant_sensitive.add(drug_lower)
                if "RESIST" in assoc_upper:
                    variant_resistant.add(drug_lower)

        # Check CGI preclinical and early-phase biomarkers
        for biomarker in self.preclinical_biomarkers + self.early_phase_biomarkers:
            if biomarker.locus_match == "variant" and biomarker.drug and biomarker.association:
                drug_lower = biomarker.drug.lower()
                assoc_upper = biomarker.association.upper()
                if "SENS" in assoc_upper or "RESPON" in assoc_upper:
                    variant_sensitive.add(drug_lower)
                if "RESIST" in assoc_upper:
                    variant_resistant.add(drug_lower)

        return variant_sensitive, variant_resistant

    def get_resistance_summary(self) -> str:
        """Get a synthesized summary of resistance evidence.

        Returns a narrative summary grouping drugs by class/generation where possible.
        Includes locus match level for each drug.

        IMPORTANT: If a drug has variant-level SENSITIVITY evidence, we exclude codon/gene-level
        resistance evidence for that drug. This prevents misleading summaries when the queried
        variant has specific sensitivity evidence.

        Example: "Resistance to 1st/2nd-gen EGFR TKIs: erlotinib (variant-level), gefitinib (variant-level)"
        """
        # Get drugs with variant-level sensitivity - we'll exclude codon/gene resistance for these
        variant_sensitive_drugs, _ = self._get_variant_level_response_drugs()

        # Collect all resistance drugs with their sources and locus match level
        # dict[str, tuple[set[str], str]] = drug -> (sources, best_locus_match)
        resistance_drugs: dict[str, tuple[set[str], str]] = {}

        def get_best_locus(existing: str | None, new: str | None) -> str:
            """Get the most specific locus match (variant > codon > gene)."""
            priority = {"variant": 3, "codon": 2, "gene": 1}
            existing_p = priority.get(existing or "", 0)
            new_p = priority.get(new or "", 0)
            if new_p > existing_p:
                return new or "gene"
            return existing or new or "gene"

        def should_include_resistance(drug_lower: str, locus: str) -> bool:
            """Check if resistance evidence should be included.

            Exclude codon/gene-level resistance if drug has variant-level sensitivity.
            """
            if locus == "variant":
                return True  # Always include variant-level evidence
            # For codon/gene level, exclude if drug has variant-level sensitivity
            return drug_lower not in variant_sensitive_drugs

        def add_drug(drug: str, source: str, locus: str | None) -> None:
            """Add or update a drug entry with source and locus."""
            drug_lower = drug.lower()
            locus_str = locus or "gene"
            if not should_include_resistance(drug_lower, locus_str):
                return  # Skip this evidence
            if drug_lower in resistance_drugs:
                sources, existing_locus = resistance_drugs[drug_lower]
                sources.add(source)
                resistance_drugs[drug_lower] = (sources, get_best_locus(existing_locus, locus_str))
            else:
                resistance_drugs[drug_lower] = ({source}, locus_str)

        # VICC resistance evidence - with locus match
        for vicc in self.vicc_evidence:
            if vicc.is_resistance and vicc.drugs:
                locus = vicc.locus_match or "gene"
                for drug in vicc.drugs:
                    add_drug(drug, "VICC", locus)

        # CIViC resistance assertions - with locus match
        for assertion in self.civic_assertions:
            if assertion.is_resistance and assertion.therapies:
                locus = assertion.locus_match or "gene"
                for therapy in assertion.therapies:
                    add_drug(therapy, "CIViC", locus)

        # CIViC resistance evidence items - with locus match
        for evidence in self.civic_evidence:
            if evidence.clinical_significance and "RESIST" in evidence.clinical_significance.upper():
                if evidence.drugs:
                    locus = evidence.locus_match or "gene"
                    for drug in evidence.drugs:
                        add_drug(drug, "CIViC", locus)

        # CGI resistance biomarkers (FDA-approved) - with locus match
        for biomarker in self.cgi_biomarkers:
            if biomarker.association and "RESIST" in biomarker.association.upper():
                if biomarker.drug:
                    locus = biomarker.locus_match or "gene"
                    add_drug(biomarker.drug, "CGI", locus)

        # CGI preclinical resistance biomarkers - with locus match
        for biomarker in self.preclinical_biomarkers:
            if biomarker.association and "RESIST" in biomarker.association.upper():
                if biomarker.drug:
                    locus = biomarker.locus_match or "gene"
                    add_drug(biomarker.drug, "CGI (preclinical)", locus)

        # CGI early-phase resistance biomarkers - with locus match
        for biomarker in self.early_phase_biomarkers:
            if biomarker.association and "RESIST" in biomarker.association.upper():
                if biomarker.drug:
                    locus = biomarker.locus_match or "gene"
                    add_drug(biomarker.drug, "CGI (early phase)", locus)

        # Literature resistance signals from LLM extraction (no locus match available)
        if self.literature_knowledge and self.literature_knowledge.resistant_to:
            for entry in self.literature_knowledge.resistant_to:
                drug = entry.drug
                if drug:
                    add_drug(drug, "Literature", None)

        # PubMed articles with resistance evidence (no locus match available)
        for article in self.pubmed_articles:
            if article.is_resistance_evidence() and article.drugs_mentioned:
                for drug in article.drugs_mentioned:
                    add_drug(drug, f"PubMed:{article.pmid}", None)

        # Filter: only include resistance for drugs that have sensitivity evidence
        # for some variant in this gene. This removes drugs like lapatinib that were
        # tested as negative controls but have no clinical relevance for this gene/tumor.
        resistance_drugs = {
            drug: (sources, locus) for drug, (sources, locus) in resistance_drugs.items()
            if self._drug_has_sensitivity_for_gene(drug)
        }

        if not resistance_drugs:
            return ""

        # Check for conflicting evidence
        conflicting_drugs: dict[str, str] = {}  # drug -> conflict description
        if self.evidence_gaps:
            from oncomind.models.extracted.evidence_gaps import GapCategory
            discordant_gaps = self.evidence_gaps.get_gaps_by_category(GapCategory.DISCORDANT)
            for gap in discordant_gaps:
                if "Conflicting drug response for" in gap.description:
                    # Extract drug name from "Conflicting drug response for X: ..."
                    parts = gap.description.split(":")
                    if parts:
                        drug_part = parts[0].replace("Conflicting drug response for", "").strip()
                        conflicting_drugs[drug_part.lower()] = gap.description

        # Helper to format drug with locus level
        def format_drug_with_locus(drug: str) -> str:
            """Format a drug name with its locus match level."""
            _, locus = resistance_drugs[drug]
            return f"{drug} ({locus}-level)"

        # Group drugs by class for synthesis
        summaries = []
        gene = self.identifiers.gene.upper()

        # EGFR TKI grouping
        egfr_1st_2nd_gen = {"erlotinib", "gefitinib", "afatinib", "lapatinib", "dacomitinib", "neratinib"}
        egfr_3rd_gen = {"osimertinib"}

        # Build synthesized summary
        if gene == "EGFR":
            resistant_1st_2nd = [d for d in resistance_drugs if d in egfr_1st_2nd_gen]
            resistant_3rd = [d for d in resistance_drugs if d in egfr_3rd_gen]
            other_drugs = [d for d in resistance_drugs if d not in egfr_1st_2nd_gen and d not in egfr_3rd_gen]

            if resistant_1st_2nd:
                drug_list = ", ".join(format_drug_with_locus(d) for d in sorted(resistant_1st_2nd))
                summaries.append(f"Resistance to 1st/2nd-gen EGFR TKIs: {drug_list}")

            for drug in resistant_3rd:
                _, locus = resistance_drugs[drug]
                if drug in conflicting_drugs:
                    summaries.append(f"Variable response to {drug} ({locus}-level, conflicting evidence)")
                else:
                    summaries.append(f"Resistance to {drug} ({locus}-level)")

            for drug in other_drugs[:3]:
                _, locus = resistance_drugs[drug]
                if drug in conflicting_drugs:
                    summaries.append(f"Variable response to {drug} ({locus}-level, conflicting evidence)")
                else:
                    summaries.append(f"{drug} resistance ({locus}-level)")
        else:
            # Generic grouping for non-EGFR genes
            for drug in list(resistance_drugs.keys())[:8]:
                _, locus = resistance_drugs[drug]
                if drug in conflicting_drugs:
                    summaries.append(f"Variable response to {drug} ({locus}-level, conflicting evidence)")
                else:
                    summaries.append(f"{drug} resistance ({locus}-level)")

        return "; ".join(summaries[:5]) if summaries else ""

    def get_sensitivity_summary(self) -> str:
        """Get a synthesized summary of sensitivity evidence.

        Returns a narrative summary grouping drugs by evidence tier (FDA > Clinical > Preclinical).
        Includes cancer specificity and locus match level for each drug.

        IMPORTANT: If a drug has variant-level RESISTANCE evidence, we exclude codon/gene-level
        sensitivity evidence for that drug. This prevents misleading summaries like
        "sunitinib sensitive (codon-level)" when sunitinib has variant-level resistance for the queried variant.

        Example: "FDA-approved for melanoma: osimertinib (variant-level); Clinical evidence: ponatinib (codon-level)"
        """
        # Get drugs with variant-level resistance - we'll exclude codon/gene sensitivity for these
        _, variant_resistant_drugs = self._get_variant_level_response_drugs()

        # Collect drugs by evidence tier with locus match level
        # Using dict[str, str] where key=drug, value=locus_match (variant/codon/gene)
        fda_matching_drugs: dict[str, str] = {}  # FDA approved for THIS tumor
        fda_other_drugs: dict[str, tuple[str, str]] = {}  # drug -> (cancer_type, locus_match)
        clinical_drugs: dict[str, str] = {}  # VICC, CIViC, CGI -> drug: locus_match
        literature_drugs: set[str] = set()

        tumor_type = self.context.tumor_type or "cancer"

        def get_best_locus(existing: str | None, new: str | None) -> str:
            """Get the most specific locus match (variant > codon > gene)."""
            priority = {"variant": 3, "codon": 2, "gene": 1}
            existing_p = priority.get(existing or "", 0)
            new_p = priority.get(new or "", 0)
            if new_p > existing_p:
                return new or "gene"
            return existing or new or "gene"

        def should_include_sensitivity(drug_lower: str, locus: str) -> bool:
            """Check if sensitivity evidence should be included.

            Exclude codon/gene-level sensitivity if drug has variant-level resistance.
            """
            if locus == "variant":
                return True  # Always include variant-level evidence
            # For codon/gene level, exclude if drug has variant-level resistance
            return drug_lower not in variant_resistant_drugs

        # FDA approvals - with cancer specificity and biomarker match level
        # Use biomarker_match_level (not locus_match) for FDA because it indicates
        # whether the approval covers the variant (e.g., "AKT1 alteration" covers E17K)
        for approval in self.fda_approvals:
            drug = approval.generic_name or approval.brand_name or approval.drug_name
            if drug:
                drug_lower = drug.lower()
                # Use biomarker_match_level for FDA approvals
                locus = approval.biomarker_match_level or "gene"
                cancer_spec = self._get_fda_cancer_specificity(approval)
                if cancer_spec == "cancer_specific" or cancer_spec == "pan_cancer":
                    if drug_lower in fda_matching_drugs:
                        fda_matching_drugs[drug_lower] = get_best_locus(fda_matching_drugs[drug_lower], locus)
                    else:
                        fda_matching_drugs[drug_lower] = locus
                else:
                    # FDA approved for a different cancer
                    if drug_lower in fda_other_drugs:
                        _, existing_locus = fda_other_drugs[drug_lower]
                        fda_other_drugs[drug_lower] = (cancer_spec, get_best_locus(existing_locus, locus))
                    else:
                        fda_other_drugs[drug_lower] = (cancer_spec, locus)

        # All FDA drugs for deduplication
        all_fda_drugs = set(fda_matching_drugs.keys()) | set(fda_other_drugs.keys())

        # VICC sensitivity evidence - with locus match
        # Level D = preclinical, Level A/B/C = clinical
        # Keep combination drugs together (e.g., "capivasertib + fulvestrant")
        vicc_preclinical_drugs: dict[str, str] = {}
        for vicc in self.vicc_evidence:
            if vicc.is_sensitivity and vicc.drugs:
                locus = vicc.locus_match or "gene"
                is_preclinical = vicc.evidence_level and vicc.evidence_level.upper() == "D"
                # Join multiple drugs as combination therapy
                drug_name = " + ".join(sorted(vicc.drugs)) if len(vicc.drugs) > 1 else vicc.drugs[0]
                drug_lower = drug_name.lower()
                if drug_lower not in all_fda_drugs and should_include_sensitivity(drug_lower, locus):
                    if is_preclinical:
                        if drug_lower in vicc_preclinical_drugs:
                            vicc_preclinical_drugs[drug_lower] = get_best_locus(vicc_preclinical_drugs[drug_lower], locus)
                        else:
                            vicc_preclinical_drugs[drug_lower] = locus
                    else:
                        if drug_lower in clinical_drugs:
                            clinical_drugs[drug_lower] = get_best_locus(clinical_drugs[drug_lower], locus)
                        else:
                            clinical_drugs[drug_lower] = locus

        # CIViC sensitivity assertions - with locus match (assertions are curated, treat as clinical)
        # Keep combination therapies together
        for assertion in self.civic_assertions:
            if assertion.is_sensitivity and assertion.therapies:
                locus = assertion.locus_match or "gene"
                # Join multiple therapies as combination
                drug_name = " + ".join(sorted(assertion.therapies)) if len(assertion.therapies) > 1 else assertion.therapies[0]
                drug_lower = drug_name.lower()
                if drug_lower not in all_fda_drugs and should_include_sensitivity(drug_lower, locus):
                    if drug_lower in clinical_drugs:
                        clinical_drugs[drug_lower] = get_best_locus(clinical_drugs[drug_lower], locus)
                    else:
                        clinical_drugs[drug_lower] = locus

        # CIViC sensitivity evidence items - with locus match
        # Level D/E = preclinical, Level A/B/C = clinical
        # Keep combination drugs together (e.g., "capivasertib + fulvestrant")
        civic_preclinical_drugs: dict[str, str] = {}
        for evidence in self.civic_evidence:
            if evidence.is_sensitivity and evidence.drugs:
                locus = evidence.locus_match or "gene"
                level = (evidence.evidence_level or "").upper()
                is_preclinical = level in ("D", "E")
                # Join multiple drugs as combination therapy
                drug_name = " + ".join(sorted(evidence.drugs)) if len(evidence.drugs) > 1 else evidence.drugs[0]
                drug_lower = drug_name.lower()
                if drug_lower not in all_fda_drugs and should_include_sensitivity(drug_lower, locus):
                    if is_preclinical:
                        if drug_lower in civic_preclinical_drugs:
                            civic_preclinical_drugs[drug_lower] = get_best_locus(civic_preclinical_drugs[drug_lower], locus)
                        else:
                            civic_preclinical_drugs[drug_lower] = locus
                    else:
                        if drug_lower in clinical_drugs:
                            clinical_drugs[drug_lower] = get_best_locus(clinical_drugs[drug_lower], locus)
                        else:
                            clinical_drugs[drug_lower] = locus

        # CGI sensitivity biomarkers (FDA-approved) - with locus match
        # Only explicit "Responsive" (not "No Responsive", "Increased Toxicity", etc.)
        for biomarker in self.cgi_biomarkers:
            if biomarker.association and biomarker.association.upper() == "RESPONSIVE":
                if biomarker.drug:
                    drug_lower = biomarker.drug.lower()
                    locus = biomarker.locus_match or "gene"
                    if drug_lower not in all_fda_drugs and should_include_sensitivity(drug_lower, locus):
                        if drug_lower in clinical_drugs:
                            clinical_drugs[drug_lower] = get_best_locus(clinical_drugs[drug_lower], locus)
                        else:
                            clinical_drugs[drug_lower] = locus

        # CGI preclinical sensitivity biomarkers - with locus match
        # Only explicit "Responsive" (not "No Responsive", "Increased Toxicity", etc.)
        cgi_preclinical_drugs: dict[str, str] = {}
        for biomarker in self.preclinical_biomarkers:
            if biomarker.association and biomarker.association.upper() == "RESPONSIVE":
                if biomarker.drug:
                    drug_lower = biomarker.drug.lower()
                    locus = biomarker.locus_match or "gene"
                    if drug_lower not in all_fda_drugs and drug_lower not in clinical_drugs and should_include_sensitivity(drug_lower, locus):
                        if drug_lower in cgi_preclinical_drugs:
                            cgi_preclinical_drugs[drug_lower] = get_best_locus(cgi_preclinical_drugs[drug_lower], locus)
                        else:
                            cgi_preclinical_drugs[drug_lower] = locus

        # CGI early-phase sensitivity biomarkers (clinical trials, case reports) - with locus match
        # Only explicit "Responsive" (not "No Responsive", "Increased Toxicity", etc.)
        early_phase_drugs: dict[str, str] = {}
        for biomarker in self.early_phase_biomarkers:
            if biomarker.association and biomarker.association.upper() == "RESPONSIVE":
                if biomarker.drug:
                    drug_lower = biomarker.drug.lower()
                    locus = biomarker.locus_match or "gene"
                    if drug_lower not in all_fda_drugs and drug_lower not in clinical_drugs and should_include_sensitivity(drug_lower, locus):
                        if drug_lower in early_phase_drugs:
                            early_phase_drugs[drug_lower] = get_best_locus(early_phase_drugs[drug_lower], locus)
                        else:
                            early_phase_drugs[drug_lower] = locus

        # Literature sensitivity signals from LLM extraction (no locus match available)
        if self.literature_knowledge and self.literature_knowledge.sensitive_to:
            for entry in self.literature_knowledge.sensitive_to:
                drug = entry.drug
                if drug:
                    drug_lower = drug.lower()
                    if drug_lower not in all_fda_drugs and drug_lower not in clinical_drugs:
                        literature_drugs.add(drug_lower)

        # PubMed articles with sensitivity evidence (no locus match available)
        for article in self.pubmed_articles:
            if article.is_sensitivity_evidence() and article.drugs_mentioned:
                for drug in article.drugs_mentioned:
                    drug_lower = drug.lower()
                    if drug_lower not in all_fda_drugs and drug_lower not in clinical_drugs:
                        literature_drugs.add(drug_lower)

        # DepMap/PRISM drug sensitivities (preclinical) - gene-level by nature
        depmap_preclinical_drugs: dict[str, str] = {}
        if self.depmap_evidence and self.depmap_evidence.drug_sensitivities:
            for ds in self.depmap_evidence.get_top_sensitive_drugs(5):
                drug_lower = ds.drug_name.lower()
                if drug_lower not in all_fda_drugs and drug_lower not in clinical_drugs:
                    depmap_preclinical_drugs[drug_lower] = "gene"  # DepMap is gene-level

        # Combine all preclinical drugs (VICC level D, CIViC level D/E, CGI preclinical, DepMap)
        all_preclinical_drugs = {**vicc_preclinical_drugs, **civic_preclinical_drugs, **cgi_preclinical_drugs, **depmap_preclinical_drugs}

        # Build synthesized summary with locus match levels
        summaries = []

        def format_drugs_with_locus(drugs: dict[str, str], limit: int = 5) -> str:
            """Format drugs with their locus match level."""
            parts = []
            for drug, locus in sorted(drugs.items())[:limit]:
                parts.append(f"{drug} ({locus}-level)")
            return ", ".join(parts)

        # FDA approved for THIS tumor type
        if fda_matching_drugs:
            drug_list = format_drugs_with_locus(fda_matching_drugs, 5)
            summaries.append(f"FDA-approved for {tumor_type}: {drug_list}")

        # FDA approved for OTHER tumor types (IMPORTANT: LLM must know this is NOT for the queried tumor)
        if fda_other_drugs:
            other_parts = []
            for drug, (cancer, locus) in sorted(fda_other_drugs.items())[:3]:
                other_parts.append(f"{drug} ({locus}-level, {cancer})")
            summaries.append(f"FDA-approved for OTHER cancers (NOT {tumor_type}): {', '.join(other_parts)}")

        if clinical_drugs:
            drug_list = format_drugs_with_locus(clinical_drugs, 5)
            summaries.append(f"Clinical evidence: {drug_list}")

        if early_phase_drugs:
            drug_list = format_drugs_with_locus(early_phase_drugs, 3)
            summaries.append(f"Early phase trials: {drug_list}")

        if all_preclinical_drugs:
            drug_list = format_drugs_with_locus(all_preclinical_drugs, 3)
            summaries.append(f"Preclinical: {drug_list}")

        if literature_drugs and not fda_matching_drugs and not fda_other_drugs and not clinical_drugs and not early_phase_drugs and not all_preclinical_drugs:
            # Only show literature if no higher-tier evidence (no locus match for literature)
            drug_list = ", ".join(sorted(literature_drugs)[:3])
            summaries.append(f"Literature signals: {drug_list}")

        return "; ".join(summaries) if summaries else ""

    def get_evidence_summary_for_llm(self) -> str:
        """Generate a compact evidence summary for LLM prompt consumption."""
        lines = []
        tumor_type = self.context.tumor_type
        gene = self.identifiers.gene

        # FDA Approvals - use get_fda_approved_therapies() which computes response_type from VICC
        # Filter out:
        # 1. Biomarker selection drugs (e.g., DATROWAY for EGFR - targets TROP2, not EGFR)
        # 2. Drugs approved for OTHER tumor types (e.g., RYDAPT approved for AML, not GIST)
        #    - cancer_specific = matches queried tumor type
        #    - pan_cancer = applies to all tumors
        #    - Other values (e.g., "AML") = approved for different tumor type, EXCLUDE
        fda_therapies = self.get_fda_approved_therapies()
        fda_from_fda = [
            t for t in fda_therapies
            if t.source == "FDA"
            and not is_biomarker_selection_drug(t.drug_name, gene)
            and t.cancer_specificity in ("cancer_specific", "pan_cancer")
        ]
        if fda_from_fda:
            fda_drugs = []
            for therapy in fda_from_fda[:4]:
                drug = therapy.drug_name

                # Build annotation parts
                # NOTE: Do NOT include locus_match level for matched drugs - it causes LLM to hedge
                # If the drug is in fda_from_fda, it's already matched (biomarker_matched=True)
                parts = []

                # Add response type (Sensitivity/Resistance) - critical for LLM
                if therapy.response_type:
                    parts.append(therapy.response_type.lower())

                # Add clinical context if available
                if therapy.clinical_context:
                    parts.append(therapy.clinical_context)

                # Add tumor match info
                if therapy.cancer_specificity == "cancer_specific":
                    parts.append(tumor_type)
                elif therapy.cancer_specificity and therapy.cancer_specificity not in ("pan_cancer", "cancer_specific"):
                    parts.append(f"approved for {therapy.cancer_specificity}")

                fda_drugs.append(f"{drug} ({', '.join(parts)})" if parts else drug)
            lines.append(f"FDA Approved: {', '.join(fda_drugs)}")

        # FDA Codon-Level Near-Misses - drugs approved for OTHER variants at same codon
        # e.g., KRAS G12C drugs (sotorasib, adagrasib) when patient has G12A
        # These are NOT matched (biomarker_matched=False) but are at codon-level
        codon_near_misses = []
        for approval in self.fda_approvals:
            # Only include if NOT matched but has codon-level match
            if not approval.biomarker_matched and approval.biomarker_match_level == "codon":
                # Use drug_name first - it may contain combination name (e.g., "capivasertib + fulvestrant")
                drug_name = approval.drug_name or approval.generic_name or approval.brand_name
                if drug_name and not is_biomarker_selection_drug(drug_name, gene):
                    # Extract what variant it's actually approved for
                    approved_variant = approval.extract_approved_variant()
                    if approved_variant:
                        codon_near_misses.append(f"{drug_name} (approved for {approved_variant})")
                    else:
                        codon_near_misses.append(drug_name)
        if codon_near_misses:
            lines.append(f"FDA Codon-Level (not for queried variant): {', '.join(codon_near_misses[:4])}")

        # CGI Biomarkers - compact (filter out biomarker selection drugs)
        if self.cgi_biomarkers:
            approved = [
                b for b in self.cgi_biomarkers
                if b.fda_approved and not is_biomarker_selection_drug(b.drug or "", gene)
            ]
            if approved:
                resistance = [b.drug for b in approved if b.association and 'RESIST' in b.association.upper()]
                # Only explicit "Responsive" (not "No Responsive", "Increased Toxicity", etc.)
                sensitivity = [b.drug for b in approved if b.association and b.association.upper() == 'RESPONSIVE']
                if resistance:
                    lines.append(f"CGI Resistance: {', '.join(resistance[:3])}")
                if sensitivity:
                    lines.append(f"CGI Sensitivity: {', '.join(sensitivity[:3])}")

        # CIViC Assertions - compact (filter out biomarker selection drugs)
        if self.civic_assertions:
            predictive = [a for a in self.civic_assertions if a.assertion_type == "PREDICTIVE"]
            if predictive:
                civic_drugs = []
                for a in predictive[:3]:
                    # Filter therapies that are biomarker selection drugs
                    filtered_therapies = [
                        t for t in (a.therapies or [])
                        if not is_biomarker_selection_drug(t, gene)
                    ]
                    therapies = ", ".join(filtered_therapies[:2]) if filtered_therapies else ""
                    sig = "sens" if a.is_sensitivity else "res" if a.is_resistance else ""
                    if therapies:
                        civic_drugs.append(f"{therapies} ({sig})")
                if civic_drugs:
                    lines.append(f"CIViC Assertions: {'; '.join(civic_drugs)}")

        # CIViC Evidence Items - compact (includes per-publication drug response data)
        # Filter out biomarker selection drugs
        if self.civic_evidence:
            civic_eid_drugs = []
            for e in self.civic_evidence[:5]:
                if e.drugs:
                    filtered_drugs = [d for d in e.drugs if not is_biomarker_selection_drug(d, gene)]
                    if not filtered_drugs:
                        continue
                    drug_str = ", ".join(filtered_drugs[:2])
                    # Use computed properties that check both clinical_significance AND evidence_direction
                    sig = ""
                    if e.is_resistance:
                        sig = "res"
                    elif e.is_sensitivity:
                        sig = "sens"
                    if sig:
                        civic_eid_drugs.append(f"{drug_str} ({sig})")
                    else:
                        civic_eid_drugs.append(drug_str)
            if civic_eid_drugs:
                lines.append(f"CIViC Evidence: {'; '.join(civic_eid_drugs)}")

        # VICC MetaKB evidence - compact (aggregated from CIViC, CGI, JAX, OncoKB, PMKB)
        # Filter out biomarker selection drugs
        if self.vicc_evidence:
            vicc_drugs = []
            for v in self.vicc_evidence[:5]:
                if v.drugs:
                    filtered_drugs = [d for d in v.drugs if not is_biomarker_selection_drug(d, gene)]
                    if not filtered_drugs:
                        continue
                    drug_str = ", ".join(filtered_drugs[:2])
                    sig = ""
                    if v.is_resistance:
                        sig = "res"
                    elif v.is_sensitivity:
                        sig = "sens"
                    source = f", {v.source}" if v.source else ""
                    if sig:
                        vicc_drugs.append(f"{drug_str} ({sig}{source})")
                    else:
                        vicc_drugs.append(f"{drug_str}{source}")
            if vicc_drugs:
                lines.append(f"VICC MetaKB: {'; '.join(vicc_drugs)}")

        # ClinVar - include all entries with significance, review status, and conditions
        if self.clinvar_entries:
            clinvar_summaries = []
            for entry in self.clinvar_entries[:10]:  # Limit to 10 entries
                parts = []
                if entry.clinical_significance:
                    parts.append(entry.clinical_significance)
                if entry.review_status:
                    parts.append(f"review: {entry.review_status}")
                if entry.conditions:
                    # Show first 2 conditions
                    conds = ", ".join(entry.conditions[:2])
                    if len(entry.conditions) > 2:
                        conds += f" (+{len(entry.conditions) - 2} more)"
                    parts.append(f"conditions: {conds}")
                if parts:
                    summary_text = " | ".join(parts)
                    # Add link if variation_id exists
                    url = entry.get_url()
                    if url:
                        clinvar_summaries.append(f"[{summary_text}]({url})")
                    else:
                        clinvar_summaries.append(summary_text)
            if clinvar_summaries:
                lines.append(f"ClinVar ({len(self.clinvar_entries)} entries):")
                for i, summary in enumerate(clinvar_summaries, 1):
                    lines.append(f"  {i}. {summary}")
        elif self.clinvar_significance:
            # Fallback to single significance if no entries but significance exists
            lines.append(f"ClinVar: {self.clinvar_significance}")

        # Clinical Trials - compact
        if self.clinical_trials:
            recruiting = [t for t in self.clinical_trials if t.status == "RECRUITING"]
            if recruiting:
                trial_drugs = []
                for trial in recruiting[:5]:
                    drugs = trial.get_drug_names()
                    if drugs:
                        phase_str = f" ({trial.phase})" if trial.phase else ""
                        trial_drugs.append(f"{drugs[0]}{phase_str}")
                if trial_drugs:
                    lines.append(f"Active Trials: {'; '.join(trial_drugs)}")
                else:
                    lines.append(f"Active Trials: {len(recruiting)} recruiting")
            else:
                lines.append(f"Clinical Trials: {len(self.clinical_trials)} found (not recruiting)")

        # PubMed Literature - compact with links
        if self.pubmed_articles:
            resistance_articles = [a for a in self.pubmed_articles if a.is_resistance_evidence()]
            if resistance_articles:
                pmid_links = [
                    f"[{a.pmid}](https://pubmed.ncbi.nlm.nih.gov/{a.pmid}/)"
                    for a in resistance_articles[:3]
                ]
                lines.append(f"Literature (resistance): PMIDs {', '.join(pmid_links)}")

        return "\n".join(lines) if lines else ""

    def get_biological_context_for_llm(self) -> str:
        """Get biological context formatted for LLM prompt.

        Combines cBioPortal and DepMap data with gene context for research-focused synthesis.

        Returns:
            Formatted string with gene role, pathway, cBioPortal, and DepMap data
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

        # cBioPortal data - only include if we have data
        # When no cBioPortal data exists (e.g., acquired resistance mutations like T790M
        # that don't appear in baseline sequencing), skip prevalence section entirely
        # rather than misleading the LLM with "No data available"
        if self.cbioportal_evidence and self.cbioportal_evidence.has_data():
            lines.append(self.cbioportal_evidence.to_prompt_context())
        else:
            # Check if this is a known acquired resistance mutation
            gene = self.identifiers.gene
            variant = self.identifiers.variant
            if gene and variant and is_acquired_resistance_mutation(gene, variant):
                lines.append(
                    f"PREVALENCE: {gene} {variant} is a known acquired resistance mutation. "
                    "Prevalence data from baseline tumor sequencing databases (cBioPortal/TCGA) is not applicable "
                    "as these mutations emerge under treatment pressure."
                )
                lines.append("")

        # DepMap data (gene dependencies and drug sensitivities)
        if self.depmap_evidence and self.depmap_evidence.has_data():
            lines.append(self.depmap_evidence.to_prompt_context())

            # Add oncogene driver context qualifier for pan-cancer dependency
            # When gene is an oncogene with modest pan-cancer dependency but tumor-matched models exist,
            # note that dependency may be higher in the relevant tumor context
            if (self.context.gene_role == "oncogene" and
                self.depmap_evidence.gene_dependency and
                not self.depmap_evidence.is_essential() and
                self.context.tumor_type):
                # Check if we have tumor-matched cell line models
                mutant_models = [cl for cl in self.depmap_evidence.cell_line_models if cl.has_mutation]
                tumor_matched_models = [
                    cl for cl in mutant_models
                    if cl.primary_disease and tumor_types_match(cl.primary_disease, self.context.tumor_type)
                ]
                if tumor_matched_models:
                    score = self.depmap_evidence.gene_dependency.mean_dependency_score
                    lines.append(
                        f"NOTE: Pan-cancer CERES ({score:.2f}) may underestimate {self.context.tumor_type}-specific "
                        f"dependency for {self.identifiers.gene}, a known oncogene driver. "
                        f"{len(tumor_matched_models)} {self.context.tumor_type} cell line model(s) available for validation."
                    )
                    lines.append("")

        # Mutation class (for oncogenes like BRAF)
        if self.context.mutation_class:
            lines.append(f"MUTATION CLASS: {self.context.mutation_class}")
            lines.append("")

        # Cancer Hotspots data (mutation recurrence across cancers)
        # Include both direct hotspot matches AND adjacent hotspots
        if self.hotspots_evidence and (
            self.hotspots_evidence.has_data() or self.hotspots_evidence.is_adjacent_to_hotspot()
        ):
            hotspot_context = self.hotspots_evidence.to_prompt_context()
            if hotspot_context:
                lines.append(hotspot_context)
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
                    pmid_link = f"[PMID {article.pmid}](https://pubmed.ncbi.nlm.nih.gov/{article.pmid}/)"
                    lines.append(f"  - {pmid_link}: {article.title[:80]}...")
                    lines.append(f"    Signal: {article.signal_type}{drugs_str}{impact_str}")
                    if summary:
                        lines.append(f"    Summary: {summary}")
                lines.append("")

            if sensitivity_articles:
                lines.append(f"SENSITIVITY LITERATURE ({len(sensitivity_articles)} articles):")
                for article in sensitivity_articles[:5]:
                    drugs_str = f" [Drugs: {', '.join(article.drugs_mentioned[:3])}]" if article.drugs_mentioned else ""
                    summary = article.get_best_summary(200)
                    pmid_link = f"[PMID {article.pmid}](https://pubmed.ncbi.nlm.nih.gov/{article.pmid}/)"
                    lines.append(f"  - {pmid_link}: {article.title[:80]}...")
                    lines.append(f"    Signal: {article.signal_type}{drugs_str}")
                    if summary:
                        lines.append(f"    Summary: {summary}")
                lines.append("")

            if other_articles:
                lines.append(f"OTHER RELEVANT LITERATURE ({len(other_articles)} articles):")
                for article in other_articles[:3]:
                    summary = article.get_best_summary(150)
                    pmid_link = f"[PMID {article.pmid}](https://pubmed.ncbi.nlm.nih.gov/{article.pmid}/)"
                    lines.append(f"  - {pmid_link}: {article.title[:80]}...")
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

    def get_gene_context_for_llm(self) -> str:
        """Get gene context (role, pathway) formatted for LLM prompt.

        Returns:
            Formatted string with gene role and pathway information
        """
        lines = []

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

        return "\n".join(lines) if lines else ""

    def get_cross_source_synthesis_for_llm(self) -> str:
        """Generate cross-source synthesis grouping evidence by drug/drug-class.

        This method identifies:
        1. Corroboration - when multiple sources agree on a drug's effect
        2. Conflicts - when sources disagree on the same drug
        3. Unique insights - evidence only found in one source

        Returns:
            Formatted string showing drug-centric evidence synthesis
        """
        from collections import defaultdict

        gene = self.identifiers.gene

        # Structure: drug_name -> {sources: [], signals: [], pmids: [], locus_levels: []}
        drug_evidence: dict[str, dict] = defaultdict(lambda: {
            "sources": [],
            "signals": [],  # "sensitivity", "resistance", "unknown"
            "pmids": [],
            "locus_levels": [],  # "variant", "codon", "gene"
            "tumor_types": [],
            "evidence_levels": [],  # "FDA", "clinical", "preclinical"
        })

        # Normalize drug names for grouping
        def normalize_drug(name: str) -> str:
            if not name:
                return ""
            # Lowercase, strip common suffixes
            n = name.lower().strip()
            # Remove parenthetical info
            if "(" in n:
                n = n.split("(")[0].strip()
            return n

        def get_signal(is_sens: bool | None, is_res: bool | None, assoc: str | None) -> str:
            if is_res or (assoc and "resist" in assoc.lower()):
                return "resistance"
            if is_sens or (assoc and ("sens" in assoc.lower() or "respon" in assoc.lower())):
                return "sensitivity"
            return "unknown"

        # 1. CGI Biomarkers (FDA-approved)
        for b in (self.cgi_biomarkers or []):
            if b.drug and not is_biomarker_selection_drug(b.drug, gene):
                key = normalize_drug(b.drug)
                if key:
                    drug_evidence[key]["sources"].append("CGI")
                    drug_evidence[key]["signals"].append(get_signal(None, None, b.association))
                    drug_evidence[key]["locus_levels"].append(b.locus_match or "gene")
                    if b.tumor_type:
                        drug_evidence[key]["tumor_types"].append(b.tumor_type)
                    drug_evidence[key]["evidence_levels"].append("FDA" if b.fda_approved else "clinical")

        # 2. CGI Preclinical/Early Phase
        for b in (self.preclinical_biomarkers or []) + (self.early_phase_biomarkers or []):
            if b.drug and not is_biomarker_selection_drug(b.drug, gene):
                key = normalize_drug(b.drug)
                if key:
                    drug_evidence[key]["sources"].append("CGI")
                    drug_evidence[key]["signals"].append(get_signal(None, None, b.association))
                    drug_evidence[key]["locus_levels"].append(b.locus_match or "gene")
                    if b.tumor_type:
                        drug_evidence[key]["tumor_types"].append(b.tumor_type)
                    level = "preclinical" if b in (self.preclinical_biomarkers or []) else "clinical"
                    drug_evidence[key]["evidence_levels"].append(level)

        # 3. CIViC Evidence
        for e in (self.civic_evidence or []):
            for drug in (e.drugs or []):
                if not is_biomarker_selection_drug(drug, gene):
                    key = normalize_drug(drug)
                    if key:
                        drug_evidence[key]["sources"].append("CIViC")
                        sig = get_signal(
                            "sens" in (e.clinical_significance or "").lower(),
                            "resist" in (e.clinical_significance or "").lower(),
                            e.clinical_significance
                        )
                        drug_evidence[key]["signals"].append(sig)
                        drug_evidence[key]["locus_levels"].append(
                            getattr(e, 'locus_variant_match', None) and
                            e.locus_variant_match.level if hasattr(e, 'locus_variant_match') and e.locus_variant_match else "variant"
                        )
                        if e.pmid:
                            drug_evidence[key]["pmids"].append(e.pmid)
                        drug_evidence[key]["evidence_levels"].append(
                            "clinical" if e.evidence_level in ("A", "B") else "preclinical"
                        )

        # 4. VICC Evidence
        for v in (self.vicc_evidence or []):
            for drug in (v.drugs or []):
                if not is_biomarker_selection_drug(drug, gene):
                    key = normalize_drug(drug)
                    if key:
                        drug_evidence[key]["sources"].append(f"VICC:{v.source}" if v.source else "VICC")
                        drug_evidence[key]["signals"].append(get_signal(v.is_sensitivity, v.is_resistance, None))
                        locus = "variant"
                        if hasattr(v, 'locus_variant_match') and v.locus_variant_match:
                            locus = v.locus_variant_match.level
                        drug_evidence[key]["locus_levels"].append(locus)
                        drug_evidence[key]["evidence_levels"].append("clinical")

        # 5. Literature (PubMed)
        for article in (self.pubmed_articles or []):
            for drug in (article.drugs_mentioned or []):
                if not is_biomarker_selection_drug(drug, gene):
                    key = normalize_drug(drug)
                    if key:
                        drug_evidence[key]["sources"].append("Literature")
                        sig = "resistance" if article.is_resistance_evidence() else (
                            "sensitivity" if article.is_sensitivity_evidence() else "unknown"
                        )
                        drug_evidence[key]["signals"].append(sig)
                        drug_evidence[key]["pmids"].append(article.pmid)
                        drug_evidence[key]["locus_levels"].append("variant")  # Literature is usually variant-specific
                        drug_evidence[key]["evidence_levels"].append("literature")

        # Now analyze and format output
        if not drug_evidence:
            return ""

        lines = []
        lines.append("CROSS-SOURCE DRUG EVIDENCE SYNTHESIS:")
        lines.append("(Drugs mentioned across multiple sources with corroboration/conflict analysis)")
        lines.append("")

        # Sort by number of sources (most corroborated first)
        sorted_drugs = sorted(
            drug_evidence.items(),
            key=lambda x: (len(set(x[1]["sources"])), len(x[1]["sources"])),
            reverse=True
        )

        corroborated = []
        conflicting = []
        single_source = []

        for drug, data in sorted_drugs:
            unique_sources = set(data["sources"])
            unique_signals = set(s for s in data["signals"] if s != "unknown")

            # Classify
            if len(unique_sources) >= 2:
                if len(unique_signals) > 1:
                    # Has both sensitivity and resistance signals
                    conflicting.append((drug, data))
                else:
                    corroborated.append((drug, data))
            else:
                single_source.append((drug, data))

        # Corroborated evidence (multiple sources agree)
        if corroborated:
            lines.append("**CORROBORATED (multiple sources agree):**")
            for drug, data in corroborated[:10]:  # Limit to top 10
                sources = sorted(set(data["sources"]))
                signal = next((s for s in data["signals"] if s != "unknown"), "unknown")
                locus = "variant" if "variant" in data["locus_levels"] else (
                    "codon" if "codon" in data["locus_levels"] else "gene"
                )
                pmids = list(set(data["pmids"]))[:3]
                pmid_str = f" (PMIDs: {', '.join(pmids)})" if pmids else ""

                lines.append(f"  • {drug.upper()}: {signal} signal")
                lines.append(f"    Sources: {', '.join(sources)} | Locus: {locus}{pmid_str}")
            lines.append("")

        # Conflicting evidence (sources disagree)
        if conflicting:
            lines.append("**CONFLICTING (sources disagree - requires investigation):**")
            for drug, data in conflicting[:5]:  # Limit to top 5
                sources = sorted(set(data["sources"]))
                signals = sorted(set(s for s in data["signals"] if s != "unknown"))
                pmids = list(set(data["pmids"]))[:3]
                pmid_str = f" (PMIDs: {', '.join(pmids)})" if pmids else ""

                lines.append(f"  • {drug.upper()}: MIXED signals ({', '.join(signals)})")
                lines.append(f"    Sources: {', '.join(sources)}{pmid_str}")
                lines.append(f"    ⚠️ Review evidence for context (resistance mutation? different tumor?)")
            lines.append("")

        # Single source (unique insights, lower confidence)
        # Only show if notable (FDA or has PMID)
        notable_single = [
            (d, data) for d, data in single_source
            if "FDA" in data["evidence_levels"] or data["pmids"]
        ]
        if notable_single:
            lines.append("**SINGLE SOURCE (lower confidence, needs validation):**")
            for drug, data in notable_single[:8]:
                source = data["sources"][0] if data["sources"] else "Unknown"
                signal = next((s for s in data["signals"] if s != "unknown"), "unknown")
                level = data["evidence_levels"][0] if data["evidence_levels"] else "unknown"
                pmids = list(set(data["pmids"]))[:2]
                pmid_str = f" (PMID: {', '.join(pmids)})" if pmids else ""

                lines.append(f"  • {drug.upper()}: {signal} ({source}, {level}){pmid_str}")
            lines.append("")

        # Summary statistics
        lines.append(f"Summary: {len(corroborated)} drugs with corroborated evidence, "
                    f"{len(conflicting)} with conflicting signals, "
                    f"{len(single_source)} from single source only")

        return "\n".join(lines)

    def _count_locus_matches(self, items: list, attr: str = 'locus_match') -> dict[str, int]:
        """Count items by locus match level.

        Args:
            items: List of evidence items with locus_match attribute
            attr: Attribute name to check for locus match level

        Returns:
            Dict with counts for variant, codon, and gene levels
        """
        counts = {'variant': 0, 'codon': 0, 'gene': 0}
        for item in items:
            level = getattr(item, attr, None) or 'gene'
            if level in counts:
                counts[level] += 1
            else:
                counts['gene'] += 1
        return counts

    def get_locus_match_summary(self) -> dict:
        """Compute summary statistics about evidence match specificity.

        Returns:
            Dict with counts and summary text for LLM prompt
        """
        variant_count = 0
        codon_count = 0
        gene_count = 0

        # FDA approvals
        fda_counts = self._count_locus_matches(self.fda_approvals)
        variant_count += fda_counts['variant']
        codon_count += fda_counts['codon']
        gene_count += fda_counts['gene']

        # VICC evidence
        vicc_counts = self._count_locus_matches(self.vicc_evidence)
        variant_count += vicc_counts['variant']
        codon_count += vicc_counts['codon']
        gene_count += vicc_counts['gene']

        # CIViC assertions
        civic_a_counts = self._count_locus_matches(self.civic_assertions)
        variant_count += civic_a_counts['variant']
        codon_count += civic_a_counts['codon']
        gene_count += civic_a_counts['gene']

        # CIViC evidence items
        civic_e_counts = self._count_locus_matches(self.civic_evidence)
        variant_count += civic_e_counts['variant']
        codon_count += civic_e_counts['codon']
        gene_count += civic_e_counts['gene']

        # CGI biomarkers
        cgi_counts = self._count_locus_matches(self.cgi_biomarkers)
        variant_count += cgi_counts['variant']
        codon_count += cgi_counts['codon']
        gene_count += cgi_counts['gene']

        # Clinical trials (use match_scope attribute)
        for trial in self.clinical_trials:
            scope = getattr(trial, 'match_scope', None)
            if scope == 'specific':
                variant_count += 1
            elif scope == 'ambiguous':
                codon_count += 1
            else:
                gene_count += 1

        # Literature knowledge - resistance and sensitivity signals
        if self.literature_knowledge:
            for entry in self.literature_knowledge.resistant_to:
                level = getattr(entry, 'locus_match', None) or 'gene'
                if level == 'variant':
                    variant_count += 1
                elif level == 'codon':
                    codon_count += 1
                else:
                    gene_count += 1
            for entry in self.literature_knowledge.sensitive_to:
                level = getattr(entry, 'locus_match', None) or 'gene'
                if level == 'variant':
                    variant_count += 1
                elif level == 'codon':
                    codon_count += 1
                else:
                    gene_count += 1

        # cBioPortal - check if we have variant-specific data
        cbioportal_variant_specific = False
        if self.cbioportal_evidence:
            if self.cbioportal_evidence.samples_with_exact_variant > 0:
                cbioportal_variant_specific = True

        total = variant_count + codon_count + gene_count

        # Build summary text
        if total == 0:
            summary_text = "No therapeutic evidence available."
        elif variant_count == total:
            summary_text = f"All {total} therapeutic evidence items are variant-specific."
        elif gene_count == total:
            summary_text = f"All {total} therapeutic evidence items are gene-level (no variant-specific data)."
        else:
            parts = []
            if variant_count > 0:
                parts.append(f"{variant_count} variant-specific")
            if codon_count > 0:
                parts.append(f"{codon_count} codon-level")
            if gene_count > 0:
                parts.append(f"{gene_count} gene-level")
            summary_text = f"Therapeutic evidence: {', '.join(parts)}."

        # Add cBioPortal note
        if self.cbioportal_evidence and self.cbioportal_evidence.has_data():
            if not cbioportal_variant_specific:
                summary_text += f" cBioPortal has gene-level data only (0 samples with exact variant)."

        return {
            "variant_count": variant_count,
            "codon_count": codon_count,
            "gene_count": gene_count,
            "total": total,
            "summary_text": summary_text,
            "has_variant_specific": variant_count > 0,
            "is_all_gene_level": gene_count == total and total > 0,
            "cbioportal_variant_specific": cbioportal_variant_specific,
        }

    def get_tumor_match_summary(self) -> dict:
        """Compute summary statistics about tumor type match specificity.

        Counts evidence items by tumor match level:
        - cancer_specific: Evidence matches the user's queried tumor type
        - pan_cancer: Evidence is tumor-agnostic (e.g., "Solid Tumor", "MSI-H")
        - other: Evidence is for a different specific tumor type
        - unknown: Evidence has no tumor type information

        Returns:
            Dict with counts, other tumor types found, and summary text for LLM prompt
        """
        cancer_specific_count = 0
        pan_cancer_count = 0
        other_count = 0
        unknown_count = 0
        other_tumor_types: set[str] = set()

        def count_tumor_match(item) -> None:
            """Count tumor match for an evidence item."""
            nonlocal cancer_specific_count, pan_cancer_count, other_count, unknown_count
            # Try tumor_match property first (bool: True = cancer_specific)
            tumor_match = getattr(item, 'tumor_match', None)
            if tumor_match is True:
                cancer_specific_count += 1
                return
            elif tumor_match is False:
                # False means either pan_cancer or other - check cancer_specificity
                pass

            # Check cancer_specificity for more detail
            specificity = getattr(item, 'cancer_specificity', None)
            if specificity == "cancer_specific":
                cancer_specific_count += 1
            elif specificity == "pan_cancer":
                pan_cancer_count += 1
            elif specificity:
                # It's a specific cancer name (e.g., "ovarian cancer")
                other_count += 1
                other_tumor_types.add(specificity)
            else:
                # No tumor info - count as unknown (NOT pan_cancer)
                unknown_count += 1

        # FDA approvals
        for approval in self.fda_approvals:
            count_tumor_match(approval)

        # VICC evidence
        for vicc in self.vicc_evidence:
            count_tumor_match(vicc)

        # CIViC assertions
        for assertion in self.civic_assertions:
            count_tumor_match(assertion)

        # CIViC evidence items
        for evidence in self.civic_evidence:
            count_tumor_match(evidence)

        # CGI biomarkers
        for biomarker in self.cgi_biomarkers:
            count_tumor_match(biomarker)

        # Clinical trials
        for trial in self.clinical_trials:
            count_tumor_match(trial)

        total = cancer_specific_count + pan_cancer_count + other_count + unknown_count
        queried_tumor = self.context.tumor_type or "unknown"

        # Build summary text
        if total == 0:
            summary_text = "No therapeutic evidence available."
        elif cancer_specific_count == total:
            summary_text = f"All {total} evidence items are specific to {queried_tumor}."
        elif other_count == total:
            others_str = ", ".join(sorted(other_tumor_types)[:3])
            summary_text = f"All {total} evidence items are from other tumor types ({others_str})."
        elif pan_cancer_count == total:
            summary_text = f"All {total} evidence items are pan-cancer/tumor-agnostic."
        elif unknown_count == total:
            summary_text = f"All {total} evidence items have unknown tumor type."
        else:
            parts = []
            if cancer_specific_count > 0:
                parts.append(f"{cancer_specific_count} {queried_tumor}-specific")
            if pan_cancer_count > 0:
                parts.append(f"{pan_cancer_count} pan-cancer")
            if other_count > 0:
                others_str = ", ".join(sorted(other_tumor_types)[:3])
                if len(other_tumor_types) > 3:
                    others_str += f" +{len(other_tumor_types) - 3} more"
                parts.append(f"{other_count} other ({others_str})")
            if unknown_count > 0:
                parts.append(f"{unknown_count} unknown")
            summary_text = f"Tumor match: {', '.join(parts)}."

        return {
            "cancer_specific_count": cancer_specific_count,
            "pan_cancer_count": pan_cancer_count,
            "other_count": other_count,
            "unknown_count": unknown_count,
            "other_tumor_types": sorted(other_tumor_types),
            "total": total,
            "summary_text": summary_text,
            "has_tumor_specific": cancer_specific_count > 0,
            "is_all_other": other_count == total and total > 0,
            "queried_tumor": queried_tumor,
        }

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

    # DepMap/CCLE (gene dependencies and drug sensitivities)
    depmap_evidence: DepMapEvidence | None = Field(
        None, description="DepMap gene dependency and drug sensitivity data"
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

                # Determine cancer_specificity based on tumor type match
                cancer_specificity = self._get_fda_cancer_specificity(approval)

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
                    match_level=approval.match_level,
                    cancer_specificity=cancer_specificity,
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

                        # Determine cancer specificity
                        cancer_specificity = self._get_cancer_specificity_from_disease(assertion.disease)

                        evidence_list.append(TherapeuticEvidence(
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
                            match_level=assertion.match_level,
                            cancer_specificity=cancer_specificity,
                        ))

        # From CGI biomarkers - FDA approved
        for biomarker in self.cgi_biomarkers:
            # Ensure drug is a non-empty string (not a list or empty)
            if biomarker.fda_approved and biomarker.drug and isinstance(biomarker.drug, str):
                drug_key = biomarker.drug.lower()
                if drug_key not in seen_drugs:
                    seen_drugs.add(drug_key)

                    response_type = None
                    if biomarker.association:
                        if "RESIST" in biomarker.association.upper():
                            response_type = "Resistance"
                        else:
                            response_type = "Sensitivity"

                    # Determine cancer specificity
                    cancer_specificity = self._get_cancer_specificity_from_disease(biomarker.tumor_type)

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
                        match_level=biomarker.match_level,
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

                        evidence_list.append(TherapeuticEvidence(
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
                            match_level=vicc.match_level,
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

                        response_type = None
                        if biomarker.association:
                            if "RESIST" in biomarker.association.upper():
                                response_type = "Resistance"
                            else:
                                response_type = "Sensitivity"

                        # Determine cancer specificity
                        cancer_specificity = self._get_cancer_specificity_from_disease(biomarker.tumor_type)

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
                            match_level=biomarker.match_level,
                            cancer_specificity=cancer_specificity,
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

    def _get_fda_cancer_specificity(self, approval) -> str:
        """Determine cancer_specificity for an FDA approval.

        Returns:
            - "cancer_specific": if the FDA indication matches the queried tumor type
            - The actual cancer type (e.g., "ovarian cancer"): if no match but we know the cancer
            - "pan_cancer": if tumor-agnostic or unknown
        """
        # First check if there's a tumor type match
        if self.context.tumor_type and approval.indication:
            parsed = approval.parse_indication_for_tumor(self.context.tumor_type)
            if parsed.get('tumor_match'):
                return "cancer_specific"

        # No match - extract what cancer the drug IS approved for
        indication_cancer = approval.extract_indication_cancer_type()
        if indication_cancer:
            # Check if it's a pan-cancer approval
            if "pan-cancer" in indication_cancer.lower():
                return "pan_cancer"
            return indication_cancer

        # Unknown
        return "pan_cancer"

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
            return "pan_cancer"

        queried_tumor = self.context.tumor_type
        if not queried_tumor:
            return "pan_cancer"

        disease_lower = disease.lower()
        queried_lower = queried_tumor.lower()

        # Check for tumor-agnostic terms
        pan_cancer_terms = [
            "solid tumor", "all solid", "pan-cancer", "any solid",
            "advanced solid", "all tumor", "tumor agnostic"
        ]
        if any(term in disease_lower for term in pan_cancer_terms):
            return "pan_cancer"

        # Build keyword mappings for flexible matching
        tumor_keywords = {
            'colorectal': ['colorectal', 'colon', 'rectal', 'crc', 'mcrc'],
            'melanoma': ['melanoma'],
            'lung': ['lung', 'nsclc', 'non-small cell'],
            'breast': ['breast'],
            'thyroid': ['thyroid', 'atc', 'anaplastic thyroid'],
            'gist': ['gist', 'gastrointestinal stromal'],
            'bladder': ['bladder', 'urothelial', 'transitional cell'],
            'cholangiocarcinoma': ['cholangiocarcinoma', 'bile duct', 'biliary'],
            'myeloproliferative': ['myelofibrosis', 'polycythemia vera', 'myeloproliferative', 'mpn'],
            'ovarian': ['ovarian', 'ovary'],
            'pancreatic': ['pancreatic', 'pancreas'],
            'prostate': ['prostate'],
            'gastric': ['gastric', 'stomach'],
            'esophageal': ['esophageal', 'esophagus'],
            'renal': ['renal', 'kidney'],
            'hepatocellular': ['hepatocellular', 'liver', 'hcc'],
            'endometrial': ['endometrial', 'uterine', 'uterus'],
            'cervical': ['cervical', 'cervix'],
            'head and neck': ['head and neck', 'hnscc', 'oral', 'laryngeal', 'pharyngeal'],
            'glioblastoma': ['glioblastoma', 'gbm', 'glioma'],
            'leukemia': ['leukemia', 'aml', 'cml', 'all', 'cll'],
            'lymphoma': ['lymphoma', 'hodgkin', 'non-hodgkin'],
            'myeloma': ['myeloma', 'multiple myeloma'],
        }

        # Find keywords for queried tumor
        queried_keywords = []
        for key, keywords in tumor_keywords.items():
            if any(kw in queried_lower for kw in keywords):
                queried_keywords = keywords
                break
        if not queried_keywords:
            queried_keywords = [queried_lower]

        # Check if disease matches queried tumor
        if any(kw in disease_lower for kw in queried_keywords):
            return "cancer_specific"

        # Also check direct substring match
        if queried_lower in disease_lower or disease_lower in queried_lower:
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

    def get_resistance_summary(self) -> str:
        """Get a synthesized summary of resistance evidence.

        Returns a narrative summary grouping drugs by class/generation where possible.
        Example: "Gatekeeper mutation conferring resistance to 1st/2nd-gen EGFR TKIs
                  (erlotinib, gefitinib, afatinib); variable response to osimertinib (conflicting evidence)"
        """
        # Collect all resistance drugs with their sources
        resistance_drugs: dict[str, set[str]] = {}  # drug -> set of sources

        # VICC resistance evidence
        for vicc in self.vicc_evidence:
            if vicc.is_resistance and vicc.drugs:
                for drug in vicc.drugs:
                    drug_lower = drug.lower()
                    if drug_lower not in resistance_drugs:
                        resistance_drugs[drug_lower] = set()
                    resistance_drugs[drug_lower].add("VICC")

        # CIViC resistance assertions
        for assertion in self.civic_assertions:
            if assertion.is_resistance and assertion.therapies:
                for therapy in assertion.therapies:
                    drug_lower = therapy.lower()
                    if drug_lower not in resistance_drugs:
                        resistance_drugs[drug_lower] = set()
                    resistance_drugs[drug_lower].add("CIViC")

        # CIViC resistance evidence items
        for evidence in self.civic_evidence:
            if evidence.clinical_significance and "RESIST" in evidence.clinical_significance.upper():
                if evidence.drugs:
                    for drug in evidence.drugs:
                        drug_lower = drug.lower()
                        if drug_lower not in resistance_drugs:
                            resistance_drugs[drug_lower] = set()
                        resistance_drugs[drug_lower].add("CIViC")

        # CGI resistance biomarkers
        for biomarker in self.cgi_biomarkers:
            if biomarker.association and "RESIST" in biomarker.association.upper():
                if biomarker.drug:
                    drug_lower = biomarker.drug.lower()
                    if drug_lower not in resistance_drugs:
                        resistance_drugs[drug_lower] = set()
                    resistance_drugs[drug_lower].add("CGI")

        # Literature resistance signals from LLM extraction
        if self.literature_knowledge and self.literature_knowledge.resistant_to:
            for entry in self.literature_knowledge.resistant_to:
                drug = entry.get("drug", "")
                if drug:
                    drug_lower = drug.lower()
                    if drug_lower not in resistance_drugs:
                        resistance_drugs[drug_lower] = set()
                    resistance_drugs[drug_lower].add("Literature")

        # PubMed articles with resistance evidence
        for article in self.pubmed_articles:
            if article.is_resistance_evidence() and article.drugs_mentioned:
                for drug in article.drugs_mentioned:
                    drug_lower = drug.lower()
                    if drug_lower not in resistance_drugs:
                        resistance_drugs[drug_lower] = set()
                    resistance_drugs[drug_lower].add(f"PubMed:{article.pmid}")

        if not resistance_drugs:
            return ""

        # Check for conflicting evidence
        conflicting_drugs: dict[str, str] = {}  # drug -> conflict description
        if self.evidence_gaps:
            from oncomind.models.evidence.evidence_gaps import GapCategory
            discordant_gaps = self.evidence_gaps.get_gaps_by_category(GapCategory.DISCORDANT)
            for gap in discordant_gaps:
                if "Conflicting drug response for" in gap.description:
                    # Extract drug name from "Conflicting drug response for X: ..."
                    parts = gap.description.split(":")
                    if parts:
                        drug_part = parts[0].replace("Conflicting drug response for", "").strip()
                        conflicting_drugs[drug_part.lower()] = gap.description

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
                drug_list = ", ".join(sorted(resistant_1st_2nd))
                summaries.append(f"Resistance to 1st/2nd-gen EGFR TKIs ({drug_list})")

            for drug in resistant_3rd:
                if drug in conflicting_drugs:
                    summaries.append(f"Variable response to {drug} (conflicting evidence)")
                else:
                    summaries.append(f"Resistance to {drug}")

            for drug in other_drugs[:3]:
                if drug in conflicting_drugs:
                    summaries.append(f"Variable response to {drug} (conflicting evidence)")
                else:
                    sources = ", ".join(sorted(resistance_drugs[drug]))
                    summaries.append(f"{drug} resistance ({sources})")
        else:
            # Generic grouping for non-EGFR genes
            # First add any conflicting drugs
            for drug in list(resistance_drugs.keys())[:8]:
                if drug in conflicting_drugs:
                    summaries.append(f"Variable response to {drug} (conflicting evidence)")
                else:
                    sources = ", ".join(sorted(resistance_drugs[drug]))
                    summaries.append(f"{drug} resistance ({sources})")

        return "; ".join(summaries[:5]) if summaries else ""

    def get_sensitivity_summary(self) -> str:
        """Get a synthesized summary of sensitivity evidence.

        Returns a narrative summary grouping drugs by evidence tier (FDA > Clinical > Preclinical).
        Now includes cancer specificity to distinguish FDA approvals that match vs don't match
        the queried tumor type.
        Example: "FDA-approved for melanoma: osimertinib; FDA-approved for OTHER cancer (ovarian cancer): avutometinib"
        """
        # Collect drugs by evidence tier and cancer specificity
        fda_matching_drugs: set[str] = set()  # FDA approved for THIS tumor
        fda_other_drugs: dict[str, str] = {}  # drug -> cancer type (for other tumors)
        clinical_drugs: set[str] = set()  # VICC, CIViC, CGI
        literature_drugs: set[str] = set()

        tumor_type = self.context.tumor_type or "cancer"

        # FDA approvals - now with cancer specificity
        for approval in self.fda_approvals:
            drug = approval.generic_name or approval.brand_name or approval.drug_name
            if drug:
                cancer_spec = self._get_fda_cancer_specificity(approval)
                if cancer_spec == "cancer_specific":
                    fda_matching_drugs.add(drug.lower())
                elif cancer_spec == "pan_cancer":
                    fda_matching_drugs.add(drug.lower())  # pan-cancer counts as matching
                else:
                    # FDA approved for a different cancer
                    fda_other_drugs[drug.lower()] = cancer_spec

        # All FDA drugs for deduplication
        all_fda_drugs = fda_matching_drugs | set(fda_other_drugs.keys())

        # VICC sensitivity evidence
        for vicc in self.vicc_evidence:
            if vicc.is_sensitivity and vicc.drugs:
                for drug in vicc.drugs:
                    drug_lower = drug.lower()
                    if drug_lower not in all_fda_drugs:
                        clinical_drugs.add(drug_lower)

        # CIViC sensitivity assertions
        for assertion in self.civic_assertions:
            if assertion.is_sensitivity and assertion.therapies:
                for therapy in assertion.therapies:
                    drug_lower = therapy.lower()
                    if drug_lower not in all_fda_drugs:
                        clinical_drugs.add(drug_lower)

        # CIViC sensitivity evidence items
        for evidence in self.civic_evidence:
            if evidence.clinical_significance:
                sig_upper = evidence.clinical_significance.upper()
                if ("SENS" in sig_upper or "RESPON" in sig_upper) and "RESIST" not in sig_upper:
                    if evidence.drugs:
                        for drug in evidence.drugs:
                            drug_lower = drug.lower()
                            if drug_lower not in all_fda_drugs:
                                clinical_drugs.add(drug_lower)

        # CGI sensitivity biomarkers
        for biomarker in self.cgi_biomarkers:
            if biomarker.association and "RESIST" not in biomarker.association.upper():
                if biomarker.drug:
                    drug_lower = biomarker.drug.lower()
                    if drug_lower not in all_fda_drugs:
                        clinical_drugs.add(drug_lower)

        # Literature sensitivity signals from LLM extraction
        if self.literature_knowledge and self.literature_knowledge.sensitive_to:
            for entry in self.literature_knowledge.sensitive_to:
                drug = entry.get("drug", "")
                if drug:
                    drug_lower = drug.lower()
                    if drug_lower not in all_fda_drugs and drug_lower not in clinical_drugs:
                        literature_drugs.add(drug_lower)

        # PubMed articles with sensitivity evidence
        for article in self.pubmed_articles:
            if article.is_sensitivity_evidence() and article.drugs_mentioned:
                for drug in article.drugs_mentioned:
                    drug_lower = drug.lower()
                    if drug_lower not in all_fda_drugs and drug_lower not in clinical_drugs:
                        literature_drugs.add(drug_lower)

        # DepMap/PRISM drug sensitivities (preclinical)
        preclinical_drugs: set[str] = set()
        if self.depmap_evidence and self.depmap_evidence.drug_sensitivities:
            for ds in self.depmap_evidence.get_top_sensitive_drugs(5):
                drug_lower = ds.drug_name.lower()
                if drug_lower not in all_fda_drugs and drug_lower not in clinical_drugs:
                    preclinical_drugs.add(drug_lower)

        # Build synthesized summary
        summaries = []

        # FDA approved for THIS tumor type
        if fda_matching_drugs:
            drug_list = ", ".join(sorted(fda_matching_drugs)[:5])
            summaries.append(f"FDA-approved for {tumor_type}: {drug_list}")

        # FDA approved for OTHER tumor types (IMPORTANT: LLM must know this is NOT for the queried tumor)
        if fda_other_drugs:
            other_parts = []
            for drug, cancer in sorted(fda_other_drugs.items())[:3]:
                other_parts.append(f"{drug} ({cancer})")
            summaries.append(f"FDA-approved for OTHER cancers (NOT {tumor_type}): {', '.join(other_parts)}")

        if clinical_drugs:
            drug_list = ", ".join(sorted(clinical_drugs)[:5])
            summaries.append(f"Clinical evidence: {drug_list}")

        if preclinical_drugs:
            drug_list = ", ".join(sorted(preclinical_drugs)[:3])
            summaries.append(f"Preclinical (DepMap): {drug_list}")

        if literature_drugs and not fda_matching_drugs and not fda_other_drugs and not clinical_drugs and not preclinical_drugs:
            # Only show literature if no higher-tier evidence
            drug_list = ", ".join(sorted(literature_drugs)[:3])
            summaries.append(f"Literature signals: {drug_list}")

        return "; ".join(summaries) if summaries else ""

    def get_evidence_summary_for_llm(self) -> str:
        """Generate a compact evidence summary for LLM prompt consumption."""
        lines = []
        tumor_type = self.context.tumor_type

        # FDA Approvals - very compact
        if self.fda_approvals:
            fda_drugs = []
            for approval in self.fda_approvals[:4]:
                drug = approval.generic_name or approval.brand_name or approval.drug_name
                if tumor_type:
                    parsed = approval.parse_indication_for_tumor(tumor_type)
                    if parsed['tumor_match']:
                        line_info = parsed['line_of_therapy']
                        fda_drugs.append(f"{drug} ({line_info})")
                    else:
                        fda_drugs.append(f"{drug} (other indication)")
                else:
                    fda_drugs.append(drug)
            lines.append(f"FDA Approved: {', '.join(fda_drugs)}")

        # CGI Biomarkers - compact
        if self.cgi_biomarkers:
            approved = [b for b in self.cgi_biomarkers if b.fda_approved]
            if approved:
                resistance = [b.drug for b in approved if b.association and 'RESIST' in b.association.upper()]
                sensitivity = [b.drug for b in approved if b.association and 'RESIST' not in b.association.upper()]
                if resistance:
                    lines.append(f"CGI Resistance: {', '.join(resistance[:3])}")
                if sensitivity:
                    lines.append(f"CGI Sensitivity: {', '.join(sensitivity[:3])}")

        # CIViC Assertions - compact
        if self.civic_assertions:
            predictive = [a for a in self.civic_assertions if a.assertion_type == "PREDICTIVE"]
            if predictive:
                civic_drugs = []
                for a in predictive[:3]:
                    therapies = ", ".join(a.therapies[:2]) if a.therapies else ""
                    sig = "sens" if a.is_sensitivity else "res" if a.is_resistance else ""
                    if therapies:
                        civic_drugs.append(f"{therapies} ({sig})")
                if civic_drugs:
                    lines.append(f"CIViC Assertions: {'; '.join(civic_drugs)}")

        # CIViC Evidence Items - compact (includes per-publication drug response data)
        if self.civic_evidence:
            civic_eid_drugs = []
            for e in self.civic_evidence[:5]:
                if e.drugs:
                    drug_str = ", ".join(e.drugs[:2])
                    sig = ""
                    if e.clinical_significance:
                        sig_upper = e.clinical_significance.upper()
                        if "RESIST" in sig_upper:
                            sig = "res"
                        elif "SENS" in sig_upper or "RESPON" in sig_upper:
                            sig = "sens"
                    if sig:
                        civic_eid_drugs.append(f"{drug_str} ({sig})")
                    else:
                        civic_eid_drugs.append(drug_str)
            if civic_eid_drugs:
                lines.append(f"CIViC Evidence: {'; '.join(civic_eid_drugs)}")

        # ClinVar - one line
        if self.clinvar_significance:
            lines.append(f"ClinVar: {self.clinvar_significance}")

        # PubMed Literature - compact
        if self.pubmed_articles:
            resistance_articles = [a for a in self.pubmed_articles if a.is_resistance_evidence()]
            if resistance_articles:
                pmids = [a.pmid for a in resistance_articles[:3]]
                lines.append(f"Literature (resistance): PMIDs {', '.join(pmids)}")

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

        # cBioPortal data
        if self.cbioportal_evidence and self.cbioportal_evidence.has_data():
            lines.append(self.cbioportal_evidence.to_prompt_context())
        else:
            lines.append("PREVALENCE: No cBioPortal data available")
            lines.append("")

        # DepMap data (gene dependencies and drug sensitivities)
        if self.depmap_evidence and self.depmap_evidence.has_data():
            lines.append(self.depmap_evidence.to_prompt_context())

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

    def _count_match_levels(self, items: list, attr: str = 'match_level') -> dict[str, int]:
        """Count items by match level.

        Args:
            items: List of evidence items with match_level attribute
            attr: Attribute name to check for match level

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

    def get_match_level_summary(self) -> dict:
        """Compute summary statistics about evidence match specificity.

        Returns:
            Dict with counts and summary text for LLM prompt
        """
        variant_count = 0
        codon_count = 0
        gene_count = 0

        # FDA approvals
        fda_counts = self._count_match_levels(self.fda_approvals)
        variant_count += fda_counts['variant']
        codon_count += fda_counts['codon']
        gene_count += fda_counts['gene']

        # VICC evidence
        vicc_counts = self._count_match_levels(self.vicc_evidence)
        variant_count += vicc_counts['variant']
        codon_count += vicc_counts['codon']
        gene_count += vicc_counts['gene']

        # CIViC assertions
        civic_a_counts = self._count_match_levels(self.civic_assertions)
        variant_count += civic_a_counts['variant']
        codon_count += civic_a_counts['codon']
        gene_count += civic_a_counts['gene']

        # CIViC evidence items
        civic_e_counts = self._count_match_levels(self.civic_evidence)
        variant_count += civic_e_counts['variant']
        codon_count += civic_e_counts['codon']
        gene_count += civic_e_counts['gene']

        # CGI biomarkers
        cgi_counts = self._count_match_levels(self.cgi_biomarkers)
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
                level = getattr(entry, 'match_level', None) or 'gene'
                if level == 'variant':
                    variant_count += 1
                elif level == 'codon':
                    codon_count += 1
                else:
                    gene_count += 1
            for entry in self.literature_knowledge.sensitive_to:
                level = getattr(entry, 'match_level', None) or 'gene'
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

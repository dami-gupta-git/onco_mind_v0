"""Evidence - Top-level per-variant evidence aggregation model.

This is the primary output of the annotation pipeline, providing a strongly-typed
container for all structured evidence gathered about a variant from databases and APIs.

For LLM narrative synthesis, see LLMInsight.
"""

from __future__ import annotations

from typing import Any
from pydantic import BaseModel, Field

from oncomind.models.evidence.civic import CIViCEvidence, CIViCAssertionEvidence
from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.cosmic import COSMICEvidence
from oncomind.models.evidence.fda import FDAApproval
from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
from oncomind.models.evidence.vicc import VICCEvidence
from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
from oncomind.models.evidence.pubmed import PubMedEvidence
from oncomind.models.evidence.literature_knowledge import LiteratureKnowledge
from oncomind.models.recommended_therapies import RecommendedTherapy


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


class KnowledgebaseEvidence(BaseModel):
    """Evidence from curated knowledgebases."""

    civic: list[CIViCEvidence] = Field(default_factory=list, description="CIViC evidence entries")
    civic_assertions: list[CIViCAssertionEvidence] = Field(
        default_factory=list, description="CIViC curated assertions"
    )
    clinvar: list[ClinVarEvidence] = Field(default_factory=list, description="ClinVar entries")
    cosmic: list[COSMICEvidence] = Field(default_factory=list, description="COSMIC entries")
    cgi_biomarkers: list[CGIBiomarkerEvidence] = Field(
        default_factory=list, description="CGI biomarker evidence"
    )
    vicc: list[VICCEvidence] = Field(default_factory=list, description="VICC MetaKB evidence")

    def has_evidence(self) -> bool:
        """Check if any knowledgebase evidence exists."""
        return bool(
            self.civic or self.civic_assertions or self.clinvar or
            self.cosmic or self.cgi_biomarkers or self.vicc
        )

    def get_evidence_sources(self) -> list[str]:
        """Get list of sources that have evidence."""
        sources = []
        if self.civic:
            sources.append("CIViC")
        if self.civic_assertions:
            sources.append("CIViC Assertions")
        if self.clinvar:
            sources.append("ClinVar")
        if self.cosmic:
            sources.append("COSMIC")
        if self.cgi_biomarkers:
            sources.append("CGI")
        if self.vicc:
            sources.append("VICC")
        return sources


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


class ClinicalContext(BaseModel):
    """Clinical context including drugs, trials, and gene context."""

    tumor_type: str | None = Field(None, description="Tumor type for context")
    tumor_type_resolved: str | None = Field(None, description="OncoTree-resolved tumor type")

    # FDA approvals
    fda_approvals: list[FDAApproval] = Field(default_factory=list, description="FDA drug approvals")

    # Clinical trials
    clinical_trials: list[ClinicalTrialEvidence] = Field(
        default_factory=list, description="Matching clinical trials"
    )

    # Gene context
    gene_role: str | None = Field(
        None, description="Gene role (oncogene, TSG, DDR, MMR, etc.)"
    )
    gene_class: str | None = Field(
        None, description="Gene functional class from context"
    )
    mutation_class: str | None = Field(
        None, description="Mutation class (e.g., BRAF Class I/II/III)"
    )
    pathway: str | None = Field(None, description="Associated pathway if relevant")

    # ClinVar significance
    clinvar_clinical_significance: str | None = Field(
        None, description="ClinVar clinical significance"
    )
    clinvar_accession: str | None = Field(None, description="ClinVar accession")

    def has_fda_approval(self) -> bool:
        """Check if any FDA approvals exist."""
        return bool(self.fda_approvals)

    def has_clinical_trials(self) -> bool:
        """Check if any clinical trials exist."""
        return bool(self.clinical_trials)

    def get_approved_drugs(self) -> list[str]:
        """Get list of FDA-approved drug names."""
        drugs = []
        for approval in self.fda_approvals:
            name = approval.brand_name or approval.generic_name or approval.drug_name
            if name and name not in drugs:
                drugs.append(name)
        return drugs


class LiteratureEvidence(BaseModel):
    """Literature and publication evidence."""

    pubmed_articles: list[PubMedEvidence] = Field(
        default_factory=list, description="PubMed articles"
    )
    literature_knowledge: LiteratureKnowledge | None = Field(
        None, description="LLM-extracted structured knowledge from literature"
    )

    # Key references
    key_pmids: list[str] = Field(default_factory=list, description="Key PubMed IDs")

    def has_literature(self) -> bool:
        """Check if any literature evidence exists."""
        return bool(self.pubmed_articles) or self.literature_knowledge is not None

    def get_resistance_drugs(self) -> list[str]:
        """Get drugs the variant is resistant to from literature."""
        if self.literature_knowledge:
            return self.literature_knowledge.get_resistance_drugs()
        return []

    def get_sensitivity_drugs(self) -> list[str]:
        """Get drugs the variant may be sensitive to from literature."""
        if self.literature_knowledge:
            return self.literature_knowledge.get_sensitivity_drugs()
        return []


class Evidence(BaseModel):
    """Top-level evidence aggregation for a single variant.

    This is the primary output of the OncoMind annotation pipeline.
    It provides strongly-typed access to all structured evidence gathered
    about a variant from databases and APIs.

    For LLM narrative synthesis, see LLMInsight.

    Structure:
    - identifiers: Core variant info and database IDs
    - kb: Knowledgebase evidence (CIViC, ClinVar, COSMIC, etc.)
    - functional: Computational predictions (AlphaMissense, CADD, etc.)
    - clinical: Clinical context (FDA, trials, gene role)
    - literature: Publications and extracted knowledge

    Example:
        >>> evidence = build_evidence("BRAF V600E", tumor_type="Melanoma")
        >>> print(evidence.identifiers.gene, evidence.identifiers.variant)
        BRAF V600E
        >>> print(evidence.clinical.get_approved_drugs())
        ['Dabrafenib', 'Vemurafenib', 'Encorafenib']
    """

    identifiers: VariantIdentifiers = Field(..., description="Variant identifiers and notation")
    kb: KnowledgebaseEvidence = Field(
        default_factory=KnowledgebaseEvidence, description="Knowledgebase evidence"
    )
    functional: FunctionalScores = Field(
        default_factory=FunctionalScores, description="Functional predictions"
    )
    clinical: ClinicalContext = Field(
        default_factory=ClinicalContext, description="Clinical context"
    )
    literature: LiteratureEvidence = Field(
        default_factory=LiteratureEvidence, description="Literature evidence"
    )

    def has_evidence(self) -> bool:
        """Check if any evidence was found."""
        return (
            self.kb.has_evidence() or
            self.clinical.has_fda_approval() or
            self.clinical.has_clinical_trials() or
            self.literature.has_literature()
        )

    def get_summary(self) -> str:
        """Generate a summary of the variant evidence.

        This is the primary human-readable summary used by both CLI and UI.

        Example:
            "BRAF V600E is a missense variant (NC_000007.14:g.140753336A>T).
             It is associated with FDA-approved therapies: Dabrafenib, Vemurafenib, Encorafenib."
        """
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
        drugs = self.clinical.get_approved_drugs()
        if drugs:
            ellipsis = "..." if len(drugs) > 3 else ""
            parts.append(f"It is associated with FDA-approved therapies: {', '.join(drugs[:3])}{ellipsis}.")

        return " ".join(parts)

    def get_summary_stats(self) -> dict[str, Any]:
        """Get summary statistics about the evidence."""
        return {
            "gene": self.identifiers.gene,
            "variant": self.identifiers.variant,
            "tumor_type": self.clinical.tumor_type,
            "evidence_sources": self.kb.get_evidence_sources(),
            "fda_approved_drugs": self.clinical.get_approved_drugs(),
            "clinical_trials_count": len(self.clinical.clinical_trials),
            "pubmed_articles_count": len(self.literature.pubmed_articles),
            "has_evidence": self.has_evidence(),
        }

    def get_recommended_therapies(self, max_results: int = 10) -> list[RecommendedTherapy]:
        """Get recommended therapies derived from FDA approvals and CGI biomarkers.

        This provides evidence-based therapy recommendations without LLM synthesis.
        For LLM-curated recommendations, see LLMInsight.recommended_therapies.

        Args:
            max_results: Maximum number of therapies to return (default 10)

        Returns:
            List of RecommendedTherapy objects, deduplicated by drug name
        """
        therapies: list[RecommendedTherapy] = []
        seen_drugs: set[str] = set()

        # From FDA approvals (highest confidence)
        for approval in self.clinical.fda_approvals:
            # Format as "Generic (Brand)" when both available
            if approval.generic_name and approval.brand_name:
                drug_name = f"{approval.generic_name} ({approval.brand_name})"
            else:
                drug_name = approval.brand_name or approval.generic_name or approval.drug_name

            # Use generic name for deduplication
            dedup_key = (approval.generic_name or approval.brand_name or approval.drug_name or "").lower()
            if drug_name and dedup_key not in seen_drugs:
                seen_drugs.add(dedup_key)
                therapies.append(RecommendedTherapy(
                    drug_name=drug_name,
                    evidence_level="A",
                    approval_status="FDA Approved",
                    clinical_context=approval.indication or "",
                ))

        # From CGI biomarkers (FDA-approved only)
        for biomarker in self.kb.cgi_biomarkers:
            if biomarker.fda_approved and biomarker.drug:
                dedup_key = biomarker.drug.lower()
                if dedup_key not in seen_drugs:
                    seen_drugs.add(dedup_key)
                    therapies.append(RecommendedTherapy(
                        drug_name=biomarker.drug,
                        evidence_level=biomarker.evidence_level or "A",
                        approval_status="FDA Approved (CGI)",
                        clinical_context=biomarker.tumor_type or "",
                    ))

        return therapies[:max_results]

    def to_flat_dict(self) -> dict[str, Any]:
        """Convert to a flat dictionary for DataFrame/tabular output.

        This is useful for creating DataFrames from batch results.
        """
        flat = {
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
            "tumor_type": self.clinical.tumor_type,
            "gene_role": self.clinical.gene_role,
            "clinvar_significance": self.clinical.clinvar_clinical_significance,
            "fda_approved_drugs": ", ".join(self.clinical.get_approved_drugs()),
            "clinical_trials_count": len(self.clinical.clinical_trials),

            # Functional
            "alphamissense_score": self.functional.alphamissense_score,
            "alphamissense_prediction": self.functional.alphamissense_prediction,
            "cadd_score": self.functional.cadd_score,
            "gnomad_af": self.functional.gnomad_exome_af,

            # Counts
            "civic_evidence_count": len(self.kb.civic),
            "vicc_evidence_count": len(self.kb.vicc),
            "pubmed_articles_count": len(self.literature.pubmed_articles),
        }
        return flat

    class Config:
        """Pydantic config."""
        json_schema_extra = {
            "example": {
                "identifiers": {
                    "variant_id": "BRAF:V600E",
                    "gene": "BRAF",
                    "variant": "V600E",
                    "cosmic_id": "COSM476",
                },
                "kb": {
                    "civic": [],
                    "clinvar": [],
                },
                "functional": {
                    "alphamissense_score": 0.98,
                    "cadd_score": 32.0,
                },
                "clinical": {
                    "tumor_type": "Melanoma",
                    "fda_approvals": [],
                },
                "literature": {
                    "pubmed_articles": [],
                },
            }
        }

    def get_evidence_summary_for_llm(self) -> str:
        """Generate a compact evidence summary for LLM prompt consumption.

        Returns a string with FDA approvals, CGI biomarkers, CIViC assertions,
        ClinVar significance, and literature - formatted for the LLM to synthesize.
        """
        lines = [f"Evidence for {self.identifiers.gene} {self.identifiers.variant}:\n"]
        tumor_type = self.clinical.tumor_type

        # FDA Approvals
        if self.clinical.fda_approvals:
            lines.append(f"FDA Approved Drugs ({len(self.clinical.fda_approvals)}):")
            for approval in self.clinical.fda_approvals[:5]:
                drug = approval.brand_name or approval.generic_name or approval.drug_name
                if tumor_type:
                    parsed = approval.parse_indication_for_tumor(tumor_type)
                    if parsed['tumor_match']:
                        line_info = parsed['line_of_therapy'].upper()
                        approval_info = parsed['approval_type'].upper()
                        lines.append(f"  • {drug} [FOR {tumor_type.upper()}]:")
                        lines.append(f"      Line of therapy: {line_info}")
                        lines.append(f"      Approval type: {approval_info}")
                        lines.append(f"      Excerpt: {parsed['indication_excerpt'][:200]}...")
                    else:
                        indication = (approval.indication or "")[:300]
                        lines.append(f"  • {drug} [DIFFERENT TUMOR TYPE]: {indication}...")
                else:
                    indication = (approval.indication or "")[:300]
                    date_str = f" (Approved: {approval.approval_date})" if approval.approval_date else ""
                    lines.append(f"  • {drug}{date_str}: {indication}...")
            lines.append("")

        # CGI Biomarkers
        if self.kb.cgi_biomarkers:
            approved = [b for b in self.kb.cgi_biomarkers if b.fda_approved]
            if approved:
                resistance = [b for b in approved if b.association and 'RESIST' in b.association.upper()]
                sensitivity = [b for b in approved if b.association and 'RESIST' not in b.association.upper()]

                if resistance:
                    lines.append(f"CGI FDA-APPROVED RESISTANCE MARKERS ({len(resistance)}):")
                    for b in resistance[:5]:
                        lines.append(f"  • {b.drug} [{b.association.upper()}] in {b.tumor_type or 'solid tumors'} - Evidence: {b.evidence_level}")
                    lines.append("")

                if sensitivity:
                    lines.append(f"CGI FDA-Approved Sensitivity Biomarkers ({len(sensitivity)}):")
                    for b in sensitivity[:5]:
                        lines.append(f"  • {b.drug} [{b.association}] in {b.tumor_type or 'solid tumors'} - Evidence: {b.evidence_level}")
                    lines.append("")

        # CIViC Assertions
        if self.kb.civic_assertions:
            predictive_tier_i = [a for a in self.kb.civic_assertions
                                  if a.amp_tier == "Tier I" and a.assertion_type == "PREDICTIVE"]
            if predictive_tier_i:
                lines.append(f"CIViC PREDICTIVE TIER I ASSERTIONS ({len(predictive_tier_i)}):")
                for a in predictive_tier_i[:5]:
                    therapies = ", ".join(a.therapies) if a.therapies else "N/A"
                    lines.append(f"  • {a.molecular_profile}: {therapies} [{a.significance}]")
                    lines.append(f"      AMP Level: {a.amp_level}, Disease: {a.disease}")
                lines.append("")

        # ClinVar
        if self.clinical.clinvar_clinical_significance:
            lines.append(f"ClinVar: {self.clinical.clinvar_clinical_significance}")
            lines.append("")

        # PubMed Literature
        if self.literature.pubmed_articles:
            resistance_articles = [a for a in self.literature.pubmed_articles if a.is_resistance_evidence()]
            if resistance_articles:
                lines.append(f"PUBMED RESISTANCE LITERATURE ({len(resistance_articles)} articles):")
                for article in resistance_articles[:3]:
                    drugs_str = f" [Drugs: {', '.join(article.drugs_mentioned[:3])}]" if article.drugs_mentioned else ""
                    lines.append(f"  • PMID {article.pmid}: {article.title[:100]}...{drugs_str}")
                lines.append("")

        return "\n".join(lines) if len(lines) > 1 else ""



"""EvidencePanel - Top-level per-variant evidence aggregation model.

This is the primary output of the annotation pipeline, providing a strongly-typed
container for all evidence gathered about a variant.
"""

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


class EvidenceMeta(BaseModel):
    """Metadata about the evidence collection."""

    sources_queried: list[str] = Field(default_factory=list, description="Data sources queried")
    sources_with_data: list[str] = Field(default_factory=list, description="Sources that returned data")
    sources_failed: list[str] = Field(default_factory=list, description="Sources that failed")

    evidence_strength: str | None = Field(
        None, description="Overall evidence strength (Strong/Moderate/Weak)"
    )

    # Experimental - not part of core API
    experimental_tier: str | None = Field(
        None, description="Experimental AMP/ASCO/CAP tier (not authoritative)"
    )

    processing_notes: list[str] = Field(
        default_factory=list, description="Notes from processing"
    )


class EvidencePanel(BaseModel):
    """Top-level evidence aggregation for a single variant.

    This is the primary output of the OncoMind annotation pipeline.
    It provides strongly-typed access to all evidence gathered about a variant.

    Structure:
    - identifiers: Core variant info and database IDs
    - kb: Knowledgebase evidence (CIViC, ClinVar, COSMIC, etc.)
    - functional: Computational predictions (AlphaMissense, CADD, etc.)
    - clinical: Clinical context (FDA, trials, gene role)
    - literature: Publications and extracted knowledge
    - meta: Processing metadata

    Example:
        >>> panel = get_insight("BRAF V600E", tumor_type="Melanoma")
        >>> print(panel.identifiers.gene, panel.identifiers.variant)
        BRAF V600E
        >>> print(panel.clinical.get_approved_drugs())
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
    meta: EvidenceMeta = Field(
        default_factory=EvidenceMeta, description="Processing metadata"
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
            parts.append(f"It is associated with FDA-approved therapies: {', '.join(drugs[:3])}.")

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
            "evidence_strength": self.meta.evidence_strength,
            "has_evidence": self.has_evidence(),
        }

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

            # Meta
            "evidence_strength": self.meta.evidence_strength,
            "sources_with_data": ", ".join(self.meta.sources_with_data),
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
                "meta": {
                    "evidence_strength": "Strong",
                },
            }
        }

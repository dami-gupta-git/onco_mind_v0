from oncomind.config.constants import BROAD_VARIANTS
from oncomind.models.evidence.base import (
    EvidenceItemBase,
    EvidenceLevel,
    VariantGeneLevel,
    Scope,
    Origin,
    CancerSpecificity,
    is_ambiguous_variant,
    extract_variant_codon,
    extract_variant_position,
)
from oncomind.models.evidence.cbioportal import CBioPortalEvidence, CoMutationEntry
from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
from oncomind.models.evidence.civic import CIViCEvidence, CIViCAssertionEvidence
from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
from oncomind.models.evidence.cosmic import COSMICEvidence
from oncomind.models.evidence.depmap import (
    DepMapEvidence,
    GeneDependency,
    DrugSensitivity,
    CellLineModel,
)
from oncomind.models.evidence.hotspots import (
    HotspotsEvidence,
    HotspotEntry,
    HotspotVariant,
    HotspotTumorType,
)
from oncomind.models.evidence.evidence import (
    # Core model
    Evidence,
    VariantIdentifiers,
    FunctionalScores,
    VariantContext,
)
from oncomind.models.extracted.evidence_gaps import (
    EvidenceGap,
    EvidenceGaps,
    GapCategory,
    GapSeverity,
)
from oncomind.models.evidence.fda import (
    FDAApproval,
    FDALabelEvidence,
    ClinicalStudyEvidence,
    MechanismEvidence,
    AdverseReactionsEvidence,
    CombinationPartner,
)
from oncomind.models.extracted.literature_knowledge import LiteratureKnowledge
from oncomind.models.evidence.pubmed import PubMedEvidence
from oncomind.models.evidence.vicc import VICCEvidence
from oncomind.models.evidence.tumor_evidence import TumorEvidenceMatch, SourceMatch

__all__ = [
    # Base class and types
    "EvidenceItemBase",
    "EvidenceLevel",
    "VariantGeneLevel",
    "Scope",
    "Origin",
    "CancerSpecificity",
    "BROAD_VARIANTS",
    # Utility functions
    "is_ambiguous_variant",
    "extract_variant_codon",
    "extract_variant_position",
    # Core model
    "Evidence",
    # Core components
    "VariantIdentifiers",
    "FunctionalScores",
    "VariantContext",
    # Evidence gaps
    "EvidenceGap",
    "EvidenceGaps",
    "GapCategory",
    "GapSeverity",
    # Individual evidence types (used in Evidence lists)
    "CBioPortalEvidence",
    "CoMutationEntry",
    "CGIBiomarkerEvidence",
    "CIViCEvidence",
    "CIViCAssertionEvidence",
    "ClinicalTrialEvidence",
    "ClinVarEvidence",
    "COSMICEvidence",
    "DepMapEvidence",
    "GeneDependency",
    "DrugSensitivity",
    "CellLineModel",
    "FDAApproval",
    "FDALabelEvidence",
    "ClinicalStudyEvidence",
    "MechanismEvidence",
    "AdverseReactionsEvidence",
    "CombinationPartner",
    "HotspotsEvidence",
    "HotspotEntry",
    "HotspotVariant",
    "HotspotTumorType",
    "LiteratureKnowledge",
    "PubMedEvidence",
    "VICCEvidence",
    # Tumor evidence aggregation
    "TumorEvidenceMatch",
    "SourceMatch",
]

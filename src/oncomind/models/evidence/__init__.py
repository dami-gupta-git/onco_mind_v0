from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
from oncomind.models.evidence.civic import CIViCEvidence, CIViCAssertionEvidence
from oncomind.models.evidence.clinvar import ClinVarEvidence
from oncomind.models.evidence.clinical_trials import ClinicalTrialEvidence
from oncomind.models.evidence.cosmic import COSMICEvidence
from oncomind.models.evidence.evidence import EvidenceForLLM
from oncomind.models.evidence.evidence_panel import (
    EvidencePanel,
    VariantIdentifiers,
    KnowledgebaseEvidence,
    FunctionalScores,
    ClinicalContext,
    LiteratureEvidence,
    EvidenceMeta,
)
from oncomind.models.evidence.fda import FDAApproval
from oncomind.models.evidence.literature_knowledge import LiteratureKnowledge
from oncomind.models.evidence.pubmed import PubMedEvidence
from oncomind.models.evidence.vicc import VICCEvidence

__all__ = [
    # Core aggregation models
    "EvidenceForLLM",
    "EvidencePanel",
    # EvidencePanel components
    "VariantIdentifiers",
    "KnowledgebaseEvidence",
    "FunctionalScores",
    "ClinicalContext",
    "LiteratureEvidence",
    "EvidenceMeta",
    # Individual evidence types
    "CGIBiomarkerEvidence",
    "CIViCEvidence",
    "CIViCAssertionEvidence",
    "ClinicalTrialEvidence",
    "ClinVarEvidence",
    "COSMICEvidence",
    "FDAApproval",
    "LiteratureKnowledge",
    "PubMedEvidence",
    "VICCEvidence",
]

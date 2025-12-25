"""Evidence data models from external databases.

This module previously contained EvidenceForLLM, which has been removed.
Evidence models are now organized in individual modules:
- cgi.py: CGIBiomarkerEvidence
- civic.py: CIViCEvidence, CIViCAssertionEvidence
- clinical_trials.py: ClinicalTrialEvidence
- clinvar.py: ClinVarEvidence
- cosmic.py: COSMICEvidence
- fda.py: FDAApproval
- pubmed.py: PubMedEvidence
- vicc.py: VICCEvidence
- literature_knowledge.py: LiteratureKnowledge
- insight.py: Insight (main aggregated evidence model)
"""

# This file is kept for backwards compatibility during migration.
# All evidence models should be imported from their individual modules
# or from oncomind.models.evidence (the package __init__.py).

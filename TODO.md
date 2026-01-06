# TODO

## Code Quality

- [ ] Check if need to use `_tumor_type_matches` everywhere for tumor matching consistency, e.g. in `evidence.py:524` (clinical trial cancer specificity matching). The curated lists in `gene_context.py` are already normalized so simple substring matching is appropriate there, but external data sources (ClinicalTrials.gov, CIViC, CGI) use varied terminology that may benefit from alias resolution.

## Data Sources

- [ ] **DGIdb Integration**: Replace the manual `BIOMARKER_SELECTION_DRUGS` mapping in `constants.py` with DGIdb (Drug-Gene Interaction Database) integration for comprehensive drug-target relationships. This will allow us to automatically distinguish between:
  - Drugs that directly target a gene (e.g., osimertinib targets EGFR)
  - Drugs where a gene mutation is a patient selection biomarker but not the drug target (e.g., DATROWAY/datopotamab deruxtecan targets TROP2, but is indicated for EGFR-mutant NSCLC patients)

  DGIdb is free and provides ~40k drug-gene interactions with mechanism annotations.
  API: https://dgidb.org/api

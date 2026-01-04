# TODO

## Code Quality

- [ ] Check if need to use `_tumor_type_matches` everywhere for tumor matching consistency, e.g. in `evidence.py:524` (clinical trial cancer specificity matching). The curated lists in `gene_context.py` are already normalized so simple substring matching is appropriate there, but external data sources (ClinicalTrials.gov, CIViC, CGI) use varied terminology that may benefit from alias resolution.

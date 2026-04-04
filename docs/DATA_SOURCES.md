# Data Sources & Caching

This document describes the external data sources used by OncoMind, their caching strategy, and data files.

## Data Directory Structure

Each data source has its own subdirectory under `data/`:

```
data/
├── cgi/                    # Cancer Genome Interpreter
│   └── cgi_biomarkers.tsv
├── depmap/                 # DepMap cell line data
│   ├── OmicsSomaticMutations.csv
│   ├── primary-screen-replicate-collapsed-logfold-change.csv
│   ├── primary-screen-replicate-collapsed-treatment-info.csv
│   └── sample_info.csv
├── fda/                    # Raw FDA downloads
│   └── Table of Pharmacogenomic Biomarkers...xlsx
├── hotspots_msk/           # MSK cancer hotspots
│   └── hotspots.txt
├── samples/                # Sample/template files
│   └── sample_fda_oncology_biomarkers.xlsx
└── processing/             # Filtered/generated/cached files
    ├── fda_oncology_biomarkers.xlsx
    └── labels/
        └── fda_labels.json
```

Path constants are defined in `src/oncomind/config/constants.py`:
- `CGI_DATA_DIR` - CGI data directory
- `FDA_DATA_DIR` - Raw FDA downloads
- `PROCESSING_DATA_DIR` - Filtered/generated/cached files
- `CGI_BIOMARKERS_FILE` - CGI biomarkers TSV
- `FDA_LABELS_FILE` - FDA drug labels JSON cache

## Auto-Downloaded Data Files

| File | Location | Source | TTL | Description |
|------|----------|--------|-----|-------------|
| `cgi_biomarkers.tsv` | `data/cgi/` | [CGI](https://www.cancergenomeinterpreter.org) | 7 days | Cancer biomarkers with drug associations |

## Data Source Details

### CGI Biomarkers (`data/cgi/cgi_biomarkers.tsv`)  
Contains FDA-approved and clinical-stage biomarkers linking gene alterations to drug response.  

**Source:** Cancer Genome Interpreter  
**URL:** https://www.cancergenomeinterpreter.org/biomarkers  
**Format:** TSV with columns: gene, alteration, drug, association, evidence_level, source, tumor_type  
**Refresh:** Auto-downloaded if missing or >7 days old  

### FDA  
#### Biomarkers (`data/fda/Table of Pharmacogenomic Biomarkers in Drug Labeling  FDA.xlsx`)  
Raw data file -  FDA biomarker table (original download from FDA).  
**URL:** https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling  
**Refresh:** Manual: Should update quarterly. It should be filtered manually to generate
fda_oncology_biomarkers.xlsx

#### Drug Labels (`data/processing/labels/fda_labels.json`)  
FDA drug label indications from OpenFDA API. Initial approval dates are fetched from the drugsfda endpoint.  
**Source:** OpenFDA Drug Labels API + DrugsFDA API    
**URLs:**  
- https://api.fda.gov/drug/label.json (label content)  
- https://api.fda.gov/drug/drugsfda.json (approval dates)  
**Refresh:** Populated on-demand during queries, or manually by prefetch script every quarter  

### MSK Hotspots (`data/hotspots_msk/hotspots.txt`)  
Curated list of statistically significant mutation hotspots in cancer.  

**Source:** Cancer Hotspots (cancerhotspots.org)  
**Refresh:** Static file, TBD if update is needed  

### DepMap (`data/depmap/`)  
Cell line dependency and drug sensitivity data.  
**Source:** DepMap Portal (Broad Institute)  
**URL:** https://depmap.org/portal/  
**Files:**  
- `OmicsSomaticMutations.csv` - Cell line mutation profiles  
- `primary-screen-replicate-collapsed-logfold-change.csv` - PRISM drug sensitivity data  
- `primary-screen-replicate-collapsed-treatment-info.csv` - PRISM treatment info  
- `sample_info.csv` - Cell line metadata  
**Refresh:** Manual: should download quarterly  


## FDA Labels Setup

To populate the FDA labels cache:

1. Download the FDA Biomarkers file from FDA website
2. Filter it to create `fda_oncology_biomarkers.xlsx` (oncology drugs only)
   - There is a sample file at `data/samples/sample_fda_oncology_biomarkers.xlsx`
3. Run the prefetch script:

```bash
python scripts/prefetch_fda_labels.py
```

**Input:** `data/processing/fda_oncology_biomarkers.xlsx`
**Output:** `data/processing/labels/fda_labels.json`

The cache is also populated on-demand when querying drugs not yet cached.

## Manually Updated Data Files

| File                                                    | Location                 | Source                                                                                                     | Cadence | Next Refresh | Description                            |
|---------------------------------------------------------|--------------------------|------------------------------------------------------------------------------------------------------------|---------|--------------|----------------------------------------|
| `Table of Pharmacogenomic Biomarkers...xlsx`            | `data/fda/`              | [FDA](https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling) | Quarterly | 2025-04-01 | Raw FDA biomarker table                |
| `fda_oncology_biomarkers.xlsx`                          | `data/processing/`       | Filtered from above                                                                                        | Quarterly | 2025-04-01 | Oncology-only subset of FDA biomarkers |
| `OmicsSomaticMutations.csv`                             | `data/depmap/`           | [DepMap](https://depmap.org/portal/)                                                                       | Quarterly | 2025-04-01 | Cell line mutation profiles            |
| `primary-screen-replicate-collapsed-logfold-change.csv` | `data/depmap/`           | [DepMap](https://depmap.org/portal/)                                                                       | Quarterly | 2025-04-01 | PRISM drug sensitivity data            |
| `primary-screen-replicate-collapsed-treatment-info.csv` | `data/depmap/`           | [DepMap](https://depmap.org/portal/)                                                                       | Quarterly | 2025-04-01 | PRISM treatment info                   |
| `sample_info.csv`                                       | `data/depmap/`           | [DepMap](https://depmap.org/portal/)                                                                       | Quarterly | 2025-04-01 | Cell line metadata                     |
| `hotspots.txt`                                          | `data/hotspots_msk/`     | [Cancer Hotspots](https://cancerhotspots.org)                                                              | TBD | TBD | MSK mutation hotspots                  |

## Processing/Cached Files

Files in `data/processing/` are filtered, generated, or cached:
- `fda_oncology_biomarkers.xlsx` - FDA biomarker table filtered to oncology drugs only
- `labels/fda_labels.json` - Cache of FDA label data (populated by prefetch script or on-demand)

## Adding New Data Sources

When adding a new data source:

1. Create directory under `data/` (e.g., `data/new_source/`)
2. Add path constants to `src/oncomind/config/constants.py`
3. Update the API client to use the constants
4. Implement `_cache_is_valid()` with appropriate TTL
5. Implement `_download_data()` for fetching
6. Add entry to this document

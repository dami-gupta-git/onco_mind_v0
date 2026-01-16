# Annotate Workflow

VCF-based variant annotation pipeline for OncoMind.

## Overview

The `mind annotate` command provides a separate workflow from the existing `mind insight` command. It takes a VCF file as input, fetches annotations from MyVariant.info API, and outputs JSON.

```bash
mind annotate variants.vcf
mind annotate variants.vcf --output results.json
```

## Input

Standard VCF format (VCFv4.2). Required columns:
- CHROM
- POS
- ID
- REF
- ALT

Required INFO fields:
- `GENE` - Gene symbol (e.g., BRAF)
- `PROTEIN` - Protein change notation (e.g., p.V600E)
- `TUMOR_TYPE` - Tumor type (e.g., NSCLC, Melanoma)

Example VCF line:
```
7	140453136	COSV54736915	T	A	.	PASS	GENE=BRAF;PROTEIN=p.V600E;TUMOR_TYPE=Melanoma
```

## Output

JSON structure with full MyVariant.info annotations:

```json
{
  "variants": [
    {
      "gene": "BRAF",
      "protein": "p.V600E",
      "tumor_type": "Melanoma",
      "myvariant": {
        "vcf": {"ref": "A", "alt": "T"},
        "cadd": {
          "phred": 32,
          "rawscore": 6.641785,
          "consequence": "NON_SYNONYMOUS",
          "gene": {
            "genename": "BRAF",
            "gene_id": "ENSG00000157764",
            "feature_id": "ENST00000288602"
          }
        },
        "clinvar": {
          "allele_id": 29000,
          "gene": {"symbol": "BRAF", "id": "673"},
          "rcv": [{"accession": "RCV000014992", "clinical_significance": "Pathogenic"}]
        },
        "cosmic": {"cosmic_id": "COSM476"},
        "dbnsfp": {
          "polyphen2": {"hdiv": {"pred": "D", "score": 0.971}},
          "alphamissense": {"pred": ["P"], "score": [0.9853]}
        },
        "dbsnp": {"rsid": "rs113488022"},
        "gnomad_exome": {"af": {"af": 3.97994e-06}}
      },
      "timings": {"myvariant": 0.234}
    }
  ]
}
```

### MyVariant Fields

| Field | Description |
|-------|-------------|
| `vcf` | Reference and alternate alleles |
| `cadd` | CADD pathogenicity scores (phred, rawscore, consequence) |
| `clinvar` | ClinVar clinical significance and RCV records |
| `cosmic` | COSMIC mutation ID |
| `dbnsfp` | Functional predictions (PolyPhen2, SIFT, AlphaMissense) |
| `dbsnp` | dbSNP rsID |
| `gnomad_exome` | gnomAD exome allele frequency |

## CLI Options

| Option | Description |
|--------|-------------|
| `--output`, `-o` | Output JSON file (default: stdout) |
| `--log-level`, `-l` | Log level: DEBUG, INFO, WARN, ERROR |

## Architecture

```
VCF File → Annotator → MyVariant API → AnnotationResult → JSON
              │
              ├── AnnInputParser (VCF parsing)
              ├── AnnotatorMyVariantClient (API client)
              └── MyVariantAnnotation (Pydantic model)
```

### Components

- **Annotator** - Main pipeline orchestrator (`annotator/annotator.py`)
- **AnnInputParser** - Parses VCF into `AnnParsedVariant` objects (`annotator/ann_input_parser.py`)
- **AnnotatorMyVariantClient** - Async HTTP client for MyVariant.info API (`annotator/clients/myvariant.py`)
- **MyVariantAnnotation** - Pydantic model with validated fields (`annotator/models/myvariant.py`)

### Pydantic Models

All MyVariant API responses are validated through Pydantic models with:
- Type validation (patterns for IDs like `ENSG*`, `COSM*`, `rs*`)
- Range validation (phred 0-100, allele frequency 0-1)
- Nested models for complex structures (CADD gene info, ClinVar RCV records)

## Programmatic Usage

```python
from oncomind.annotator import Annotator

async with Annotator() as annotator:
    result = await annotator.run("variants.vcf")

    for variant in result.variants:
        print(f"{variant['gene']} {variant['protein']}")
        print(f"  ClinVar: {variant['myvariant']['clinvar']}")
        print(f"  CADD phred: {variant['myvariant']['cadd']['phred']}")
```

## Tests

```bash
# Unit tests (mocked API)
pytest tests/unit/annotate/ -v

# Integration tests (real API calls)
pytest tests/integration/annotate/ -v -m integration
```

### Test Fixtures

VCF fixtures in `tests/integration/annotate/fixtures/`:

| File | Variants | Description |
|------|----------|-------------|
| `cancer_panel.vcf` | 6 | BRAF, KRAS, TP53, EGFR, PTEN, NRAS |
| `single_variant.vcf` | 1 | BRAF V600E |
| `empty.vcf` | 0 | Headers only |

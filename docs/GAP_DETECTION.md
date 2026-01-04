# Evidence Gap Detection

## The Core Insight

Most variant annotation tools answer the question: *"What do we know about this variant?"*

OncoMind answers a different question: ***"What don't we know yet?"***

For well-characterized variants like BRAF V600E or EGFR L858R, the clinical answers are already in databases—FDA approvals, clinical guidelines, established therapies. The interesting research questions live in the gaps: variants with conflicting evidence, understudied tumor contexts, unknown resistance mechanisms, or missing preclinical validation.

**Evidence gap detection is OncoMind's core differentiator.** By systematically identifying what's missing, where sources conflict, and what research questions remain open, we help researchers prioritize investigations and identify opportunities for novel contributions.

---

## How It Works

OncoMind analyzes evidence across 10 gap categories with context-aware severity ratings:

| Category | What We Check | Gap Example |
|----------|---------------|-------------|
| **FUNCTIONAL** | AlphaMissense, CADD, PolyPhen2, gene mechanism | "Functional impact of R248W on TP53 protein unknown" |
| **CLINICAL** | CIViC assertions, FDA approvals | "No curated clinical evidence for MAP2K1 K57N" |
| **TUMOR_TYPE** | Evidence specific to queried tumor | "No evidence specific to cholangiocarcinoma for IDH1 R132H" |
| **DRUG_RESPONSE** | CGI, VICC, DepMap drug sensitivity | "No drug sensitivity/resistance data for ARID1A Q1328*" |
| **RESISTANCE** | Resistance literature, CGI resistance markers | "Resistance mechanisms for ALK G1202R not well characterized" |
| **PREVALENCE** | cBioPortal mutation frequency | "Prevalence of FGFR2 S252W in bladder cancer unknown" |
| **PRECLINICAL** | DepMap cell line models | "No cell line models identified for NF1 R1513*" |
| **PROGNOSTIC** | Survival/outcome data | "Prognostic impact unknown" |
| **DISCORDANT** | Conflicting evidence between sources | "Conflicting drug response for imatinib: sensitive (CIViC) vs resistant (CGI)" |
| **VALIDATION** | Strong oncogenic signal but limited therapeutic evidence | "Strong oncogenic signal but limited therapeutic validation" |

---

## Severity Levels

| Severity | Meaning | Weighted Score |
|----------|---------|----------------|
| **CRITICAL** | No data at all in a key area | ×3.0 |
| **SIGNIFICANT** | Limited data, needs more research | ×2.0 |
| **MINOR** | Some data exists but could be deeper | ×1.0 |

Severity is now **context-aware**:
- Pathogenic signal (AlphaMissense, CADD ≥20, ClinVar pathogenic) increases gap severity
- Known cancer genes get higher severity ratings than unknown genes
- Tumor-type gaps are CRITICAL for cancer genes with pathogenic variants but no clinical evidence

---

## Research-Oriented Weighted Scoring

Gap categories are weighted by research value (biological gaps > clinical gaps):

| Category | Weight | Rationale |
|----------|--------|-----------|
| **VALIDATION** | 3.5 | Strong signal + no validation = prime research target |
| **FUNCTIONAL** | 3.0 | Mechanism unknown = high research value |
| **PRECLINICAL** | 2.5 | No models to test hypotheses |
| **RESISTANCE** | 2.0 | Resistance mechanisms unknown |
| **DISCORDANT** | 2.0 | Conflicting evidence needs resolution |
| **DRUG_RESPONSE** | 1.5 | Drug sensitivity unknown |
| **TUMOR_TYPE** | 1.5 | Not studied in this tumor |
| **PREVALENCE** | 1.0 | Epidemiology unknown |
| **CLINICAL** | 1.0 | Lower weight for research context |
| **PROGNOSTIC** | 1.0 | Prognostic impact unknown |

**Overall Quality** uses net scoring (gap penalty minus well-characterized credit):

```
net_score = gap_score - (well_characterized_count × 1.5)
```

| Net Score | Quality |
|-----------|---------|
| < 0 | `comprehensive` |
| 0–6 | `moderate` |
| 6–12 | `limited` |
| ≥ 12 | `minimal` |

---

## Research Priority Levels

| Priority | Criteria | Icon |
|----------|----------|------|
| **very_high** | Strong oncogenic signal (pathogenic + essential gene) + biological gaps; OR hotspot-adjacent + pathogenic + biological gaps | 🔥 |
| **high** | Cancer gene with critical gaps; OR hotspot-adjacent in cancer gene | 🔴 |
| **medium** | Any critical gaps; OR cancer gene with significant gaps | 🟡 |
| **low** | Comprehensive evidence quality AND no critical/significant gaps | 🟢 |

**Biological gaps** = gaps in VALIDATION, FUNCTIONAL, or PRECLINICAL categories.

---

## Hotspot Context Detection

OncoMind detects variant proximity to known cancer hotspots (source: cancerhotspots.org, COSMIC, OncoKB):

### Known Hotspots
- **BRAF**: 600, 469, 601, 594, 597
- **KRAS**: 12, 13, 61, 117, 146
- **EGFR**: 719, 746, 790, 858, 861
- **PIK3CA**: 545, 542, 1047, 1049, 420
- **TP53**: 175, 245, 248, 249, 273, 282
- *...and 25+ more genes*

### Hotspot vs Adjacent

| Variant | Status | Well-Characterized Tag | Research Priority |
|---------|--------|------------------------|-------------------|
| BRAF V600E | Hotspot | "known cancer hotspot" | Depends on other factors |
| BRAF V598E | Adjacent (within 5 codons of 600) | "near hotspot codon 600 — structural hypothesis likely" | Boosted to high/very_high |
| BRAF V500E | Neither | — | Normal calculation |

**Why hotspot-adjacent matters**: Rare variants near activating hotspots are research gold—structural similarity suggests similar functional impact, but they lack the extensive characterization of the hotspot itself.

---

## Discordant Evidence Detection

OncoMind detects conflicts between data sources:

### Drug Response Conflicts
```
Conflicting drug response for imatinib: sensitive (CIViC) vs resistant (CGI)
```

Checks across: CIViC, CGI biomarkers, VICC MetaKB

### ClinVar Conflicts
```
ClinVar has conflicting interpretations: both pathogenic and benign submissions
```

---

## Validation Gap (Oncogenicity Potential)

A special gap category for high-potential research targets:

**Triggers when:**
1. Variant has pathogenic signal (AlphaMissense P, CADD ≥20, ClinVar pathogenic, or truncating)
2. Gene is essential in cancer cells (DepMap)
3. BUT no therapeutic validation (no CIViC assertions, FDA approvals, or VICC evidence)

**Why it matters**: These are variants with strong biological driver potential but no established therapeutic relevance—prime candidates for drug sensitivity screening or functional validation.

---

## Cross-Histology Preclinical Models

When cell line models exist with the mutation but NOT in the queried tumor type:

```
Models with V600E exist but none in Cholangiocarcinoma — cross-histology testing possible
```

**Suggested studies:**
- Test in cholangiocarcinoma-derived organoids
- Compare drug response vs other histologies
- Generate isogenic model in cholangiocarcinoma background

---

## LLM Research Hypothesis Generation

Gap detection feeds into LLM-powered hypothesis generation:

```python
# In LLM output
{
    "research_hypotheses": [
        "Given the lack of functional data for JAK1 V657F despite its recurrence in T-ALL, isogenic knock-in models could determine whether this variant causes gain- or loss-of-function signaling.",
        "The absence of preclinical drug sensitivity data for this variant, combined with its structural similarity to JAK2 V617F, suggests testing JAK inhibitor panels in cell lines harboring this mutation."
    ]
}
```

Hypotheses are:
- Specific and testable (not vague)
- Connect multiple evidence elements (gap + existing data)
- Focus on biological mechanism, preclinical testing, co-mutation effects
- Avoid clinical treatment recommendations

---

## Example: Hotspot-Adjacent Rare Variant

```
$ mind KRAS G14D --tumor NSCLC
```

**Evidence Quality: limited**

**Well Characterized:**
- near hotspot codon 12 — structural hypothesis likely
- computational pathogenicity (CADD: 28)
- gene role (oncogene)

**Gaps Detected:**

| Gap | Severity | Description |
|-----|----------|-------------|
| Functional | SIGNIFICANT | Rare variant near known hotspot (codon 12) — functional characterization needed |
| Clinical | CRITICAL | No curated clinical evidence for KRAS G14D |
| Drug Response | SIGNIFICANT | No drug sensitivity/resistance data |

**Suggested Studies:**
- Compare to nearby hotspot KRAS codon 12
- Structural modeling to assess activation mechanism
- Functional assay (transformation, signaling)

**Research Priority: very_high** 🔥

This is a rare variant adjacent to the most common KRAS hotspot. Structural similarity to G12 mutations suggests potential activating function, but it lacks characterization.

---

## API Usage

```python
from oncomind.insight_builder.gap_detector import detect_evidence_gaps
from oncomind.models.gene_context import is_hotspot_variant, is_hotspot_adjacent

# After building evidence
gaps = detect_evidence_gaps(evidence)

# Overall assessment
print(gaps.overall_evidence_quality)  # "limited"
print(gaps.research_priority)          # "very_high"

# What's known vs unknown (simple strings)
print(gaps.well_characterized)         # ["Near Hotspot Codon 12", "Drug Response", ...]
print(gaps.poorly_characterized)       # ["clinical evidence", ...]

# Detailed well-characterized with basis, match levels, and tumor tracking
for aspect in gaps.well_characterized_detailed:
    print(f"{aspect.aspect}: {aspect.basis}")
    print(f"  Locus levels: {aspect.matches_on}")  # "3 variant, 0 codon, 2 gene"
    print(f"  Tumor match: {aspect.tumor_match}")  # "4 tumor, 1 other"

# Individual gaps with actionable recommendations
for gap in gaps.gaps:
    print(f"{gap.severity}: {gap.description}")
    print(f"  Suggested: {gap.suggested_studies}")
    print(f"  Data sources: {gap.addressable_with}")

# Helper methods
gaps.has_critical_gaps()                           # True/False
gaps.get_gaps_by_severity(GapSeverity.CRITICAL)    # Filter by severity
gaps.get_gaps_by_category(GapCategory.DISCORDANT)  # Filter by category
gaps.top_gaps(n=3)                                 # Get top N gaps by severity
gaps.to_summary()                                  # Human-readable text
gaps.to_dict_for_llm()                             # Optimized for LLM prompts

# Hotspot detection
is_hotspot_variant("BRAF", "V600E")                # True
is_hotspot_adjacent("BRAF", "V598E", window=5)     # (True, 600)
```

---

## Gap Categories Reference

| Category | Checks For | Sources Used | Weight |
|----------|------------|--------------|--------|
| `FUNCTIONAL` | Pathogenicity predictions, protein impact, hotspot context | MyVariant, VEP, DepMap, Hotspot DB | 3.0 |
| `CLINICAL` | Clinical assertions, FDA approvals | CIViC, FDA, CGI | 1.0 |
| `TUMOR_TYPE` | Tumor-specific evidence | CIViC, VICC, CGI (all tiers), FDA | 1.5 |
| `DRUG_RESPONSE` | Sensitivity/resistance data | CGI (FDA tier), VICC, FDA | 1.5 |
| `RESISTANCE` | Known resistance mechanisms | PubMed, CGI, CIViC, VICC | 2.0 |
| `PRECLINICAL` | Cell line models (tumor-specific) + CGI preclinical/early phase | DepMap, CGI | 2.5 |
| `PREVALENCE` | Mutation frequency | cBioPortal | 1.0 |
| `PROGNOSTIC` | Survival/outcome data | CIViC, Literature | 1.0 |
| `DISCORDANT` | Conflicting drug response, ClinVar conflicts | Cross-source comparison | 2.0 |
| `VALIDATION` | Strong oncogenic signal, limited therapeutic validation | DepMap essentiality + pathogenicity | 3.5 |

---

## Design Philosophy

1. **Gaps over facts** — Prioritize surfacing what's unknown
2. **Severity is context-aware** — Pathogenic signals and cancer genes increase severity
3. **Research-weighted scoring** — Biological gaps weighted higher than clinical gaps
4. **Hotspot context matters** — Rare variants near hotspots are research opportunities
5. **Conflicts are valuable** — Discordant evidence highlights unresolved questions
6. **Actionable recommendations** — Every gap includes suggested studies and data sources
7. **Honest uncertainty** — If we don't know, we say so explicitly
8. **Research-first** — Optimize for hypothesis generation, not clinical decisions

---

## Match Level and Tumor Match Tracking

Each `CharacterizedAspect` in `well_characterized_detailed` includes granular tracking fields for transparency:

### CharacterizedAspect Fields

```python
class CharacterizedAspect(BaseModel):
    aspect: str           # "drug response", "resistance mechanisms", etc.
    basis: str            # "2 CGI + 3 VICC + 1 FDA"
    category: GapCategory # For grouping (DRUG_RESPONSE, RESISTANCE, etc.)
    matches_on: str | None    # "3 variant, 0 codon, 2 gene"
    tumor_match: str | None   # "4 tumor, 1 other"
    cancer_mismatch: str | None  # "Melanoma" (if FDA approval is for different cancer)
```

### matches_on: Locus Level Tracking

Shows how precisely each evidence item matches the query variant:

| Level | Meaning | Example |
|-------|---------|---------|
| `variant` | Exact variant match | BRAF V600E → evidence for V600E specifically |
| `codon` | Same position, different AA | BRAF V600K → evidence for "V600 mutations" |
| `gene` | Gene-level only | BRAF V600E → evidence for "BRAF mutations" |

**Format**: Always shows all three levels, even when 0:
```
"3 variant, 0 codon, 2 gene"
```

### tumor_match: Tumor Type Matching

Shows how many evidence items match the queried tumor type:

| Field | Meaning |
|-------|---------|
| `X tumor` | Items matching the queried tumor type (substring match) |
| `Y other` | Items for different tumor types |

**Format**: Always shows tumor count (even 0), shows "other" only if > 0:
```
"4 tumor, 1 other"  # 4 match queried tumor, 1 for different tumor
"0 tumor, 5 other"  # No matches for queried tumor
"3 tumor"           # All 3 match queried tumor
```

**Tumor matching logic**: Uses substring matching (`tumor_lower in disease.lower()`), so "Lung" matches "Lung Cancer", "Lung adenocarcinoma", etc.

---

## Data Sources by Well-Characterized Row

### Drug Response
**Sources counted**: CGI biomarkers (FDA-approved tier only), VICC evidence, FDA approvals

| Source | Field Checked |
|--------|---------------|
| CGI | `cgi_biomarkers` |
| VICC | `vicc_evidence` |
| FDA | `fda_approvals` |

**Excludes**: `preclinical_biomarkers`, `early_phase_biomarkers` (tracked separately)

### Preclinical/Early Phase Biomarkers
**Sources counted**: CGI preclinical and early phase tiers

| Source | Field Checked |
|--------|---------------|
| Preclinical | `preclinical_biomarkers` |
| Early Phase | `early_phase_biomarkers` |

### Resistance Mechanisms
**Sources counted**: Literature, CGI, CIViC, VICC

| Source | Resistance Check |
|--------|-----------------|
| PubMed | `article.is_resistance_evidence()` |
| CGI | `biomarker.association == "Resistance"` |
| CIViC | `assertion.is_resistance` |
| VICC | `evidence.response_type == "RESISTANCE"` |

### Evidence In [Tumor]
**Sources counted**: All tumor-specific evidence across sources

| Source | Evidence Type |
|--------|---------------|
| CIViC Assertions | `civic_assertions` with matching disease |
| CIViC Evidence | `civic_evidence` with matching disease |
| FDA Approvals | `fda_approvals` with matching indication |
| VICC | `vicc_evidence` with matching disease |
| CGI (all tiers) | `cgi_biomarkers` + `preclinical_biomarkers` + `early_phase_biomarkers` with matching tumor_type |

**Note**: "Evidence In [Tumor]" includes ALL CGI tiers (FDA-approved, preclinical, and early phase) because all represent tumor-specific evidence, even if at different validation levels.

### Clinical Actionability
**Sources counted**: FDA approvals, CIViC assertions

Tracks `tumor_match` and `cancer_mismatch` to show when FDA approval exists but for a different cancer than queried.

---

## Example Output

For EGFR L858R in NSCLC:

```python
well_characterized_detailed = [
    CharacterizedAspect(
        aspect="Drug Response",
        basis="5 CGI + 8 VICC + 3 FDA",
        category=GapCategory.DRUG_RESPONSE,
        matches_on="12 variant, 2 codon, 2 gene",
        tumor_match="14 tumor, 2 other"
    ),
    CharacterizedAspect(
        aspect="Resistance Mechanisms",
        basis="2 PubMed + 3 CGI + 1 CIViC + 4 VICC",
        category=GapCategory.RESISTANCE,
        matches_on="8 variant, 0 codon, 2 gene",
        tumor_match="7 tumor, 3 other"
    ),
    CharacterizedAspect(
        aspect="Evidence In NSCLC",
        basis="3 CIViC Assertions, 2 FDA, 5 VICC, 2 CGI",
        category=GapCategory.TUMOR_TYPE,
        matches_on="10 variant, 1 codon, 1 gene",
        tumor_match="12 tumor"
    ),
]
```

---

## Code References

- **Gap detector**: [gap_detector.py](../src/oncomind/insight_builder/gap_detector.py)
- **EvidenceGaps model**: [evidence_gaps.py](../src/oncomind/models/evidence/evidence_gaps.py)
- **CharacterizedAspect model**: [evidence_gaps.py](../src/oncomind/models/evidence/evidence_gaps.py) — `CharacterizedAspect` class
- **TumorEvidenceMatch model**: [tumor_evidence.py](../src/oncomind/models/evidence/tumor_evidence.py) — aggregates tumor-specific evidence
- **Hotspot detection**: [gene_context.py](../src/oncomind/models/gene_context.py) — `is_hotspot_variant()`, `is_hotspot_adjacent()`
- **LLM integration**: [prompts.py](../src/oncomind/llm/prompts.py) — gaps and hypotheses in LLM context
- **LLMInsight model**: [llm_insight.py](../src/oncomind/models/llm_insight.py) — `research_hypotheses` field

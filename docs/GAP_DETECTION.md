# Evidence Gap Detection

## The Core Insight

Most variant annotation tools answer the question: *"What do we know about this variant?"*

OncoMind answers a different question: ***"What don't we know yet?"***

For well-characterized variants like BRAF V600E or EGFR L858R, the clinical answers are already in databases—FDA approvals, clinical guidelines, established therapies. The interesting research questions live in the gaps: variants with conflicting evidence, understudied tumor contexts, unknown resistance mechanisms, or missing preclinical validation.

**Evidence gap detection is OncoMind's core differentiator.** By systematically identifying what's missing, where sources conflict, and what research questions remain open, we help researchers prioritize investigations and identify opportunities for novel contributions.

---

## Gap Severity Levels

OncoMind uses a 6-level severity scale to classify evidence gaps, from most urgent to least:

| Level | Multiplier | Meaning | Icon |
|-------|-----------|---------|------|
| **CRITICAL** | 4.0 | No data for actionable/clinically relevant variant | 🔴 |
| **HIGH** | 3.0 | Conflicting data, or data exists but unreliable/outdated | 🟠 |
| **SIGNIFICANT** | 2.0 | Limited data, clear research opportunity | 🟡 |
| **MODERATE** | 1.5 | Some data but gaps in specific contexts (tumor type, extrapolation) | 🔵 |
| **MINOR** | 1.0 | Well-characterized overall, minor enrichment possible | ⚪ |
| **INFORMATIONAL** | 0.5 | Not a true gap, just noting a limitation | ℹ️ |

---

## Severity Assignment Logic

### CRITICAL (🔴)

| Gap Type | Condition | Description |
|----------|-----------|-------------|
| No clinical evidence | No CIViC assertions, no FDA approvals | "No curated clinical evidence for {gene} {variant}" |
| No tumor-specific evidence | Cancer gene + no clinical + pathogenic signal + no tumor data | "No evidence specific to {tumor_type}" |
| No literature | Cancer gene + literature searched + 0 papers | "No published literature found" |
| No therapeutic validation | Cancer gene + strong oncogenic signal + no CIViC/FDA/VICC | "Strong oncogenic signal but limited therapeutic validation" |

### HIGH (🟠)

| Gap Type | Condition | Description |
|----------|-----------|-------------|
| Cross-source drug conflict (variant-level) | Different sources disagree on drug response AT VARIANT LEVEL | "Conflicting drug response for {drug} at variant level: sensitive (Source1) vs resistant (Source2)" |
| ClinVar conflicts (variant-level) | Both pathogenic and benign submissions at variant level | "ClinVar has conflicting interpretations at variant level: both pathogenic and benign submissions for this exact variant" |
| ClinVar conflicts (codon-level) | Both pathogenic and benign submissions at codon level | "ClinVar has conflicting interpretations at codon level: both pathogenic and benign submissions for variants at this position" |
| ClinVar vs CIViC conflict (variant-level) | ClinVar benign at variant level + CIViC actionable at variant level | "ClinVar classifies as benign at variant level but {source} suggests clinical relevance for this exact variant" |

**Important notes for HIGH severity:**
- Gene-level ClinVar "conflicts" are NOT flagged — different variants in the same gene can legitimately have different pathogenicity classifications
- Cross-source drug conflicts only fire when BOTH sensitive AND resistant sources have variant-level evidence (gene/codon-level conflicts could be for different variants)
- The ClinVar vs CIViC conflict only fires when BOTH sources have variant-level evidence for the exact same variant

### SIGNIFICANT (🟡)

| Gap Type | Condition | Description |
|----------|-----------|-------------|
| **FDA vs VICC conflict (variant-level)** | FDA approves drug (sensitive) but VICC reports resistance, both at variant level | "FDA approves {drug} (sensitive) at variant level but VICC reports resistance at variant level — requires validation" |
| **No variant-specific drug approvals (gene-level)** | FDA therapies exist for gene but NO variant-level match | "No variant-specific drug approvals for {gene} {variant} (gene-level evidence only)" |
| **No variant-specific drug approvals (codon-level)** | FDA therapies exist for codon but not this exact variant | "No variant-specific drug approvals for {gene} {variant} (codon N evidence only)" |
| FDA approval in other cancers only | FDA exists but not tumor-matched or pan-cancer | "FDA-approved therapies exist only in other cancers" |
| Near-hotspot | Variant within N codons of known hotspot | "Variant near known hotspot — functional characterization needed" |
| No pathogenicity predictions | No AlphaMissense, CADD, PolyPhen2, SIFT | "No computational pathogenicity predictions" |
| No functional mechanism | No gene role + no DepMap essentiality | "Functional impact unknown" |
| No drug sensitivity data | No CGI + no VICC + no FDA drug data | "No drug sensitivity/resistance data" |
| Resistance not characterized | Has clinical/drug data but no resistance signals | "Resistance mechanisms not well characterized" |
| Prevalence unknown (clinical cancer gene) | Cancer gene + has clinical + gene in cBioPortal but variant not seen | "Prevalence unknown" |
| No tumor-specific evidence (cancer gene) | Cancer gene OR pathogenic signal, but NO tumor data at any level | "No evidence specific to {tumor_type}" |
| Preclinical models in wrong tumor | Models exist but none match queried tumor type | "Models exist but none in {tumor_type}" |
| No literature (non-cancer gene) | Non-cancer gene + literature searched + 0 papers | "No published literature found" |
| No therapeutic validation (non-cancer gene) | Non-cancer gene + strong oncogenic signal + no validation | "Strong oncogenic signal but limited therapeutic validation" |

**Important:** Gene/codon-level drug approvals are NOT marked as "well-characterized for clinical actionability" — only variant-level matches are. This ensures users clearly see when drug evidence is extrapolated.

### MODERATE (🔵)

| Gap Type | Condition | Description |
|----------|-----------|-------------|
| Drug response in other tumors | Drug data exists but only in non-matching tumor types | "Drug response data exists only in other cancers, not {tumor_type}" |
| Gene-level drug extrapolation | Drug data (CGI/VICC) exists for gene but not variant-specific | "Drug response data exists for {gene} but not specifically for {variant}" |
| Codon-level drug extrapolation | Drug data (CGI/VICC) exists for codon but not this exact variant | "Drug response data exists for codon {N} variants but not specifically for {variant}" |
| **Tumor evidence at gene-level only** | Tumor evidence exists but only at gene level (not variant) | "Evidence in {tumor_type} exists for {gene} but not specifically for {variant} (gene-level only)" |
| **Tumor evidence at codon-level only** | Tumor evidence exists but only at codon level (not variant) | "Evidence in {tumor_type} exists for {gene} codon {N} but not specifically for {variant}" |
| Preclinical biomarkers in other tumors | Preclinical data exists but only in non-matching tumors | "Preclinical biomarker data exists only in other cancers" |
| Resistance data in other tumors | Resistance signals exist but only in non-matching tumors | "Resistance data exists only in other cancers" |

**Note:** FDA/clinical drug approvals at gene/codon level are SIGNIFICANT (above), not MODERATE. MODERATE is for CGI/VICC drug response data that's extrapolated — lower clinical relevance than FDA approvals.

### MINOR (⚪)

| Gap Type | Condition | Description |
|----------|-----------|-------------|
| Prevalence unknown (cancer gene, no clinical) | Cancer gene but no clinical evidence | "Prevalence unknown" |
| No cell line models | Has drug/clinical data but no DepMap cell lines with mutation | "No cell line models identified" or "DepMap has N cell lines but none with {variant}" |
| No tumor-specific evidence (non-cancer gene) | Non-cancer gene + no pathogenic signal + no tumor data | "No evidence specific to {tumor_type}" |

### INFORMATIONAL (ℹ️)

| Gap Type | Condition | Description |
|----------|-----------|-------------|
| Prevalence unknown (non-cancer gene) | Non-cancer gene + gene in cBioPortal but variant not seen | "Prevalence unknown" |
| No active clinical trials | Has clinical/drug data but no trials found | "No active clinical trials" |
| Literature not in databases | Literature mentions drug signals not curated in CIViC/CGI/VICC | "Literature mentions drug signals not in curated databases" |

---

## Decision Tree Summary

```
Is there NO data at all for a key area?
├─ Yes → Is it a cancer gene or actionable?
│        ├─ Yes → CRITICAL
│        └─ No  → SIGNIFICANT
└─ No  → Is the data CONFLICTING?
         ├─ Yes → HIGH
         └─ No  → Is the data LIMITED (needs more research)?
                  ├─ Yes → SIGNIFICANT
                  └─ No  → Is data in WRONG CONTEXT (other tumor, gene-level)?
                           ├─ Yes → MODERATE
                           └─ No  → Is it a MINOR enrichment opportunity?
                                    ├─ Yes → MINOR
                                    └─ No  → INFORMATIONAL (just noting)
```

---

## Gap Categories

Each gap has a category (orthogonal to severity):

| Category | Description | Weight |
|----------|-------------|--------|
| `VALIDATION` | Strong oncogenic signal but lacks therapeutic validation | 3.5 |
| `FUNCTIONAL` | Mechanism unknown | 3.0 |
| `PRECLINICAL` | No cell line/model data | 2.5 |
| `RESISTANCE` | Resistance mechanisms unknown | 2.0 |
| `DISCORDANT` | Conflicting evidence between sources | 2.0 |
| `DRUG_RESPONSE` | No drug sensitivity data | 1.5 |
| `TUMOR_TYPE` | Not studied in this tumor type | 1.5 |
| `CLINICAL` | No clinical trials/outcomes | 1.0 |
| `PREVALENCE` | Frequency unknown | 1.0 |
| `PROGNOSTIC` | Survival impact unknown | 1.0 |

---

## How Severity Affects Scoring

### Overall Evidence Quality

Quality is computed in two steps:

**Step 1: Base quality from net score**

`net_score = gap_score - (well_characterized_count × 1.5)`

Where `gap_score = Σ (category_weight × severity_multiplier)` for each gap.

| Net Score | Base Quality |
|-----------|--------------|
| ≥ 12.0 | minimal |
| ≥ 6.0 | limited |
| ≥ 0.0 | moderate |
| < 0.0 | comprehensive |

**Step 2: Apply floor caps (severity-based limits)**

High-severity gaps prevent "comprehensive" classification regardless of well-characterized count:

| Gap Severity Present | Maximum Quality |
|---------------------|-----------------|
| Any CRITICAL gap | limited (cannot be comprehensive or moderate) |
| Any HIGH gap | moderate (cannot be comprehensive) |
| Any SIGNIFICANT gap | moderate (cannot be comprehensive) |
| Only MODERATE/MINOR/INFORMATIONAL | comprehensive (no cap) |

**Rationale:** A variant with SIGNIFICANT knowledge gaps (e.g., "No variant-specific drug approvals") should not be classified as "comprehensive" even if it has many well-characterized aspects. The severity caps ensure quality ratings honestly reflect the presence of important gaps.

### Research Priority

| Priority | Condition | Icon |
|----------|-----------|------|
| `very_high` | Strong oncogenic signal + biological gaps, OR hotspot-adjacent + pathogenic | 🔥 |
| `high` | Cancer gene with CRITICAL gaps, OR hotspot-adjacent cancer gene, OR any HIGH gaps | 🔴 |
| `medium` | Any CRITICAL gaps, OR cancer gene with SIGNIFICANT gaps | 🟡 |
| `low` | Everything else (comprehensive quality with no significant gaps) | 🟢 |

**Biological gaps** = gaps in VALIDATION, FUNCTIONAL, or PRECLINICAL categories.

---

## Locus Match Levels

Evidence is tracked at three levels of specificity:

| Level | Meaning | Example |
|-------|---------|---------|
| `variant` | Exact variant match | BRAF V600E → evidence specifically for V600E |
| `codon` | Same position, different AA | BRAF V600K → evidence for "V600 mutations" |
| `gene` | Gene-level only | BRAF V600E → evidence for "BRAF mutations" |

**Why this matters:**
- Cross-source conflicts are only flagged at **variant level** — gene/codon conflicts could be for different variants
- Drug approvals at gene/codon level are flagged as **SIGNIFICANT gaps** (extrapolation needed)
- Only **variant-level** evidence qualifies as "well-characterized" for clinical actionability

---

## Discordant Evidence Detection

OncoMind detects conflicts between data sources with locus-level awareness:

### Drug Response Conflicts (HIGH)

Only flagged when **both** sensitive and resistant sources have **variant-level** evidence:

```
Conflicting drug response for imatinib at variant level:
sensitive (CIViC) vs resistant (CGI)
```

Sources checked: CIViC, CGI, VICC MetaKB, FDA

### FDA vs VICC Conflicts (SIGNIFICANT)

When FDA approves a drug but VICC reports resistance, both at variant level:

```
FDA approves imatinib (sensitive) at variant level but VICC reports resistance
at variant level — requires validation
```

VICC data is less reliable than FDA, so this is SIGNIFICANT (not HIGH).

### ClinVar Internal Conflicts (HIGH)

When both pathogenic and benign submissions exist at the **same locus level**:

```
ClinVar has conflicting interpretations at variant level: both pathogenic
and benign submissions for this exact variant
```

**Note:** Gene-level ClinVar "conflicts" are NOT flagged — different variants in the same gene can legitimately have different pathogenicity.

### ClinVar vs CIViC Conflicts (HIGH)

When ClinVar classifies as benign but CIViC suggests actionability, **both at variant level**:

```
ClinVar classifies as benign at variant level but CIViC assertion suggests
clinical relevance for this exact variant
```

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

## Validation Gap (Oncogenicity Potential)

A special gap category for high-potential research targets:

**Triggers when:**
1. Variant has pathogenic signal (AlphaMissense P, CADD ≥20, ClinVar pathogenic, or truncating)
2. Gene is essential in cancer cells (DepMap)
3. BUT no therapeutic validation (no CIViC assertions, FDA approvals, or VICC evidence)

**Why it matters**: These are variants with strong biological driver potential but no established therapeutic relevance—prime candidates for drug sensitivity screening or functional validation.

---

## Other Cancer Evidence Gaps

When evidence exists for a variant but only in **other cancer types** (not the queried tumor and not pan-cancer):

### When This Gap Fires

The gap is created when **ALL** of the following are true:
1. Evidence exists for the variant (FDA approvals, drug response data, preclinical biomarkers, or resistance data)
2. Zero evidence matches the queried tumor type
3. Zero pan-cancer evidence exists

### Example: Gap Fires

```
Query: PIK3CA H1047R in Cholangiocarcinoma

Evidence found:
- FDA approval for Breast Cancer
- FDA approval for Ovarian Cancer
- CGI biomarker for Endometrial Cancer

Result: SIGNIFICANT gap created
"FDA-approved therapies for PIK3CA H1047R exist only in other cancers
(Breast Cancer, Ovarian Cancer, Endometrial Cancer), not Cholangiocarcinoma"
```

### Example: Gap Does NOT Fire

```
Query: BRAF V600E in Thyroid Cancer

Evidence found:
- FDA approval for Melanoma (other)
- FDA approval for Solid Tumor (pan-cancer) ← Pan-cancer exists!
- FDA approval for NSCLC (other)

Result: No "other cancer" gap created
Because pan-cancer evidence exists, the gap doesn't fire.
```

---

## Match Level and Tumor Match Tracking

Each `CharacterizedAspect` in `well_characterized_detailed` includes granular tracking:

```python
class CharacterizedAspect(BaseModel):
    aspect: str           # "drug response", "resistance mechanisms", etc.
    basis: str            # "2 CGI + 3 VICC + 1 FDA"
    category: GapCategory # For grouping (DRUG_RESPONSE, RESISTANCE, etc.)
    matches_on: str | None    # "3 variant, 0 codon, 2 gene"
    tumor_match: str | None   # "4 tumor, 1 other"
```

### matches_on: Locus Level Tracking

Shows how precisely each evidence item matches the query variant:

**Format**: Always shows all three levels, even when 0:
```
"3 variant, 0 codon, 2 gene"
```

### tumor_match: Tumor Type Matching

Shows how many evidence items match the queried tumor type:

**Format**: Always shows tumor count (even 0), shows "other" only if > 0:
```
"4 tumor, 1 other"  # 4 match queried tumor, 1 for different tumor
"0 tumor, 5 other"  # No matches for queried tumor
"3 tumor"           # All 3 match queried tumor
```

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
| Functional | SIGNIFICANT | Variant near known hotspot (codon 12) — functional characterization needed |
| Clinical | CRITICAL | No curated clinical evidence for KRAS G14D |
| Drug Response | SIGNIFICANT | No drug sensitivity/resistance data |

**Suggested Studies:**
- Compare to nearby hotspot KRAS codon 12
- Structural modeling to assess activation mechanism
- Functional assay (transformation, signaling)

**Research Priority: very_high** 🔥

---

## Design Philosophy

1. **Gaps over facts** — Prioritize surfacing what's unknown
2. **Severity is context-aware** — Pathogenic signals and cancer genes increase severity
3. **Locus-level awareness** — Only flag conflicts when evidence is at the same specificity level
4. **Research-weighted scoring** — Biological gaps weighted higher than clinical gaps
5. **Hotspot context matters** — Rare variants near hotspots are research opportunities
6. **Conflicts are valuable** — Discordant evidence highlights unresolved questions
7. **Actionable recommendations** — Every gap includes suggested studies and data sources
8. **Honest uncertainty** — If we don't know, we say so explicitly
9. **Research-first** — Optimize for hypothesis generation, not clinical decisions

---

## Code References

- **Gap detector**: [gap_detector.py](../src/oncomind/insight_builder/gap_detector.py)
- **EvidenceGaps model**: [evidence_gaps.py](../src/oncomind/models/extracted/evidence_gaps.py)
- **CharacterizedAspect model**: [evidence_gaps.py](../src/oncomind/models/extracted/evidence_gaps.py) — `CharacterizedAspect` class
- **Hotspot detection**: [gene_context.py](../src/oncomind/models/gene_context.py) — `is_hotspot_variant()`, `is_hotspot_adjacent()`
- **LLM integration**: [prompts.py](../src/oncomind/llm/prompts.py) — gaps and hypotheses in LLM context

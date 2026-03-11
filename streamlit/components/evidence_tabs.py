"""Evidence tab rendering components for OncoMind Streamlit UI."""

import re
import pandas as pd
import streamlit as st

from oncomind.config.constants import (
    UI_MAX_CIVIC_EVIDENCE_ROWS,
    CADD_DELETERIOUS_THRESHOLD,
    GNOMAD_RARE_THRESHOLD,
    GNOMAD_UNCOMMON_THRESHOLD,
    GNOMAD_LOW_FREQ_THRESHOLD,
)
from oncomind.api.fda_drugs import get_best_fda_url, get_cached_fda_label
from .utils import scrollable_table


# ==============================================
# SHARED LEGEND COMPONENTS
# ==============================================

LOCUS_TUMOR_LEGEND = """<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 85px;'>
<b>Locus Match:</b><br/>
🎯 Variant (exact locus match)<br/>
📍 Codon (other variants in this codon)<br/>
🧬 Gene (other variants in this gene)<br/><br/>
<b>Tumor Match:</b><br/>
✅ Yes (match on specified tumor)<br/>
🔸 Other (match on other tumor)<br/>
🌐 Pan-cancer
</div>"""

LOCUS_TUMOR_LEGEND_TALL = """<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 180px;'>
<b>Locus Match:</b><br/>
🎯 Variant (exact locus match)<br/>
📍 Codon (other variants in this codon)<br/>
🧬 Gene (other variants in this gene)<br/><br/>
<b>Tumor Match:</b><br/>
✅ Yes (match on specified tumor)<br/>
🔸 Other (match on other tumor)<br/>
🌐 Pan-cancer
</div>"""

TRIALS_LEGEND = """<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 85px;'>
<b>Locus Match:</b><br/>
🎯 Variant (exact locus match)<br/>
📌 Broad (match on related variants)<br/>
🧬 Gene (other variants in this gene)<br/><br/>
<b>Tumor Match:</b><br/>
✅ Yes (match on specified tumor)<br/>
🔸 Other (match on other tumor)<br/>
🌐 Pan-cancer
</div>"""

CGI_LEGEND = """<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 85px;'>
<b>Locus Match:</b><br/>
🎯 Variant (exact locus match)<br/>
📍 Codon (other variants in this codon)<br/>
🧬 Gene (other variants in this gene)<br/><br/>
<b>Association:</b><br/>
<span style="color: #28a745">Responsive</span> - drug works<br/>
<span style="color: #dc3545">Resistant</span> - drug doesn't work
</div>"""


# ==============================================
# FUNCTIONAL SCORES TAB
# ==============================================

def render_functional_tab(
    annotations: dict,
    hgvs: dict,
    gene_display: str,
    get_pdb_image_url,
    get_pdb_page_url,
):
    """Render the Functional Scores tab."""
    # Two-column layout: scores on left, protein structure on right
    pdb_info = get_pdb_image_url(gene_display)
    func_col1, func_col2 = st.columns([2, 1]) if pdb_info else (st.container(), None)

    with func_col1:
        rows = []
        if annotations.get('alphamissense_score') is not None:
            pred = annotations.get('alphamissense_prediction') or '-'
            rows.append(f"| AlphaMissense | {annotations['alphamissense_score']:.3f} | {pred} |")
        if annotations.get('cadd_score') is not None:
            rows.append(f"| CADD | {annotations['cadd_score']:.1f} | {'Deleterious' if annotations['cadd_score'] > CADD_DELETERIOUS_THRESHOLD else 'Benign'} |")
        if annotations.get('polyphen2_prediction'):
            rows.append(f"| PolyPhen2 | - | {annotations['polyphen2_prediction']} |")
        if annotations.get('gnomad_exome_af') is not None:
            af = annotations['gnomad_exome_af']
            freq = f"{af:.2e}" if af < GNOMAD_UNCOMMON_THRESHOLD else f"{af:.4f}"
            if af < GNOMAD_RARE_THRESHOLD:
                af_label = 'Rare'
            elif af < GNOMAD_UNCOMMON_THRESHOLD:
                af_label = 'Uncommon'
            elif af < GNOMAD_LOW_FREQ_THRESHOLD:
                af_label = 'Low-frequency'
            else:
                af_label = 'Common'
            rows.append(f"| gnomAD AF | {freq} | {af_label} |")
        if annotations.get('snpeff_effect'):
            rows.append(f"| SnpEff | - | {annotations['snpeff_effect']} |")
        if rows:
            st.markdown("| Score | Value | Prediction |\n|-------|-------|------------|" + "\n" + "\n".join(rows))
        else:
            st.info("No functional scores available")

        hgvs_genomic = hgvs.get('genomic')
        if hgvs_genomic:
            import urllib.parse
            encoded_id = urllib.parse.quote(hgvs_genomic, safe='')
            myvariant_url = f"https://myvariant.info/v1/variant/{encoded_id}"
            st.markdown(f"*Source: [MyVariant.info]({myvariant_url})*")

    if pdb_info and func_col2:
        with func_col2:
            pdb_page_url = get_pdb_page_url(gene_display)
            # Image on left, legend/link on right
            img_col, legend_col = st.columns([1, 1])
            with img_col:
                st.image(pdb_info['url'], width=220)
            with legend_col:
                st.markdown(f"**{gene_display}**<br>{pdb_info['description']}<br>[PDB: {pdb_info['pdb_id']}]({pdb_page_url})", unsafe_allow_html=True)


# ==============================================
# CIViC TAB
# ==============================================

def render_civic_tab(
    civic_assertions: list,
    civic_evidence: list,
    tumor_display: str | None = None,
):
    """Render the CIViC Evidence tab."""
    # Use columns: tables on left, legend on right
    table_col, legend_col = st.columns([4, 1])

    # Render legend first, with padding to align with Level A table
    legend_col.markdown(LOCUS_TUMOR_LEGEND_TALL, unsafe_allow_html=True)

    with table_col:
        st.caption("Curated clinical evidence from CIViC (Clinical Interpretation of Variants in Cancer)")
        # Build match summary for CIViC
        with st.expander("📖 CIViC Evidence Level Guide", expanded=False):
            st.markdown("""
**CIViC Evidence Levels:**
- **Level A (FDA/Guideline)**: FDA-approved therapy or included in professional guidelines
- **Level B (Clinical)**: Well-powered clinical studies with expert consensus
- **Level C (Case Study)**: Case studies, small series, or early-phase trials
- **Level D (Preclinical)**: Preclinical data or inferential evidence
- **Level E (Inferential)**: Indirect evidence from related variants/pathways

**AMP/ASCO/CAP Tiers:**
- **Tier I**: Strong clinical significance (FDA-approved or guideline-recommended)
- **Tier II**: Potential clinical significance (clinical trials, case studies)
- **Tier III**: Unknown clinical significance
- **Tier IV**: Benign or likely benign variants
""")
        # Group evidence by level FIRST so we can render Level A before Assertions
        level_groups = {
            "A": [],  # FDA/Guideline
            "B": [],  # Clinical
            "C": [],  # Case Study
            "D": [],  # Preclinical
            "E": [],  # Inferential
        }
        if civic_evidence:
            for e in civic_evidence:
                level = e.get('evidence_level', '')
                if level in level_groups:
                    level_groups[level].append(e)
                else:
                    level_groups.setdefault('other', []).append(e)

        # Helper to render a CIViC evidence table
        def render_civic_table(evidence_items: list, show_fda_cols: bool = False):
            if show_fda_cols:
                # Level A FDA table: ID, Locus Match, Tumor Match, Drugs, Significance, Regulatory Status, Indication, Links
                rows = ["| ID | Locus Match | Tumor Match | Drugs | Significance | Regulatory Status | Indication | Links |",
                        "|----|-------------|-------------|-------|--------------|-------------------|------------|-------|"]
            else:
                rows = ["| ID | Locus Match   | Tumor Match | Drugs | Significance | Disease | Type |",
                        "|----|---------------|-------------|-------|--------------|---------|------|"]
            for e in evidence_items[:UI_MAX_CIVIC_EVIDENCE_ROWS]:
                drugs_list = e.get('drugs', [])
                drugs_str = ", ".join(drugs_list) or "N/A"
                drugs_str = drugs_str[:25] if len(drugs_str) > 25 else drugs_str
                eid = e.get('eid') or ''
                url = e.get('civic_url', '')
                id_link = f"[{eid}]({url})" if url else eid
                disease = (e.get('disease', '') or '')[:20]
                sig_raw = e.get('clinical_significance', 'Unknown') or 'Unknown'
                # Abbreviate long significance values
                sig_map = {
                    'SENSITIVITY/RESPONSE': 'Sensitive',
                    'SENSITIVITYRESPONSE': 'Sensitive',
                    'RESISTANCE': 'Resistant',
                    'REDUCED SENSITIVITY': 'Reduced',
                    'ADVERSE RESPONSE': 'Adverse',
                }
                sig_abbrev = sig_map.get(sig_raw.upper(), sig_raw.title())
                # Add "NOT" prefix if direction is DOES_NOT_SUPPORT (except for N/A)
                direction = e.get('evidence_direction', '') or ''
                if direction.upper() == 'DOES_NOT_SUPPORT' and sig_raw.upper() != 'N/A':
                    sig = f"NOT {sig_abbrev}"
                else:
                    sig = sig_abbrev
                # Match level indicator
                match = e.get('locus_match', '')
                match_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "")
                # Tumor match indicator
                tumor_match = e.get('tumor_match')
                tumor_match_cell = "✅ Yes" if tumor_match else "🔸 Other"

                if show_fda_cols:
                    # Get FDA label info for first drug
                    first_drug = drugs_list[0] if drugs_list else None
                    fda_info = get_cached_fda_label(first_drug) if first_drug else None

                    # Regulatory Status
                    if fda_info and fda_info.approval_date:
                        reg_status = f"FDA-approved ({fda_info.approval_date})"
                    else:
                        reg_status = "FDA-approved"

                    # Indication - use full indications_and_usage
                    if fda_info and fda_info.indications_and_usage:
                        indication = fda_info.indications_and_usage
                    elif fda_info and fda_info.biomarker_text_exact:
                        indication = fda_info.biomarker_text_exact
                    else:
                        indication = disease or "-"

                    # Links: FDA + PubMed
                    links = []
                    if fda_info and fda_info.drugs_at_fda_url:
                        links.append(f"[FDA]({fda_info.drugs_at_fda_url})")
                    elif first_drug:
                        fda_url = get_best_fda_url(first_drug)
                        if fda_url:
                            links.append(f"[FDA]({fda_url})")
                    if fda_info and fda_info.pubmed_url:
                        links.append(f"[PubMed]({fda_info.pubmed_url})")
                    links_str = " · ".join(links) if links else "-"

                    rows.append(f"| {id_link} | {match_display} | {tumor_match_cell} | {drugs_str} | {sig} | {reg_status} | {indication} | {links_str} |")
                else:
                    etype = e.get('evidence_type', '')
                    rows.append(f"| {id_link} | {match_display} | {tumor_match_cell} | {drugs_str} | {sig} | {disease} | {etype} |")
            scrollable_table("\n".join(rows))

        # Track sections rendered for dividers
        sections_rendered = 0

        # 1. Render Level A FIRST (highest clinical significance)
        if level_groups["A"]:
            st.markdown(f"**✅ Level A - FDA/Guideline ({len(level_groups['A'])})**")
            st.caption("FDA-approved therapies - see FDA tab for label details")
            render_civic_table(level_groups["A"])
            sections_rendered += 1

        # 2. Render Curated Assertions second
        if civic_assertions:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown("**Curated Assertions:**")
            rows = ["| ID | Locus Match | Tumor Match | Therapies | Significance | Disease | AMP Level |",
                    "|-----|-------------|-------------|-----------|--------------|---------|-----------|"]
            for a in civic_assertions:
                therapies_str = ", ".join(a.get('therapies', [])) or "N/A"
                aid = a.get('aid') or a.get('id', '')
                url = a.get('civic_url', '')
                id_link = f"[{aid}]({url})" if url else aid
                disease = (a.get('disease', '') or '')[:35]
                sig_raw = a.get('significance', 'Unknown') or 'Unknown'
                # Abbreviate long significance values
                sig_map = {
                    'SENSITIVITY/RESPONSE': 'Sensitive',
                    'SENSITIVITYRESPONSE': 'Sensitive',
                    'RESISTANCE': 'Resistant',
                    'REDUCED SENSITIVITY': 'Reduced',
                    'ADVERSE RESPONSE': 'Adverse',
                }
                sig = sig_map.get(sig_raw.upper(), sig_raw.title())
                amp = a.get('amp_level', '')
                # Match level indicator with label
                match = a.get('locus_match', '')
                match_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "")
                # Tumor match indicator
                tumor_match = a.get('tumor_match')
                tumor_match_cell = "✅ Yes" if tumor_match else "🔸 Other"
                rows.append(f"| {id_link} | {match_display} | {tumor_match_cell} | {therapies_str} | {sig} | {disease} | {amp} |")
            scrollable_table("\n".join(rows))
            sections_rendered += 1

        # 3. Render remaining evidence levels (B, C, D, E)
        if level_groups["B"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**🔬 Level B - Clinical ({len(level_groups['B'])})**")
            st.caption("Well-powered clinical studies with expert consensus")
            render_civic_table(level_groups["B"])
            sections_rendered += 1

        if level_groups["C"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**📋 Level C - Case Study ({len(level_groups['C'])})**")
            st.caption("Case studies, small series, or early-phase trials")
            render_civic_table(level_groups["C"])
            sections_rendered += 1

        if level_groups["D"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**🧪 Level D - Preclinical ({len(level_groups['D'])})**")
            st.caption("Preclinical data or inferential evidence")
            render_civic_table(level_groups["D"])
            sections_rendered += 1

        if level_groups["E"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**💭 Level E - Inferential ({len(level_groups['E'])})**")
            st.caption("Indirect evidence from related variants/pathways")
            render_civic_table(level_groups["E"])


# ==============================================
# VICC TAB
# ==============================================

def render_vicc_tab(vicc: list, tumor_display: str | None = None):
    """Render the VICC MetaKB Evidence tab."""
    # Use columns: tables on left, legend on right
    table_col, legend_col = st.columns([4, 1])

    # Render legend first, with padding to align with first data row
    legend_col.markdown(LOCUS_TUMOR_LEGEND, unsafe_allow_html=True)

    with table_col:
        st.caption("Aggregated evidence from VICC MetaKB (CIViC, OncoKB, MOAlmanac, etc.)")
        if tumor_display:
            st.caption(f"Filtered to **{tumor_display}**")
        with st.expander("📖 VICC Evidence Level Guide", expanded=False):
            st.markdown("""
**VICC Evidence Levels (harmonized across sources):**
- **1/A (FDA/Guideline)**: FDA-approved or standard of care therapy
- **2/B (Clinical)**: Clinical trial evidence or expert consensus
- **3/C (Case Study)**: Case reports or limited clinical evidence
- **4/D (Preclinical)**: Preclinical or computational evidence
- **R1/R2 (Resistance)**: Resistance evidence (strong/emerging)
""")

        # Group VICC evidence by level
        level_groups = {
            "fda": [],      # 1, A, amp_1, tier_1
            "clinical": [], # 2, B, amp_2, tier_2
            "case": [],     # 3, C
            "preclin": [],  # 4, D
            "resist": [],   # R1, R2
            "other": [],
        }

        for v in vicc:
            level = str(v.get('evidence_level', '')).upper()
            # Extract tier from combined formats like "2C", "3A", etc.
            tier = level[0] if level and level[0].isdigit() else None
            if level in ['1', 'A', 'AMP_1', 'TIER_1', 'TIER I', 'TIER_I'] or tier == '1':
                level_groups["fda"].append(v)
            elif level in ['2', 'B', 'AMP_2', 'TIER_2', 'TIER II', 'TIER_II'] or tier == '2':
                level_groups["clinical"].append(v)
            elif level in ['3', 'C', 'AMP_3', 'TIER_3', 'TIER III', 'TIER_III'] or tier == '3':
                level_groups["case"].append(v)
            elif level in ['4', 'D', 'AMP_4', 'TIER_4', 'TIER IV', 'TIER_IV'] or tier == '4':
                level_groups["preclin"].append(v)
            elif level.startswith('R'):
                level_groups["resist"].append(v)
            else:
                level_groups["other"].append(v)

        # Helper to render a VICC table
        def render_vicc_table(evidence_items: list, show_fda_cols: bool = False):
            if show_fda_cols:
                rows = ["| Source | Locus Match | Tumor Match | Drugs | Response | Regulatory Status | Indication | Links |",
                        "|--------|-------------|-------------|-------|----------|-------------------|------------|-------|"]
            else:
                rows = ["| Source | Locus Match | Tumor Match | Drugs | Response | Disease |",
                        "|--------|-------------|-------------|-------|----------|---------|"]
            for v in evidence_items:
                source = (v.get('source') or 'vicc').upper()
                pub_url = v.get('publication_url')
                if isinstance(pub_url, list) and pub_url:
                    pub_url = pub_url[0]
                # Handle gated/missing URLs - NCCN requires login, some sources have no URL
                # Provide DailyMed fallback for the first drug if available
                drugs_list = v.get('drugs', [])
                if pub_url and 'nccn.org' in pub_url:
                    # NCCN URLs require professional membership - use DailyMed instead
                    first_drug = drugs_list[0] if drugs_list else None
                    if first_drug:
                        dailymed_url = f"https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query={first_drug}"
                        source_link = f"[{source}]({dailymed_url})"
                    else:
                        source_link = f"{source} (NCCN)"
                elif pub_url:
                    source_link = f"[{source}]({pub_url})"
                else:
                    # No URL - try DailyMed fallback
                    first_drug = drugs_list[0] if drugs_list else None
                    if first_drug:
                        dailymed_url = f"https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query={first_drug}"
                        source_link = f"[{source}]({dailymed_url})"
                    else:
                        source_link = source
                locus_match = v.get('locus_match', '')
                match_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(locus_match, "")
                tumor_match = v.get('tumor_match')
                tumor_match_cell = "✅ Yes" if tumor_match else "🔸 Other"
                drugs = ", ".join(drugs_list) or "N/A"
                drugs = drugs[:30] if len(drugs) > 30 else drugs
                # Handle response_type
                response = v.get('response_type', 'Unknown')
                if response and re.match(r'^[1234][A-D]?$', response.upper()):
                    # This is an AMP tier code, not a response type
                    desc = (v.get('description') or '').lower()
                    if 'sensitivity' in desc or 'sensitive' in desc or 'response' in desc:
                        response = 'sensitive'
                    elif 'resistance' in desc or 'resistant' in desc:
                        response = 'resistant'
                    else:
                        response = '-'
                disease = (v.get('disease', '') or '')[:25]

                if show_fda_cols:
                    first_drug = drugs_list[0] if drugs_list else None
                    fda_info = get_cached_fda_label(first_drug) if first_drug else None

                    if fda_info and fda_info.approval_date:
                        reg_status = f"FDA-approved ({fda_info.approval_date})"
                    else:
                        reg_status = "FDA-approved"

                    if fda_info and fda_info.indications_and_usage:
                        indication = fda_info.indications_and_usage[:60]
                        if len(fda_info.indications_and_usage) > 60:
                            indication += "..."
                    else:
                        indication = disease or "-"

                    links = []
                    if fda_info and fda_info.drugs_at_fda_url:
                        links.append(f"[FDA]({fda_info.drugs_at_fda_url})")
                    elif first_drug:
                        fda_url = get_best_fda_url(first_drug)
                        if fda_url:
                            links.append(f"[FDA]({fda_url})")
                    if fda_info and fda_info.pubmed_url:
                        links.append(f"[PubMed]({fda_info.pubmed_url})")
                    links_str = " · ".join(links) if links else "-"

                    rows.append(f"| {source_link} | {match_display} | {tumor_match_cell} | {drugs} | {response} | {reg_status} | {indication} | {links_str} |")
                else:
                    rows.append(f"| {source_link} | {match_display} | {tumor_match_cell} | {drugs} | {response} | {disease} |")
            scrollable_table("\n".join(rows))

        # Render each level group as a separate section
        sections_rendered = 0

        if level_groups["fda"]:
            st.markdown(f"**✅ Level 1/A - FDA/Guideline ({len(level_groups['fda'])})**")
            st.caption("FDA-approved therapies - see FDA tab for label details")
            render_vicc_table(level_groups["fda"])
            sections_rendered += 1

        if level_groups["clinical"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**🔬 Level 2/B - Clinical ({len(level_groups['clinical'])})**")
            st.caption("Clinical trial evidence or expert consensus")
            render_vicc_table(level_groups["clinical"])
            sections_rendered += 1

        if level_groups["case"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**📋 Level 3/C - Case Study ({len(level_groups['case'])})**")
            st.caption("Case reports or limited clinical evidence")
            render_vicc_table(level_groups["case"])
            sections_rendered += 1

        if level_groups["preclin"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**🧪 Level 4/D - Preclinical ({len(level_groups['preclin'])})**")
            st.caption("Preclinical or computational evidence")
            render_vicc_table(level_groups["preclin"])
            sections_rendered += 1

        if level_groups["resist"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**🔴 Resistance Evidence ({len(level_groups['resist'])})**")
            st.caption("Drug resistance associations (R1=strong, R2=emerging)")
            render_vicc_table(level_groups["resist"])
            sections_rendered += 1

        if level_groups["other"]:
            if sections_rendered > 0:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown(f"**Other ({len(level_groups['other'])})**")
            render_vicc_table(level_groups["other"])


# ==============================================
# CGI TAB
# ==============================================

def render_cgi_tab(
    cgi_biomarkers: list,
    early_phase: list,
    preclinical: list,
):
    """Render the CGI Biomarkers tab."""
    # Use columns: tables on left, legend on right
    table_col, legend_col = st.columns([4, 1])

    # Render legend
    legend_col.markdown(CGI_LEGEND, unsafe_allow_html=True)

    with table_col:
        # Sort FDA-approved: Responsive first, Resistant last
        def assoc_sort_key(b):
            assoc = (b.get('association') or '').lower()
            if 'respons' in assoc or 'sensitiv' in assoc:
                return 0
            elif 'resist' in assoc:
                return 2
            return 1

        # Helper to render a CGI table
        def render_cgi_table(biomarkers: list, show_level: bool = False, show_fda_cols: bool = False):
            if show_fda_cols:
                rows = ["| Drug | Locus Match | Tumor Match | Association | Regulatory Status | Indication | Links |",
                        "|------|-------------|-------------|-------------|-------------------|------------|-------|"]
            elif show_level:
                rows = ["| Drug | Locus Match | Association | Tumor Type | Level |",
                        "|------|-------------|-------------|------------|-------|"]
            else:
                rows = ["| Drug | Locus Match | Association | Tumor Type |",
                        "|------|-------------|-------------|------------|"]
            for b in biomarkers:
                drug = b.get('drug', 'Unknown')
                # Link to DailyMed for specific drug names, skip link for drug classes
                drug_words = drug.lower().split()
                is_drug_class = any(w in drug_words for w in ['inhibitor', 'inhibitors', 'agonist', 'antagonist', 'blocker'])
                if is_drug_class:
                    drug_link = drug
                else:
                    dailymed_url = f"https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query={drug}"
                    drug_link = f"[{drug}]({dailymed_url})"
                # Locus match
                match = b.get('locus_match', '')
                locus_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "-")
                # Tumor match
                tumor_match = b.get('cancer_type_match')
                tumor_match_cell = "✅ Yes" if tumor_match else "🔸 Other"
                # Association with color
                assoc_raw = b.get('association', 'Unknown')
                if 'respons' in assoc_raw.lower() or 'sensitiv' in assoc_raw.lower():
                    association = f'<span style="color: #28a745">{assoc_raw}</span>'
                elif 'resist' in assoc_raw.lower():
                    association = f'<span style="color: #dc3545">{assoc_raw}</span>'
                else:
                    association = assoc_raw
                tumor = (b.get('tumor_type', '') or '')[:25]

                if show_fda_cols:
                    fda_info = get_cached_fda_label(drug) if drug else None

                    if fda_info and fda_info.approval_date:
                        reg_status = f"FDA-approved ({fda_info.approval_date})"
                    else:
                        reg_status = "FDA-approved"

                    if fda_info and fda_info.indications_and_usage:
                        indication = fda_info.indications_and_usage
                    elif fda_info and fda_info.biomarker_text_exact:
                        indication = fda_info.biomarker_text_exact
                    else:
                        indication = tumor or "-"

                    links = []
                    if fda_info and fda_info.drugs_at_fda_url:
                        links.append(f"[FDA]({fda_info.drugs_at_fda_url})")
                    else:
                        fda_url = get_best_fda_url(drug) if drug else None
                        if fda_url:
                            links.append(f"[FDA]({fda_url})")
                    if fda_info and fda_info.pubmed_url:
                        links.append(f"[PubMed]({fda_info.pubmed_url})")
                    links_str = " · ".join(links) if links else "-"

                    rows.append(f"| {drug_link} | {locus_display} | {tumor_match_cell} | {association} | {reg_status} | {indication} | {links_str} |")
                elif show_level:
                    level = b.get('evidence_level', '')
                    rows.append(f"| {drug_link} | {locus_display} | {association} | {tumor} | {level} |")
                else:
                    rows.append(f"| {drug_link} | {locus_display} | {association} | {tumor} |")
            scrollable_table("\n".join(rows))

        # 1. FDA-Approved
        if cgi_biomarkers:
            cgi_fda_sorted = sorted(cgi_biomarkers, key=assoc_sort_key)
            st.markdown("**✅ FDA Approved (CGI)**")
            st.caption("FDA-approved biomarker-drug associations - see FDA tab for label details")
            render_cgi_table(cgi_fda_sorted)

        # 2. Early Phase
        if early_phase:
            if cgi_biomarkers:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown("**🔬 Clinical Evidence (CGI)**")
            st.caption("Late trials, early trials, case reports")
            render_cgi_table(early_phase, show_level=True)

        # 3. Preclinical
        if preclinical:
            if cgi_biomarkers or early_phase:
                st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
            st.markdown("**🧪 Preclinical (CGI)**")
            st.caption("Cell line and preclinical studies")
            render_cgi_table(preclinical, show_level=True)


# ==============================================
# FDA TAB
# ==============================================

def render_fda_tab(fda_biomarker_evidence: list):
    """Render the FDA Approvals tab."""
    st.caption("FDA biomarker-drug indications parsed from drug labels including near hits")
    st.warning("Auto-parsed from FDA labels - verify before clinical use")

    rows = ["| Drug | Gene | Match | Label Level | Variants | Tumors |",
            "|------|------|-------|-------------|----------|--------|"]
    for ev in fda_biomarker_evidence:
        drug = ev.get('drug_name', 'Unknown')
        brand = ev.get('brand_name', '')
        fda_url = ev.get('fda_label_url', '')
        # Include combination partners if present
        combination_partners = ev.get('combination_partners', [])
        if combination_partners:
            partners_str = " + ".join(combination_partners)
            drug_with_partners = f"{drug} + {partners_str}"
        else:
            drug_with_partners = drug
        drug_text = f"{drug_with_partners} ({brand})" if brand else drug_with_partners
        drug_display = f"[{drug_text}]({fda_url})" if fda_url else drug_text
        gene = ev.get('gene', '')
        # Match: did the query variant match this indication?
        match_type = ev.get('variant_match_result', '')
        if match_type in ('exact', 'codon', 'gene'):
            match_display = '✅ Yes'
        else:
            match_display = '❌ No'
        # Label Level
        specificity = ev.get('specificity', '')
        spec_icons = {'variant': '🎯 Variant', 'codon': '📍 Codon', 'gene': '🧬 Gene', 'pathway': '🔀 Pathway'}
        specificity_display = spec_icons.get(specificity, specificity or '—')
        # Show specified variants
        variants = ev.get('specified_variants', [])
        if variants:
            variants_display = ', '.join(variants)[:40]
        elif specificity == 'gene':
            variants_display = 'Any mutation'
        else:
            variants_display = '—'
        tumors_list = ev.get('tumor_types', [])
        tumor_stage = ev.get('tumor_stage', '')
        if tumor_stage and tumors_list:
            tumors = f"{tumor_stage} {', '.join(tumors_list)}"
        elif tumors_list:
            tumors = ', '.join(tumors_list)
        else:
            tumors = '—'
        tumors = tumors[:40]
        rows.append(f"| {drug_display} | {gene} | {match_display} | {specificity_display} | {variants_display} | {tumors} |")

    scrollable_table("\n".join(rows))


# ==============================================
# CLINVAR TAB
# ==============================================

def render_clinvar_tab(clinvar_entries: list, clinvar_sig: str | None):
    """Render the ClinVar tab."""
    if clinvar_sig:
        st.markdown(f"**Clinical Significance**")
    if clinvar_entries:
        rows = ["| Variation ID | Significance | Conditions | Review Status |",
                "|--------------|--------------|------------|---------------|"]
        for entry in clinvar_entries:
            review = entry.get('review_status', '')
            if review.lower() in (
                'no assertion criteria provided',
                'no classification provided',
                'no classification for the individual variant',
            ):
                continue
            var_id = entry.get('variation_id', '')
            var_url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_id}/" if var_id else ''
            var_link = f"[{var_id}]({var_url})" if var_id and var_url else (var_id or '-')
            sig = entry.get('clinical_significance', 'Unknown')
            conds = entry.get('conditions', [])
            conds_str = ', '.join(conds)[:30] if conds else 'N/A'
            rows.append(f"| {var_link} | {sig} | {conds_str} | {review} |")
        scrollable_table("\n".join(rows))
    elif not clinvar_sig:
        st.info("No ClinVar entries found")


# ==============================================
# COSMIC TAB
# ==============================================

def render_cosmic_tab(cosmic_id: str):
    """Render the COSMIC tab."""
    cosmic_num = cosmic_id.replace('COSM', '').replace('COSV', '')
    st.markdown(f"**COSMIC ID:** [{cosmic_id}](https://cancer.sanger.ac.uk/cosmic/mutation/overview?id={cosmic_num})")
    st.caption("Click link to view mutation details in COSMIC database")


# ==============================================
# CLINICAL TRIALS TAB
# ==============================================

def render_trials_tab(trials: list):
    """Render the Clinical Trials tab."""
    # Use columns: tables on left, legend on right
    table_col, legend_col = st.columns([4, 1])

    # Render legend first
    legend_col.markdown(TRIALS_LEGEND, unsafe_allow_html=True)

    with table_col:
        rows = ["| Locus Match | Tumor Match | NCT ID | Phase | Status | Title |",
                "|-------------|-------------|--------|-------|--------|-------|"]
        for t in trials:
            locus_match = t.get('locus_match', '')
            match_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(locus_match, "🧬 Gene")

            tumor_match = t.get('tumor_match')
            if tumor_match is True:
                tumor_match_display = "✅ Yes"
            elif tumor_match is False:
                tumor_match_display = "🔸 Other"
            else:
                tumor_match_display = "-"

            nct_id = t.get('nct_id', '')
            nct_url = t.get('url') or f"https://clinicaltrials.gov/study/{nct_id}" if nct_id else ''
            nct_link = f"[{nct_id}]({nct_url})" if nct_id and nct_url else nct_id
            phase = t.get('phase', 'N/A')
            status = t.get('status', '')
            title = t.get('title', '') or ''
            title_display = title[:150] + "..." if len(title) > 150 else title
            rows.append(f"| {match_display} | {tumor_match_display} | {nct_link} | {phase} | {status} | {title_display} |")
        scrollable_table("\n".join(rows))


# ==============================================
# LITERATURE TAB
# ==============================================

def render_literature_tab(articles: list):
    """Render the Literature tab."""
    rows = ["| PMID | Year | Journal | Signal | Title |",
            "|------|------|---------|--------|-------|"]
    for a in articles:
        pmid = a.get('pmid', '')
        url = a.get('url') or f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ''
        pmid_link = f"[{pmid}]({url})" if pmid and url else pmid
        year = a.get('year', '')
        journal = (a.get('journal', '') or '')[:20]
        signal = a.get('signal_type', '') or '-'
        title = (a.get('title', '') or '')[:80] + "..."
        rows.append(f"| {pmid_link} | {year} | {journal} | {signal} | {title} |")
    scrollable_table("\n".join(rows))


# ==============================================
# CBIOPORTAL TAB
# ==============================================

def render_cbioportal_tab(cbioportal: dict):
    """Render the cBioPortal tab."""
    study_name = cbioportal.get('study_name', 'N/A')
    study_id = cbioportal.get('study_id', '')
    total = cbioportal.get('total_samples', 0)

    if study_id:
        study_url = f"https://www.cbioportal.org/study/summary?id={study_id}"
        st.markdown(f"**Study:** [{study_name}]({study_url}) — cohort of {total:,} samples")
    else:
        st.markdown(f"**Study:** {study_name} — cohort of {total:,} samples")

    gene_pct = cbioportal.get('gene_prevalence_pct', 0)
    variant_pct = cbioportal.get('variant_prevalence_pct', 0)
    gene_count = cbioportal.get('samples_with_gene_mutation', 0)
    variant_count = cbioportal.get('samples_with_exact_variant', 0)
    gene_symbol = cbioportal.get('gene', 'Gene')
    variant_name = cbioportal.get('variant', 'variant')

    prev_cols = st.columns(2)
    with prev_cols[0]:
        st.metric(f"Any {gene_symbol} alteration", f"{gene_pct:.1f}%")
        st.caption(f"{gene_count}/{total}")
    with prev_cols[1]:
        st.metric(f"Exact {variant_name}", f"{variant_pct:.1f}%")
        st.caption(f"{variant_count}/{total}")

    co_occurring = cbioportal.get('co_occurring', [])
    mutually_exclusive = cbioportal.get('mutually_exclusive', [])

    if co_occurring or mutually_exclusive:
        st.markdown("---")
        co_col, me_col = st.columns(2)

        with co_col:
            if co_occurring:
                st.markdown(f"**Co-occurring ({len(co_occurring)}):**")
                st.caption("_Odds > 1 — possible functional interaction_")
                co_rows = []
                for c in co_occurring:
                    odds = c.get('odds_ratio')
                    odds_str = f"{odds:.2f}" if odds else "N/A"
                    co_rows.append({
                        "Gene": c.get('gene', ''),
                        "Count": c.get('count', 0),
                        "Freq": f"{c.get('pct', 0):.1f}%",
                        "OR": odds_str,
                    })
                st.dataframe(pd.DataFrame(co_rows), hide_index=True, width="stretch", height=min(300, 35 * (len(co_rows) + 1)))

        with me_col:
            if mutually_exclusive:
                st.markdown(f"**Mutually Exclusive ({len(mutually_exclusive)}):**")
                st.caption("_Odds < 1 — likely redundant drivers_")
                me_rows = []
                for m in mutually_exclusive:
                    odds = m.get('odds_ratio')
                    odds_str = f"{odds:.2f}" if odds else "N/A"
                    me_rows.append({
                        "Gene": m.get('gene', ''),
                        "Count": m.get('count', 0),
                        "Freq": f"{m.get('pct', 0):.1f}%",
                        "OR": odds_str,
                    })
                st.dataframe(pd.DataFrame(me_rows), hide_index=True, width="stretch", height=min(300, 35 * (len(me_rows) + 1)))


# ==============================================
# DEPMAP TAB
# ==============================================

def render_depmap_tab(depmap: dict, gene_display: str, variant_display: str, tumor_display: str | None = None):
    """Render the DepMap tab."""
    gene_dep = depmap.get('gene_dependency')
    drug_sens = depmap.get('drug_sensitivities', [])
    cell_lines = depmap.get('cell_line_models', [])
    is_essential = depmap.get('is_essential', False)

    if gene_dep:
        st.markdown("### Gene Essentiality")
        score = gene_dep.get('mean_dependency_score')
        dep_pct = gene_dep.get('dependency_pct', 0)
        n_total = gene_dep.get('n_total_lines', 0)

        dep_cols = st.columns(2)
        with dep_cols[0]:
            if is_essential:
                st.error(f"🔴 **{gene_display} is ESSENTIAL**")
            else:
                st.info(f"⚪ {gene_display} is not essential")
        with dep_cols[1]:
            st.metric("Gene CERES Score", f"{score:.2f}" if score else "N/A")
            st.caption(f"{dep_pct:.0f}% of cell lines depend on {gene_display}")
        st.caption(f"Based on CRISPR screens in {n_total} cancer cell lines. CERES < -0.5 indicates essentiality.")

    if drug_sens:
        st.markdown("---")
        st.markdown("### Drug Sensitivities")
        drug_rows = []
        for ds in drug_sens:
            log2fc = ds.get('mean_log2fc')
            log2fc_str = f"{log2fc:.2f}" if log2fc is not None else "N/A"
            drug_rows.append({
                "Drug": ds.get('drug_name', ''),
                "Log2FC": log2fc_str,
                "Cell Lines": ds.get('n_cell_lines', 0),
            })
        st.dataframe(pd.DataFrame(drug_rows), width="stretch", hide_index=True, height=min(300, 35 * (len(drug_rows) + 1)))

    if cell_lines:
        st.markdown("---")
        st.markdown("### Model Cell Lines")
        mutant_lines = [cl for cl in cell_lines if cl.get('has_mutation')]
        if mutant_lines:
            st.success(f"✅ {len(mutant_lines)} cell lines with {gene_display} {variant_display} mutation")
            cl_rows = []
            for cl in mutant_lines:
                name = cl.get('name', '')
                depmap_id = cl.get('depmap_id')
                url = f"https://depmap.org/portal/cell_line/{depmap_id}" if depmap_id else None
                # Tumor match
                cl_disease = (cl.get('primary_disease') or '').lower()
                tumor_lower = (tumor_display or '').lower() if tumor_display else ''
                if tumor_lower and tumor_lower in cl_disease:
                    tumor_match_cell = "✅ Yes"
                elif tumor_lower and cl_disease:
                    tumor_match_cell = "🔸 Other"
                else:
                    tumor_match_cell = "-"
                cl_rows.append({
                    "Cell Line": name,
                    "DepMap": url,
                    "Tumor Match": tumor_match_cell,
                    "Disease": cl.get('primary_disease', ''),
                    "Subtype": cl.get('subtype', ''),
                    "Mutation": cl.get('mutation_details', variant_display),
                })
            st.dataframe(
                pd.DataFrame(cl_rows),
                width="stretch",
                hide_index=True,
                height=min(300, 35 * (len(cl_rows) + 1)),
                column_config={"DepMap": st.column_config.LinkColumn("DepMap", display_text="🔗")},
            )
        else:
            st.info(f"{len(cell_lines)} cell lines available (mutation status unknown)")

    st.markdown(f"[🔗 Explore on DepMap Portal](https://depmap.org/portal/gene/{gene_display})")


# ==============================================
# HOTSPOTS TAB
# ==============================================

def render_hotspots_tab(hotspots: dict, gene_display: str):
    """Render the Hotspots tab."""
    is_at_hotspot = hotspots.get('is_hotspot', False)
    is_adjacent = hotspots.get('is_adjacent_to_hotspot', False)

    def hotspot_summary(residue: str, hotspot_data: dict, total_samples: int) -> str:
        tumor_types = hotspot_data.get('tumor_type_composition', [])
        sorted_tumors = sorted(tumor_types, key=lambda x: x.get('count', 0), reverse=True)[:3]
        if sorted_tumors:
            tumor_strs = [t['tumor_type'].replace('_', ' ').title() for t in sorted_tumors]
            return f"The {residue} hotspot has {total_samples:,} samples, most common in {', '.join(tumor_strs)}."
        return f"The {residue} hotspot has {total_samples:,} samples."

    if is_at_hotspot:
        hotspot_data = hotspots.get('hotspot', {})
        residue = hotspot_data.get('residue', '')
        total_samples = hotspot_data.get('total_samples', 0)

        st.markdown(f"**This variant is at known cancer hotspot {gene_display} {residue}.**")
        st.caption(hotspot_summary(residue, hotspot_data, total_samples))

    elif is_adjacent:
        adjacent_data = hotspots.get('adjacent_hotspot', {})
        adjacent_distance = hotspots.get('adjacent_distance', 0)
        residue = adjacent_data.get('residue', '')
        total_samples = adjacent_data.get('total_samples', 0)

        st.markdown(f"**This variant is ±{adjacent_distance} codons from known hotspot {gene_display} {residue}.**")
        st.caption(hotspot_summary(residue, adjacent_data, total_samples))

    st.caption("[Cancer Hotspots](https://www.cancerhotspots.org/) (Chang et al., Cancer Discovery 2017)")

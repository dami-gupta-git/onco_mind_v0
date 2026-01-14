"""OncoMind Streamlit Application - Variant insight and evidence synthesis tool."""
import asyncio
import pandas as pd
import re
import json
import streamlit as st
from backend import get_variant_insight, batch_get_variant_insights
from pdb_images import get_pdb_image_url, get_pdb_page_url
from oncomind.api.fda_drugs import get_best_fda_url, get_cached_fda_label
from oncomind.config.constants import (
    LLM_DEFAULT_TEMPERATURE,
    UI_MAX_CIVIC_EVIDENCE_ROWS,
    CADD_DELETERIOUS_THRESHOLD,
    GNOMAD_RARE_THRESHOLD,
)

# Import UI components
from components import (
    apply_styles,
    scrollable_table,
    result_to_markdown,
    render_batch_tab,
    render_llm_synthesis,
    render_cross_source_analysis,
)

# Initialize logging from environment variable (ONCOMIND_LOG_LEVEL=DEBUG|INFO|WARN|ERROR)
# This import triggers the logger initialization which reads from env
from oncomind.config.debug import get_logger
logger = get_logger(__name__)

st.set_page_config(page_title="OncoMind", page_icon="🧬", layout="wide")

# Apply custom styling
apply_styles()
# scrollable_table function moved to components/utils.py

st.markdown("<h1 style='margin-bottom: 0;'><span style='font-size: 0.85em;'>🧬</span> OncoMind: Variant Insight</h1>", unsafe_allow_html=True)
st.markdown(
    "<div style='background: linear-gradient(90deg, #ff6b6b, #ffa500); color: white; padding: 8px 16px; "
    "border-radius: 4px; font-weight: 600; text-align: center; margin: 8px 0;'>"
    "Proof of  Concept - What DON'T you know about your variant?</div>",
    unsafe_allow_html=True
)
st.caption("**Note:** This tool is for research purposes only. Clinical decisions should always be made by qualified healthcare professionals.")

MODELS = {
    "Anthropic Claude Sonnet 4": "claude-sonnet-4-20250514",
    "Anthropic Claude 3.5 Haiku": "claude-3-5-haiku-20241022",
    "OpenAI GPT-4o-mini": "gpt-4o-mini",
    "OpenAI GPT-4o": "gpt-4o",
    "Google Gemini 1.5 Pro": "gemini/gemini-1.5-pro",
    "xAI Grok 3": "xai/grok-3",
}

# Pre-populated example variants (Gene, Variant, Tumor Type)
# EXAMPLE_VARIANTS = {
#     "-- Select an example --": ("", "", ""),
#     "EGFR C797S - NSCLC": ("EGFR", "C797S", "NSCLC"),
#     "BRAF V600E - Melanoma": ("BRAF", "V600E", "Melanoma"),
#     "BRAF K601E - Brain (near-hotspot)": ("BRAF", "K601E", "Brain"),
#     "EGFR L858R - NSCLC": ("EGFR", "L858R", "NSCLC"),
#     "EGFR T790M - NSCLC": ("EGFR", "T790M", "NSCLC"),
#     "KRAS G12A - Lung Cancer": ("KRAS", "G12A", "Lung Cancer"),
#     "KRAS G12D - NSCLC": ("KRAS", "G12D", "NSCLC"),
#     "KRAS G12C - NSCLC": ("KRAS", "G12C", "NSCLC"),
#     "EGFR N771H - NSCLC": ("EGFR", "N771H", "NSCLC"),
#     "PIK3CA H1047R  - Thyroid Cancer": ("PIK3CA", "H1047R", "Thyroid Cancer"),
#     "AKT1 E17K - Breast Cancer": ("AKT1", "E17K", "Breast Cancer"),
#     "GNAQ Q209L - Uveal Melanoma": ("GNAQ", "Q209L", "Uveal Melanoma"),
#     "KIT D816V - GIST": ("KIT", "D816V", "GIST"),
#     "ALK F1174L - Neuroblastoma": ("ALK", "F1174L", "Neuroblastoma"),
#     "ERBB2 S310F - Bladder Cancer": ("ERBB2", "S310F", "Bladder Cancer"),
#     "IDH1 R132H - Glioma": ("IDH1", "R132H", "Glioma"),
#     "NRAS Q61R - Melanoma": ("NRAS", "Q61R", "Melanoma"),
#     "CHEK2 p.Ile157Thr - Colorectal": ("CHEK2", "p.Ile157Thr", "Colorectal Cancer"),
# }
EXAMPLE_VARIANTS = {
    "IDH1 R132H - Glioma": ("IDH1", "R132H", "Glioma"),
    "CHEK2 p.Ile157Thr - Colorectal": ("CHEK2", "p.Ile157Thr", "Colorectal Cancer"),
    "ALK F1174L - Neuroblastoma": ("ALK", "F1174L", "Neuroblastoma"),
}

# Initialize session state for persisting results
if "single_result" not in st.session_state:
    st.session_state.single_result = None
if "single_gene" not in st.session_state:
    st.session_state.single_gene = None
if "single_variant" not in st.session_state:
    st.session_state.single_variant = None
if "batch_results" not in st.session_state:
    st.session_state.batch_results = None
if "batch_df" not in st.session_state:
    st.session_state.batch_df = None

tab1, tab2 = st.tabs(["🔬 Single Variant", "📊 Batch Upload"])

# TAB 1: Single Variant
with tab1:
    # ==============================================
    # COMPACT INPUT ROW
    # ==============================================

    # Callback for example dropdown - updates the text input keys directly
    def on_example_change():
        selected = st.session_state.get("example_selector")
        if selected and selected != "-- Select an example --":
            gene, variant, tumor = EXAMPLE_VARIANTS[selected]
            st.session_state.gene_input = gene
            st.session_state.variant_input = variant
            st.session_state.tumor_input = tumor

    input_cols = st.columns([1.5, 1.5, 2, 1.2, 1.2, 1])
    with input_cols[0]:
        gene = st.text_input("Gene", value="IDH1", placeholder="e.g. IDH1", key="gene_input")
    with input_cols[1]:
        variant = st.text_input("Variant", value="R132H", placeholder="e.g. R132H", key="variant_input")
    with input_cols[2]:
        tumor = st.text_input("Tumor Type", value="Glioma", placeholder="e.g. Glioma", key="tumor_input")
    with input_cols[3]:
        st.markdown("<div style='height: 28px'></div>", unsafe_allow_html=True)  # Spacer to align with labels
        enable_literature = st.toggle(
            "Literature",
            value=True,
            help="Search recent literature via Semantic Scholar (with citations). Falls back to PubMed if rate limited."
        )
        # Semantic Scholar by default, PubMed is automatic fallback on rate limit
        literature_source = "semantic_scholar" if enable_literature else "none"
    with input_cols[4]:
        st.markdown("<div style='height: 28px'></div>", unsafe_allow_html=True)  # Spacer to align with labels
        enable_llm = st.toggle(
            "LLM Mode",
            value=True,
            help="LLM mode includes AI-powered synthesis (~15s). Without LLM, you get fast annotation (~5s)."
        )
    with input_cols[5]:
        st.markdown("<div style='height: 28px'></div>", unsafe_allow_html=True)  # Spacer to align with labels
        insight_btn = st.button("🔍 Go", type="primary", width="stretch")

    # Example variants dropdown (experimental - below the input row)
    st.selectbox(
        "Try an example",
        options=list(EXAMPLE_VARIANTS.keys()),
        key="example_selector",
        on_change=on_example_change,
        label_visibility="collapsed"
    )

    # LLM settings expander (only if LLM enabled)
    if enable_llm:
        with st.expander("⚙️ LLM Settings", expanded=False):
            llm_cols = st.columns(2)
            with llm_cols[0]:
                model_name = st.selectbox("Model", list(MODELS.keys()))
            with llm_cols[1]:
                temperature = st.slider("Temperature", 0.0, 1.0, LLM_DEFAULT_TEMPERATURE, 0.05)
    else:
        model_name = list(MODELS.keys())[0]
        temperature = LLM_DEFAULT_TEMPERATURE

    # ==============================================
    # PROCESS REQUEST
    # ==============================================
    if insight_btn:
        if not gene or not variant:
            st.error("Gene and variant are required")
        else:
            from oncomind.utils.variant_normalization import normalize_variant, VariantNormalizer
            normalized = normalize_variant(gene, variant)
            variant_type = normalized['variant_type']

            if variant_type not in VariantNormalizer.ALLOWED_VARIANT_TYPES:
                st.error(
                    f"❌ Unsupported variant type: **{variant_type}**\n\n"
                    f"Only **SNPs and small indels** are supported:\n"
                    f"- Missense mutations (e.g., V600E)\n"
                    f"- Nonsense mutations (e.g., R172*)\n"
                    f"- Small insertions/deletions (e.g., ins, del)\n"
                    f"- Frameshift mutations (e.g., fs)\n\n"
                    f"Your variant '{variant}' is classified as '{variant_type}'."
                )
            else:
                mode_label = "LLM" if enable_llm else "Annotation"
                with st.spinner(f"🔬 Getting insight for {gene} {variant} ({mode_label} mode)..."):
                    # Check for timing flag via environment variable
                    import os
                    enable_timing = os.environ.get("ONCOMIND_TIMING", "").lower() in ("1", "true", "yes")
                    result = asyncio.run(get_variant_insight(
                        gene, variant, tumor or None,
                        enable_llm=enable_llm,
                        enable_literature=enable_literature,
                        literature_source=literature_source,
                        model=MODELS[model_name],
                        temperature=temperature,
                        enable_timing=enable_timing,
                    ))
                    if "error" in result:
                        st.error(result["error"])
                    else:
                        st.session_state.single_result = result
                        st.session_state.single_gene = gene
                        st.session_state.single_variant = variant

    # ==============================================
    # DISPLAY RESULTS
    # ==============================================
    if st.session_state.single_result is not None:
        result = st.session_state.single_result
        gene_display = st.session_state.single_gene
        variant_display = st.session_state.single_variant
        tumor_display = result.get('variant', {}).get('tumor_type')
        ids = result.get('identifiers', {})
        hgvs = result.get('hgvs', {})
        transcript = result.get('transcript', {})

        # ==============================================
        # HEADER: Variant name + Download button (no border)
        # ==============================================
        header_text = f"**{gene_display} {variant_display}**"
        if tumor_display:
            header_text += f" in {tumor_display}"

        # Header row with download button next to it
        header_col1, header_col2 = st.columns([2, 1])
        with header_col1:
            st.markdown(f"<span style='font-size: 1.1rem; font-weight: 600;'>✅ {header_text}</span>", unsafe_allow_html=True)
        with header_col2:
            st.download_button(
                "📥 Download Report (Markdown)",
                result_to_markdown(result),
                f"{gene_display}_{variant_display}_report.md",
                "text/markdown",
                key="download_top_md"
            )

        # Compact linked identifiers row
        id_links = []
        if ids.get('cosmic_id'):
            cosmic_id = ids['cosmic_id']
            cosmic_num = cosmic_id.replace('COSM', '').replace('COSV', '')
            id_links.append(f"[COSMIC:{cosmic_id}](https://cancer.sanger.ac.uk/cosmic/mutation/overview?id={cosmic_num})")
        if ids.get('dbsnp_id'):
            id_links.append(f"[dbSNP:{ids['dbsnp_id']}](https://www.ncbi.nlm.nih.gov/snp/{ids['dbsnp_id']})")
        if ids.get('clinvar_id'):
            id_links.append(f"[ClinVar:{ids['clinvar_id']}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{ids['clinvar_id']}/)")
        if ids.get('ncbi_gene_id'):
            id_links.append(f"[NCBI:{ids['ncbi_gene_id']}](https://www.ncbi.nlm.nih.gov/gene/{ids['ncbi_gene_id']})")
        # Ensembl transcript link
        transcript_id = transcript.get('id')
        if transcript_id:
            # Extract base transcript ID (remove version number for Ensembl link)
            ensembl_id = transcript_id.split('.')[0] if '.' in transcript_id else transcript_id
            id_links.append(f"[Ensembl:{transcript_id}](https://ensembl.org/Homo_sapiens/Transcript/Summary?t={ensembl_id})")
        # HGVS genomic notation link (links to Ensembl VEP)
        hgvs_g = hgvs.get('genomic')
        if hgvs_g:
            # Display shortened version, link to VEP
            hgvs_short = hgvs_g if len(hgvs_g) <= 25 else hgvs_g[:22] + "..."
            id_links.append(f"[HGVS:{hgvs_short}](https://ensembl.org/Homo_sapiens/Tools/VEP/Results?hgvs={hgvs_g})")
        # Always show these search links
        id_links.append(f"[PubMed](https://pubmed.ncbi.nlm.nih.gov/?term={gene_display}+{variant_display})")
        id_links.append(f"[CIViC](https://civicdb.org/variants?geneSearch={gene_display})")
        id_links.append(f"[OncoKB](https://www.oncokb.org/gene/{gene_display})")
        id_links.append(f"[DepMap](https://depmap.org/portal/gene/{gene_display})")

        st.markdown(" &nbsp;|&nbsp; ".join(id_links))

        # ==============================================
        # EVIDENCE SOURCES (tabs) - in bordered card
        # ==============================================
        with st.container(border=True):
            st.caption("Curated clinical and research evidence from cancer knowledge bases, clinical trials, and literature.")
            st.markdown("<span style='font-size: 1.5rem; font-weight: 600;'>📚 Evidence Sources</span>", unsafe_allow_html=True)

            # Collect available sources
            fda_approvals = result.get('fda_approvals', [])
            civic_assertions = result.get('civic_assertions', [])
            civic_evidence = result.get('civic_evidence', [])
            vicc = result.get('vicc_evidence', [])
            cgi_biomarkers = result.get('cgi_biomarkers', [])
            clinvar_entries = result.get('clinvar_entries', [])
            clinvar_sig = result.get('clinvar', {}).get('clinical_significance')
            cosmic_id = ids.get('cosmic_id')
            trials = result.get('clinical_trials', [])
            articles = result.get('pubmed_articles', [])
            preclinical = result.get('preclinical_biomarkers', [])
            early_phase = result.get('early_phase_biomarkers', [])
            annotations = result.get('annotations', {})
            cbioportal = result.get('cbioportal_evidence')
            depmap = result.get('depmap_evidence')
            hotspots = result.get('hotspots_evidence')
            therapies = result.get('recommended_therapies', [])
            fda_labels = result.get('fda_labels', [])

            # Build tab names
            tab_names = []
            has_functional = any([
                annotations.get('alphamissense_score'),
                annotations.get('alphamissense_prediction'),
                annotations.get('cadd_score'),
                annotations.get('polyphen2_prediction'),
                annotations.get('gnomad_exome_af'),
                annotations.get('snpeff_effect'),
            ])
            if has_functional:
                tab_names.append("Functional")
            if civic_assertions or civic_evidence:
                tab_names.append(f"CIViC ({len(civic_assertions) + len(civic_evidence)})")
            if vicc:
                tab_names.append(f"VICC ({len(vicc)})")
            # CGI tab - shows ALL CGI biomarkers organized by evidence level
            all_cgi_count = len(cgi_biomarkers) + len(preclinical) + len(early_phase)
            if all_cgi_count > 0:
                tab_names.append(f"CGI ({all_cgi_count})")
            # FDA tab - shows FDA label data (clinical studies, mechanism, toxicity)
            if fda_labels:
                tab_names.append(f"💊 FDA ({len(fda_labels)})")
            if clinvar_entries or clinvar_sig:
                tab_names.append("ClinVar")
            if cosmic_id:
                tab_names.append("COSMIC")
            if trials:
                tab_names.append(f"Trials ({len(trials)})")
            if articles:
                tab_names.append(f"Literature ({len(articles)})")
            # Research tab removed - CGI preclinical/early_phase data now shown in CGI Therapies tab
            if cbioportal:
                tab_names.append("cBioPortal")
            if depmap:
                tab_names.append("🧬 DepMap")
            if hotspots and (hotspots.get('is_hotspot') or hotspots.get('is_adjacent_to_hotspot')):
                tab_names.append("🔥 Hotspots")

            if tab_names:
                tabs = st.tabs(tab_names)
                tab_idx = 0

                # Functional tab
                if has_functional:
                    with tabs[tab_idx]:
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
                                freq = f"{af:.2e}" if af < GNOMAD_RARE_THRESHOLD else f"{af:.4f}"
                                rows.append(f"| gnomAD AF | {freq} | {'Rare' if af < GNOMAD_RARE_THRESHOLD else 'Common'} |")
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
                    tab_idx += 1

                # CIViC tab
                if civic_assertions or civic_evidence:
                    with tabs[tab_idx]:
                        # Use columns: tables on left, legend on right
                        table_col, legend_col = st.columns([4, 1])

                        # Render legend first, with padding to align with Level A table
                        legend_col.markdown("""<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 180px;'>
<b>Locus Match:</b><br/>
🎯 Variant (exact locus match)<br/>
📍 Codon (other variants in this codon)<br/>
🧬 Gene (other variants in this gene)<br/><br/>
<b>Tumor Match:</b><br/>
✅ Yes (match on specified tumor)<br/>
🔸 Other (match on other tumor)<br/>
🌐 Pan-cancer
</div>""", unsafe_allow_html=True)

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
                                    # Note: Drug Specificity removed - OpenFDA returns outdated label text
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

                                        # Biomarker Specificity: gene vs variant vs phenotype level
                                        specificity = "-"
                                        if fda_info and fda_info.biomarker_specificity:
                                            spec = fda_info.biomarker_specificity
                                            if spec == 'variant':
                                                variants = fda_info.specific_variants
                                                if variants:
                                                    specificity = f"🎯 {', '.join(variants[:2])}"
                                                else:
                                                    specificity = "🎯 Variant"
                                            elif spec == 'gene':
                                                specificity = "🧬 Gene"
                                            elif spec == 'phenotype':
                                                variants = fda_info.specific_variants
                                                if variants:
                                                    specificity = f"🔬 {', '.join(variants[:2])}"
                                                else:
                                                    specificity = "🔬 Phenotype"

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
                    tab_idx += 1

                # VICC tab
                if vicc:
                    with tabs[tab_idx]:
                        # Use columns: tables on left, legend on right
                        table_col, legend_col = st.columns([4, 1])

                        # Render legend first, with padding to align with first data row
                        legend_col.markdown("""<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 85px;'>
<b>Locus Match:</b><br/>
🎯 Variant (exact locus match)<br/>
📍 Codon (other variants in this codon)<br/>
🧬 Gene (other variants in this gene)<br/><br/>
<b>Tumor Match:</b><br/>
✅ Yes (match on specified tumor)<br/>
🔸 Other (match on other tumor)<br/>
🌐 Pan-cancer
</div>""", unsafe_allow_html=True)

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
                            # VICC levels can be: 1, 2, 3, 4, A, B, C, D, R1, R2, etc.
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
                                    # FDA level table: Source, Locus Match, Tumor Match, Drugs, Response, Regulatory Status, Indication, Links
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
                                    source_link = f"[{source}]({pub_url})" if pub_url else source
                                    locus_match = v.get('locus_match', '')
                                    match_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(locus_match, "")
                                    tumor_match = v.get('tumor_match')
                                    tumor_match_cell = "✅ Yes" if tumor_match else "🔸 Other"
                                    drugs_list = v.get('drugs', [])
                                    drugs = ", ".join(drugs_list) or "N/A"
                                    drugs = drugs[:30] if len(drugs) > 30 else drugs
                                    # Handle response_type - some sources (molecularmatch) return AMP level codes
                                    # like "2C" instead of clinical response types like "Sensitive"
                                    response = v.get('response_type', 'Unknown')
                                    import re
                                    if response and re.match(r'^[1234][A-D]?$', response.upper()):
                                        # This is an AMP tier code, not a response type
                                        # Try to infer from description or show as unknown
                                        desc = (v.get('description') or '').lower()
                                        if 'sensitivity' in desc or 'sensitive' in desc or 'response' in desc:
                                            response = 'sensitive'
                                        elif 'resistance' in desc or 'resistant' in desc:
                                            response = 'resistant'
                                        else:
                                            response = '-'
                                    disease = (v.get('disease', '') or '')[:25]

                                    if show_fda_cols:
                                        # Get FDA label info for first drug
                                        first_drug = drugs_list[0] if drugs_list else None
                                        fda_info = get_cached_fda_label(first_drug) if first_drug else None

                                        # Regulatory Status
                                        if fda_info and fda_info.approval_date:
                                            reg_status = f"FDA-approved ({fda_info.approval_date})"
                                        else:
                                            reg_status = "FDA-approved"

                                        # Indication - use full indications_and_usage, not biomarker_text_exact
                                        # biomarker_text_exact is a substring that often starts mid-sentence
                                        if fda_info and fda_info.indications_and_usage:
                                            indication = fda_info.indications_and_usage[:60]
                                            if len(fda_info.indications_and_usage) > 60:
                                                indication += "..."
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
                    tab_idx += 1

                # CGI tab - shows ALL CGI data organized by evidence level
                if all_cgi_count > 0:
                    with tabs[tab_idx]:
                        # Use columns: tables on left, legend on right
                        table_col, legend_col = st.columns([4, 1])

                        # Render legend
                        legend_col.markdown("""<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 85px;'>
<b>Locus Match:</b><br/>
🎯 Variant (exact locus match)<br/>
📍 Codon (other variants in this codon)<br/>
🧬 Gene (other variants in this gene)<br/><br/>
<b>Association:</b><br/>
<span style="color: #28a745">Responsive</span> - drug works<br/>
<span style="color: #dc3545">Resistant</span> - drug doesn't work
</div>""", unsafe_allow_html=True)

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
                                # Build header based on options
                                if show_fda_cols:
                                    # FDA-approved table: Drug, Locus Match, Tumor Match, Association, Regulatory Status, Indication, Links
                                    # Note: Drug Specificity removed - OpenFDA returns outdated label text
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
                                    # Drug classes are things like "allosteric AKT inhibitor" which won't have FDA labels
                                    drug_words = drug.lower().split()
                                    is_drug_class = any(w in drug_words for w in ['inhibitor', 'inhibitors', 'agonist', 'antagonist', 'blocker'])
                                    if is_drug_class:
                                        drug_link = drug  # Plain text for drug classes
                                    else:
                                        # Link to DailyMed search for specific drugs
                                        dailymed_url = f"https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query={drug}"
                                        drug_link = f"[{drug}]({dailymed_url})"
                                    # Locus match
                                    match = b.get('locus_match', '')
                                    locus_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "-")
                                    # Tumor match
                                    tumor_match = b.get('cancer_type_match')
                                    tumor_match_cell = "✅ Yes" if tumor_match else "🔸 Other"
                                    # Association with color (using HTML spans for HTML table rendering)
                                    assoc_raw = b.get('association', 'Unknown')
                                    if 'respons' in assoc_raw.lower() or 'sensitiv' in assoc_raw.lower():
                                        association = f'<span style="color: #28a745">{assoc_raw}</span>'
                                    elif 'resist' in assoc_raw.lower():
                                        association = f'<span style="color: #dc3545">{assoc_raw}</span>'
                                    else:
                                        association = assoc_raw
                                    tumor = (b.get('tumor_type', '') or '')[:25]

                                    if show_fda_cols:
                                        # Get FDA label info from cache
                                        fda_info = get_cached_fda_label(drug) if drug else None

                                        # Regulatory Status: "FDA-approved (date)"
                                        if fda_info and fda_info.approval_date:
                                            reg_status = f"FDA-approved ({fda_info.approval_date})"
                                        else:
                                            reg_status = "FDA-approved"

                                        # Indication: full indications_and_usage from FDA label
                                        if fda_info and fda_info.indications_and_usage:
                                            indication = fda_info.indications_and_usage
                                        elif fda_info and fda_info.biomarker_text_exact:
                                            # Fallback to biomarker-specific extract
                                            indication = fda_info.biomarker_text_exact
                                        else:
                                            indication = tumor or "-"

                                        # Links: FDA + PubMed
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

                                        # Biomarker Specificity: gene vs variant vs phenotype level
                                        specificity = "-"
                                        if fda_info and fda_info.biomarker_specificity:
                                            spec = fda_info.biomarker_specificity
                                            if spec == 'variant':
                                                # Show specific variants if available
                                                variants = fda_info.specific_variants
                                                if variants:
                                                    specificity = f"🎯 {', '.join(variants[:2])}"
                                                else:
                                                    specificity = "🎯 Variant"
                                            elif spec == 'gene':
                                                specificity = "🧬 Gene"
                                            elif spec == 'phenotype':
                                                # Show phenotype markers
                                                variants = fda_info.specific_variants
                                                if variants:
                                                    specificity = f"🔬 {', '.join(variants[:2])}"
                                                else:
                                                    specificity = "🔬 Phenotype"

                                        rows.append(f"| {drug_link} | {locus_display} | {tumor_match_cell} | {association} | {reg_status} | {indication} | {links_str} |")
                                    elif show_level:
                                        level = b.get('evidence_level', '')
                                        rows.append(f"| {drug_link} | {locus_display} | {association} | {tumor} | {level} |")
                                    else:
                                        rows.append(f"| {drug_link} | {locus_display} | {association} | {tumor} |")
                                scrollable_table("\n".join(rows))

                            # 1. FDA-Approved (from cgi_biomarkers which only contains fda_approved=True)
                            if cgi_biomarkers:
                                cgi_fda_sorted = sorted(cgi_biomarkers, key=assoc_sort_key)
                                st.markdown("**✅ FDA Approved (CGI)**")
                                st.caption("FDA-approved biomarker-drug associations - see FDA tab for label details")
                                render_cgi_table(cgi_fda_sorted)

                            # 2. Early Phase (clinical trials, late trials, case reports)
                            if early_phase:
                                if cgi_biomarkers:
                                    st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
                                st.markdown("**🔬 Clinical Evidence (CGI)**")
                                st.caption("Late trials, early trials, case reports")
                                render_cgi_table(early_phase, show_level=True)

                            # 3. Preclinical (cell line, pre-clinical)
                            if preclinical:
                                if cgi_biomarkers or early_phase:
                                    st.markdown('<hr class="evidence-section-divider">', unsafe_allow_html=True)
                                st.markdown("**🧪 Preclinical (CGI)**")
                                st.caption("Cell line and preclinical studies")
                                render_cgi_table(preclinical, show_level=True)
                    tab_idx += 1

                # FDA Labels tab - shows detailed drug label information
                if fda_labels:
                    with tabs[tab_idx]:
                        st.caption("FDA drug label data from OpenFDA • Click to expand")

                        # Use columns: drugs on left, legend on right
                        fda_col, legend_col = st.columns([4, 1])

                        # Render legend (consistent with other tabs)
                        legend_col.markdown("""<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 10px;'>
<b>FDA Approval Match:</b><br/>
🟢 Matched<br/>
🔴 Not Matched<br/>
<br/>
<b>Locus Level:</b><br/>
🎯 Variant (exact)<br/>
◐ Codon (same position)<br/>
🧬 Gene (any alteration)<br/>
<br/>
<b>Tumor Match:</b><br/>
✅ Yes<br/>
🌐 Pan-cancer<br/>
🔸 Other
</div>""", unsafe_allow_html=True)

                        # Scrollable container for FDA labels
                        with fda_col:
                            fda_container = st.container(height=400)
                            with fda_container:
                                for label in fda_labels:
                                    drug_name = label.get('drug', 'Unknown')
                                    brand = label.get('brand_name', '')

                                    # Clinical studies data
                                    cs = label.get('clinical_studies') or {}
                                    nct_id = cs.get('nct_id', '')
                                    hr = cs.get('hazard_ratio')
                                    hr_ci = cs.get('hazard_ratio_ci')

                                    # Format HR with 95% CI
                                    if hr:
                                        hr_str = f"{hr:.2f}"
                                        if hr_ci:
                                            hr_str += f" ({hr_ci[0]:.2f}-{hr_ci[1]:.2f})"
                                    else:
                                        hr_str = None

                                    # Mechanism targets
                                    moa = label.get('mechanism_of_action') or {}
                                    targets = moa.get('targets', [])

                                    # Biomarker match icons:
                                    # - 🟢/🔴 for whether variant is covered by FDA approval
                                    # - 🎯/◐/🧬 for match level (variant/codon/gene)
                                    # - ✅/🌐/🔸 for tumor match
                                    biomarker_match = label.get('biomarker_match') or {}
                                    is_matched = biomarker_match.get('matched', False)
                                    match_icon = "🟢" if is_matched else "🔴"
                                    match_level = biomarker_match.get('match_level')
                                    level_icon = {"variant": "🎯", "codon": "◐", "gene": "🧬"}.get(match_level, "")

                                    # Tumor match icon (consistent with other tabs)
                                    # tumor=✅, pan_cancer=🌐, other=🔸
                                    tumor_matched = biomarker_match.get('tumor_matched')
                                    tumor_match_type = biomarker_match.get('tumor_match_type')
                                    if tumor_matched is True:
                                        tumor_icon = "🌐" if tumor_match_type == "pan_cancer" else "✅"
                                    elif tumor_matched is False:
                                        tumor_icon = "🔸"
                                    else:
                                        tumor_icon = ""  # None means no tumor type was queried

                                    # Build expander header with biomarker icons at the beginning, tumor icon at end
                                    icons = f"{match_icon} {level_icon}".strip()
                                    drug_display = f"{icons} **{drug_name}** ({brand})" if brand else f"{icons} **{drug_name}**"
                                    header_parts = [drug_display]
                                    if targets:
                                        header_parts.append(f"Targets: {', '.join(targets)}")
                                    if nct_id:
                                        header_parts.append(f"NCT: {nct_id}")
                                    if hr_str:
                                        header_parts.append(f"HR: {hr_str}")
                                    # Add tumor icon at the end of the header
                                    if tumor_icon:
                                        header_parts.append(tumor_icon)
                                    header = " | ".join(header_parts)

                                    search_term = brand if brand else drug_name
                                    fda_url = f"https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&varq={search_term}"

                                    with st.expander(header):
                                        approved_indications = label.get('approved_indications', [])
                                        indications = label.get('indications_and_usage', '') or '-'

                                        st.markdown(f"[View FDA Label]({fda_url})")
                                        if approved_indications:
                                            st.markdown(f"**Approved for:** {', '.join(approved_indications)}")
                                        st.markdown(f"**Indications & Usage:** {indications}")

                                        # Add a divider before match status
                                        st.divider()

                                        # Show match summary at the end (consistent with other tabs)
                                        match_parts = []
                                        if is_matched:
                                            level_text = {"variant": "exact variant", "codon": "codon", "gene": "gene-level"}.get(match_level, "")
                                            match_parts.append(f"✓ FDA Approval ({level_text})")
                                        else:
                                            match_parts.append("✗ FDA Approval not matched")

                                        if tumor_matched is True:
                                            if tumor_match_type == "pan_cancer":
                                                match_parts.append("🌐 Tumor (pan-cancer)")
                                            else:
                                                match_parts.append("✅ Tumor")
                                        elif tumor_matched is False:
                                            match_parts.append("🔸 Tumor (other)")

                                        if match_parts:
                                            st.markdown(f"**Match Status:** {' | '.join(match_parts)}")

                    tab_idx += 1

                # ClinVar tab
                if clinvar_entries or clinvar_sig:
                    with tabs[tab_idx]:
                        if clinvar_sig:
                            st.markdown(f"**Clinical Significance:** {clinvar_sig}")
                        if clinvar_entries:
                            # Use markdown table for clickable variation IDs
                            rows = ["| Variation ID | Significance | Conditions | Review Status |",
                                    "|--------------|--------------|------------|---------------|"]
                            for entry in clinvar_entries:
                                var_id = entry.get('variation_id', '')
                                var_url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_id}/" if var_id else ''
                                var_link = f"[{var_id}]({var_url})" if var_id and var_url else (var_id or '-')
                                sig = entry.get('clinical_significance', 'Unknown')
                                conds = entry.get('conditions', [])
                                conds_str = ', '.join(conds)[:30] if conds else 'N/A'
                                review = entry.get('review_status', '')
                                rows.append(f"| {var_link} | {sig} | {conds_str} | {review} |")
                            scrollable_table("\n".join(rows))
                        elif not clinvar_sig:
                            st.info("No ClinVar entries found")
                    tab_idx += 1

                # COSMIC tab
                if cosmic_id:
                    with tabs[tab_idx]:
                        cosmic_num = cosmic_id.replace('COSM', '').replace('COSV', '')
                        st.markdown(f"**COSMIC ID:** [{cosmic_id}](https://cancer.sanger.ac.uk/cosmic/mutation/overview?id={cosmic_num})")
                        st.caption("Click link to view mutation details in COSMIC database")
                    tab_idx += 1

                # Trials tab
                if trials:
                    with tabs[tab_idx]:
                        # Use columns: tables on left, legend on right
                        table_col, legend_col = st.columns([4, 1])

                        # Render legend first, with padding to align with first data row
                        legend_col.markdown("""<div style='font-size: 0.85rem; line-height: 1.8; padding-top: 85px;'>
<b>Locus Match:</b><br/>
🎯 Variant (exact locus match)<br/>
📌 Broad (match on related variants)<br/>
🧬 Gene (other variants in this gene)<br/><br/>
<b>Tumor Match:</b><br/>
✅ Yes (match on specified tumor)<br/>
🔸 Other (match on other tumor)<br/>
🌐 Pan-cancer
</div>""", unsafe_allow_html=True)

                        with table_col:
                            # Count matches by type
                            specific_trials = [t for t in trials if t.get('match_scope') == 'specific']
                            ambiguous_trials = [t for t in trials if t.get('match_scope') == 'ambiguous']
                            gene_only_trials = [t for t in trials if not t.get('variant_specific', False)]

                            # Use markdown table for clickable NCT IDs
                            rows = ["| Locus Match | Tumor Match | NCT ID | Phase | Status | Title |",
                                    "|-------------|-------------|--------|-------|--------|-------|"]
                            for t in trials:
                                locus_match = t.get('locus_match', '')
                                # Use consistent format with other tabs: icon + level name only
                                match_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(locus_match, "🧬 Gene")

                                # Tumor match display
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
                                # Show more of the title (truncation handled by CSS)
                                title_display = title[:150] + "..." if len(title) > 150 else title
                                rows.append(f"| {match_display} | {tumor_match_display} | {nct_link} | {phase} | {status} | {title_display} |")
                            scrollable_table("\n".join(rows))
                    tab_idx += 1

                # Literature tab
                if articles:
                    with tabs[tab_idx]:
                        # Use markdown table for clickable PMIDs
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
                    tab_idx += 1

                # cBioPortal tab
                if cbioportal:
                    with tabs[tab_idx]:
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
                    tab_idx += 1

                # DepMap tab
                if depmap:
                    with tabs[tab_idx]:
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
                                    # Tumor match: compare cell line disease to user's tumor type
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
                    tab_idx += 1

                # Hotspots tab
                if hotspots and (hotspots.get('is_hotspot') or hotspots.get('is_adjacent_to_hotspot')):
                    with tabs[tab_idx]:
                        is_at_hotspot = hotspots.get('is_hotspot', False)
                        is_adjacent = hotspots.get('is_adjacent_to_hotspot', False)

                        # Helper to build one-line hotspot summary
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
                    tab_idx += 1

            else:
                st.info("No evidence found from any source")

        # ==============================================
        # EVIDENCE ASSESSMENT (after Evidence by Source) - in bordered card
        # ==============================================
        # Note about match levels - displayed above both Gap Analysis and LLM Synthesis
        st.markdown(
            "<p style='color: #5a5a5a; font-size: 0.95rem; margin: 1rem 0 0.5rem 0;'>"
            "⚠️ <strong>Note:</strong> Some of the evidence may only be a gene- or codon-level match, "
            "but that is not accounted for in the quality ranking or LLM Research Synthesis.</p>",
            unsafe_allow_html=True
        )

        with st.container(border=True):
            evidence_gaps = result.get('evidence_gaps', {})

            # Build badge HTML
            evidence_quality = evidence_gaps.get('overall_quality', 'unknown')
            quality_colors = {"comprehensive": "🟢", "moderate": "🟡", "limited": "🟠", "minimal": "🔴"}
            badge = quality_colors.get(evidence_quality.lower(), "⚪")

            research_priority = evidence_gaps.get('research_priority', 'unknown')
            priority_colors = {"very_high": "🔥", "high": "🔴", "medium": "🟡", "low": "🟢"}
            priority_badge = priority_colors.get(research_priority.lower(), "⚪")
            display_priority = research_priority.replace("_", " ").title()

            # Description above title, then title with badges on the right
            st.caption("What's known vs. unknown about this variant — identifying opportunities for further research.")
            st.markdown(
                f"<div style='display: flex; justify-content: space-between; align-items: flex-start;'>"
                f"<span style='font-size: 1.5rem; font-weight: 600;'>🔍 Gap Analysis</span>"
                f"<span style='font-size: 0.9rem; text-align: right;'>"
                f"<strong>Evidence Quality:</strong> {badge} {evidence_quality.capitalize()} &nbsp;&nbsp; "
                f"<strong>Research Priority:</strong> {priority_badge} {display_priority}"
                f"</span></div>",
                unsafe_allow_html=True
            )

            # Two tables side by side (Well Characterized ~20% wider than Evidence Gaps)
            table_cols = st.columns([7, 4])

            with table_cols[0]:
                st.markdown("**✅ Well Characterized** — _what we know_")

                well_characterized_detailed = evidence_gaps.get('well_characterized_detailed', [])

                # Compute trial match breakdown once using match_scope
                specific_count = len([t for t in trials if t.get('match_scope') == 'specific']) if trials else 0
                ambiguous_count = len([t for t in trials if t.get('match_scope') == 'ambiguous']) if trials else 0
                gene_only_count = len([t for t in trials if t.get('match_scope') not in ('specific', 'ambiguous')]) if trials else 0

                # Build match string for trials
                match_parts = []
                if specific_count > 0:
                    match_parts.append(f"🎯 {specific_count} variant")
                if ambiguous_count > 0:
                    match_parts.append(f"⚠️ {ambiguous_count} broad")
                if gene_only_count > 0:
                    match_parts.append(f"🧬 {gene_only_count} gene")
                trial_match_str = ", ".join(match_parts) if match_parts else ""

                # Compute drug response match breakdown from VICC and CGI
                vicc_variant = len([v for v in vicc if v.get('locus_match') == 'variant']) if vicc else 0
                vicc_codon = len([v for v in vicc if v.get('locus_match') == 'codon']) if vicc else 0
                vicc_gene = len([v for v in vicc if v.get('locus_match') == 'gene']) if vicc else 0
                cgi_variant = len([b for b in cgi_biomarkers if b.get('locus_match') == 'variant']) if cgi_biomarkers else 0
                cgi_codon = len([b for b in cgi_biomarkers if b.get('locus_match') == 'codon']) if cgi_biomarkers else 0
                cgi_gene = len([b for b in cgi_biomarkers if b.get('locus_match') == 'gene']) if cgi_biomarkers else 0

                drug_variant = vicc_variant + cgi_variant
                drug_codon = vicc_codon + cgi_codon
                drug_gene = vicc_gene + cgi_gene

                drug_match_parts = []
                if drug_variant > 0:
                    drug_match_parts.append(f"🎯 {drug_variant} variant")
                if drug_codon > 0:
                    drug_match_parts.append(f"📍 {drug_codon} codon")
                if drug_gene > 0:
                    drug_match_parts.append(f"🧬 {drug_gene} gene")
                drug_match_str = ", ".join(drug_match_parts) if drug_match_parts else ""

                # Compute FDA approval match breakdown
                fda_variant = len([a for a in fda_approvals if a.get('locus_match') == 'variant']) if fda_approvals else 0
                fda_codon = len([a for a in fda_approvals if a.get('locus_match') == 'codon']) if fda_approvals else 0
                fda_gene = len([a for a in fda_approvals if a.get('locus_match') == 'gene']) if fda_approvals else 0

                fda_match_parts = []
                if fda_variant > 0:
                    fda_match_parts.append(f"🎯 {fda_variant} variant")
                if fda_codon > 0:
                    fda_match_parts.append(f"📍 {fda_codon} codon")
                if fda_gene > 0:
                    fda_match_parts.append(f"🧬 {fda_gene} gene")
                fda_match_str = ", ".join(fda_match_parts) if fda_match_parts else ""

                # Helper to convert matches_on value to icon string
                def format_locus_match(matches_on: str) -> str:
                    """Convert matches_on value like 'gene', 'variant', 'codon' or '1 codon, 2 gene' to icon format."""
                    if not matches_on:
                        return ""

                    # Check if it contains commas (multiple match types like "1 codon, 2 gene")
                    if "," in matches_on:
                        formatted_parts = []
                        for part in matches_on.split(","):
                            part = part.strip()
                            parts = part.split()
                            if len(parts) == 2 and parts[0].isdigit():
                                count = parts[0]
                                match_type = parts[1].lower()
                                icon = {"variant": "🎯", "codon": "📍", "gene": "🧬"}.get(match_type, "")
                                formatted_parts.append(f"{icon} {count} {match_type}" if icon else part)
                            else:
                                formatted_parts.append(part)
                        return ", ".join(formatted_parts) if formatted_parts else matches_on

                    # Check for format like "5 variant" (count + type)
                    parts = matches_on.strip().split()
                    if len(parts) == 2 and parts[0].isdigit():
                        count = parts[0]
                        match_type = parts[1].lower()
                        icon = {"variant": "🎯", "codon": "📍", "gene": "🧬"}.get(match_type, "")
                        return f"{icon} {count} {match_type}" if icon else matches_on

                    # Check for simple type like "gene", "variant", "codon"
                    match_type = matches_on.strip().lower()
                    icon = {"variant": "🎯", "codon": "📍", "gene": "🧬"}.get(match_type, "")
                    if icon:
                        return f"{icon} {match_type}"

                    # Already formatted or unknown - return as-is
                    return matches_on

                def format_tumor_match(tumor_match: str) -> str:
                    """Convert tumor_match value like '2 tumor, 1 pan_cancer, 1 other' to icon format."""
                    if not tumor_match:
                        return ""

                    # Handle special "Yes" value (e.g., from cBioPortal prevalence data)
                    if tumor_match.strip().lower() == "yes":
                        return "✅ Yes"

                    # Parse format like "2 tumor, 1 pan_cancer, 1 other" or "1 tumor"
                    formatted_parts = []
                    for part in tumor_match.split(","):
                        part = part.strip()
                        parts = part.split()
                        if len(parts) == 2 and parts[0].isdigit():
                            count = parts[0]
                            match_type = parts[1].lower()
                            # Match icons from Therapies tab: tumor=✅, pan_cancer=🌐, other=⚠️
                            icon = {"tumor": "✅", "pan_cancer": "🌐", "other": "⚠️"}.get(match_type, "")
                            # Display "pan_cancer" as "pan-cancer" for readability
                            display_type = "pan-cancer" if match_type == "pan_cancer" else match_type
                            formatted_parts.append(f"{icon} {count} {display_type}" if icon else part)
                        else:
                            formatted_parts.append(part)

                    return ", ".join(formatted_parts) if formatted_parts else tumor_match

                # Build rows from well_characterized_detailed
                wc_rows = []
                if well_characterized_detailed:
                    for item in well_characterized_detailed:
                        aspect = item.get('aspect', '')
                        basis = item.get('basis', '').lower()

                        # First check if matches_on is already set in the data
                        locus_str = item.get('matches_on', '') or ''

                        # Fall back to computed match string for backwards compatibility
                        if not locus_str:
                            is_trial_row = 'trial' in aspect.lower()
                            is_drug_row = any(term in aspect.lower() for term in ['drug', 'therapy', 'therapeutic', 'treatment', 'response', 'resistance'])
                            is_fda_row = 'fda' in basis or 'actionability' in aspect.lower()

                            if is_trial_row:
                                locus_str = trial_match_str
                            elif is_fda_row:
                                locus_str = fda_match_str
                            elif is_drug_row:
                                locus_str = drug_match_str

                        # Convert simple match types to icon format
                        locus_str = format_locus_match(locus_str)

                        # Tumor match column - first check if data has tumor_match field
                        tumor_match_data = item.get('tumor_match', '')
                        if tumor_match_data:
                            # Use the new structured tumor_match field (e.g., "2 tumor, 1 other")
                            tumor_str = format_tumor_match(tumor_match_data)
                        else:
                            # Fall back to old logic for rows that don't have tumor_match
                            # Show tumor match info if available
                            tumor_match = item.get('tumor_match', '')
                            if tumor_match:
                                # Parse tumor match string like "2 tumor, 1 other"
                                if 'tumor' in tumor_match:
                                    tumor_str = f"✅ {tumor_match}"
                                else:
                                    tumor_str = f"⚠️ {tumor_match}"
                            else:
                                tumor_str = ""  # No tumor match indicator

                        wc_rows.append({
                            "Category": (item.get('category') or '').replace('_', ' ').title(),
                            "Aspect": aspect,
                            "Basis": item.get('basis', ''),
                            "Locus Match": locus_str,
                            "Tumor Match": tumor_str,
                        })

                if wc_rows:
                    # Use HTML table for column width control
                    html_rows = []
                    for row in wc_rows:
                        aspect = row.get("Aspect", "")
                        basis = row.get("Basis", "").replace("<", "&lt;").replace(">", "&gt;")
                        locus = row.get("Locus Match", "")
                        tumor = row.get("Tumor Match", "")
                        html_rows.append(f"<tr><td>{aspect}</td><td>{basis}</td><td>{locus}</td><td>{tumor}</td></tr>")

                    html_table = f"""
                    <style>
                        .wc-table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
                        .wc-table th, .wc-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                        .wc-table th {{ background-color: #f8f9fa; font-weight: 600; }}
                        .wc-table td {{ word-wrap: break-word; }}
                        .wc-table th:nth-child(1), .wc-table td:nth-child(1) {{ width: 25%; }}
                        .wc-table th:nth-child(2), .wc-table td:nth-child(2) {{ width: 25%; }}
                        .wc-table th:nth-child(3), .wc-table td:nth-child(3) {{ width: 25%; }}
                        .wc-table th:nth-child(4), .wc-table td:nth-child(4) {{ width: 25%; }}
                    </style>
                    <table class="wc-table">
                        <thead><tr><th>Aspect</th><th>Basis</th><th>Locus Match</th><th>Tumor Match</th></tr></thead>
                        <tbody>{"".join(html_rows)}</tbody>
                    </table>
                    """
                    st.markdown(html_table, unsafe_allow_html=True)
                else:
                    well_characterized = evidence_gaps.get('well_characterized', [])
                    if well_characterized:
                        wc_df = pd.DataFrame({"Aspect": well_characterized})
                        st.dataframe(wc_df, width="stretch", hide_index=True, height=min(300, 35 * (len(well_characterized) + 1)))

                # Legend below the table
                st.markdown("""<div style='font-size: 0.85rem; margin-top: 0.5rem;'>
<b>Locus Match:</b> 🎯 Variant (exact locus match) · 📍 Codon (other variants in this codon) · 🧬 Gene (other variants in this gene)<br/>
<b>Tumor Match:</b> ✅ Yes (match on specified tumor) · 🔸 Other (match on other tumor) · 🌐 Pan-cancer
</div>""", unsafe_allow_html=True)

            with table_cols[1]:
                gaps = evidence_gaps.get('gaps', [])

                st.markdown("**❓ Evidence Gaps** — _what we don't know_")

                # Build gaps rows
                gaps_data = []
                if gaps:
                    for gap in gaps:
                        severity = gap.get('severity', 'unknown')
                        severity_icon = {
                            "critical": "🔴",
                            "high": "🟠",
                            "significant": "🟡",
                            "moderate": "🔵",
                            "minor": "⚪",
                            "informational": "ℹ️",
                        }.get(severity, "⚪")
                        desc = gap.get('description', '')
                        desc = re.sub(r'\s+for\s+\w+\s+\S+$', '', desc)
                        desc = re.sub(r'\s+of\s+\w+\s+\S+\s+in\s+\S+$', '', desc)
                        desc = re.sub(r'\s+of\s+\w+\s+\S+$', '', desc)
                        category = gap.get('category', '').replace('_', ' ').title()
                        # Check if this is a trial-related gap
                        is_trial_gap = 'trial' in desc.lower() or 'trial' in category.lower()
                        gaps_data.append({
                            "Severity": f"{severity_icon} {severity.capitalize()}",
                            "Category": category,
                            "Description": desc,
                            "Matches On": trial_match_str if is_trial_gap else "",
                        })

                if gaps_data:
                    # Use HTML table for column width control
                    html_rows = []
                    for row in gaps_data:
                        severity = row.get("Severity", "")
                        category = row.get("Category", "")
                        desc = row.get("Description", "").replace("<", "&lt;").replace(">", "&gt;")
                        html_rows.append(f"<tr><td>{severity}</td><td>{category}</td><td>{desc}</td></tr>")

                    html_table = f"""
                    <style>
                        .gaps-table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
                        .gaps-table th, .gaps-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                        .gaps-table th {{ background-color: #f8f9fa; font-weight: 600; }}
                        .gaps-table td {{ word-wrap: break-word; }}
                        .gaps-table th:nth-child(1), .gaps-table td:nth-child(1) {{ width: 25%; white-space: nowrap; }}
                        .gaps-table th:nth-child(2), .gaps-table td:nth-child(2) {{ width: 25%; }}
                        .gaps-table th:nth-child(3), .gaps-table td:nth-child(3) {{ width: 50%; }}
                    </style>
                    <table class="gaps-table">
                        <thead><tr><th>Severity</th><th>Category</th><th>Description</th></tr></thead>
                        <tbody>{"".join(html_rows)}</tbody>
                    </table>
                    """
                    st.markdown(html_table, unsafe_allow_html=True)
                else:
                    st.success("✅ No significant evidence gaps identified for this well-characterized variant.")

            # Divider before Suggested Studies
            st.markdown("<hr style='margin: 1rem 0; border: none; border-top: 1px solid #ddd;'>", unsafe_allow_html=True)

            # Suggested studies (collapsible) - full width below tables
            gaps = evidence_gaps.get('gaps', [])
            if gaps:
                all_suggested = []
                for gap in gaps:
                    all_suggested.extend(gap.get('suggested_studies', []))
                if all_suggested:
                    with st.expander("📋 Suggested Studies"):
                        for study in list(set(all_suggested)):
                            st.markdown(f"- {study}")

        # ==============================================
        # LLM RESEARCH SYNTHESIS (at bottom, if enabled) - in bordered card
        # ==============================================
        render_llm_synthesis(result)

        # ==============================================
        # LLM CROSS-SOURCE DRUG ANALYSIS (after LLM Research Synthesis)
        # ==============================================
        render_cross_source_analysis(result)

        # ==============================================
        # FOOTER: Download & Clear
        # ==============================================
        st.markdown("---")
        footer_cols = st.columns([2, 2, 2, 1])
        with footer_cols[0]:
            st.download_button(
                "📥 Download JSON",
                json.dumps(result, indent=2),
                f"{gene_display}_{variant_display}_insight.json",
                "application/json",
                key="download_single_json"
            )
        with footer_cols[1]:
            st.download_button(
                "📥 Download Markdown",
                result_to_markdown(result),
                f"{gene_display}_{variant_display}_report.md",
                "text/markdown",
                key="download_single_md"
            )
        with footer_cols[2]:
            with st.expander("🔧 Raw JSON"):
                st.json(result)
        with footer_cols[3]:
            if st.button("🗑️ Clear", key="clear_single"):
                st.session_state.single_result = None
                st.session_state.single_gene = None
                st.session_state.single_variant = None
                st.rerun()

# TAB 2: Batch Upload
with tab2:
    render_batch_tab(batch_get_variant_insights, MODELS)

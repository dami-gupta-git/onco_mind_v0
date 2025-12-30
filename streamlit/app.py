"""OncoMind Streamlit Application - Variant insight and evidence synthesis tool."""
import streamlit as st
import pandas as pd
import asyncio
import json
import re
import os
from backend import get_variant_insight, batch_get_variant_insights
from pdb_images import get_pdb_image_url, get_pdb_page_url

# Initialize logging from environment variable (ONCOMIND_LOG_LEVEL=DEBUG|INFO|WARN|ERROR)
# This import triggers the logger initialization which reads from env
from oncomind.config.debug import get_logger
logger = get_logger(__name__)

st.set_page_config(page_title="OncoMind", page_icon="🧬", layout="wide")

# Custom styling
st.markdown("""
<style>
    .block-container {
        padding-top: 1rem;
    }
    /* Make text inputs clearly look like input fields */
    .stTextInput > div > div > input {
        border: 1.5px solid #d0d0d0 !important;
        border-radius: 8px !important;
        padding: 0.5rem 0.75rem !important;
        background-color: #fafafa !important;
        transition: border-color 0.2s, box-shadow 0.2s !important;
    }
    .stTextInput > div > div > input:focus {
        border-color: #ff4b4b !important;
        box-shadow: 0 0 0 2px rgba(255, 75, 75, 0.2) !important;
        background-color: #fff !important;
    }
    .stTextInput > div > div > input:hover {
        border-color: #a0a0a0 !important;
    }
    /* Placeholder text styling */
    .stTextInput > div > div > input::placeholder {
        color: #999 !important;
        font-style: italic !important;
    }
    /* Smaller text in dataframes/tables */
    .stDataFrame [data-testid="stDataFrameResizable"] {
        font-size: 0.8rem !important;
    }
    .stDataFrame table {
        font-size: 0.8rem !important;
    }
    .stDataFrame th {
        font-size: 0.75rem !important;
    }
    .stDataFrame td {
        font-size: 0.8rem !important;
        padding: 4px 8px !important;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("<h1 style='margin-bottom: 0;'><span style='font-size: 0.85em;'>🧬</span> OncoMind: Variant Insight</h1>", unsafe_allow_html=True)
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
EXAMPLE_VARIANTS = {
    "-- Select an example --": ("", "", ""),
    "BRAF V600E - Melanoma": ("BRAF", "V600E", "Melanoma"),
    "EGFR L858R - NSCLC": ("EGFR", "L858R", "NSCLC"),
    "EGFR T790M - NSCLC": ("EGFR", "T790M", "NSCLC"),
    "KRAS G12A - Lung Cancer": ("KRAS", "G12A", "Lung Cancer"),
    "KRAS G12D - Pancreatic": ("KRAS", "G12D", "Pancreatic Adenocarcinoma"),
    "PIK3CA H1047R - Breast Cancer": ("PIK3CA", "H1047R", "Breast Cancer"),
    "AKT1 E17K - Breast Cancer": ("AKT1", "E17K", "Breast Cancer"),
    "TP53 R248W - Multiple": ("TP53", "R248W", ""),
    "TP53 R273H - Multiple": ("TP53", "R273H", ""),
    "GNAQ Q209L - Uveal Melanoma": ("GNAQ", "Q209L", "Uveal Melanoma"),
    "KIT D816V - GIST": ("KIT", "D816V", "GIST"),
    "ALK F1174L - Neuroblastoma": ("ALK", "F1174L", "Neuroblastoma"),
    "ERBB2 S310F - Bladder Cancer": ("ERBB2", "S310F", "Bladder Cancer"),
    "IDH1 R132H - Glioma": ("IDH1", "R132H", "Glioma"),
    "NRAS Q61R - Melanoma": ("NRAS", "Q61R", "Melanoma"),
    "MET exon 14 skip - NSCLC": ("MET", "exon14del", "NSCLC"),
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
        gene = st.text_input("Gene", value="AKT1", placeholder="e.g. AKT1", key="gene_input")
    with input_cols[1]:
        variant = st.text_input("Variant", value="E17K", placeholder="e.g. E17K", key="variant_input")
    with input_cols[2]:
        tumor = st.text_input("Tumor Type", value="Breast Cancer", placeholder="e.g. Breast Cancer", key="tumor_input")
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
            value=False,
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
                temperature = st.slider("Temperature", 0.0, 1.0, 0.1, 0.05)
    else:
        model_name = list(MODELS.keys())[0]
        temperature = 0.1

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
                    result = asyncio.run(get_variant_insight(
                        gene, variant, tumor or None,
                        enable_llm=enable_llm,
                        enable_literature=enable_literature,
                        literature_source=literature_source,
                        model=MODELS[model_name],
                        temperature=temperature
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
        # HEADER: Variant name + Quick Links (no border)
        # ==============================================
        header_text = f"**{gene_display} {variant_display}**"
        if tumor_display:
            header_text += f" in {tumor_display}"

        # Header panel (just variant info, no summary)
        st.markdown(f"<span style='font-size: 1.1rem; font-weight: 600;'>✅ {header_text}</span>", unsafe_allow_html=True)

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
            st.markdown("<span style='font-size: 1.5rem; font-weight: 600;'>📚 Evidence Sources</span>", unsafe_allow_html=True)
            st.caption("Curated clinical and research evidence from cancer knowledge bases, clinical trials, and literature.")

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
            therapies = result.get('recommended_therapies', [])

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
            if cgi_biomarkers:
                tab_names.append(f"CGI ({len(cgi_biomarkers)})")
            if clinvar_entries or clinvar_sig:
                tab_names.append("ClinVar")
            if cosmic_id:
                tab_names.append("COSMIC")
            if trials:
                tab_names.append(f"Trials ({len(trials)})")
            if articles:
                tab_names.append(f"Literature ({len(articles)})")
            if preclinical or early_phase:
                tab_names.append(f"Research ({len(preclinical) + len(early_phase)})")
            if cbioportal:
                tab_names.append("cBioPortal")
            if depmap:
                tab_names.append("🧬 DepMap")
            if therapies:
                tab_names.append(f"Therapies ({len(therapies)})")

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
                                rows.append(f"| CADD | {annotations['cadd_score']:.1f} | {'Deleterious' if annotations['cadd_score'] > 20 else 'Benign'} |")
                            if annotations.get('polyphen2_prediction'):
                                rows.append(f"| PolyPhen2 | - | {annotations['polyphen2_prediction']} |")
                            if annotations.get('gnomad_exome_af') is not None:
                                af = annotations['gnomad_exome_af']
                                freq = f"{af:.2e}" if af < 0.01 else f"{af:.4f}"
                                rows.append(f"| gnomAD AF | {freq} | {'Rare' if af < 0.01 else 'Common'} |")
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
                        # Build match summary for CIViC
                        all_civic = civic_assertions + civic_evidence
                        civic_variant = len([c for c in all_civic if c.get('match_level') == 'variant'])
                        civic_codon = len([c for c in all_civic if c.get('match_level') == 'codon'])
                        civic_gene = len([c for c in all_civic if c.get('match_level') == 'gene'])
                        civic_match_parts = []
                        if civic_variant > 0:
                            civic_match_parts.append(f"🎯 **{civic_variant}** variant")
                        if civic_codon > 0:
                            civic_match_parts.append(f"📍 **{civic_codon}** codon")
                        if civic_gene > 0:
                            civic_match_parts.append(f"🧬 **{civic_gene}** gene")
                        tumor_filter_note = f"filtered by **{tumor_display}**" if tumor_display else ""
                        if civic_match_parts:
                            civic_match_parts.append(tumor_filter_note) if tumor_filter_note else None
                            st.caption(" &nbsp;|&nbsp; ".join([p for p in civic_match_parts if p]))

                        with st.expander("📖 Evidence Level Guide", expanded=False):
                            st.markdown("""
    **AMP/ASCO/CAP Tiers:**
    - **Tier I**: Variants with strong clinical significance (FDA-approved or guideline-recommended)
    - **Tier II**: Variants with potential clinical significance (clinical trials, case studies)
    - **Tier III**: Variants of unknown clinical significance
    - **Tier IV**: Benign or likely benign variants

    **Evidence Levels (A-D):**
    - **A**: FDA-approved therapy or included in professional guidelines
    - **B**: Well-powered studies with consensus
    - **C**: Case studies or small studies
    - **D**: Preclinical or inferential evidence
    """)
                        if civic_assertions:
                            st.markdown("**Curated Assertions:**")
                            rows = ["| ID | Match | Therapies | Significance | Disease | AMP Level |",
                                    "|-----|-------|-----------|--------------|---------|-----------|"]
                            for a in civic_assertions:
                                therapies_str = ", ".join(a.get('therapies', [])) or "N/A"
                                aid = a.get('aid') or a.get('id', '')
                                url = a.get('civic_url', '')
                                id_link = f"[{aid}]({url})" if url else aid
                                disease = (a.get('disease', '') or '')[:35]
                                sig = a.get('significance', 'Unknown')
                                amp = a.get('amp_level', '')
                                # Match level indicator with label
                                match = a.get('match_level', '')
                                match_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "")
                                rows.append(f"| {id_link} | {match_display} | {therapies_str} | {sig} | {disease} | {amp} |")
                            st.markdown("\n".join(rows))

                        if civic_evidence:
                            if civic_assertions:
                                st.markdown("---")
                            st.markdown(f"**Evidence Items ({len(civic_evidence)}):**")
                            # Use markdown table so IDs are clickable
                            rows = ["| ID | Match | Drugs | Significance | Disease | Level | Type |",
                                    "|----|-------|-------|--------------|---------|-------|------|"]
                            for e in civic_evidence[:15]:  # Limit to 15 rows
                                drugs_str = ", ".join(e.get('drugs', [])) or "N/A"
                                drugs_str = drugs_str[:25] if len(drugs_str) > 25 else drugs_str
                                eid = e.get('eid') or ''
                                url = e.get('civic_url', '')
                                id_link = f"[{eid}]({url})" if url else eid
                                disease = (e.get('disease', '') or '')[:20]
                                sig = e.get('clinical_significance', 'Unknown')
                                level = e.get('evidence_level', '')
                                etype = e.get('evidence_type', '')
                                # Match level indicator with label
                                match = e.get('match_level', '')
                                match_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "")
                                rows.append(f"| {id_link} | {match_display} | {drugs_str} | {sig} | {disease} | {level} | {etype} |")
                            st.markdown("\n".join(rows))
                    tab_idx += 1

                # VICC tab
                if vicc:
                    with tabs[tab_idx]:
                        with st.expander("📖 Evidence Level Guide", expanded=False):
                            st.markdown("""
    **Evidence Levels:**
    - **1/A**: FDA-approved or standard of care
    - **2/B**: Clinical trial evidence or expert consensus
    - **3/C**: Case reports or limited evidence
    - **4/D**: Preclinical or computational evidence
    - **R1/R2**: Resistance evidence (strong/emerging)
    """)
                        # Use markdown table for clickable source links
                        rows = ["| Source | Drugs | Response | Disease | Level |",
                                "|--------|-------|----------|---------|-------|"]
                        for v in vicc:
                            source = (v.get('source') or 'vicc').upper()
                            pub_url = v.get('publication_url')
                            # Handle publication_url being a list or string
                            if isinstance(pub_url, list) and pub_url:
                                pub_url = pub_url[0]
                            source_link = f"[{source}]({pub_url})" if pub_url else source
                            drugs = ", ".join(v.get('drugs', [])) or "N/A"
                            drugs = drugs[:30] if len(drugs) > 30 else drugs
                            response = v.get('response_type', 'Unknown')
                            disease = (v.get('disease', '') or '')[:25]
                            level = v.get('evidence_level', '')
                            rows.append(f"| {source_link} | {drugs} | {response} | {disease} | {level} |")
                        st.markdown("\n".join(rows))
                    tab_idx += 1

                # CGI tab
                if cgi_biomarkers:
                    with tabs[tab_idx]:
                        # Use markdown table with links to CGI database
                        rows = ["| Drug | Association | Tumor Type | Level |",
                                "|------|-------------|------------|-------|"]
                        for b in cgi_biomarkers:
                            drug = b.get('drug', 'Unknown')
                            gene = b.get('gene', gene_display)
                            # Link to CGI biomarkers search for this gene
                            cgi_url = f"https://www.cancergenomeinterpreter.org/biomarkers?gene={gene}"
                            drug_link = f"[{drug}]({cgi_url})"
                            association = b.get('association', 'Unknown')
                            tumor = (b.get('tumor_type', '') or '')[:25]
                            level = b.get('evidence_level', '')
                            rows.append(f"| {drug_link} | {association} | {tumor} | {level} |")
                        st.markdown("\n".join(rows))
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
                            st.markdown("\n".join(rows))
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
                        # Count matches by type
                        specific_trials = [t for t in trials if t.get('match_scope') == 'specific']
                        ambiguous_trials = [t for t in trials if t.get('match_scope') == 'ambiguous']
                        gene_only_trials = [t for t in trials if not t.get('variant_specific', False)]

                        # Show counts and tumor filter status
                        tumor_filter_note = f"filtered by **{tumor_display}**" if tumor_display else "across **all cancer types**"
                        count_parts = []
                        if len(specific_trials) > 0:
                            count_parts.append(f"🎯 **{len(specific_trials)}** variant")
                        if len(ambiguous_trials) > 0:
                            count_parts.append(f"⚠️ **{len(ambiguous_trials)}** broad")
                        if len(gene_only_trials) > 0:
                            count_parts.append(f"🧬 **{len(gene_only_trials)}** gene")
                        count_parts.append(tumor_filter_note)
                        st.caption(" &nbsp;|&nbsp; ".join(count_parts))

                        # Use markdown table for clickable NCT IDs
                        rows = ["| Match | NCT ID | Phase | Status | Title |",
                                "|-------|--------|-------|--------|-------|"]
                        for t in trials:
                            matched_biomarker = t.get('matched_biomarker', '')
                            match_scope = t.get('match_scope')

                            # Build display string with icons
                            if match_scope == 'ambiguous':
                                match_display = f"⚠️ {matched_biomarker}" if matched_biomarker else "⚠️ Broad"
                            elif match_scope == 'specific':
                                match_display = f"🎯 {matched_biomarker}" if matched_biomarker else "🎯 Variant"
                            elif t.get('variant_specific', False):
                                match_display = f"🎯 {matched_biomarker}" if matched_biomarker else "🎯 Variant"
                            else:
                                match_display = f"🧬 {matched_biomarker}" if matched_biomarker else "🧬 Gene"

                            nct_id = t.get('nct_id', '')
                            nct_url = t.get('url') or f"https://clinicaltrials.gov/study/{nct_id}" if nct_id else ''
                            nct_link = f"[{nct_id}]({nct_url})" if nct_id and nct_url else nct_id
                            phase = t.get('phase', 'N/A')
                            status = t.get('status', '')
                            title = (t.get('title', '') or '')[:50] + "..."
                            rows.append(f"| {match_display} | {nct_link} | {phase} | {status} | {title} |")
                        st.markdown("\n".join(rows))
                    tab_idx += 1

                # Literature tab
                if articles:
                    with tabs[tab_idx]:
                        # Use markdown table for clickable PMIDs
                        rows = ["| PMID | Year | Journal | Title |",
                                "|------|------|---------|-------|"]
                        for a in articles:
                            pmid = a.get('pmid', '')
                            url = a.get('url') or f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ''
                            pmid_link = f"[{pmid}]({url})" if pmid and url else pmid
                            year = a.get('year', '')
                            journal = (a.get('journal', '') or '')[:20]
                            title = (a.get('title', '') or '')[:50] + "..."
                            rows.append(f"| {pmid_link} | {year} | {journal} | {title} |")
                        st.markdown("\n".join(rows))
                    tab_idx += 1

                # Research tab
                if preclinical or early_phase:
                    with tabs[tab_idx]:
                        st.warning("⚠️ Preclinical/investigational data - not FDA-approved")
                        if preclinical:
                            st.markdown("**Preclinical (Cell Line/Animal Models):**")
                            for b in preclinical:
                                drug = b.get('drug', 'Unknown')
                                assoc = b.get('association', 'Unknown')
                                st.markdown(f"- {drug}: {assoc}")
                        if early_phase:
                            st.markdown("**Early Phase/Case Reports:**")
                            for b in early_phase:
                                drug = b.get('drug', 'Unknown')
                                assoc = b.get('association', 'Unknown')
                                st.markdown(f"- {drug}: {assoc}")
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
                            st.metric(f"Any {gene_symbol} alteration", f"{gene_pct:.1f}%", f"{gene_count}/{total}")
                        with prev_cols[1]:
                            st.metric(f"Exact {variant_name}", f"{variant_pct:.1f}%", f"{variant_count}/{total}")

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
                                st.metric("CERES Score", f"{score:.2f}" if score else "N/A",
                                          f"{dep_pct:.0f}% of cell lines depend on {gene_display}")
                            st.caption(f"Based on CRISPR screens in {n_total} cancer cell lines. CERES < -0.5 indicates essentiality.")

                        if drug_sens:
                            st.markdown("---")
                            st.markdown("### Drug Sensitivities")
                            drug_rows = []
                            for ds in drug_sens:
                                ic50 = ds.get('ic50_nm')
                                ic50_str = f"{ic50:.0f} nM" if ic50 else "N/A"
                                drug_rows.append({
                                    "Drug": ds.get('drug_name', ''),
                                    "IC50": ic50_str,
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
                                    cl_rows.append({
                                        "Cell Line": cl.get('name', ''),
                                        "Disease": cl.get('primary_disease', ''),
                                        "Subtype": cl.get('subtype', ''),
                                        "Mutation": cl.get('mutation_details', variant_display),
                                    })
                                st.dataframe(pd.DataFrame(cl_rows), width="stretch", hide_index=True, height=min(300, 35 * (len(cl_rows) + 1)))
                            else:
                                st.info(f"{len(cell_lines)} cell lines available (mutation status unknown)")

                        st.markdown(f"[🔗 Explore on DepMap Portal](https://depmap.org/portal/gene/{gene_display})")
                    tab_idx += 1

                # Therapies tab
                if therapies:
                    with tabs[tab_idx]:
                        # Build match summary
                        therapy_variant = len([t for t in therapies if t.get('match_level') == 'variant'])
                        therapy_codon = len([t for t in therapies if t.get('match_level') == 'codon'])
                        therapy_gene = len([t for t in therapies if t.get('match_level') == 'gene'])
                        therapy_match_parts = []
                        if therapy_variant > 0:
                            therapy_match_parts.append(f"🎯 **{therapy_variant}** variant")
                        if therapy_codon > 0:
                            therapy_match_parts.append(f"📍 **{therapy_codon}** codon")
                        if therapy_gene > 0:
                            therapy_match_parts.append(f"🧬 **{therapy_gene}** gene")
                        tumor_filter_note = f"filtered by **{tumor_display}**" if tumor_display else ""
                        if therapy_match_parts:
                            therapy_match_parts.append(tumor_filter_note) if tumor_filter_note else None
                            st.caption(" &nbsp;|&nbsp; ".join([p for p in therapy_match_parts if p]))

                        fda_approved = [t for t in therapies if t.get('evidence_level', '').lower() == 'fda-approved']
                        clinical = [t for t in therapies if t.get('evidence_level', '').lower() in ('phase 3', 'phase 2', 'phase 1', 'case report')]
                        preclinical_therapies = [t for t in therapies if t.get('evidence_level', '').lower() in ('preclinical', 'in vitro')]

                        if fda_approved:
                            st.markdown("**✅ FDA-Approved:**")
                            # Use markdown table for clickable links
                            rows = ["| Drug | Locus Match | Tumor Match | Response | Source |",
                                    "|------|-------------|-------------|----------|--------|"]
                            for t in fda_approved:
                                drug = t.get('drug_name', 'Unknown')
                                source_url = t.get('source_url', '')
                                drug_display = f"[{drug}]({source_url})" if source_url else drug
                                response = t.get('response_type', '') or "Sensitivity"
                                source = t.get('source', '')
                                match = t.get('match_level', '')
                                locus_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "-")

                                # Tumor match column
                                cancer_spec = t.get('cancer_specificity', '')
                                if cancer_spec == 'cancer_specific':
                                    tumor_display = "✅ Yes"
                                elif cancer_spec == 'pan_cancer':
                                    tumor_display = "🌐 Pan-cancer"
                                elif cancer_spec:
                                    tumor_display = "⚠️ Other"
                                else:
                                    tumor_display = "-"

                                rows.append(f"| {drug_display} | {locus_display} | {tumor_display} | {response} | {source} |")
                            st.markdown("\n".join(rows))

                        if clinical:
                            if fda_approved:
                                st.markdown("---")
                            st.markdown("**🔬 Clinical Evidence:**")
                            rows = ["| Drug | Locus Match | Tumor Match | Level | Response | Source |",
                                    "|------|-------------|-------------|-------|----------|--------|"]
                            for t in clinical:
                                drug = t.get('drug_name', 'Unknown')
                                source_url = t.get('source_url', '')
                                drug_display = f"[{drug}]({source_url})" if source_url else drug
                                level = t.get('evidence_level', '')
                                response = t.get('response_type', '') or "-"
                                source = t.get('source', '')
                                match = t.get('match_level', '')
                                locus_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "-")

                                # Tumor match column
                                cancer_spec = t.get('cancer_specificity', '')
                                if cancer_spec == 'cancer_specific':
                                    tumor_display = "✅ Yes"
                                elif cancer_spec == 'pan_cancer':
                                    tumor_display = "🌐 Pan-cancer"
                                elif cancer_spec:
                                    tumor_display = "⚠️ Other"
                                else:
                                    tumor_display = "-"

                                rows.append(f"| {drug_display} | {locus_display} | {tumor_display} | {level} | {response} | {source} |")
                            st.markdown("\n".join(rows))

                        if preclinical_therapies:
                            if fda_approved or clinical:
                                st.markdown("---")
                            st.markdown("**🧪 Preclinical:**")
                            st.warning("⚠️ Preclinical data - not validated in humans")
                            rows = ["| Drug | Locus Match | Tumor Match | Response | Source |",
                                    "|------|-------------|-------------|----------|--------|"]
                            for t in preclinical_therapies:
                                drug = t.get('drug_name', 'Unknown')
                                source_url = t.get('source_url', '')
                                drug_display = f"[{drug}]({source_url})" if source_url else drug
                                response = t.get('response_type', '') or "-"
                                source = t.get('source', '')
                                match = t.get('match_level', '')
                                locus_display = {"variant": "🎯 Variant", "codon": "📍 Codon", "gene": "🧬 Gene"}.get(match, "-")

                                # Tumor match column
                                cancer_spec = t.get('cancer_specificity', '')
                                if cancer_spec == 'cancer_specific':
                                    tumor_display = "✅ Yes"
                                elif cancer_spec == 'pan_cancer':
                                    tumor_display = "🌐 Pan-cancer"
                                elif cancer_spec:
                                    tumor_display = "⚠️ Other"
                                else:
                                    tumor_display = "-"

                                rows.append(f"| {drug_display} | {locus_display} | {tumor_display} | {response} | {source} |")
                            st.markdown("\n".join(rows))
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

            # Title with badges on the right
            st.markdown(
                f"<div style='display: flex; justify-content: space-between; align-items: flex-start;'>"
                f"<span style='font-size: 1.5rem; font-weight: 600;'>🔍 Gap Analysis</span>"
                f"<span style='font-size: 0.9rem; text-align: right;'>"
                f"<strong>Evidence Quality:</strong> {badge} {evidence_quality.capitalize()} &nbsp;&nbsp; "
                f"<strong>Research Priority:</strong> {priority_badge} {display_priority}"
                f"</span></div>"
                f"<p style='color: rgba(49, 51, 63, 0.6); font-size: 0.875rem; margin-top: 0.25rem;'>What's known vs. unknown about this variant — identifying opportunities for further research.</p>",
                unsafe_allow_html=True
            )

            # Two tables side by side
            table_cols = st.columns(2)

            with table_cols[0]:
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
                vicc_variant = len([v for v in vicc if v.get('match_level') == 'variant']) if vicc else 0
                vicc_codon = len([v for v in vicc if v.get('match_level') == 'codon']) if vicc else 0
                vicc_gene = len([v for v in vicc if v.get('match_level') == 'gene']) if vicc else 0
                cgi_variant = len([b for b in cgi_biomarkers if b.get('match_level') == 'variant']) if cgi_biomarkers else 0
                cgi_codon = len([b for b in cgi_biomarkers if b.get('match_level') == 'codon']) if cgi_biomarkers else 0
                cgi_gene = len([b for b in cgi_biomarkers if b.get('match_level') == 'gene']) if cgi_biomarkers else 0

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
                fda_variant = len([a for a in fda_approvals if a.get('match_level') == 'variant']) if fda_approvals else 0
                fda_codon = len([a for a in fda_approvals if a.get('match_level') == 'codon']) if fda_approvals else 0
                fda_gene = len([a for a in fda_approvals if a.get('match_level') == 'gene']) if fda_approvals else 0

                fda_match_parts = []
                if fda_variant > 0:
                    fda_match_parts.append(f"🎯 {fda_variant} variant")
                if fda_codon > 0:
                    fda_match_parts.append(f"📍 {fda_codon} codon")
                if fda_gene > 0:
                    fda_match_parts.append(f"🧬 {fda_gene} gene")
                fda_match_str = ", ".join(fda_match_parts) if fda_match_parts else ""

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

                        # Tumor match column - check for cancer mismatch
                        cancer_mismatch = item.get('cancer_mismatch', '')
                        if cancer_mismatch:
                            tumor_str = "⚠️ Other"
                        else:
                            tumor_str = "✅ Yes"

                        wc_rows.append({
                            "Category": (item.get('category') or '').replace('_', ' ').title(),
                            "Aspect": aspect,
                            "Basis": item.get('basis', ''),
                            "Locus Match": locus_str,
                            "Tumor Match": tumor_str,
                        })

                if wc_rows:
                    st.markdown("**✅ Well Characterized** — _what we know_")
                    # Use markdown table for text wrapping
                    md_rows = ["| Aspect | Basis | Locus Match | Tumor Match |",
                               "|--------|-------|-------------|-------------|"]
                    for row in wc_rows:
                        aspect = row.get("Aspect", "")
                        basis = row.get("Basis", "")
                        locus = row.get("Locus Match", "")
                        tumor = row.get("Tumor Match", "")
                        md_rows.append(f"| {aspect} | {basis} | {locus} | {tumor} |")
                    st.markdown("\n".join(md_rows))
                else:
                    well_characterized = evidence_gaps.get('well_characterized', [])
                    if well_characterized:
                        st.markdown("**✅ Well Characterized** — _what we know_")
                        wc_df = pd.DataFrame({"Aspect": well_characterized})
                        st.dataframe(wc_df, width="stretch", hide_index=True, height=min(300, 35 * (len(well_characterized) + 1)))

            with table_cols[1]:
                gaps = evidence_gaps.get('gaps', [])

                # Build gaps rows
                gaps_data = []
                if gaps:
                    for gap in gaps:
                        severity = gap.get('severity', 'unknown')
                        severity_icon = {"critical": "🔴", "significant": "🟠", "minor": "🟡"}.get(severity, "⚪")
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
                    st.markdown("**❓ Evidence Gaps** — _what we don't know_")
                    # Use markdown table for text wrapping
                    md_rows = ["| Severity | Category | Description |",
                               "|----------|----------|-------------|"]
                    for row in gaps_data:
                        severity = row.get("Severity", "")
                        category = row.get("Category", "")
                        desc = row.get("Description", "")
                        md_rows.append(f"| {severity} | {category} | {desc} |")
                    st.markdown("\n".join(md_rows))

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
        # Check for LLM error first
        llm_rationale = result['insight'].get('rationale', '')
        if llm_rationale and llm_rationale.startswith("LLM narrative generation failed:"):
            error_msg = llm_rationale.replace("LLM narrative generation failed: ", "")
            st.warning(f"⚠️ **LLM synthesis failed:** {error_msg}\n\nEvidence-only results are shown above.")

        llm_narrative = result['insight'].get('llm_narrative')
        if llm_narrative and not (llm_rationale and llm_rationale.startswith("LLM narrative generation failed:")):
            with st.container(border=True):
                st.markdown("### 🤖 LLM Research Synthesis")

                functional_summary = result['insight'].get('functional_summary')
                biological_context = result['insight'].get('biological_context')
                therapeutic_landscape = result['insight'].get('therapeutic_landscape')

                if functional_summary:
                    st.markdown(f"**Functional Impact:** {functional_summary}")
                if biological_context:
                    st.markdown(f"**Biological Context:** {biological_context}")
                if therapeutic_landscape:
                    tl_parts = []
                    if therapeutic_landscape.get("fda_approved"):
                        tl_parts.append(f"FDA-approved: {', '.join(therapeutic_landscape['fda_approved'])}")
                    if therapeutic_landscape.get("clinical_evidence"):
                        tl_parts.append(f"Clinical evidence: {', '.join(therapeutic_landscape['clinical_evidence'])}")
                    if therapeutic_landscape.get("preclinical"):
                        tl_parts.append(f"Preclinical: {', '.join(therapeutic_landscape['preclinical'])}")
                    if therapeutic_landscape.get("resistance_mechanisms"):
                        tl_parts.append(f"Resistance: {', '.join(therapeutic_landscape['resistance_mechanisms'])}")
                    if tl_parts:
                        st.markdown(f"**Therapeutic Landscape:** {'; '.join(tl_parts)}")
                    # NOTE: match_level_note removed from LLM output
                    # Match level info is shown in Evidence Specificity panel instead

                if not any([functional_summary, biological_context, therapeutic_landscape]):
                    st.markdown(llm_narrative)

                conflicting_evidence = result['insight'].get('conflicting_evidence', [])
                if conflicting_evidence:
                    st.markdown(f"**⚠️ Conflicting Evidence:** {' · '.join(conflicting_evidence)}")

                evidence_tags = [tag.title() for tag in result['insight'].get('evidence_tags', [])]
                if evidence_tags:
                    st.markdown(f"**🏷️ Evidence Types:** {' · '.join(evidence_tags)}")

                research_implications = result['insight'].get('research_implications')
                if research_implications and research_implications != result['insight'].get('rationale'):
                    st.markdown(f"**🔬 Research Implications:** {research_implications}")

                research_hypotheses = result['insight'].get('research_hypotheses', [])
                if research_hypotheses:
                    st.markdown("**💡 Emerging Research Hypotheses:**")
                    for i, hypothesis in enumerate(research_hypotheses[:3], 1):
                        st.markdown(f"  {i}. {hypothesis}")

                references = result['insight'].get('references', [])
                if references:
                    clickable_refs = []
                    for ref in references[:5]:
                        ref_str = str(ref).strip()
                        if ref_str.startswith("PMID"):
                            pmid = ref_str.replace("PMID", "").replace(":", "").strip()
                            clickable_refs.append(f"[PMID {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
                        elif ref_str.isdigit():
                            clickable_refs.append(f"[PMID {ref_str}](https://pubmed.ncbi.nlm.nih.gov/{ref_str}/)")
                        elif "cBioPortal" in ref_str or "cbioportal" in ref_str.lower():
                            if ":" in ref_str:
                                study_id = ref_str.split(":")[-1].strip()
                                clickable_refs.append(f"[{ref_str}](https://www.cbioportal.org/study/summary?id={study_id})")
                            else:
                                clickable_refs.append(ref_str)
                        elif ref_str.startswith("NCT"):
                            clickable_refs.append(f"[{ref_str}](https://clinicaltrials.gov/study/{ref_str})")
                        else:
                            clickable_refs.append(ref_str)
                    st.markdown(f"**📚 Key References:** {', '.join(clickable_refs)}")

                st.caption("_Synthesis incorporates established domain knowledge beyond queried databases._")

        # ==============================================
        # FOOTER: Download & Clear
        # ==============================================
        st.markdown("---")
        footer_cols = st.columns([2, 2, 1])
        with footer_cols[0]:
            st.download_button(
                "📥 Download JSON",
                json.dumps(result, indent=2),
                f"{gene_display}_{variant_display}_insight.json",
                "application/json",
                key="download_single"
            )
        with footer_cols[1]:
            with st.expander("🔧 Raw JSON"):
                st.json(result)
        with footer_cols[2]:
            if st.button("🗑️ Clear", key="clear_single"):
                st.session_state.single_result = None
                st.session_state.single_gene = None
                st.session_state.single_variant = None
                st.rerun()

# TAB 2: Batch Upload
with tab2:
    st.subheader("Batch Variant Insight")
    st.markdown("**CSV Format:** Must contain `gene`, `variant`, and optionally `tumor_type` columns")

    batch_cols = st.columns([1.5, 1.5, 1])
    with batch_cols[0]:
        enable_literature_batch = st.toggle(
            "Literature",
            value=True,
            help="Search recent literature via Semantic Scholar (with citations). Falls back to PubMed if rate limited.",
            key="batch_literature"
        )
        literature_source_value_batch = "semantic_scholar" if enable_literature_batch else "none"
    with batch_cols[1]:
        enable_llm_batch = st.toggle(
            "LLM Mode",
            value=False,
            help="LLM mode includes AI synthesis (~25s/variant). Without LLM: fast annotation (~7s/variant).",
            key="batch_llm"
        )
    with batch_cols[2]:
        if enable_llm_batch:
            model_name_batch = st.selectbox("Model", list(MODELS.keys()), key="batch_model")
        else:
            model_name_batch = list(MODELS.keys())[0]

    uploaded_file = st.file_uploader("Upload CSV", type=['csv'])
    if uploaded_file:
        df = pd.read_csv(uploaded_file)
        st.dataframe(df.head(), width="stretch")
        if st.button("🚀 Get Batch Insights", type="primary"):
            if 'gene' not in df.columns or 'variant' not in df.columns:
                st.error("CSV must contain 'gene' and 'variant' columns")
            else:
                progress_bar = st.progress(0)
                status_text = st.empty()
                variants = [{"gene": row.get('gene'), "variant": row.get('variant'),
                           "tumor_type": row.get('tumor_type', None)} for _, row in df.iterrows()]
                results = asyncio.run(batch_get_variant_insights(
                    variants=variants,
                    enable_llm=enable_llm_batch,
                    enable_literature=enable_literature_batch,
                    literature_source=literature_source_value_batch,
                    model=MODELS[model_name_batch],
                    temperature=0.1,
                    progress_callback=lambda i, t: (progress_bar.progress(i/t), status_text.text(f"Processing {i}/{t}..."))
                ))
                status_text.text("✅ Batch processing complete!")
                progress_bar.progress(1.0)

                st.session_state.batch_results = results

                # Check for LLM errors in batch results
                llm_errors = []
                for r in results:
                    if 'error' not in r:
                        rationale = r.get('insight', {}).get('rationale', '')
                        if rationale and rationale.startswith("LLM narrative generation failed:"):
                            error_msg = rationale.replace("LLM narrative generation failed: ", "")
                            llm_errors.append(f"{r['variant']['gene']} {r['variant']['variant']}: {error_msg}")

                if llm_errors:
                    st.warning(f"⚠️ **LLM errors ({len(llm_errors)}):**\n" + "\n".join([f"- {e}" for e in llm_errors]))

                results_df = pd.DataFrame([{"Gene": r['variant']['gene'], "Variant": r['variant']['variant'],
                    "Tumor": r['variant'].get('tumor_type', 'N/A'),
                    "Therapies": len(r.get('recommended_therapies', []))} for r in results if 'error' not in r])
                st.session_state.batch_df = results_df

    if st.session_state.batch_results is not None:
        st.dataframe(st.session_state.batch_df, width="stretch")
        batch_footer_cols = st.columns(3)
        with batch_footer_cols[0]:
            st.download_button(
                "📥 Download Results CSV",
                st.session_state.batch_df.to_csv(index=False),
                "batch_results.csv",
                "text/csv",
                key="download_batch_csv"
            )
        with batch_footer_cols[1]:
            st.download_button(
                "📥 Download Full JSON",
                json.dumps(st.session_state.batch_results, indent=2),
                "batch_results.json",
                "application/json",
                key="download_batch_json"
            )
        with batch_footer_cols[2]:
            if st.button("🗑️ Clear Results", key="clear_batch"):
                st.session_state.batch_results = None
                st.session_state.batch_df = None
                st.rerun()

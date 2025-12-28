"""OncoMind Streamlit Application - Variant insight and evidence synthesis tool."""
import streamlit as st
import pandas as pd
import asyncio
import json
import re
from backend import get_variant_insight, batch_get_variant_insights
from pdb_images import get_pdb_image_url, get_pdb_page_url

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
    "OpenAI GPT-4o-mini": "gpt-4o-mini",
    "OpenAI GPT-4o": "gpt-4o",
    "Anthropic Claude 3 Haiku": "claude-3-haiku-20240307",
    "Google Gemini 1.5 Pro": "gemini/gemini-1.5-pro",
    "Groq Llama 3.1 70B": "groq/llama-3.1-70b-versatile"
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
    input_cols = st.columns([1.5, 1.5, 2, 2, 1.5])
    with input_cols[0]:
        gene = st.text_input("Gene", value="AKT1", placeholder="e.g. AKT1", key="gene_input")
    with input_cols[1]:
        variant = st.text_input("Variant", value="E17K", placeholder="e.g. E17K", key="variant_input")
    with input_cols[2]:
        tumor = st.text_input("Tumor Type", value="Breast Cancer", placeholder="e.g. Breast Cancer", key="tumor_input")
    with input_cols[3]:
        st.markdown("<div style='height: 28px'></div>", unsafe_allow_html=True)  # Spacer to align with labels
        enable_llm = st.toggle(
            "LLM Mode",
            value=False,
            help="LLM mode includes literature search and AI-powered synthesis (~25s). Without LLM, you get fast annotation (~7s)."
        )
        enable_literature = enable_llm
    with input_cols[4]:
        st.markdown("<div style='height: 28px'></div>", unsafe_allow_html=True)  # Spacer to align with labels
        insight_btn = st.button("🔍 Get Insight", type="primary", use_container_width=True)

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
                            rows = ["| ID | Therapies | Significance | Disease | AMP Level |",
                                    "|-----|-----------|--------------|---------|-----------|"]
                            for a in civic_assertions:
                                therapies_str = ", ".join(a.get('therapies', [])) or "N/A"
                                aid = a.get('aid') or a.get('id', '')
                                url = a.get('civic_url', '')
                                id_link = f"[{aid}]({url})" if url else aid
                                disease = (a.get('disease', '') or '')[:40]
                                sig = a.get('significance', 'Unknown')
                                amp = a.get('amp_level', '')
                                rows.append(f"| {id_link} | {therapies_str} | {sig} | {disease} | {amp} |")
                            st.markdown("\n".join(rows))

                        if civic_evidence:
                            if civic_assertions:
                                st.markdown("---")
                            st.markdown(f"**Evidence Items ({len(civic_evidence)}):**")
                            evidence_rows = []
                            for e in civic_evidence:
                                drugs_str = ", ".join(e.get('drugs', [])) or "N/A"
                                drugs_str = drugs_str[:30] if len(drugs_str) > 30 else drugs_str
                                eid = e.get('eid') or ''
                                url = e.get('civic_url', '')
                                evidence_rows.append({
                                    "ID": eid,
                                    "Drugs": drugs_str,
                                    "Significance": e.get('clinical_significance', 'Unknown'),
                                    "Disease": (e.get('disease', '') or '')[:25],
                                    "Level": e.get('evidence_level', ''),
                                    "Type": e.get('evidence_type', ''),
                                    "Rating": e.get('trust_rating') or e.get('rating') or '',
                                    "_url": url,
                                })
                            evidence_df = pd.DataFrame(evidence_rows)
                            # Make ID clickable via column config
                            st.dataframe(
                                evidence_df.drop(columns=['_url']),
                                width="stretch",
                                hide_index=True,
                                height=min(300, 35 * (len(evidence_rows) + 1))  # ~8 rows visible
                            )
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
                        vicc_rows = []
                        for v in vicc:
                            drugs = ", ".join(v.get('drugs', [])) or "N/A"
                            drugs = drugs[:30] if len(drugs) > 30 else drugs
                            vicc_rows.append({
                                "Source": (v.get('source') or 'vicc').upper(),
                                "Drugs": drugs,
                                "Response": v.get('response_type', 'Unknown'),
                                "Disease": (v.get('disease', '') or '')[:25],
                                "Level": v.get('evidence_level', ''),
                            })
                        st.dataframe(
                            pd.DataFrame(vicc_rows),
                            width="stretch",
                            hide_index=True,
                            height=min(300, 35 * (len(vicc_rows) + 1))
                        )
                    tab_idx += 1

                # CGI tab
                if cgi_biomarkers:
                    with tabs[tab_idx]:
                        cgi_rows = []
                        for b in cgi_biomarkers:
                            cgi_rows.append({
                                "Drug": b.get('drug', 'Unknown'),
                                "Association": b.get('association', 'Unknown'),
                                "Tumor Type": (b.get('tumor_type', '') or '')[:25],
                                "Level": b.get('evidence_level', ''),
                            })
                        if cgi_rows:
                            st.dataframe(
                                pd.DataFrame(cgi_rows),
                                width="stretch",
                                hide_index=True,
                                height=min(300, 35 * (len(cgi_rows) + 1))
                            )
                    tab_idx += 1

                # ClinVar tab
                if clinvar_entries or clinvar_sig:
                    with tabs[tab_idx]:
                        if clinvar_sig:
                            st.markdown(f"**Clinical Significance:** {clinvar_sig}")
                        for entry in clinvar_entries:
                            sig = entry.get('clinical_significance', 'Unknown')
                            conds = entry.get('conditions', [])
                            review = entry.get('review_status', '')
                            st.markdown(f"- {sig}: {', '.join(conds) if conds else 'N/A'} ({review})")
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
                        # Count matches - all trials match on gene, some also match on variant
                        variant_specific_trials = [t for t in trials if t.get('variant_specific', False)]
                        gene_only_trials = [t for t in trials if not t.get('variant_specific', False)]

                        trial_rows = []
                        for t in trials:
                            is_variant_specific = t.get('variant_specific', False)
                            match_type = "🎯 Variant" if is_variant_specific else "🧬 Gene"
                            trial_rows.append({
                                "Matches On": match_type,
                                "NCT ID": t.get('nct_id', ''),
                                "Phase": t.get('phase', 'N/A'),
                                "Status": t.get('status', ''),
                                "Title": (t.get('title', '') or '')[:50] + "...",
                            })

                        # Show counts and tumor filter status
                        tumor_filter_note = f"filtered by **{tumor_display}**" if tumor_display else "across **all cancer types**"
                        count_parts = [f"🎯 **{len(variant_specific_trials)}** variant"]
                        if len(gene_only_trials) > 0:
                            count_parts.append(f"🧬 **{len(gene_only_trials)}** gene")
                        count_parts.append(tumor_filter_note)
                        st.caption(" &nbsp;|&nbsp; ".join(count_parts))

                        st.dataframe(
                            pd.DataFrame(trial_rows),
                            width="stretch",
                            hide_index=True,
                            height=min(300, 35 * (len(trial_rows) + 1))
                        )
                    tab_idx += 1

                # Literature tab
                if articles:
                    with tabs[tab_idx]:
                        article_rows = []
                        for a in articles:
                            article_rows.append({
                                "PMID": a.get('pmid', ''),
                                "Year": a.get('year', ''),
                                "Journal": (a.get('journal', '') or '')[:15],
                                "Title": (a.get('title', '') or '')[:50] + "...",
                            })
                        st.dataframe(
                            pd.DataFrame(article_rows),
                            width="stretch",
                            hide_index=True,
                            height=min(300, 35 * (len(article_rows) + 1))
                        )
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
                                st.success(f"✅ {len(mutant_lines)} cell lines with {variant_display} mutation")
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
                        fda_approved = [t for t in therapies if t.get('evidence_level', '').lower() == 'fda-approved']
                        clinical = [t for t in therapies if t.get('evidence_level', '').lower() in ('phase 3', 'phase 2', 'phase 1', 'case report')]
                        preclinical_therapies = [t for t in therapies if t.get('evidence_level', '').lower() in ('preclinical', 'in vitro')]

                        if fda_approved:
                            st.markdown("**✅ FDA-Approved:**")
                            fda_rows = []
                            for t in fda_approved:
                                drug = t.get('drug_name', 'Unknown')
                                response = t.get('response_type', '')
                                source = t.get('source', '')
                                context = t.get('clinical_context', '') or ''
                                context = context[:40] + '...' if len(context) > 40 else context
                                fda_rows.append({
                                    "Drug": drug,
                                    "Response": response or "Sensitivity",
                                    "Context": context,
                                    "Source": source,
                                })
                            st.dataframe(pd.DataFrame(fda_rows), width="stretch", hide_index=True, height=min(300, 35 * (len(fda_rows) + 1)))

                        if clinical:
                            if fda_approved:
                                st.markdown("---")
                            st.markdown("**🔬 Clinical Evidence:**")
                            clinical_rows = []
                            for t in clinical:
                                clinical_rows.append({
                                    "Drug": t.get('drug_name', 'Unknown'),
                                    "Level": t.get('evidence_level', ''),
                                    "Response": t.get('response_type', '') or "-",
                                    "Source": t.get('source', ''),
                                })
                            st.dataframe(pd.DataFrame(clinical_rows), width="stretch", hide_index=True, height=min(300, 35 * (len(clinical_rows) + 1)))

                        if preclinical_therapies:
                            if fda_approved or clinical:
                                st.markdown("---")
                            st.markdown("**🧪 Preclinical:**")
                            st.warning("⚠️ Preclinical data - not validated in humans")
                            preclin_rows = []
                            for t in preclinical_therapies:
                                preclin_rows.append({
                                    "Drug": t.get('drug_name', 'Unknown'),
                                    "Response": t.get('response_type', '') or "-",
                                    "Source": t.get('source', ''),
                                })
                            st.dataframe(pd.DataFrame(preclin_rows), width="stretch", hide_index=True, height=min(300, 35 * (len(preclin_rows) + 1)))
                    tab_idx += 1
            else:
                st.info("No evidence found from any source")

        # ==============================================
        # EVIDENCE ASSESSMENT (after Evidence by Source) - in bordered card
        # ==============================================
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
                f"<div style='display: flex; justify-content: space-between; align-items: center;'>"
                f"<span style='font-size: 1.5rem; font-weight: 600;'>🔍 Gap Analysis</span>"
                f"<span style='font-size: 0.9rem;'>"
                f"<strong>Evidence Quality:</strong> {badge} {evidence_quality.capitalize()} &nbsp;&nbsp; "
                f"<strong>Research Priority:</strong> {priority_badge} {display_priority}"
                f"</span></div>",
                unsafe_allow_html=True
            )
            st.caption("What's known vs. unknown about this variant — identifying opportunities for further research.")

            # Two tables side by side
            table_cols = st.columns(2)

            with table_cols[0]:
                well_characterized_detailed = evidence_gaps.get('well_characterized_detailed', [])

                # Compute trial match breakdown once
                variant_specific_count = len([t for t in trials if t.get('variant_specific', False)]) if trials else 0
                gene_only_count = len([t for t in trials if not t.get('variant_specific', False)]) if trials else 0

                # Build match string for trials
                match_parts = []
                if variant_specific_count > 0:
                    match_parts.append(f"🎯 {variant_specific_count} variant")
                if gene_only_count > 0:
                    match_parts.append(f"🧬 {gene_only_count} gene")
                trial_match_str = ", ".join(match_parts) if match_parts else ""

                # Build rows from well_characterized_detailed
                wc_rows = []
                if well_characterized_detailed:
                    for item in well_characterized_detailed:
                        aspect = item.get('aspect', '')
                        # Check if this row is about clinical trials
                        is_trial_row = 'trial' in aspect.lower()
                        wc_rows.append({
                            "Category": (item.get('category') or '').replace('_', ' ').title(),
                            "Aspect": aspect,
                            "Basis": item.get('basis', ''),
                            "Matches On": trial_match_str if is_trial_row else "",
                        })

                if wc_rows:
                    st.markdown("**✅ Well Characterized** — _what we know_")
                    wc_df = pd.DataFrame(wc_rows)
                    st.dataframe(wc_df, use_container_width=True, hide_index=True, height=min(300, 35 * (len(wc_rows) + 1)))
                else:
                    well_characterized = evidence_gaps.get('well_characterized', [])
                    if well_characterized:
                        st.markdown("**✅ Well Characterized** — _what we know_")
                        wc_df = pd.DataFrame({"Aspect": well_characterized})
                        st.dataframe(wc_df, use_container_width=True, hide_index=True, height=min(300, 35 * (len(well_characterized) + 1)))

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
                    st.dataframe(pd.DataFrame(gaps_data), use_container_width=True, hide_index=True, height=min(300, 35 * (len(gaps_data) + 1)))

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
        llm_narrative = result['insight'].get('llm_narrative')
        if llm_narrative:
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

    batch_cols = st.columns([2, 1])
    with batch_cols[0]:
        enable_llm_batch = st.toggle(
            "🤖 Enable LLM Mode",
            value=False,
            help="LLM mode includes literature search + AI synthesis (~25s/variant). Without LLM: fast annotation (~7s/variant).",
            key="batch_llm"
        )
        enable_literature_batch = enable_llm_batch
    with batch_cols[1]:
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
                    variants, MODELS[model_name_batch], 0.1,
                    lambda i, t: (progress_bar.progress(i/t), status_text.text(f"Processing {i}/{t}...")),
                    enable_llm=enable_llm_batch,
                    enable_literature=enable_literature_batch
                ))
                status_text.text("✅ Batch processing complete!")
                progress_bar.progress(1.0)

                st.session_state.batch_results = results
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

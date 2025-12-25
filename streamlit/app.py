"""OncoMind Streamlit Application - Variant insight and evidence synthesis tool."""
import streamlit as st
import pandas as pd
import asyncio
import json
from backend import get_variant_insight, batch_get_variant_insights

st.set_page_config(page_title="OncoMind", page_icon="🧬", layout="wide")
st.title("🧬 OncoMind: Variant Insight")
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
    col1, col2 = st.columns([1, 2])
    with col1:
        st.subheader("Input")

        gene = st.text_input("Gene Symbol", value="BRAF", placeholder="e.g., BRAF, EGFR, TP53", key="gene_input")
        variant = st.text_input("Variant", value="V600E", placeholder="e.g., V600E, L858R", key="variant_input")
        tumor = st.text_input("Tumor Type (optional)", value="Melanoma", placeholder="e.g., Melanoma, NSCLC", key="tumor_input")

        st.subheader("Mode")
        mode = st.radio(
            "Analysis Mode",
            options=["Default", "Lite", "Full"],
            index=0,
            horizontal=True,
            help="Lite: ~7s, no LLM | Default: ~12s, with LLM | Full: ~25s, with literature"
        )

        # Mode descriptions
        mode_info = {
            "Lite": "⚡ **Lite** (~7s): Structured evidence only, no LLM narrative",
            "Default": "🎯 **Default** (~12s): Structured evidence + LLM clinical summary",
            "Full": "📚 **Full** (~25s): + Literature search + enhanced narrative"
        }
        st.caption(mode_info[mode])

        # Derive settings from mode
        enable_llm = mode in ["Default", "Full"]
        enable_literature = mode == "Full"

        # Only show LLM options if LLM is enabled
        if enable_llm:
            with st.expander("LLM Settings"):
                model_name = st.selectbox("LLM Model", list(MODELS.keys()))
                temperature = st.slider("Temperature", 0.0, 1.0, 0.1, 0.05)
        else:
            model_name = list(MODELS.keys())[0]
            temperature = 0.1

        insight_btn = st.button("🔍 Get Insight", type="primary", use_container_width=True)

    with col2:
        # Run analysis if button clicked
        if insight_btn:
            if not gene or not variant:
                st.error("Gene and variant are required")
            else:
                # Validate variant type before processing
                from oncomind.utils.variant_normalization import normalize_variant, VariantNormalizer
                normalized = normalize_variant(gene, variant)
                variant_type = normalized['variant_type']

                if variant_type not in VariantNormalizer.ALLOWED_VARIANT_TYPES:
                    st.error(
                        f"❌ Unsupported variant type: **{variant_type}**\n\n"
                        f"Only **SNPs and small indels** are supported:\n"
                        f"- Missense mutations (e.g., V600E)\n"
                        f"- Nonsense mutations (e.g., R172*)\n"
                        f"- Small insertions (e.g., ins)\n"
                        f"- Small deletions (e.g., del)\n"
                        f"- Frameshift mutations (e.g., fs)\n\n"
                        f"Your variant '{variant}' is classified as '{variant_type}'."
                    )
                else:
                    with st.spinner(f"🔬 Getting insight for {gene} {variant} ({mode} mode)..."):
                        # Use fast annotation API or legacy insight engine
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
                            # Store in session state to persist across reruns
                            st.session_state.single_result = result
                            st.session_state.single_gene = gene
                            st.session_state.single_variant = variant

        # Display results from session state (persists across reruns)
        if st.session_state.single_result is not None:
            result = st.session_state.single_result
            gene_display = st.session_state.single_gene
            variant_display = st.session_state.single_variant

            tumor_display = result.get('variant', {}).get('tumor_type')
            header_text = f"**{gene_display} {variant_display}**"
            if tumor_display:
                header_text += f" in {tumor_display}"

            st.success(f"✅ Insight Ready: {header_text}")

            # Card 1: Summary (always expanded) - 1-line summary
            with st.expander("📋 Summary", expanded=True):
                st.info(result['insight'].get('summary', 'No summary available'))

            # Card 2: Identifiers & Links
            with st.expander("🔗 Identifiers & External Links", expanded=False):
                ids = result.get('identifiers', {})
                hgvs = result.get('hgvs', {})

                col1, col2 = st.columns(2)
                with col1:
                    st.markdown("**Database IDs:**")
                    if ids.get('cosmic_id'):
                        cosmic_id = ids['cosmic_id']
                        st.markdown(f"- COSMIC: [{cosmic_id}](https://cancer.sanger.ac.uk/cosmic/mutation/overview?id={cosmic_id.replace('COSM', '').replace('COSV', '')})")
                    if ids.get('dbsnp_id'):
                        dbsnp = ids['dbsnp_id']
                        st.markdown(f"- dbSNP: [{dbsnp}](https://www.ncbi.nlm.nih.gov/snp/{dbsnp})")
                    if ids.get('clinvar_id'):
                        clinvar_id = ids['clinvar_id']
                        st.markdown(f"- ClinVar: [{clinvar_id}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/)")
                    if ids.get('ncbi_gene_id'):
                        ncbi_id = ids['ncbi_gene_id']
                        st.markdown(f"- NCBI Gene: [{ncbi_id}](https://www.ncbi.nlm.nih.gov/gene/{ncbi_id})")

                with col2:
                    st.markdown("**HGVS Notation:**")
                    if hgvs.get('protein'):
                        st.markdown(f"- Protein: `{hgvs['protein']}`")
                    if hgvs.get('transcript'):
                        st.markdown(f"- Transcript: `{hgvs['transcript']}`")
                    if hgvs.get('genomic'):
                        st.markdown(f"- Genomic: `{hgvs['genomic']}`")

                # External search links
                st.markdown("**Search:**")
                search_cols = st.columns(4)
                search_cols[0].markdown(f"[PubMed](https://pubmed.ncbi.nlm.nih.gov/?term={gene_display}+{variant_display})")
                search_cols[1].markdown(f"[CIViC](https://civicdb.org/variants?geneSearch={gene_display})")
                search_cols[2].markdown(f"[OncoKB](https://www.oncokb.org/gene/{gene_display})")
                search_cols[3].markdown(f"[cBioPortal](https://www.cbioportal.org/results?cancer_study_list=msk_impact_2017&gene_list={gene_display})")

            # ==============================================
            # SOURCE CARDS - Tabbed structure (flat lists from backend)
            # ==============================================
            st.markdown("---")
            st.markdown("### 📊 Evidence by Source")

            # Get data from flat lists
            ids = result.get('identifiers', {})

            # Collect available sources with counts (flat list structure)
            fda_approvals = result.get('fda_approvals', [])
            civic_assertions = result.get('civic_assertions', [])
            vicc = result.get('vicc_evidence', [])
            cgi_biomarkers = result.get('cgi_biomarkers', [])
            clinvar_entries = result.get('clinvar_entries', [])
            clinvar_sig = result.get('clinvar', {}).get('clinical_significance')
            cosmic_id = ids.get('cosmic_id')
            trials = result.get('clinical_trials', [])
            articles = result.get('pubmed_articles', [])
            preclinical = result.get('preclinical_biomarkers', [])
            early_phase = result.get('early_phase_biomarkers', [])

            # Build tab names with counts
            tab_names = []
            if fda_approvals:
                tab_names.append(f"FDA ({len(fda_approvals)})")
            if civic_assertions:
                tab_names.append(f"CIViC ({len(civic_assertions)})")
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

            if tab_names:
                tabs = st.tabs(tab_names)
                tab_idx = 0

                # FDA tab
                if fda_approvals:
                    with tabs[tab_idx]:
                        rows = []
                        for a in fda_approvals:
                            drug = a.get('brand_name') or a.get('drug_name', 'Unknown')
                            indication = (a.get('indication') or '')[:100]
                            rows.append(f"| {drug} | {indication}... |")
                        st.markdown("| Drug | Indication |\n|------|------------|" + "\n" + "\n".join(rows))
                    tab_idx += 1

                # CIViC tab
                if civic_assertions:
                    with tabs[tab_idx]:
                        for a in civic_assertions:
                            therapies_str = ", ".join(a.get('therapies', [])) or "N/A"
                            sig = a.get('significance', 'Unknown')
                            amp = a.get('amp_level', '')
                            disease = a.get('disease', '')
                            aid = a.get('id', '')
                            link = f"[AID:{aid}](https://civicdb.org/assertions/{aid})" if aid else ""
                            st.markdown(f"- {link} **{therapies_str}** → {sig} ({disease}) [{amp}]")
                    tab_idx += 1

                # VICC tab
                if vicc:
                    with tabs[tab_idx]:
                        sources = {}
                        for v in vicc:
                            src = v.get('source', 'unknown')
                            if src not in sources:
                                sources[src] = []
                            sources[src].append(v)
                        for source, entries in sources.items():
                            st.markdown(f"**{source.upper()}** ({len(entries)})")
                            for v in entries[:3]:
                                drugs = ", ".join(v.get('drugs', [])) or "N/A"
                                response = v.get('response_type', 'Unknown')
                                disease = v.get('disease', '')[:30]
                                st.markdown(f"- {drugs} → {response} ({disease})")
                    tab_idx += 1

                # CGI tab
                if cgi_biomarkers:
                    with tabs[tab_idx]:
                        for b in cgi_biomarkers:
                            drug = b.get('drug', 'Unknown')
                            assoc = b.get('association', 'Unknown')
                            tumor = b.get('tumor_type', '')
                            level = b.get('evidence_level', '')
                            st.markdown(f"- **{drug}**: {assoc} in {tumor} ({level})")
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

                # Clinical Trials tab
                if trials:
                    with tabs[tab_idx]:
                        rows = []
                        for t in trials:
                            nct = t.get('nct_id', '')
                            phase = t.get('phase', 'N/A')
                            status = t.get('status', '')
                            title = (t.get('title', '') or '')[:50]
                            link = f"[{nct}](https://clinicaltrials.gov/study/{nct})"
                            rows.append(f"| {link} | {phase} | {status} | {title}... |")
                        st.markdown("| NCT ID | Phase | Status | Title |\n|--------|-------|--------|-------|" + "\n" + "\n".join(rows))
                    tab_idx += 1

                # Literature tab
                if articles:
                    with tabs[tab_idx]:
                        rows = []
                        for a in articles:
                            pmid = a.get('pmid', '')
                            year = a.get('year', '')
                            title = (a.get('title', '') or '')[:50]
                            journal = (a.get('journal', '') or '')[:15]
                            link = f"[{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
                            rows.append(f"| {link} | {year} | {journal} | {title}... |")
                        st.markdown("| PMID | Year | Journal | Title |\n|------|------|---------|-------|" + "\n" + "\n".join(rows))
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
            else:
                st.info("No evidence found from any source")

            # ==============================================
            # LLM CARDS - After evidence cards
            # ==============================================
            llm_narrative = result['insight'].get('llm_narrative')
            therapies = result.get('recommended_therapies', [])

            if llm_narrative or therapies:
                st.markdown("---")
                st.markdown("### 🤖 LLM Analysis")

                # LLM Insight card
                if llm_narrative:
                    with st.expander("🤖 LLM Insight", expanded=True):
                        st.markdown(llm_narrative)

                # Recommended Therapies card
                if therapies:
                    with st.expander(f"💊 Recommended Therapies ({len(therapies)})", expanded=True):
                        table_rows = []
                        for t in therapies:
                            drug = t.get('drug_name', 'Unknown')
                            drug_link = f"[{drug}](https://www.drugs.com/search.php?searchterm={drug.replace(' ', '+')})"
                            table_rows.append(f"| {drug_link} | {t.get('evidence_level', 'N/A')} | {t.get('approval_status', '')} |")
                        table_header = "| Drug | Level | Status |\n|------|-------|--------|"
                        st.markdown(table_header + "\n" + "\n".join(table_rows))

            # Download and Clear buttons
            st.markdown("---")
            col1, col2, col3 = st.columns([2, 2, 1])
            with col1:
                st.download_button(
                    "📥 Download JSON",
                    json.dumps(result, indent=2),
                    f"{gene_display}_{variant_display}_insight.json",
                    "application/json",
                    key="download_single"
                )
            with col2:
                # Raw JSON expander
                with st.expander("🔧 Raw JSON"):
                    st.json(result)
            with col3:
                if st.button("🗑️ Clear", key="clear_single"):
                    st.session_state.single_result = None
                    st.session_state.single_gene = None
                    st.session_state.single_variant = None
                    st.rerun()

# TAB 2: Batch Upload
with tab2:
    st.subheader("Batch Variant Insight")
    st.markdown("**CSV Format:** Must contain `gene`, `variant`, and optionally `tumor_type` columns")

    col1, col2 = st.columns([2, 1])
    with col1:
        mode_batch = st.radio(
            "Analysis Mode",
            options=["Default", "Lite", "Full"],
            index=0,
            horizontal=True,
            help="Lite: ~7s/variant | Default: ~12s/variant | Full: ~25s/variant",
            key="batch_mode"
        )
        # Derive settings from mode
        enable_llm_batch = mode_batch in ["Default", "Full"]
        enable_literature_batch = mode_batch == "Full"
    with col2:
        if enable_llm_batch:
            model_name_batch = st.selectbox("Model", list(MODELS.keys()), key="batch_model")
        else:
            model_name_batch = list(MODELS.keys())[0]

    uploaded_file = st.file_uploader("Upload CSV", type=['csv'])
    if uploaded_file:
        df = pd.read_csv(uploaded_file)
        st.dataframe(df.head(), use_container_width=True)
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

                # Store in session state
                st.session_state.batch_results = results
                results_df = pd.DataFrame([{"Gene": r['variant']['gene'], "Variant": r['variant']['variant'],
                    "Tumor": r['variant'].get('tumor_type', 'N/A'),
                    "Therapies": len(r.get('recommended_therapies', []))} for r in results if 'error' not in r])
                st.session_state.batch_df = results_df

    # Display batch results from session state
    if st.session_state.batch_results is not None:
        st.dataframe(st.session_state.batch_df, use_container_width=True)
        col1, col2, col3 = st.columns(3)
        with col1:
            st.download_button(
                "📥 Download Results CSV",
                st.session_state.batch_df.to_csv(index=False),
                "batch_results.csv",
                "text/csv",
                key="download_batch_csv"
            )
        with col2:
            st.download_button(
                "📥 Download Full JSON",
                json.dumps(st.session_state.batch_results, indent=2),
                "batch_results.json",
                "application/json",
                key="download_batch_json"
            )
        with col3:
            if st.button("🗑️ Clear Results", key="clear_batch"):
                st.session_state.batch_results = None
                st.session_state.batch_df = None
                st.rerun()

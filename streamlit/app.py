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
            panel = result.get('evidence_panel', {})

            st.success(f"✅ Insight Ready: **{gene_display} {variant_display}**")

            # Top metrics row
            metrics_col = st.columns(4)
            metrics_col[0].metric("Evidence Strength", result['insight'].get('evidence_strength', 'N/A'))
            metrics_col[1].metric("Therapies", len(result.get('recommended_therapies', [])))
            clinvar_sig = result.get('clinvar', {}).get('clinical_significance', 'N/A') or 'N/A'
            metrics_col[2].metric("ClinVar", clinvar_sig)
            am_score = result.get('annotations', {}).get('alphamissense_score')
            metrics_col[3].metric("AlphaMissense", f"{am_score:.2f}" if am_score else 'N/A')

            # Card 1: Summary (always expanded)
            with st.expander("📋 Summary", expanded=True):
                st.markdown(result['insight'].get('summary', 'No summary available'))

            # Card 2: Identifiers & Links
            with st.expander("🔗 Identifiers & External Links", expanded=True):
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

            # Card 3: Recommended Therapies
            therapies = result.get('recommended_therapies', [])
            if therapies:
                with st.expander(f"💊 Recommended Therapies ({len(therapies)})", expanded=True):
                    for t in therapies:
                        drug = t.get('drug_name', 'Unknown')
                        level = t.get('evidence_level', 'N/A')
                        status = t.get('approval_status', '')
                        context = t.get('clinical_context', '')

                        # Create clickable drug link
                        drug_link = f"[{drug}](https://www.drugs.com/search.php?searchterm={drug.replace(' ', '+')})"
                        st.markdown(f"- **{drug_link}** (Level {level}) - {status}")
                        if context:
                            st.caption(f"   {context[:150]}...")

            # Card 4: Functional Predictions
            annotations = result.get('annotations', {})
            with st.expander("🧬 Functional Predictions", expanded=False):
                col1, col2, col3 = st.columns(3)

                with col1:
                    st.markdown("**Pathogenicity:**")
                    am_pred = annotations.get('alphamissense_prediction')
                    am_score_val = annotations.get('alphamissense_score')
                    if am_score_val:
                        pred_label = {'P': 'Pathogenic', 'B': 'Benign', 'A': 'Ambiguous'}.get(am_pred, am_pred)
                        st.metric("AlphaMissense", f"{am_score_val:.3f}", pred_label)

                    pp2 = annotations.get('polyphen2_prediction')
                    if pp2:
                        st.metric("PolyPhen2", pp2)

                with col2:
                    st.markdown("**Deleteriousness:**")
                    cadd = annotations.get('cadd_score')
                    if cadd:
                        st.metric("CADD", f"{cadd:.1f}", ">20 = top 1%")

                    effect = annotations.get('snpeff_effect')
                    if effect:
                        st.write(f"**Effect:** {effect}")

                with col3:
                    st.markdown("**Population Frequency:**")
                    gnomad = annotations.get('gnomad_exome_af')
                    if gnomad:
                        st.metric("gnomAD AF", f"{gnomad:.2e}")
                    else:
                        st.write("Not found in gnomAD (rare)")

            # Card 5: Clinical Trials
            if panel and panel.get('clinical', {}).get('clinical_trials'):
                trials = panel['clinical']['clinical_trials']
                with st.expander(f"🏥 Clinical Trials ({len(trials)})", expanded=False):
                    for trial in trials[:10]:
                        nct_id = trial.get('nct_id', 'Unknown')
                        title = trial.get('title', 'No title')[:80]
                        status = trial.get('status', 'Unknown')
                        phase = trial.get('phase', '')

                        # Clickable NCT link
                        nct_link = f"[{nct_id}](https://clinicaltrials.gov/study/{nct_id})"
                        variant_tag = " 🎯" if trial.get('variant_specific') else ""
                        st.markdown(f"- {nct_link}{variant_tag} - {title}...")
                        st.caption(f"   {status} | {phase}")

            # Card 6: Knowledge Base Evidence
            if panel and panel.get('kb'):
                kb = panel['kb']
                has_kb_data = kb.get('civic_assertions') or kb.get('cgi_biomarkers') or kb.get('vicc')
                if has_kb_data:
                    with st.expander("📚 Knowledge Base Evidence", expanded=False):
                        # CIViC Assertions
                        if kb.get('civic_assertions'):
                            st.markdown("**CIViC Assertions:**")
                            for a in kb['civic_assertions'][:5]:
                                therapies_str = ", ".join(a.get('therapies', [])) or "N/A"
                                sig = a.get('significance', 'Unknown')
                                disease = a.get('disease', '')
                                aid = a.get('id', '')
                                link = f"[AID:{aid}](https://civicdb.org/assertions/{aid})" if aid else f"AID:{aid}"
                                st.markdown(f"- {link}: {therapies_str} → {sig} ({disease})")

                        # CGI Biomarkers
                        if kb.get('cgi_biomarkers'):
                            st.markdown("**CGI Biomarkers:**")
                            for b in kb['cgi_biomarkers'][:5]:
                                drug = b.get('drug', 'Unknown')
                                assoc = b.get('association', 'Unknown')
                                fda = " ✓ FDA" if b.get('fda_approved') else ""
                                st.markdown(f"- {drug}: {assoc}{fda}")

                        # VICC MetaKB
                        if kb.get('vicc'):
                            st.markdown("**VICC MetaKB:**")
                            for v in kb['vicc'][:5]:
                                source = v.get('source', 'Unknown')
                                drugs = ", ".join(v.get('drugs', [])) or "N/A"
                                sig = v.get('clinical_significance', 'Unknown')
                                st.markdown(f"- {source}: {drugs} → {sig}")

            # Card 7: Literature (Full mode only)
            if panel and panel.get('literature', {}).get('pubmed_articles'):
                articles = panel['literature']['pubmed_articles']
                with st.expander(f"📄 Literature ({len(articles)} articles)", expanded=False):
                    for article in articles[:10]:
                        pmid = article.get('pmid', '')
                        title = article.get('title', 'No title')[:70]
                        year = article.get('year', '')
                        journal = article.get('journal', '')

                        # Clickable PubMed link
                        pmid_link = f"[PMID:{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
                        st.markdown(f"- {pmid_link} - {title}...")
                        st.caption(f"   {journal} ({year})")

                    # Literature knowledge summary if available
                    lit_knowledge = panel['literature'].get('literature_knowledge')
                    if lit_knowledge:
                        st.markdown("---")
                        st.markdown("**LLM-Synthesized Insights:**")
                        st.info(lit_knowledge)

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
                    "Evidence": r['insight'].get('evidence_strength', 'N/A'),
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

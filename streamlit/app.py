"""OncoMind Streamlit Application - Variant insight and evidence synthesis tool."""
import streamlit as st
import pandas as pd
import asyncio
import json
from backend import get_variant_insight, get_variant_annotation, batch_get_variant_insights

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

        gene = st.text_input("Gene Symbol", placeholder="e.g., BRAF, EGFR, TP53", key="gene_input")
        variant = st.text_input("Variant", placeholder="e.g., V600E, L858R", key="variant_input")
        tumor = st.text_input("Tumor Type (optional)", placeholder="e.g., Melanoma, NSCLC", key="tumor_input")

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
                        result = asyncio.run(get_variant_annotation(
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

            st.success(f"✅ Insight Ready")
            metrics_col = st.columns(3)
            metrics_col[0].metric("Evidence Strength", result['insight'].get('evidence_strength', 'N/A'))
            metrics_col[1].metric("Therapies", len(result.get('recommended_therapies', [])))
            metrics_col[2].metric("ClinVar", result.get('clinvar', {}).get('clinical_significance', 'N/A') or 'N/A')

            st.subheader("Summary")
            st.markdown(result['insight'].get('summary', 'No summary available'))

            st.subheader("Complete Insight")
            st.json(result)
            st.download_button(
                "📥 Download JSON",
                json.dumps(result, indent=2),
                f"{gene_display}_{variant_display}_insight.json",
                "application/json",
                key="download_single"
            )
            # Future features placeholders
            with st.expander("🧬 Protein Structure (Coming Soon)"):
                st.info("ESMFold visualization will be added here")
            with st.expander("🤖 Agent Workflow (Coming Soon)"):
                st.info("Multi-agent analysis pipeline will be added here")

            # Clear button
            if st.button("🗑️ Clear Results", key="clear_single"):
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

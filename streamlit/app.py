"""OncoMind Streamlit Application - Variant annotation and evidence synthesis tool."""
import streamlit as st
import pandas as pd
import asyncio
import json
from backend import get_variant_insight, batch_get_variant_insights

st.set_page_config(page_title="OncoMind", page_icon="🧬", layout="wide")
st.title("🧬 OncoMind: Variant Annotation")
st.caption("**Note:** This tool is for research purposes only. Clinical decisions should always be made by qualified healthcare professionals.")

MODELS = {
    "OpenAI GPT-4o-mini": "gpt-4o-mini",
    "OpenAI GPT-4o": "gpt-4o",
    "Anthropic Claude 3 Haiku": "claude-3-haiku-20240307",
    "Google Gemini 1.5 Pro": "gemini/gemini-1.5-pro",
    "Groq Llama 3.1 70B": "groq/llama-3.1-70b-versatile"
}

tab1, tab2 = st.tabs(["🔬 Single Variant", "📊 Batch Upload"])

# TAB 1: Single Variant
with tab1:
    col1, col2 = st.columns([1, 2])
    with col1:
        st.subheader("Input")

        gene = st.text_input("Gene Symbol", placeholder="e.g., BRAF, EGFR, TP53", key="gene_input")
        variant = st.text_input("Variant", placeholder="e.g., V600E, L858R", key="variant_input")
        tumor = st.text_input("Tumor Type (optional)", placeholder="e.g., Melanoma, NSCLC", key="tumor_input")

        model_name = st.selectbox("LLM Model", list(MODELS.keys()))
        temperature = st.slider("Temperature", 0.0, 1.0, 0.1, 0.05)

        insight_btn = st.button("🔍 Process Variant", type="primary", use_container_width=True)

    with col2:
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
                    with st.spinner(f"🔬 Analyzing {gene} {variant}... Fetching evidence from CIViC, ClinVar, and COSMIC databases"):
                        result = asyncio.run(get_variant_insight(gene, variant, tumor or None, MODELS[model_name], temperature))
                        if "error" in result:
                            st.error(result["error"])
                        else:
                            st.success(f"✅ Annotation Complete")
                            metrics_col = st.columns(3)
                            metrics_col[0].metric("Evidence Strength", result['insight'].get('evidence_strength', 'N/A'))
                            metrics_col[1].metric("Therapies", len(result.get('recommended_therapies', [])))
                            metrics_col[2].metric("ClinVar", result.get('clinvar', {}).get('clinical_significance', 'N/A') or 'N/A')

                            st.subheader("Summary")
                            st.markdown(result['insight'].get('summary', 'No summary available'))

                            st.subheader("Complete Annotation")
                            st.json(result)
                            st.download_button("📥 Download JSON", json.dumps(result, indent=2),
                                             f"{gene}_{variant}_annotation.json", "application/json")
                            # Future features placeholders
                            with st.expander("🧬 Protein Structure (Coming Soon)"):
                                st.info("ESMFold visualization will be added here")
                            with st.expander("🤖 Agent Workflow (Coming Soon)"):
                                st.info("Multi-agent analysis pipeline will be added here")

# TAB 2: Batch Upload
with tab2:
    st.subheader("Batch Variant Annotation")
    st.markdown("**CSV Format:** Must contain `gene`, `variant`, and optionally `tumor_type` columns")
    col1, col2 = st.columns([1, 1])
    with col1:
        model_name_batch = st.selectbox("LLM Model", list(MODELS.keys()), key="batch_model")
    with col2:
        temperature_batch = st.slider("Temperature", 0.0, 1.0, 0.1, 0.05, key="batch_temp")

    uploaded_file = st.file_uploader("Upload CSV", type=['csv'])
    if uploaded_file:
        df = pd.read_csv(uploaded_file)
        st.dataframe(df.head(), use_container_width=True)
        if st.button("🚀 Process Batch", type="primary"):
            if 'gene' not in df.columns or 'variant' not in df.columns:
                st.error("CSV must contain 'gene' and 'variant' columns")
            else:
                progress_bar = st.progress(0)
                status_text = st.empty()
                variants = [{"gene": row.get('gene'), "variant": row.get('variant'),
                           "tumor_type": row.get('tumor_type', None)} for _, row in df.iterrows()]
                results = asyncio.run(batch_get_variant_insights(variants, MODELS[model_name_batch],
                          temperature_batch, lambda i, t: (progress_bar.progress(i/t),
                          status_text.text(f"Processing {i}/{t}..."))))
                status_text.text("✅ Batch processing complete!")
                progress_bar.progress(1.0)
                results_df = pd.DataFrame([{"Gene": r['variant']['gene'], "Variant": r['variant']['variant'],
                    "Tumor": r['variant'].get('tumor_type', 'N/A'),
                    "Evidence": r['insight'].get('evidence_strength', 'N/A'),
                    "Therapies": len(r.get('recommended_therapies', []))} for r in results if 'error' not in r])
                st.dataframe(results_df, use_container_width=True)
                st.download_button("📥 Download Results CSV", results_df.to_csv(index=False),
                                 "batch_results.csv", "text/csv")
                st.download_button("📥 Download Full JSON", json.dumps(results, indent=2),
                                 "batch_results.json", "application/json")

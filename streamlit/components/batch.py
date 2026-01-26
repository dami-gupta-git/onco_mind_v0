"""Batch upload component for OncoMind."""

import asyncio
import io
import json
import zipfile
import pandas as pd
import streamlit as st

from oncomind.config.constants import LLM_DEFAULT_TEMPERATURE
from .utils import result_to_markdown


def render_batch_tab(batch_get_variant_insights, models: dict[str, str]) -> None:
    """Render the batch upload tab.

    Args:
        batch_get_variant_insights: Async function to get batch insights
        models: Dictionary mapping model display names to model IDs
    """
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
            model_name_batch = st.selectbox("Model", list(models.keys()), key="batch_model")
        else:
            model_name_batch = list(models.keys())[0]

    uploaded_file = st.file_uploader("Upload CSV", type=['csv'], help="CSV files only. Must contain 'gene' and 'variant' columns.")
    if uploaded_file:
        df = pd.read_csv(uploaded_file)
        st.dataframe(df.head(), width="stretch")
        if st.button("🚀 Get Batch Insights", type="primary"):
            if 'gene' not in df.columns or 'variant' not in df.columns:
                st.error("CSV must contain 'gene' and 'variant' columns")
            else:
                _process_batch(
                    df=df,
                    batch_get_variant_insights=batch_get_variant_insights,
                    models=models,
                    model_name=model_name_batch,
                    enable_llm=enable_llm_batch,
                    enable_literature=enable_literature_batch,
                    literature_source=literature_source_value_batch,
                )

    _render_batch_results()


def _process_batch(
    df: pd.DataFrame,
    batch_get_variant_insights,
    models: dict[str, str],
    model_name: str,
    enable_llm: bool,
    enable_literature: bool,
    literature_source: str,
) -> None:
    """Process a batch of variants."""
    progress_bar = st.progress(0)
    status_text = st.empty()

    variants = [
        {
            "gene": row.get('gene'),
            "variant": row.get('variant'),
            "tumor_type": row.get('tumor_type', None)
        }
        for _, row in df.iterrows()
    ]

    results = asyncio.run(batch_get_variant_insights(
        variants=variants,
        enable_llm=enable_llm,
        enable_literature=enable_literature,
        literature_source=literature_source,
        model=models[model_name],
        temperature=LLM_DEFAULT_TEMPERATURE,
        progress_callback=lambda i, t: (
            progress_bar.progress(i/t),
            status_text.text(f"Processing {i}/{t}...")
        )
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

    results_df = pd.DataFrame([
        {
            "Gene": r['variant']['gene'],
            "Variant": r['variant']['variant'],
            "Tumor": r['variant'].get('tumor_type', 'N/A'),
            "Therapies": len(r.get('recommended_therapies', []))
        }
        for r in results if 'error' not in r
    ])
    st.session_state.batch_df = results_df


def _render_batch_results() -> None:
    """Render batch results if available."""
    if st.session_state.get('batch_results') is None:
        return

    st.dataframe(st.session_state.batch_df, width="stretch")

    batch_footer_cols = st.columns(3)
    with batch_footer_cols[0]:
        # Create ZIP file with markdown reports
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
            for r in st.session_state.batch_results:
                if 'error' not in r:
                    gene = r['variant']['gene']
                    variant = r['variant']['variant']
                    filename = f"{gene}_{variant}_report.md"
                    markdown_content = result_to_markdown(r)
                    zf.writestr(filename, markdown_content)
        zip_buffer.seek(0)

        st.download_button(
            "📥 Download Reports (ZIP)",
            zip_buffer.getvalue(),
            "batch_reports.zip",
            "application/zip",
            key="download_batch_zip"
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

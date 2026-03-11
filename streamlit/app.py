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
    # render_cross_source_analysis,  # Disabled
    # Evidence tab renderers
    render_functional_tab,
    render_civic_tab,
    render_vicc_tab,
    render_cgi_tab,
    render_fda_tab,
    render_clinvar_tab,
    render_cosmic_tab,
    render_trials_tab,
    render_literature_tab,
    render_cbioportal_tab,
    render_depmap_tab,
    render_hotspots_tab,
)

# Initialize logging from environment variable (ONCOMIND_LOG_LEVEL=DEBUG|INFO|WARN|ERROR)
# This import triggers the logger initialization which reads from env
from oncomind.config.debug import get_logger
logger = get_logger(__name__)

st.set_page_config(page_title="OncoMind", page_icon="🧬", layout="wide")

# Apply custom styling
apply_styles()
# scrollable_table function moved to components/utils.py

st.markdown(
    "<div class='header-container'>"
    "<h1 class='header-title'>🧬 OncoMind: Find the gaps in your Cancer Research</h1>"
    "<span class='tagline-text'>What's missing — and what you can do next.</span>"
    "</div>",
    unsafe_allow_html=True
)
st.markdown(
    "<div style='background: #d4a5a5; color: white; padding: 8px 16px; "
    "border-radius: 4px; font-weight: 600; text-align: center; margin: 8px 0;'>"
    "PROOF OF CONCEPT | Results are representative </div>",
    unsafe_allow_html=True
)
#st.markdown("<p style='text-align: center; margin-top: -4px; color: #666; font-size: 0.90em;'>What DON'T we know about your variant?</p>", unsafe_allow_html=True)
# st.caption("**Note:** This tool is for demo purposes only. The results are representative.")


MODELS = {
    "Google Gemini 2.0 Flash": "gemini/gemini-2.0-flash",
    "Google Gemini 1.5 Pro": "gemini/gemini-1.5-pro",
    "Anthropic Claude Sonnet 4": "claude-sonnet-4-20250514",
    "Anthropic Claude 3.5 Haiku": "claude-3-5-haiku-20241022",
    "OpenAI GPT-4o-mini": "gpt-4o-mini",
    "OpenAI GPT-4o": "gpt-4o",
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
    "BRAF V600E - Melanoma": ("BRAF", "V600E", "Melanoma"),
    "PIK3CA H1047R  - Thyroid Cancer": ("PIK3CA", "H1047R", "Thyroid Cancer"),
    "BRAF V600K - Colorectal": ("BRAF", "V600E", "Colorectal"),
    "EGFR L858R - NSCLC": ("EGFR", "L858R", "NSCLC"),
    "KRAS G12C - NSCLC": ("KRAS", "G12C", "NSCLC"),
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


st.subheader("Pick a Variant")

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

    # Initialize default values in session state if not present
    if "gene_input" not in st.session_state:
        st.session_state.gene_input = "IDH1"
    if "variant_input" not in st.session_state:
        st.session_state.variant_input = "R132H"
    if "tumor_input" not in st.session_state:
        st.session_state.tumor_input = "Glioma"

    input_cols = st.columns([1.5, 1.5, 2, 1.2, 1.2, 1])
    with input_cols[0]:
        gene = st.text_input("Gene", placeholder="e.g. IDH1", key="gene_input")
    with input_cols[1]:
        variant = st.text_input("Variant", placeholder="e.g. R132H", key="variant_input")
    with input_cols[2]:
        tumor = st.text_input("Tumor Type", placeholder="e.g. Glioma", key="tumor_input")
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

    st.caption("Or try an example:")

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
        # Get HGVS protein notation if available (e.g., p.Gly12Cys)
        hgvs_protein = hgvs.get('protein')

        header_text = f"**{gene_display} {variant_display}"
        if hgvs_protein:
            header_text += f" ({hgvs_protein})"
        header_text += "**"
        if tumor_display:
            header_text += f" in {tumor_display}"
        header_text += " | GRCh38"

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
            # Note: fda_approvals has been removed - using fda_biomarker_evidence exclusively
            civic_assertions = result.get('civic_assertions', [])
            civic_evidence = result.get('civic_evidence', [])
            vicc = result.get('vicc_evidence', [])
            cgi_biomarkers = result.get('cgi_biomarkers', [])
            _clinvar_entries_raw = result.get('clinvar_entries', [])
            clinvar_entries = [
                e for e in _clinvar_entries_raw
                if (e.get('review_status') or '').lower() not in (
                    'no assertion criteria provided',
                    'no classification provided',
                    'no classification for the individual variant',
                )
            ]
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
            fda_biomarker_evidence = result.get('fda_biomarker_evidence', [])

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
            # FDA tab - FDA biomarker-drug indications with negation detection
            if fda_biomarker_evidence:
                tab_names.append(f"💊 FDA ({len(fda_biomarker_evidence)})")
            if clinvar_entries or clinvar_sig:
                tab_names.append(f"ClinVar ({len(clinvar_entries)})")
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
                        render_functional_tab(
                            annotations=annotations,
                            hgvs=hgvs,
                            gene_display=gene_display,
                            get_pdb_image_url=get_pdb_image_url,
                            get_pdb_page_url=get_pdb_page_url,
                        )
                    tab_idx += 1

                # CIViC tab
                if civic_assertions or civic_evidence:
                    with tabs[tab_idx]:
                        render_civic_tab(
                            civic_assertions=civic_assertions,
                            civic_evidence=civic_evidence,
                            tumor_display=tumor_display,
                        )
                    tab_idx += 1

                # VICC tab
                if vicc:
                    with tabs[tab_idx]:
                        render_vicc_tab(vicc=vicc, tumor_display=tumor_display)
                    tab_idx += 1

                # CGI tab
                if all_cgi_count > 0:
                    with tabs[tab_idx]:
                        render_cgi_tab(
                            cgi_biomarkers=cgi_biomarkers,
                            early_phase=early_phase,
                            preclinical=preclinical,
                        )
                    tab_idx += 1

                # FDA tab
                if fda_biomarker_evidence:
                    with tabs[tab_idx]:
                        render_fda_tab(fda_biomarker_evidence=fda_biomarker_evidence)
                    tab_idx += 1

                # ClinVar tab
                if clinvar_entries or clinvar_sig:
                    with tabs[tab_idx]:
                        render_clinvar_tab(clinvar_entries=clinvar_entries, clinvar_sig=clinvar_sig)
                    tab_idx += 1

                # COSMIC tab
                if cosmic_id:
                    with tabs[tab_idx]:
                        render_cosmic_tab(cosmic_id=cosmic_id)
                    tab_idx += 1

                # Trials tab
                if trials:
                    with tabs[tab_idx]:
                        render_trials_tab(trials=trials)
                    tab_idx += 1

                # Literature tab
                if articles:
                    with tabs[tab_idx]:
                        render_literature_tab(articles=articles)
                    tab_idx += 1

                # cBioPortal tab
                if cbioportal:
                    with tabs[tab_idx]:
                        render_cbioportal_tab(cbioportal=cbioportal)
                    tab_idx += 1

                # DepMap tab
                if depmap:
                    with tabs[tab_idx]:
                        render_depmap_tab(
                            depmap=depmap,
                            gene_display=gene_display,
                            variant_display=variant_display,
                            tumor_display=tumor_display,
                        )
                    tab_idx += 1

                # Hotspots tab
                if hotspots and (hotspots.get('is_hotspot') or hotspots.get('is_adjacent_to_hotspot')):
                    with tabs[tab_idx]:
                        render_hotspots_tab(hotspots=hotspots, gene_display=gene_display)
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
                f"<div class='gap-analysis-header'>"
                f"<h2 style='margin: 0; font-size: 1.75rem !important;'>🔍 Gap Analysis</h2>"
                f"<span class='gap-analysis-badges'>"
                f"<strong>Evidence Quality:</strong> {badge} {evidence_quality.capitalize()} &nbsp;&nbsp; "
                f"<strong>Research Priority:</strong> {priority_badge} {display_priority}"
                f"</span></div>",
                unsafe_allow_html=True
            )

            # Two tables side by side (Well Characterized ~20% wider than Evidence Gaps)
            table_cols = st.columns([7, 4])

            with table_cols[0]:
                st.markdown("<span style='font-weight: 600;'>✅ Well Characterized</span> — <em>what we know</em>", unsafe_allow_html=True)

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

                # FDA match breakdown now comes from fda_biomarker_evidence
                # This is handled in the FDA Biomarker tab directly
                fda_match_str = ""

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
        # LLM CROSS-SOURCE DRUG ANALYSIS (disabled)
        # ==============================================
        # render_cross_source_analysis(result)

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

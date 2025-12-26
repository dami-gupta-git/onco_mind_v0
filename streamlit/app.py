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

            # Build tab names with counts
            tab_names = []
            # Check if we have any functional scores
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
            # FDA tab removed - now shown in Therapeutic Evidence section below
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

            # Therapeutic Evidence tab - consolidated drug/therapy evidence
            therapies = result.get('recommended_therapies', [])
            if therapies:
                tab_names.append(f"Therapeutic Evidence ({len(therapies)})")

            if tab_names:
                tabs = st.tabs(tab_names)
                tab_idx = 0

                # Functional tab
                if has_functional:
                    with tabs[tab_idx]:
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
                    tab_idx += 1

                # CIViC tab
                if civic_assertions or civic_evidence:
                    with tabs[tab_idx]:
                        # Assertions (curated AMP/ASCO/CAP) - show as table
                        if civic_assertions:
                            st.markdown("**Curated Assertions:**")
                            assertion_rows = []
                            for a in civic_assertions:
                                therapies_str = ", ".join(a.get('therapies', [])) or "N/A"
                                assertion_rows.append({
                                    "ID": a.get('id', ''),
                                    "Therapies": therapies_str,
                                    "Significance": a.get('significance', 'Unknown'),
                                    "Disease": a.get('disease', '')[:40] if a.get('disease') else '',
                                    "AMP Level": a.get('amp_level', ''),
                                })
                            if assertion_rows:
                                st.dataframe(pd.DataFrame(assertion_rows), use_container_width=True, hide_index=True)

                        # Evidence items - scrollable table
                        if civic_evidence:
                            if civic_assertions:
                                st.markdown("---")
                            st.markdown(f"**Evidence Items ({len(civic_evidence)}):**")
                            evidence_rows = []
                            for e in civic_evidence:
                                drugs_str = ", ".join(e.get('drugs', [])) or "N/A"
                                evidence_rows.append({
                                    "Drugs": drugs_str[:30] if len(drugs_str) > 30 else drugs_str,
                                    "Significance": e.get('clinical_significance', 'Unknown'),
                                    "Disease": (e.get('disease', '') or '')[:25],
                                    "Level": e.get('evidence_level', ''),
                                    "Type": e.get('evidence_type', ''),
                                    "Rating": e.get('trust_rating') or e.get('rating') or '',
                                })
                            if evidence_rows:
                                st.dataframe(
                                    pd.DataFrame(evidence_rows),
                                    use_container_width=True,
                                    hide_index=True,
                                    height=300,  # Scrollable with fixed height
                                )
                    tab_idx += 1

                # VICC tab - scrollable table
                if vicc:
                    with tabs[tab_idx]:
                        # Build rows for all VICC entries
                        vicc_rows = []
                        for v in vicc:
                            drugs = ", ".join(v.get('drugs', [])) or "N/A"
                            vicc_rows.append({
                                "Source": (v.get('source') or 'vicc').upper(),
                                "Drugs": drugs[:30] if len(drugs) > 30 else drugs,
                                "Response": v.get('response_type', 'Unknown'),
                                "Disease": (v.get('disease', '') or '')[:25],
                                "Level": v.get('evidence_level', ''),
                            })
                        if vicc_rows:
                            st.dataframe(
                                pd.DataFrame(vicc_rows),
                                use_container_width=True,
                                hide_index=True,
                                height=300,  # Scrollable with fixed height
                            )
                    tab_idx += 1

                # CGI tab - scrollable table
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
                            df_kwargs = {"use_container_width": True, "hide_index": True}
                            if len(cgi_rows) > 8:
                                df_kwargs["height"] = 300
                            st.dataframe(pd.DataFrame(cgi_rows), **df_kwargs)
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
                    tab_idx += 1

                # cBioPortal tab - co-mutation and prevalence data
                if cbioportal:
                    with tabs[tab_idx]:
                        # Prevalence stats
                        st.markdown("**Prevalence:**")
                        total = cbioportal.get('total_samples', 0)
                        gene_pct = cbioportal.get('gene_prevalence_pct', 0)
                        variant_pct = cbioportal.get('variant_prevalence_pct', 0)
                        gene_count = cbioportal.get('samples_with_gene_mutation', 0)
                        variant_count = cbioportal.get('samples_with_exact_variant', 0)
                        study = cbioportal.get('study_id', 'N/A')

                        col1, col2 = st.columns(2)
                        with col1:
                            st.metric("Gene Mutation", f"{gene_pct:.1f}%", delta=f"{gene_count}/{total} samples")
                        with col2:
                            st.metric("Exact Variant", f"{variant_pct:.1f}%", delta=f"{variant_count}/{total} samples")
                        st.caption(f"Study: {study}")

                        # Co-occurring mutations
                        co_occurring = cbioportal.get('co_occurring', [])
                        if co_occurring:
                            st.markdown("---")
                            st.markdown(f"**Co-occurring Mutations ({len(co_occurring)}):**")
                            co_rows = []
                            for c in co_occurring:
                                odds = c.get('odds_ratio')
                                odds_str = f"{odds:.2f}" if odds else "N/A"
                                co_rows.append({
                                    "Gene": c.get('gene', ''),
                                    "Count": c.get('count', 0),
                                    "Frequency": f"{c.get('pct', 0):.1f}%",
                                    "Odds Ratio": odds_str,
                                })
                            df_kwargs = {"use_container_width": True, "hide_index": True}
                            if len(co_rows) > 8:
                                df_kwargs["height"] = 250
                            st.dataframe(pd.DataFrame(co_rows), **df_kwargs)

                        # Mutually exclusive mutations
                        mutually_exclusive = cbioportal.get('mutually_exclusive', [])
                        if mutually_exclusive:
                            st.markdown("---")
                            st.markdown(f"**Mutually Exclusive ({len(mutually_exclusive)}):**")
                            me_rows = []
                            for m in mutually_exclusive:
                                odds = m.get('odds_ratio')
                                odds_str = f"{odds:.2f}" if odds else "N/A"
                                me_rows.append({
                                    "Gene": m.get('gene', ''),
                                    "Count": m.get('count', 0),
                                    "Frequency": f"{m.get('pct', 0):.1f}%",
                                    "Odds Ratio": odds_str,
                                })
                            df_kwargs = {"use_container_width": True, "hide_index": True}
                            if len(me_rows) > 8:
                                df_kwargs["height"] = 250
                            st.dataframe(pd.DataFrame(me_rows), **df_kwargs)
                    tab_idx += 1

                # Therapeutics tab - consolidated drug/therapy evidence
                if therapies:
                    with tabs[tab_idx]:
                        # Group therapies by evidence tier
                        fda_approved = [t for t in therapies if t.get('evidence_level', '').lower() == 'fda-approved']
                        clinical = [t for t in therapies if t.get('evidence_level', '').lower() in ('phase 3', 'phase 2', 'phase 1', 'case report')]
                        preclinical_therapies = [t for t in therapies if t.get('evidence_level', '').lower() in ('preclinical', 'in vitro')]

                        # FDA-Approved therapies
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
                            st.dataframe(pd.DataFrame(fda_rows), use_container_width=True, hide_index=True)

                        # Clinical trial evidence
                        if clinical:
                            if fda_approved:
                                st.markdown("---")
                            st.markdown("**🔬 Clinical Evidence:**")
                            clinical_rows = []
                            for t in clinical:
                                drug = t.get('drug_name', 'Unknown')
                                level = t.get('evidence_level', '')
                                response = t.get('response_type', '')
                                source = t.get('source', '')
                                clinical_rows.append({
                                    "Drug": drug,
                                    "Level": level,
                                    "Response": response or "-",
                                    "Source": source,
                                })
                            st.dataframe(pd.DataFrame(clinical_rows), use_container_width=True, hide_index=True)

                        # Preclinical evidence
                        if preclinical_therapies:
                            if fda_approved or clinical:
                                st.markdown("---")
                            st.markdown("**🧪 Preclinical:**")
                            st.warning("⚠️ Preclinical data - not validated in humans")
                            preclin_rows = []
                            for t in preclinical_therapies:
                                drug = t.get('drug_name', 'Unknown')
                                response = t.get('response_type', '')
                                source = t.get('source', '')
                                preclin_rows.append({
                                    "Drug": drug,
                                    "Response": response or "-",
                                    "Source": source,
                                })
                            st.dataframe(pd.DataFrame(preclin_rows), use_container_width=True, hide_index=True)
                    tab_idx += 1
            else:
                st.info("No evidence found from any source")

            # ==============================================
            # LLM RESEARCH INSIGHT - After evidence tabs
            # ==============================================
            llm_narrative = result['insight'].get('llm_narrative')
            if llm_narrative:
                st.markdown("---")
                st.markdown("### 🤖 LLM Research Synthesis")

                # Main narrative
                st.markdown(llm_narrative)

                # Evidence assessment section
                st.markdown("#### Evidence Assessment")

                # Evidence quality badge
                evidence_quality = result['insight'].get('evidence_quality')
                if evidence_quality:
                    quality_colors = {
                        "comprehensive": "🟢",
                        "moderate": "🟡",
                        "limited": "🟠",
                        "minimal": "🔴",
                    }
                    badge = quality_colors.get(evidence_quality.lower(), "⚪")
                    st.markdown(f"**Overall Quality:** {badge} {evidence_quality.capitalize()}")

                # Well-characterized aspects (single line)
                well_characterized = result['insight'].get('well_characterized', [])
                if well_characterized:
                    st.markdown(f"**✅ Well Characterized:** {' · '.join(well_characterized)}")

                # Knowledge gaps (single line)
                knowledge_gaps = result['insight'].get('knowledge_gaps', [])
                if knowledge_gaps:
                    st.markdown(f"**❓ Knowledge Gaps:** {' · '.join(knowledge_gaps)}")

                # Conflicting evidence (single line)
                conflicting_evidence = result['insight'].get('conflicting_evidence', [])
                if conflicting_evidence:
                    st.markdown(f"**⚠️ Conflicting Evidence:** {' · '.join(conflicting_evidence)}")

                # Evidence tags (single line)
                evidence_tags = result['insight'].get('evidence_tags', [])
                if evidence_tags:
                    st.markdown(f"**🏷️ Evidence Types:** {' · '.join(evidence_tags)}")

                # Research implications
                research_implications = result['insight'].get('research_implications')
                if research_implications and research_implications != result['insight'].get('rationale'):
                    st.markdown(f"**🔬 Research Implications:** {research_implications}")

                # References (make clickable)
                references = result['insight'].get('references', [])
                if references:
                    clickable_refs = []
                    for ref in references[:5]:
                        ref_str = str(ref).strip()
                        # Check if it's a PMID
                        if ref_str.startswith("PMID"):
                            pmid = ref_str.replace("PMID", "").replace(":", "").strip()
                            clickable_refs.append(f"[PMID {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
                        elif ref_str.isdigit():
                            # Just a number, assume it's a PMID
                            clickable_refs.append(f"[PMID {ref_str}](https://pubmed.ncbi.nlm.nih.gov/{ref_str}/)")
                        elif "cBioPortal" in ref_str or "cbioportal" in ref_str.lower():
                            # Extract study ID if present
                            if ":" in ref_str:
                                study_id = ref_str.split(":")[-1].strip()
                                clickable_refs.append(f"[{ref_str}](https://www.cbioportal.org/study/summary?id={study_id})")
                            else:
                                clickable_refs.append(ref_str)
                        elif ref_str.startswith("NCT"):
                            # Clinical trial ID
                            clickable_refs.append(f"[{ref_str}](https://clinicaltrials.gov/study/{ref_str})")
                        else:
                            clickable_refs.append(ref_str)
                    st.markdown(f"**📚 Key References:** {', '.join(clickable_refs)}")

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

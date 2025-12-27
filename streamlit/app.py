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
        enable_llm = st.toggle(
            "🤖 Enable LLM Mode",
            value=False,
            help="LLM mode includes literature search and AI-powered synthesis (~25s). Without LLM, you get fast annotation (~7s)."
        )

        # Mode description
        if enable_llm:
            st.caption("📚 **LLM Mode** (~25s): Literature search + AI synthesis")
        else:
            st.caption("⚡ **Annotation and Gap Analysis Mode** (~7s): Fast structured evidence and gap analysis")

        # LLM mode also enables literature
        enable_literature = enable_llm

        # Only show LLM options if LLM is enabled
        if enable_llm:
            with st.expander("LLM Settings"):
                model_name = st.selectbox("LLM Model", list(MODELS.keys()))
                temperature = st.slider("Temperature", 0.0, 1.0, 0.1, 0.05)
        else:
            model_name = list(MODELS.keys())[0]
            temperature = 0.1

        insight_btn = st.button("🔍 Get Insight", type="primary", width="stretch")

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
                    mode_label = "LLM" if enable_llm else "Annotation"
                    with st.spinner(f"🔬 Getting insight for {gene} {variant} ({mode_label} mode)..."):
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
            depmap = result.get('depmap_evidence')

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
            if depmap:
                tab_names.append("🧬 DepMap")

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

                        # Add MyVariant.info source link
                        hgvs_genomic = hgvs.get('genomic')
                        if hgvs_genomic:
                            import urllib.parse
                            encoded_id = urllib.parse.quote(hgvs_genomic, safe='')
                            myvariant_url = f"https://myvariant.info/v1/variant/{encoded_id}"
                            st.markdown(f"*Source: [MyVariant.info]({myvariant_url})*")
                    tab_idx += 1

                # CIViC tab
                if civic_assertions or civic_evidence:
                    with tabs[tab_idx]:
                        # Evidence level legend
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
                        # Assertions (curated AMP/ASCO/CAP) - show as markdown table with clickable IDs
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

                        # Evidence items - markdown table with clickable IDs
                        if civic_evidence:
                            if civic_assertions:
                                st.markdown("---")
                            st.markdown(f"**Evidence Items ({len(civic_evidence)}):**")
                            rows = ["| ID | Drugs | Significance | Disease | Level | Type | Rating |",
                                    "|----|-------|--------------|---------|-------|------|--------|"]
                            for e in civic_evidence:
                                drugs_str = ", ".join(e.get('drugs', [])) or "N/A"
                                drugs_str = drugs_str[:30] if len(drugs_str) > 30 else drugs_str
                                eid = e.get('eid') or ''
                                url = e.get('civic_url', '')
                                id_link = f"[{eid}]({url})" if url else eid
                                sig = e.get('clinical_significance', 'Unknown')
                                disease = (e.get('disease', '') or '')[:25]
                                level = e.get('evidence_level', '')
                                etype = e.get('evidence_type', '')
                                rating = e.get('trust_rating') or e.get('rating') or ''
                                rows.append(f"| {id_link} | {drugs_str} | {sig} | {disease} | {level} | {etype} | {rating} |")
                            st.markdown("\n".join(rows))
                    tab_idx += 1

                # VICC tab - markdown table with clickable sources
                if vicc:
                    with tabs[tab_idx]:
                        # Evidence level legend
                        with st.expander("📖 Evidence Level Guide", expanded=False):
                            st.markdown("""
**Evidence Levels:**
- **1/A**: FDA-approved or standard of care
- **2/B**: Clinical trial evidence or expert consensus
- **3/C**: Case reports or limited evidence
- **4/D**: Preclinical or computational evidence
- **R1/R2**: Resistance evidence (strong/emerging)
""")
                        rows = ["| Source | Drugs | Response | Disease | Level |",
                                "|--------|-------|----------|---------|-------|"]
                        for v in vicc:
                            drugs = ", ".join(v.get('drugs', [])) or "N/A"
                            drugs = drugs[:30] if len(drugs) > 30 else drugs
                            source = (v.get('source') or 'vicc').upper()
                            url = v.get('publication_url', '')
                            source_link = f"[{source}]({url})" if url else source
                            response = v.get('response_type', 'Unknown')
                            disease = (v.get('disease', '') or '')[:25]
                            level = v.get('evidence_level', '')
                            rows.append(f"| {source_link} | {drugs} | {response} | {disease} | {level} |")
                        st.markdown("\n".join(rows))
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
                            df_kwargs = {"width": "stretch", "hide_index": True}
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
                        # Study header with cohort size
                        study_name = cbioportal.get('study_name', 'N/A')
                        study_id = cbioportal.get('study_id', '')
                        total = cbioportal.get('total_samples', 0)

                        if study_id:
                            study_url = f"https://www.cbioportal.org/study/summary?id={study_id}"
                            st.markdown(f"**Study:** [{study_name}]({study_url}) — cohort of {total:,} samples")
                        else:
                            st.markdown(f"**Study:** {study_name} — cohort of {total:,} samples")

                        # Prevalence stats
                        gene_pct = cbioportal.get('gene_prevalence_pct', 0)
                        variant_pct = cbioportal.get('variant_prevalence_pct', 0)
                        gene_count = cbioportal.get('samples_with_gene_mutation', 0)
                        variant_count = cbioportal.get('samples_with_exact_variant', 0)
                        gene_symbol = cbioportal.get('gene', 'Gene')
                        variant_name = cbioportal.get('variant', 'variant')

                        col1, col2 = st.columns(2)
                        with col1:
                            st.markdown(f"**% samples with any {gene_symbol} alteration**")
                            st.markdown(f"### {gene_pct:.1f}%")
                            st.caption(f"{gene_count}/{total} samples")
                        with col2:
                            st.markdown(f"**% samples with exact {variant_name} variant**")
                            st.markdown(f"### {variant_pct:.1f}%")
                            st.caption(f"{variant_count}/{total} samples")

                        # Co-occurring and Mutually exclusive mutations side by side
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
                                    st.dataframe(pd.DataFrame(co_rows), hide_index=True, width="stretch")

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
                                    st.dataframe(pd.DataFrame(me_rows), hide_index=True, width="stretch")
                    tab_idx += 1

                # DepMap tab - gene essentiality and drug sensitivity from cell lines
                if depmap:
                    with tabs[tab_idx]:
                        gene_dep = depmap.get('gene_dependency')
                        drug_sens = depmap.get('drug_sensitivities', [])
                        cell_lines = depmap.get('cell_line_models', [])
                        is_essential = depmap.get('is_essential', False)

                        # Gene Essentiality - the key insight
                        if gene_dep:
                            st.markdown("### Gene Essentiality")
                            score = gene_dep.get('mean_dependency_score')
                            dep_pct = gene_dep.get('dependency_pct', 0)
                            n_dep = gene_dep.get('n_dependent_lines', 0)
                            n_total = gene_dep.get('n_total_lines', 0)

                            col1, col2 = st.columns(2)
                            with col1:
                                if is_essential:
                                    st.error(f"🔴 **{gene_display} is ESSENTIAL**")
                                else:
                                    st.info(f"⚪ {gene_display} is not essential")
                            with col2:
                                st.metric(
                                    "CERES Score",
                                    f"{score:.2f}" if score else "N/A",
                                    delta=f"{dep_pct:.0f}% of cell lines depend on {gene_display}",
                                    delta_color="off"
                                )

                            st.caption(f"Based on CRISPR screens in {n_total} cancer cell lines. CERES < -0.5 indicates essentiality.")

                        # Drug Sensitivities
                        if drug_sens:
                            st.markdown("---")
                            st.markdown("### Drug Sensitivities")
                            st.caption("IC50 values from PRISM drug screens in cell lines")

                            drug_rows = []
                            for ds in drug_sens:
                                ic50 = ds.get('ic50_nm')
                                ic50_str = f"{ic50:.0f} nM" if ic50 else "N/A"
                                drug_rows.append({
                                    "Drug": ds.get('drug_name', ''),
                                    "IC50": ic50_str,
                                    "Cell Lines Tested": ds.get('n_cell_lines', 0),
                                })
                            st.dataframe(pd.DataFrame(drug_rows), width="stretch", hide_index=True)

                        # Cell Line Models
                        if cell_lines:
                            st.markdown("---")
                            st.markdown("### Model Cell Lines")
                            mutant_lines = [cl for cl in cell_lines if cl.get('has_mutation')]
                            if mutant_lines:
                                st.success(f"✅ {len(mutant_lines)} cell lines with {variant_display} mutation available for research")
                                cl_rows = []
                                for cl in mutant_lines:
                                    cl_rows.append({
                                        "Cell Line": cl.get('name', ''),
                                        "Disease": cl.get('primary_disease', ''),
                                        "Subtype": cl.get('subtype', ''),
                                        "Mutation": cl.get('mutation_details', variant_display),
                                    })
                                st.dataframe(pd.DataFrame(cl_rows), width="stretch", hide_index=True)
                            else:
                                st.info(f"{len(cell_lines)} cell lines available (mutation status unknown)")

                        # Link to DepMap
                        st.markdown("---")
                        st.markdown(f"[🔗 Explore {gene_display} on DepMap Portal](https://depmap.org/portal/gene/{gene_display})")
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
                            st.dataframe(pd.DataFrame(fda_rows), width="stretch", hide_index=True)

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
                            st.dataframe(pd.DataFrame(clinical_rows), width="stretch", hide_index=True)

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
                            st.dataframe(pd.DataFrame(preclin_rows), width="stretch", hide_index=True)
                    tab_idx += 1
            else:
                st.info("No evidence found from any source")

            # ==============================================
            # EVIDENCE ASSESSMENT - Always shown (deterministic gap detection)
            # ==============================================
            st.markdown("---")
            st.markdown("### 📊 Evidence Assessment")

            # Get structured evidence gaps (deterministic analysis)
            evidence_gaps = result.get('evidence_gaps', {})

            # Evidence quality and research priority badges
            col_quality, col_priority = st.columns(2)
            with col_quality:
                evidence_quality = evidence_gaps.get('overall_quality', 'unknown')
                quality_colors = {
                    "comprehensive": "🟢",
                    "moderate": "🟡",
                    "limited": "🟠",
                    "minimal": "🔴",
                }
                badge = quality_colors.get(evidence_quality.lower(), "⚪")
                st.markdown(f"**Evidence Quality:** {badge} {evidence_quality.capitalize()}")

            with col_priority:
                research_priority = evidence_gaps.get('research_priority', 'unknown')
                priority_colors = {
                    "very_high": "🔥",  # Prime research target
                    "high": "🔴",
                    "medium": "🟡",
                    "low": "🟢",
                }
                priority_badge = priority_colors.get(research_priority.lower(), "⚪")
                # Format "very_high" nicely
                display_priority = research_priority.replace("_", " ").title()
                st.markdown(f"**Research Priority:** {priority_badge} {display_priority}")

            # Well-characterized aspects as table (with basis/reasoning)
            well_characterized_detailed = evidence_gaps.get('well_characterized_detailed', [])
            if well_characterized_detailed:
                st.markdown("**✅ Well Characterized:**")
                wc_df = pd.DataFrame([
                    {
                        "Category": (item.get('category') or '').replace('_', ' ').title(),
                        "Aspect": item.get('aspect', ''),
                        "Basis": item.get('basis', ''),
                    }
                    for item in well_characterized_detailed
                ])
                st.dataframe(wc_df, width="stretch", hide_index=True)
            else:
                # Fallback to legacy format if detailed not available
                well_characterized = evidence_gaps.get('well_characterized', [])
                if well_characterized:
                    st.markdown("**✅ Well Characterized:**")
                    wc_df = pd.DataFrame({"Aspect": well_characterized})
                    st.dataframe(wc_df, width="stretch", hide_index=True)

            # Evidence gaps table
            gaps = evidence_gaps.get('gaps', [])
            if gaps:
                st.markdown("**❓ Evidence Gaps:**")

                # Build gaps table
                gaps_data = []
                for gap in gaps:
                    severity = gap.get('severity', 'unknown')
                    severity_icon = {"critical": "🔴", "significant": "🟠", "minor": "🟡"}.get(severity, "⚪")
                    # Strip variant name from description (user already knows)
                    desc = gap.get('description', '')
                    # Remove "for GENE VARIANT" or "for GENE VARIANT" pattern at end
                    import re
                    desc = re.sub(r'\s+for\s+\w+\s+\S+$', '', desc)
                    # Also handle "of GENE VARIANT" pattern
                    desc = re.sub(r'\s+of\s+\w+\s+\S+\s+in\s+\S+$', '', desc)
                    desc = re.sub(r'\s+of\s+\w+\s+\S+$', '', desc)
                    gaps_data.append({
                        "Severity": f"{severity_icon} {severity.capitalize()}",
                        "Category": gap.get('category', '').replace('_', ' ').title(),
                        "Description": desc,
                    })

                if gaps_data:
                    gaps_df = pd.DataFrame(gaps_data)
                    st.dataframe(gaps_df, width="stretch", hide_index=True)

                # Suggested studies (collapsible)
                all_suggested = []
                for gap in gaps:
                    all_suggested.extend(gap.get('suggested_studies', []))
                if all_suggested:
                    with st.expander("📋 Suggested Studies"):
                        for study in list(set(all_suggested)):  # Deduplicate
                            st.markdown(f"  - {study}")

            # ==============================================
            # LLM RESEARCH INSIGHT - After evidence assessment
            # ==============================================
            llm_narrative = result['insight'].get('llm_narrative')
            if llm_narrative:
                st.markdown("---")
                st.markdown("### 🤖 LLM Research Synthesis")

                # Build formatted narrative from raw components
                functional_summary = result['insight'].get('functional_summary')
                biological_context = result['insight'].get('biological_context')
                therapeutic_landscape = result['insight'].get('therapeutic_landscape')

                # Display structured components with formatting
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

                # Fall back to plain narrative if no structured components
                if not any([functional_summary, biological_context, therapeutic_landscape]):
                    st.markdown(llm_narrative)

                # LLM-specific: Conflicting evidence, evidence tags, research hypotheses, references
                conflicting_evidence = result['insight'].get('conflicting_evidence', [])
                if conflicting_evidence:
                    st.markdown(f"**⚠️ Conflicting Evidence:** {' · '.join(conflicting_evidence)}")

                # Evidence tags (single line)
                evidence_tags = [tag.title() for tag in result['insight'].get('evidence_tags', [])]
                if evidence_tags:
                    st.markdown(f"**🏷️ Evidence Types:** {' · '.join(evidence_tags)}")

                # Research implications
                research_implications = result['insight'].get('research_implications')
                if research_implications and research_implications != result['insight'].get('rationale'):
                    st.markdown(f"**🔬 Research Implications:** {research_implications}")

                # Research hypotheses (LLM-generated testable hypotheses)
                research_hypotheses = result['insight'].get('research_hypotheses', [])
                if research_hypotheses:
                    st.markdown("**💡 Emerging Research Hypotheses:**")
                    for i, hypothesis in enumerate(research_hypotheses[:3], 1):
                        st.markdown(f"  {i}. {hypothesis}")

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
        enable_llm_batch = st.toggle(
            "🤖 Enable LLM Mode",
            value=False,
            help="LLM mode includes literature search + AI synthesis (~25s/variant). Without LLM: fast annotation (~7s/variant).",
            key="batch_llm"
        )
        # LLM mode also enables literature
        enable_literature_batch = enable_llm_batch
    with col2:
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

                # Store in session state
                st.session_state.batch_results = results
                results_df = pd.DataFrame([{"Gene": r['variant']['gene'], "Variant": r['variant']['variant'],
                    "Tumor": r['variant'].get('tumor_type', 'N/A'),
                    "Therapies": len(r.get('recommended_therapies', []))} for r in results if 'error' not in r])
                st.session_state.batch_df = results_df

    # Display batch results from session state
    if st.session_state.batch_results is not None:
        st.dataframe(st.session_state.batch_df, width="stretch")
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

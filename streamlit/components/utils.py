"""Utility functions for OncoMind Streamlit components."""

import re
import html as html_module
import streamlit as st

from oncomind.config.constants import (
    format_civic_significance,
    CADD_DELETERIOUS_THRESHOLD,
    GNOMAD_RARE_THRESHOLD,
)


def scrollable_table(markdown_content: str) -> None:
    """Render a markdown table inside a scrollable container.

    Converts markdown table to HTML for proper rendering with CSS styles.
    Long cells are truncated with ellipsis and show full content on hover (title)
    or click (expands the cell).
    """
    # Parse markdown table to HTML
    lines = [line.strip() for line in markdown_content.strip().split('\n') if line.strip()]
    if len(lines) < 2:
        st.markdown(markdown_content)
        return

    # Parse header
    header_cells = [cell.strip() for cell in lines[0].split('|') if cell.strip()]

    # Map header names to CSS classes for column widths
    col_class_map = {
        'source': 'col-source',
        'id': 'col-source',
        'fda': 'col-fda',
        'locus': 'col-locus',
        'locus match': 'col-locus',
        'tumor': 'col-tumor',
        'tumor match': 'col-tumor',
        'drugs': 'col-drugs',
        'drug': 'col-drugs',
        'response': 'col-response',
        'association': 'col-response',
        'significance': 'col-response',
        'match': 'col-match',
        'specificity': 'col-specificity',
        'drug specificity': 'col-specificity',
        'regulatory status': 'col-status',
        'status': 'col-status',
        'indication': 'col-indication',
        'links': 'col-links',
        'nct id': 'col-nct',
        'phase': 'col-phase',
        'title': 'col-title',
    }

    # Skip separator line (line with dashes)
    data_lines = [l for l in lines[2:] if l and not all(c in '|-: ' for c in l)]

    # Build HTML table
    html_parts = ['<table>']

    # Header row with column classes
    html_parts.append('<thead><tr>')
    header_classes = []
    for cell in header_cells:
        col_class = col_class_map.get(cell.lower(), '')
        header_classes.append(col_class)
        class_attr = f' class="{col_class}"' if col_class else ''
        # Use line breaks for two-word headers to save space
        display_cell = cell
        if cell.lower() == 'locus match':
            display_cell = 'Locus<br>Match'
        elif cell.lower() == 'tumor match':
            display_cell = 'Tumor<br>Match'
        elif cell.lower() == 'drug specificity':
            display_cell = 'Drug<br>Specificity'
        html_parts.append(f'<th{class_attr}>{display_cell}</th>')
    html_parts.append('</tr></thead>')

    # Body rows
    html_parts.append('<tbody>')
    for line in data_lines:
        cells = [cell.strip() for cell in line.split('|') if cell.strip() or cell == '']
        # Filter out empty strings from split but keep actual empty cells
        cells = [c for i, c in enumerate(line.split('|')) if i > 0 and i < len(line.split('|')) - 1]
        cells = [c.strip() for c in cells]

        html_parts.append('<tr>')
        for idx, cell in enumerate(cells):
            # Convert markdown links to HTML links
            cell_html = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<a href="\2" target="_blank">\1</a>', cell)

            # Get column class
            col_class = header_classes[idx] if idx < len(header_classes) else ''

            # For long text, truncate with ellipsis and add tooltip
            # Strip HTML tags to get plain text length
            plain_text = re.sub(r'<[^>]+>', '', cell_html)
            is_long = len(plain_text) > 35

            if is_long:
                # Escape for title attribute (replace quotes and special chars)
                escaped_full = html_module.escape(plain_text).replace('"', '&quot;')
                # Build class list
                classes = ['truncated']
                if col_class:
                    classes.append(col_class)
                class_str = ' '.join(classes)
                # Use title for hover tooltip (works in Streamlit iframe)
                # Use inline onclick for click-to-expand (works without external JS)
                # Wrap content in div for proper truncation
                onclick = "this.classList.toggle('expanded')"
                html_parts.append(f'<td class="{class_str}" title="{escaped_full}" onclick="{onclick}"><div class="cell-content">{cell_html}</div></td>')
            else:
                class_attr = f'class="{col_class}"' if col_class else ''
                class_str = f' {class_attr}' if class_attr else ''
                html_parts.append(f'<td{class_str}>{cell_html}</td>')
        html_parts.append('</tr>')
    html_parts.append('</tbody>')
    html_parts.append('</table>')

    html_table = '\n'.join(html_parts)

    st.markdown(
        f'<div class="scrollable-table">{html_table}</div>',
        unsafe_allow_html=True
    )


def format_pmid_link(pmid: str) -> str:
    """Format a PMID as a clickable link."""
    pmid_clean = str(pmid).strip()
    return f"[PMID {pmid_clean}](https://pubmed.ncbi.nlm.nih.gov/{pmid_clean}/)"


def format_nct_link(nct_id: str) -> str:
    """Format an NCT ID as a clickable link."""
    nct_clean = str(nct_id).strip()
    return f"[{nct_clean}](https://clinicaltrials.gov/study/{nct_clean})"


def result_to_markdown(result: dict) -> str:
    """Convert a variant insight result to markdown format for export."""
    lines = []

    # Header
    variant = result.get('variant', {})
    gene = variant.get('gene', 'Unknown')
    var = variant.get('variant', 'Unknown')
    tumor = variant.get('tumor_type', 'Not specified')

    lines.append(f"# {gene} {var} Variant Report")
    lines.append(f"\n**Tumor Type:** {tumor}")
    lines.append(f"\n*Generated by OncoMind*\n")

    # Identifiers
    ids = result.get('identifiers', {})
    hgvs = result.get('hgvs', {})
    if ids or hgvs:
        lines.append("## Identifiers\n")
        if ids.get('dbsnp'):
            lines.append(f"- **dbSNP:** {ids['dbsnp']}")
        if ids.get('cosmic'):
            lines.append(f"- **COSMIC:** {ids['cosmic']}")
        if hgvs.get('genomic'):
            lines.append(f"- **HGVS Genomic:** {hgvs['genomic']}")
        if hgvs.get('coding'):
            lines.append(f"- **HGVS Coding:** {hgvs['coding']}")
        if hgvs.get('protein'):
            lines.append(f"- **HGVS Protein:** {hgvs['protein']}")
        lines.append("")

    # Functional Scores (from annotations)
    annotations = result.get('annotations', {})
    func_rows = []
    if annotations.get('alphamissense_score') is not None:
        pred = annotations.get('alphamissense_prediction') or '-'
        func_rows.append(f"| AlphaMissense | {annotations['alphamissense_score']:.3f} | {pred} |")
    if annotations.get('cadd_score') is not None:
        cadd = annotations['cadd_score']
        pred = 'Deleterious' if cadd > CADD_DELETERIOUS_THRESHOLD else 'Benign'
        func_rows.append(f"| CADD | {cadd:.1f} | {pred} |")
    if annotations.get('polyphen2_prediction'):
        func_rows.append(f"| PolyPhen2 | - | {annotations['polyphen2_prediction']} |")
    if annotations.get('gnomad_exome_af') is not None:
        af = annotations['gnomad_exome_af']
        freq = f"{af:.2e}" if af < GNOMAD_RARE_THRESHOLD else f"{af:.4f}"
        pred = 'Rare' if af < GNOMAD_RARE_THRESHOLD else 'Common'
        func_rows.append(f"| gnomAD AF | {freq} | {pred} |")
    if annotations.get('snpeff_effect'):
        func_rows.append(f"| SnpEff | - | {annotations['snpeff_effect']} |")
    if func_rows:
        lines.append("## Functional Scores\n")
        lines.append("| Score | Value | Prediction |")
        lines.append("|-------|-------|------------|")
        lines.extend(func_rows)
        lines.append("")

    # CIViC Evidence - show all Level A entries, limit others to 15
    civic = result.get('civic_evidence', [])
    if civic:
        lines.append("## CIViC Evidence\n")
        lines.append("| EID | Drug | Disease | Evidence Level | Significance | Locus Match | Tumor Match |")
        lines.append("|-----|------|---------|----------------|--------------|-------------|-------------|")
        seen_entries = set()
        # Separate Level A from others
        level_a_items = [item for item in civic if (item.get('evidence_level') or '').upper() == 'A']
        other_items = [item for item in civic if (item.get('evidence_level') or '').upper() != 'A']

        # Process all Level A items (no limit, deduplicate by EID)
        seen_eids = set()
        for item in level_a_items:
            evidence_id = item.get('evidence_id')
            if evidence_id in seen_eids:
                continue
            seen_eids.add(evidence_id)
            eid = f"EID{evidence_id}" if evidence_id else ""
            drugs_list = item.get('drugs', [])
            drug = ", ".join(drugs_list) if drugs_list else "N/A"
            disease = item.get('disease', '')
            level = "A"
            sig_raw = item.get('clinical_significance', '')
            sig = format_civic_significance(sig_raw)
            locus = item.get('locus_match', '')
            tumor = "Yes" if item.get('tumor_match') else ("No" if item.get('tumor_match') is False else "")
            lines.append(f"| {eid} | {drug} | {disease} | {level} | {sig} | {locus} | {tumor} |")

        # Process other items (limit to 15)
        other_count = 0
        for item in other_items:
            if other_count >= 15:
                break
            evidence_id = item.get('evidence_id')
            eid = f"EID{evidence_id}" if evidence_id else ""
            drugs_list = item.get('drugs', [])
            drug = ", ".join(drugs_list) if drugs_list else "N/A"
            disease = item.get('disease', '')
            level = item.get('evidence_level', '')
            sig_raw = item.get('clinical_significance', '')
            sig = format_civic_significance(sig_raw)
            locus = item.get('locus_match', '')
            tumor = "Yes" if item.get('tumor_match') else ("No" if item.get('tumor_match') is False else "")
            # Deduplicate by drug+disease+significance+level
            entry_key = (
                drug.lower() if drug else '',
                disease.lower() if disease else '',
                sig_raw.lower() if sig_raw else '',
                level.lower() if level else ''
            )
            if entry_key in seen_entries:
                continue
            seen_entries.add(entry_key)
            lines.append(f"| {eid} | {drug} | {disease} | {level} | {sig} | {locus} | {tumor} |")
            other_count += 1
        lines.append("")

    # VICC Evidence
    vicc = result.get('vicc_evidence', [])
    if vicc:
        lines.append("## VICC MetaKB Evidence\n")
        lines.append("| Source | Drug | Disease | Response | Locus Match | Tumor Match |")
        lines.append("|--------|------|---------|----------|-------------|-------------|")
        for item in vicc[:15]:
            source = item.get('source', '')
            # drugs is a list, join them
            drugs_list = item.get('drugs', [])
            drug = ', '.join(drugs_list) if drugs_list else ''
            disease = item.get('disease', '')
            response = item.get('response_type', '')
            locus = item.get('locus_match', '')
            tumor = "Yes" if item.get('tumor_match') else ("No" if item.get('tumor_match') is False else "")
            lines.append(f"| {source} | {drug} | {disease} | {response} | {locus} | {tumor} |")
        lines.append("")

    # CGI Biomarkers
    cgi = result.get('cgi_biomarkers', [])
    if cgi:
        lines.append("## CGI Biomarkers\n")
        lines.append("| Drug | Association | Evidence Level | Locus Match | Tumor Match |")
        lines.append("|------|-------------|----------------|-------------|-------------|")
        for item in cgi[:15]:
            drug = item.get('drug', '')
            assoc = item.get('association', '')
            level = item.get('evidence_level', '')
            locus = item.get('locus_match', '')
            tumor = "Yes" if item.get('tumor_match') else ("No" if item.get('tumor_match') is False else "")
            lines.append(f"| {drug} | {assoc} | {level} | {locus} | {tumor} |")
        lines.append("")

    # FDA Approvals (from fda_biomarker_evidence)
    # Always show this section - it's important to know if there are no FDA approvals
    # Matches UI: Drug | Gene | Match | Label Level | Variants | Tumors
    fda = result.get('fda_biomarker_evidence', [])
    lines.append("## FDA Approvals\n")
    if fda:
        lines.append("| Drug | Gene | Match | Label Level | Variants | Tumors |")
        lines.append("|------|------|-------|-------------|----------|--------|")
        for item in fda:
            drug = item.get('drug_name', '')
            brand = item.get('brand_name', '')
            if brand:
                drug = f"{drug} ({brand})"
            # Include combination partners if present
            combination_partners = item.get('combination_partners', [])
            if combination_partners:
                partners_str = " + ".join(combination_partners)
                drug = f"{drug} + {partners_str}"
            gene = item.get('gene', '')
            # Match type from variant_match_result
            match_type = item.get('variant_match_result', '')
            if match_type == 'exact':
                match = '✅ exact'
            elif match_type == 'codon':
                match = '🔸 codon'
            elif match_type == 'gene':
                match = '🧬 gene'
            else:
                match = match_type or ''
            # Label level from specificity
            specificity = item.get('specificity', '')
            # Variants specified in the label
            specified_variants = item.get('specified_variants', [])
            variants = ", ".join(specified_variants[:3]) if specified_variants else '-'
            # Tumor types
            tumor_types = item.get('tumor_types', [])
            tumors = ", ".join(tumor_types[:2]) if tumor_types else ''
            lines.append(f"| {drug} | {gene} | {match} | {specificity} | {variants} | {tumors} |")
    else:
        lines.append("None")
    lines.append("")

    # ClinVar Evidence
    clinvar_entries = result.get('clinvar_entries', [])
    clinvar_sig = result.get('clinvar', {}).get('clinical_significance')
    if clinvar_entries or clinvar_sig:
        lines.append("## ClinVar\n")
        if clinvar_sig:
            lines.append(f"**Clinical Significance:** {clinvar_sig}\n")
        if clinvar_entries:
            lines.append("| Variation ID | Clinical Significance | Conditions | Review Status |")
            lines.append("|--------------|----------------------|------------|---------------|")
            for entry in clinvar_entries[:15]:
                var_id = entry.get('variation_id', '')
                sig = entry.get('clinical_significance', '')
                conditions = entry.get('conditions', [])
                conditions_str = ", ".join(conditions[:2]) if conditions else '-'
                review = entry.get('review_status', '')
                var_link = f"[{var_id}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_id}/)" if var_id else '-'
                lines.append(f"| {var_link} | {sig} | {conditions_str} | {review} |")
            lines.append("")

    # Clinical Trials
    trials = result.get('clinical_trials', [])
    if trials:
        lines.append("## Clinical Trials\n")
        lines.append("| NCT ID | Phase | Title | Locus Match | Tumor Match |")
        lines.append("|--------|-------|-------|-------------|-------------|")
        for trial in trials[:15]:
            nct = trial.get('nct_id', '')
            title = trial.get('title', '')
            phase = trial.get('phase', '')
            locus = trial.get('locus_match', '')
            tumor = "Yes" if trial.get('tumor_match') else ("No" if trial.get('tumor_match') is False else "")
            # Truncate long titles for table display
            title_display = title[:60] + "..." if len(title) > 60 else title
            nct_link = f"[{nct}](https://clinicaltrials.gov/study/{nct})"
            lines.append(f"| {nct_link} | {phase} | {title_display} | {locus} | {tumor} |")
        lines.append("")

    # cBioPortal Prevalence
    cbio = result.get('cbioportal_evidence')
    if cbio:
        lines.append("## cBioPortal Prevalence\n")
        study_name = cbio.get('study_name', '')
        tumor_type = cbio.get('tumor_type', '')
        total = cbio.get('total_samples', 0)
        gene_prev = cbio.get('gene_prevalence_pct')
        var_prev = cbio.get('variant_prevalence_pct')
        gene_samples = cbio.get('samples_with_gene_mutation', 0)
        var_samples = cbio.get('samples_with_exact_variant', 0)

        lines.append("| Metric | Value |")
        lines.append("|--------|-------|")
        if study_name:
            lines.append(f"| Study | {study_name} |")
        if tumor_type:
            lines.append(f"| Tumor Type | {tumor_type} |")
        if total:
            lines.append(f"| Total Samples | {total:,} |")
        if gene_prev is not None:
            lines.append(f"| Gene Mutation Prevalence | {gene_prev:.1f}% ({gene_samples:,} samples) |")
        if var_prev is not None:
            lines.append(f"| Exact Variant Prevalence | {var_prev:.1f}% ({var_samples:,} samples) |")

        # Co-occurring mutations (first 3)
        co_occurring = cbio.get('co_occurring', [])
        if co_occurring:
            co_str = ", ".join([f"{c.get('gene', '')} ({c.get('percentage', 0):.0f}%)" for c in co_occurring[:3]])
            lines.append(f"| Co-occurring Mutations | {co_str} |")

        # Mutually exclusive (first 3)
        mutually_exclusive = cbio.get('mutually_exclusive', [])
        if mutually_exclusive:
            me_str = ", ".join([f"{m.get('gene', '')}" for m in mutually_exclusive[:3]])
            lines.append(f"| Mutually Exclusive | {me_str} |")

        lines.append("")

    # DepMap Evidence
    depmap = result.get('depmap_evidence')
    if depmap:
        lines.append("## DepMap\n")
        lines.append("| Metric | Value |")
        lines.append("|--------|-------|")
        gene_dep = depmap.get('gene_dependency')
        if gene_dep:
            mean_score = gene_dep.get('mean_dependency_score')
            if mean_score is not None:
                lines.append(f"| Mean Dependency Score | {mean_score:.3f} |")
            dep_pct = gene_dep.get('dependency_pct')
            if dep_pct is not None:
                lines.append(f"| Dependency % | {dep_pct:.1f}% |")
            n_dep = gene_dep.get('n_dependent_lines', 0)
            n_total = gene_dep.get('n_total_lines', 0)
            if n_total > 0:
                lines.append(f"| Dependent Cell Lines | {n_dep} / {n_total} |")
        is_essential = depmap.get('is_essential')
        if is_essential is not None:
            lines.append(f"| Essential Gene | {'Yes' if is_essential else 'No'} |")
        n_screened = depmap.get('n_cell_lines_screened', 0)
        if n_screened > 0:
            lines.append(f"| Cell Lines Screened | {n_screened} |")
        lines.append("")

    # Hotspots Evidence
    hotspots = result.get('hotspots_evidence')
    if hotspots:
        is_hotspot = hotspots.get('is_hotspot', False)
        is_adjacent = hotspots.get('is_adjacent_to_hotspot', False)
        if is_hotspot or is_adjacent:
            lines.append("## Cancer Hotspots\n")
            lines.append("| Metric | Value |")
            lines.append("|--------|-------|")
            if is_hotspot:
                lines.append("| Hotspot Status | **Yes** |")
                hotspot_data = hotspots.get('hotspot', {})
                if hotspot_data:
                    residue = hotspot_data.get('residue', '')
                    total_count = hotspot_data.get('total_count', 0)
                    if residue:
                        lines.append(f"| Residue | {residue} |")
                    if total_count > 0:
                        lines.append(f"| Total Samples | {total_count:,} |")
            elif is_adjacent:
                adj_dist = hotspots.get('adjacent_distance')
                lines.append(f"| Hotspot Status | Adjacent ({adj_dist} codons away) |")
            lines.append("")

    # Literature
    articles = result.get('pubmed_articles', [])
    if articles:
        lines.append("## Literature\n")
        for article in articles[:10]:
            pmid = article.get('pmid', '')
            title = article.get('title', '')
            year = article.get('year', '')
            lines.append(f"- [{title}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/) ({year})")
        lines.append("")

    # Gap Analysis (just above LLM Synthesis)
    evidence_gaps = result.get('evidence_gaps', {})
    if evidence_gaps:
        lines.append("## Evidence Gap Analysis\n")
        lines.append("*What's known vs. unknown about this variant — identifying opportunities for further research.*\n")

        overall_quality = evidence_gaps.get('overall_quality', 'unknown')
        research_priority = evidence_gaps.get('research_priority', 'unknown')
        lines.append(f"**Evidence Quality:** {overall_quality.upper()}")
        lines.append(f"**Research Priority:** {research_priority.upper()}\n")

        # Well-characterized aspects
        well_characterized_detailed = evidence_gaps.get('well_characterized_detailed', [])
        if well_characterized_detailed:
            lines.append("### Well Characterized\n")
            for wc in well_characterized_detailed:
                aspect = wc.get('aspect', '')
                basis = wc.get('basis', '')
                matches = wc.get('matches_on', '')
                tumor = wc.get('tumor_match', '')
                match_info = []
                if matches:
                    match_info.append(matches)
                if tumor:
                    match_info.append(tumor)
                match_str = f" ({', '.join(match_info)})" if match_info else ""
                lines.append(f"- **{aspect}**: {basis}{match_str}")
            lines.append("")

        # Gaps by severity
        gaps = evidence_gaps.get('gaps', [])
        if gaps:
            # Group gaps by severity
            severity_order = ['critical', 'high', 'significant', 'moderate', 'minor', 'informational']
            gaps_by_severity = {}
            for gap in gaps:
                sev = gap.get('severity', 'unknown')
                if sev not in gaps_by_severity:
                    gaps_by_severity[sev] = []
                gaps_by_severity[sev].append(gap)

            lines.append("### Evidence Gaps\n")
            for sev in severity_order:
                if sev in gaps_by_severity:
                    severity_label = sev.upper()
                    for gap in gaps_by_severity[sev]:
                        desc = gap.get('description', '')
                        category = gap.get('category', '')
                        lines.append(f"- **[{severity_label}]** ({category}) {desc}")
            lines.append("")

    # LLM Synthesis (at the bottom, like in the UI)
    insight = result.get('insight', {})
    functional_summary = insight.get('functional_summary')
    biological_context = insight.get('biological_context')
    therapeutic_summary = insight.get('therapeutic_summary')
    llm_narrative = insight.get('llm_narrative')
    research_implications = insight.get('research_implications')

    has_llm_content = any([functional_summary, biological_context, therapeutic_summary, llm_narrative])
    if has_llm_content:
        lines.append("## Experimental - LLM Research Synthesis\n")
        lines.append("*AI-generated — may contain inaccuracies. Verify against primary sources.*\n")
        lines.append("*Locus Match: variant-level = evidence match at exact locus · codon-level = evidence match at other variants in this codon · gene-level = evidence match at other variants in this gene*\n")
        if functional_summary:
            lines.append(f"### Functional Summary\n{functional_summary}\n")
        if biological_context:
            lines.append(f"### Biological Context\n{biological_context}\n")
        if therapeutic_summary:
            lines.append(f"### Therapeutic Summary\n{therapeutic_summary}\n")
        if research_implications:
            lines.append(f"### Research Implications\n{research_implications}\n")
        # Fallback to llm_narrative if no structured fields
        if not any([functional_summary, biological_context, therapeutic_summary]) and llm_narrative:
            lines.append(f"{llm_narrative}\n")

    # Cross-source analysis (disabled)
    # cross_source = result.get('cross_source_analysis', {})
    # if cross_source:
    #     lines.append("## Cross-Source Drug Analysis\n")
    #     if cross_source.get('summary'):
    #         lines.append(f"{cross_source['summary']}\n")
    #
    #     strongest = cross_source.get('strongest_evidence', [])
    #     if strongest:
    #         lines.append("### Strongest Evidence\n")
    #         for item in strongest:
    #             drug = item.get('drug', 'Unknown')
    #             signal = item.get('signal', '')
    #             sources = ", ".join(item.get('sources', []))
    #             level = item.get('evidence_level', '')
    #             lines.append(f"- **{drug}** ({signal}, {level}) - Sources: {sources}")
    #         lines.append("")
    #
    #     conflicts = cross_source.get('conflicting_signals', [])
    #     if conflicts:
    #         lines.append("### Conflicting Signals\n")
    #         for item in conflicts:
    #             drug = item.get('drug', 'Unknown')
    #             conflict = item.get('conflict', '')
    #             lines.append(f"- **{drug}**: {conflict}")
    #         lines.append("")

    return "\n".join(lines)

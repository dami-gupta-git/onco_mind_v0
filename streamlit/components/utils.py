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

    # Gap Analysis
    evidence_gaps = result.get('evidence_gaps', {})
    if evidence_gaps:
        lines.append("## Evidence Gap Analysis\n")

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

        # Poorly characterized aspects
        poorly_characterized = evidence_gaps.get('poorly_characterized', [])
        if poorly_characterized:
            lines.append("### Needs More Research\n")
            for pc in poorly_characterized:
                lines.append(f"- {pc}")
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

    # FDA Biomarker Evidence (from fda_biomarker_evidence)
    # Always show this section - it's important to know if there are no FDA approvals
    # Matches UI: Drug | Gene | Match | Label Level | Variants | Tumors
    fda = result.get('fda_biomarker_evidence', [])
    lines.append("## FDA Biomarker Evidence\n")
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

    # CIViC Evidence - show all entries including those without drugs (oncogenic, functional, etc.)
    civic = result.get('civic_evidence', [])
    if civic:
        lines.append("## CIViC Evidence\n")
        lines.append("| Drug | Disease | Evidence Level | Significance | Locus Match | Tumor Match |")
        lines.append("|------|---------|----------------|--------------|-------------|-------------|")
        seen_entries = set()
        count = 0
        for item in civic:
            if count >= 10:
                break
            drugs_list = item.get('drugs', [])
            drug = ", ".join(drugs_list) if drugs_list else "N/A"
            disease = item.get('disease', '')
            level = item.get('evidence_level', '')
            sig_raw = item.get('clinical_significance', '')
            sig = format_civic_significance(sig_raw)
            locus = item.get('locus_match', '')
            tumor = "Yes" if item.get('tumor_match') else ("No" if item.get('tumor_match') is False else "")
            # Deduplicate by drug+disease+significance+level (include level to show different evidence strengths)
            entry_key = (
                drug.lower() if drug else '',
                disease.lower() if disease else '',
                sig_raw.lower() if sig_raw else '',
                level.lower() if level else ''
            )
            if entry_key in seen_entries:
                continue
            seen_entries.add(entry_key)
            lines.append(f"| {drug} | {disease} | {level} | {sig} | {locus} | {tumor} |")
            count += 1
        lines.append("")

    # VICC Evidence
    vicc = result.get('vicc_evidence', [])
    if vicc:
        lines.append("## VICC MetaKB Evidence\n")
        lines.append("| Source | Drug | Disease | Response | Locus Match | Tumor Match |")
        lines.append("|--------|------|---------|----------|-------------|-------------|")
        for item in vicc[:10]:
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
        for item in cgi[:10]:
            drug = item.get('drug', '')
            assoc = item.get('association', '')
            level = item.get('evidence_level', '')
            locus = item.get('locus_match', '')
            tumor = "Yes" if item.get('tumor_match') else ("No" if item.get('tumor_match') is False else "")
            lines.append(f"| {drug} | {assoc} | {level} | {locus} | {tumor} |")
        lines.append("")

    # Clinical Trials
    trials = result.get('clinical_trials', [])
    if trials:
        lines.append("## Clinical Trials\n")
        lines.append("| NCT ID | Phase | Title | Locus Match | Tumor Match |")
        lines.append("|--------|-------|-------|-------------|-------------|")
        for trial in trials[:10]:
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

"""Utility functions for OncoMind Streamlit components."""

import re
import html as html_module
import streamlit as st

from oncomind.config.constants import format_civic_significance


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

    # FDA Biomarker Evidence (from fda_biomarker_evidence)
    fda = result.get('fda_biomarker_evidence', [])
    if fda:
        lines.append("## FDA Biomarker Evidence\n")
        lines.append("| Drug | Tumor Types | Match Type |")
        lines.append("|------|-------------|------------|")
        for item in fda:
            drug = item.get('drug_name', '')
            # Include combination partners if present
            combination_partners = item.get('combination_partners', [])
            if combination_partners:
                partners_str = " + ".join(combination_partners)
                drug = f"{drug} + {partners_str}"
            tumor_types = item.get('tumor_types', [])
            tumor_stage = item.get('tumor_stage', '')
            # Combine tumor stage with tumor types
            if tumor_stage and tumor_types:
                tumors = f"{tumor_stage} {', '.join(tumor_types[:2])}"
            elif tumor_types:
                tumors = ", ".join(tumor_types[:2])
            else:
                tumors = ''
            match_type = item.get('match_type', '')
            lines.append(f"| {drug} | {tumors} | {match_type} |")
        lines.append("")

    # CIViC Evidence
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
            # Skip entries without drugs (prognostic/diagnostic evidence)
            if not drugs_list:
                continue
            drug = ", ".join(drugs_list)
            disease = item.get('disease', '')
            level = item.get('evidence_level', '')
            sig_raw = item.get('clinical_significance', '')
            sig = format_civic_significance(sig_raw)
            locus = item.get('locus_match', '')
            tumor = "Yes" if item.get('tumor_match') else ("No" if item.get('tumor_match') is False else "")
            # Deduplicate by drug+disease+significance
            entry_key = (drug.lower(), disease.lower(), sig_raw.lower() if sig_raw else '')
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
        for trial in trials[:5]:
            nct = trial.get('nct_id', '')
            title = trial.get('title', '')
            phase = trial.get('phase', '')
            status = trial.get('status', '')
            lines.append(f"### [{nct}](https://clinicaltrials.gov/study/{nct})")
            lines.append(f"**{title}**\n")
            lines.append(f"- Phase: {phase}")
            lines.append(f"- Status: {status}\n")

    # Literature
    articles = result.get('pubmed_articles', [])
    if articles:
        lines.append("## Literature\n")
        for article in articles[:5]:
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

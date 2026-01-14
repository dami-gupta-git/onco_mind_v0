"""Utility functions for OncoMind Streamlit components."""

import re
import html as html_module
import streamlit as st


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

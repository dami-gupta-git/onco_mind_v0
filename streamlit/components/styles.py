"""CSS styles for OncoMind Streamlit application."""

import streamlit as st


def apply_styles() -> None:
    """Apply custom CSS styles to the Streamlit app."""
    st.markdown("""
<style>
    /* Global font size reduction */
    html, body, [class*="css"] {
        font-size: 13px !important;
    }
    .stMarkdown, .stText, p, span, div {
        font-size: 0.85rem !important;
    }
    h1 { font-size: 1.5rem !important; }
    h2 { font-size: 1.25rem !important; }
    h3 { font-size: 1.05rem !important; }
    h4 { font-size: 0.95rem !important; }
    /* Tagline text on the right of title */
    .tagline-text {
        font-size: 16px !important;
        color: #666 !important;
    }
    .block-container {
        padding-top: 2.5rem;
    }
    /* Make text inputs clearly look like input fields */
    .stTextInput > div > div > input {
        border: 1.5px solid #d0d0d0 !important;
        border-radius: 8px !important;
        padding: 0.5rem 0.75rem !important;
        background-color: #fafafa !important;
        transition: border-color 0.2s, box-shadow 0.2s !important;
    }
    .stTextInput > div > div > input:focus {
        border-color: #ff4b4b !important;
        box-shadow: 0 0 0 2px rgba(255, 75, 75, 0.2) !important;
        background-color: #fff !important;
    }
    .stTextInput > div > div > input:hover {
        border-color: #a0a0a0 !important;
    }
    /* Placeholder text styling */
    .stTextInput > div > div > input::placeholder {
        color: #999 !important;
        font-style: italic !important;
    }
    /* Smaller text in dataframes/tables */
    .stDataFrame [data-testid="stDataFrameResizable"] {
        font-size: 0.8rem !important;
    }
    .stDataFrame table {
        font-size: 0.8rem !important;
    }
    .stDataFrame th {
        font-size: 0.75rem !important;
    }
    .stDataFrame td {
        font-size: 0.8rem !important;
        padding: 4px 8px !important;
    }
    /* Scrollable table container for evidence tables */
    .scrollable-table {
        display: block;
        max-height: 400px;
        overflow-y: auto;
        overflow-x: hidden;
        border: 1px solid #e0e0e0;
        border-radius: 4px;
        padding: 0.5rem;
        margin-bottom: 0.25rem;
    }
    .scrollable-table table {
        font-size: 0.85rem !important;
        width: 100%;
        table-layout: fixed;
        border-collapse: collapse;
    }
    .scrollable-table th, .scrollable-table td {
        padding: 6px 10px !important;
        border-bottom: 1px solid #e0e0e0;
        text-align: left;
        vertical-align: top;
    }
    .scrollable-table th {
        background-color: #f5f5f5;
        font-weight: 600;
        position: sticky;
        top: 0;
        z-index: 1;
    }
    .scrollable-table a {
        color: #1a73e8;
        text-decoration: none;
    }
    .scrollable-table a:hover {
        text-decoration: underline;
    }
    /* Column widths */
    .scrollable-table .col-source { width: 100px; }
    .scrollable-table .col-fda { width: 55px !important; max-width: 55px !important; min-width: 55px !important; text-align: center; }
    .scrollable-table .col-locus { width: 100px; }
    .scrollable-table .col-tumor { width: 70px; }
    .scrollable-table .col-drugs { width: 180px; }
    .scrollable-table .col-response { width: 110px; }
    .scrollable-table .col-specificity { width: 100px; }
    .scrollable-table .col-status { width: 150px; }
    .scrollable-table .col-indication { width: 250px; }
    .scrollable-table .col-links { width: 80px; }
    .scrollable-table .col-nct { width: 120px; }
    .scrollable-table .col-phase { width: 100px; }
    .scrollable-table .col-title { min-width: 300px; }
    /* Truncated text cells - inner div handles the truncation */
    .scrollable-table td.truncated {
        cursor: pointer;
        padding: 0 !important;
    }
    .scrollable-table td.truncated .cell-content {
        max-width: 200px;
        width: 200px;
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
        direction: ltr !important;
        text-align: left !important;
        unicode-bidi: isolate;
        display: block;
        padding: 6px 10px;
        /* Force text to start from left edge */
        text-indent: 0;
        margin-left: 0;
        padding-left: 10px;
    }
    /* Click to expand - toggle class via inline onclick */
    .scrollable-table td.expanded .cell-content {
        white-space: normal !important;
        word-wrap: break-word;
        max-width: none !important;
        overflow: visible !important;
        text-overflow: clip !important;
    }
    .scrollable-table td.expanded {
        background-color: #fffde7 !important;
    }
    /* Visual hint that cell is clickable */
    .scrollable-table td.truncated:hover {
        background-color: #f0f7ff;
    }
    /* Title column gets more width for readability */
    .scrollable-table td.col-title .cell-content {
        max-width: 400px;
        width: 400px;
    }
    /* Tighter section dividers for evidence groupings */
    .evidence-section-divider {
        margin: 0.5rem 0 !important;
        border: none;
        border-top: 1px solid #e0e0e0;
    }
</style>
""", unsafe_allow_html=True)

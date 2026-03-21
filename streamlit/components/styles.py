"""CSS styles for OncoMind Streamlit application."""

import streamlit as st


def apply_styles() -> None:
    """Apply custom CSS styles to the Streamlit app."""
    st.markdown("""
<style>

    /* System font stack for reliable rendering across environments */
    html, body, [class*="css"] {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif !important;
        font-size: 15px !important;
    }
    .stMarkdown, .stText, p, span, div {
        font-size: 0.95rem !important;
    }
    h1, .stMarkdown h1, [data-testid="stMarkdown"] h1 { font-size: 2.2rem !important; }
    h2 { font-size: 1.5rem !important; }
    h3 { font-size: 1.25rem !important; }
    h4 { font-size: 1.1rem !important; }
    /* Header container - responsive layout */
    .header-container {
        display: flex;
        justify-content: space-between;
        align-items: center;
        gap: 1rem;
        flex-wrap: wrap;
    }
    .header-title {
        margin-bottom: 0;
        font-size: 28px !important;
    }
    /* Tagline text on the right of title */
    .tagline-text {
        font-size: 16px !important;
        color: #666 !important;
    }
    /* Mobile responsive header */
    @media (max-width: 768px) {
        .header-container {
            flex-direction: column;
            align-items: flex-start;
            gap: 0;
        }
        .header-title {
            font-size: 22px !important;
            line-height: 1.3;
            margin-bottom: 0 !important;
        }
        .tagline-text {
            font-size: 14px !important;
            margin-top: 0 !important;
        }
    }
    .block-container {
        padding-top: 1rem !important;
    }
    /* Gap Analysis header - responsive layout */
    .gap-analysis-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        flex-wrap: wrap;
        gap: 0.5rem;
    }
    .gap-analysis-title {
        font-size: 1.75rem;
        font-weight: 600;
    }
    .gap-analysis-badges {
        font-size: 0.9rem;
        text-align: right;
    }
    @media (max-width: 768px) {
        .gap-analysis-header {
            flex-direction: column;
            align-items: flex-start;
        }
        .gap-analysis-title {
            font-size: 1.25rem;
        }
        .gap-analysis-badges {
            font-size: 0.8rem;
            text-align: left;
        }
        /* Hide Locus Match and Tumor Match columns on mobile */
        .wc-table th:nth-child(3),
        .wc-table td:nth-child(3),
        .wc-table th:nth-child(4),
        .wc-table td:nth-child(4) {
            display: none !important;
        }
        /* Make Evidence Sources tabs scrollable on mobile */
        [data-testid="stTabs"] > [data-baseweb="tab-list"] {
            overflow-x: auto;
            -webkit-overflow-scrolling: touch;
            flex-wrap: nowrap;
            padding-bottom: 8px;
        }
        /* Visible scrollbar for tabs */
        [data-testid="stTabs"] > [data-baseweb="tab-list"]::-webkit-scrollbar {
            height: 4px;
        }
        [data-testid="stTabs"] > [data-baseweb="tab-list"]::-webkit-scrollbar-track {
            background: #f0f0f0;
            border-radius: 2px;
        }
        [data-testid="stTabs"] > [data-baseweb="tab-list"]::-webkit-scrollbar-thumb {
            background: #ccc;
            border-radius: 2px;
        }
    }
    /* Make text inputs clearly look like input fields */
    .stTextInput > div > div > input {
        border: 1.5px solid #d0d0d0 !important;
        border-radius: 8px !important;
        padding: 0.5rem 0.75rem !important;
        background-color: #fafafa !important;
    }
    /* Placeholder text styling */
    .stTextInput > div > div > input::placeholder {
        color: #999 !important;
        font-style: italic !important;
    }
    /* Dataframe tables - slightly larger for readability */
    .stDataFrame [data-testid="stDataFrameResizable"] {
        font-size: 0.9rem !important;
    }
    .stDataFrame table {
        font-size: 0.9rem !important;
    }
    .stDataFrame th {
        font-size: 0.85rem !important;
    }
    .stDataFrame td {
        font-size: 0.9rem !important;
        padding: 5px 10px !important;
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
        font-size: 0.75rem !important;
        width: 100%;
        table-layout: fixed;
        border-collapse: collapse;
    }
    .scrollable-table th, .scrollable-table td {
        padding: 6px 10px !important;
        border-bottom: 1px solid #e0e0e0;
        text-align: left;
        vertical-align: top;
        font-size: 0.75rem !important;
    }
    .scrollable-table th {
        background-color: #f5f5f5;
        font-weight: 600;
        position: sticky;
        top: 0;
        z-index: 1;
    }
    .scrollable-table td * {
        font-size: inherit !important;
    }
    .scrollable-table a {
        color: #1a73e8;
        text-decoration: none;
        font-size: inherit !important;
    }
    .scrollable-table a:hover {
        text-decoration: underline;
    }
    /* Column widths */
    .scrollable-table .col-source { width: 100px; }
    .scrollable-table .col-fda { width: 55px !important; max-width: 55px !important; min-width: 55px !important; text-align: center; }
    .scrollable-table .col-locus { width: 100px; }
    .scrollable-table .col-tumor { width: 70px; }
    .scrollable-table .col-drugs { width: 500px; min-width: 300px; max-width: 500px; }
    .scrollable-table .col-response { width: 110px; }
    .scrollable-table .col-match { width: 60px; text-align: center; }
    .scrollable-table .col-specificity { width: 90px; }
    .scrollable-table .col-status { width: 150px; }
    .scrollable-table .col-indication { width: 250px; }
    .scrollable-table .col-links { width: 80px; }
    .scrollable-table .col-nct { width: 120px; }
    .scrollable-table .col-phase { width: 100px; }
    .scrollable-table .col-title { width: 450px; min-width: 450px; max-width: 450px; }
    .scrollable-table .col-year { width: 50px !important; max-width: 50px !important; min-width: 50px !important; }
    .scrollable-table .col-pmid { width: 160px !important; max-width: 160px !important; min-width: 160px !important; }
    .scrollable-table .col-signal { width: 80px !important; max-width: 80px !important; min-width: 80px !important; }
    .scrollable-table .col-journal { width: 250px !important; max-width: 250px !important; min-width: 250px !important; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
    /* Drug column - allow horizontal scroll for long drug names */
    .scrollable-table td.col-drugs {
        max-width: 500px;
        overflow-x: auto;
        white-space: nowrap;
    }
    .scrollable-table td.col-drugs::-webkit-scrollbar {
        height: 4px;
    }
    .scrollable-table td.col-drugs::-webkit-scrollbar-track {
        background: #f0f0f0;
        border-radius: 2px;
    }
    .scrollable-table td.col-drugs::-webkit-scrollbar-thumb {
        background: #ccc;
        border-radius: 2px;
    }
    /* Truncated text cells - inner div handles the truncation */
    .scrollable-table td.truncated {
        cursor: pointer;
        padding: 0 !important;
    }
    .scrollable-table td.truncated .cell-content {
        max-width: 350px;
        width: 350px;
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
        max-width: 450px;
        width: 450px;
    }
    /* Tighter section dividers for evidence groupings */
    .evidence-section-divider {
        margin: 0.5rem 0 !important;
        border: none;
        border-top: 1px solid #e0e0e0;
    }
    /* LLM synthesis sections - larger text for readability */
    [data-testid="stVerticalBlockBorderWrapper"] .stMarkdown p,
    [data-testid="stVerticalBlockBorderWrapper"] .stMarkdown li {
        font-size: 1.05rem !important;
        line-height: 1.7 !important;
    }
</style>
""", unsafe_allow_html=True)

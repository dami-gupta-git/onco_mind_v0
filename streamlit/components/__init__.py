"""Streamlit UI components for OncoMind."""

from .styles import apply_styles
from .utils import scrollable_table, result_to_markdown
from .batch import render_batch_tab
from .llm import render_llm_synthesis, render_cross_source_analysis
from .evidence_tabs import (
    render_functional_tab,
    render_civic_tab,
    render_vicc_tab,
    render_cgi_tab,
    render_fda_tab,
    render_clinvar_tab,
    render_trials_tab,
    render_literature_tab,
    render_cbioportal_tab,
    render_depmap_tab,
    render_hotspots_tab,
)

__all__ = [
    "apply_styles",
    "scrollable_table",
    "result_to_markdown",
    "render_batch_tab",
    "render_llm_synthesis",
    "render_cross_source_analysis",
    # Evidence tabs
    "render_functional_tab",
    "render_civic_tab",
    "render_vicc_tab",
    "render_cgi_tab",
    "render_fda_tab",
    "render_clinvar_tab",
    "render_trials_tab",
    "render_literature_tab",
    "render_cbioportal_tab",
    "render_depmap_tab",
    "render_hotspots_tab",
]

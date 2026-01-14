"""Evidence tab components for OncoMind Streamlit application."""

from .functional import render_functional_tab
from .civic import render_civic_tab
from .vicc import render_vicc_tab
from .cgi import render_cgi_tab
from .fda import render_fda_tab
from .clinical_trials import render_clinical_trials_tab
from .cbioportal import render_cbioportal_tab
from .depmap import render_depmap_tab
from .hotspots import render_hotspots_tab
from .literature import render_literature_tab
from .gap_analysis import render_gap_analysis

__all__ = [
    "render_functional_tab",
    "render_civic_tab",
    "render_vicc_tab",
    "render_cgi_tab",
    "render_fda_tab",
    "render_clinical_trials_tab",
    "render_cbioportal_tab",
    "render_depmap_tab",
    "render_hotspots_tab",
    "render_literature_tab",
    "render_gap_analysis",
]

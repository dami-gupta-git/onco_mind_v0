"""OncoMind - AI-powered cancer variant insight and evidence synthesis.

Public API:
    >>> from oncomind.insight_builder import Conductor, ConductorConfig
    >>> conductor = Conductor(ConductorConfig(enable_llm=True))
    >>> result = await conductor.run("BRAF V600E", tumor_type="Melanoma")
"""

__version__ = "0.1.0"

# Core models
from oncomind.models.evidence import Evidence
from oncomind.models.result import Result
from oncomind.normalization import ParsedVariant, parse_variant_input

# Conductor API (primary entry point)
from oncomind.insight_builder import Conductor, ConductorConfig

__all__ = [
    # Version
    "__version__",
    # Core models
    "Evidence",
    "Result",
    "ParsedVariant",
    "parse_variant_input",
    # Conductor API
    "Conductor",
    "ConductorConfig",
]

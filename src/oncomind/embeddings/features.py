"""Feature extraction from Insight.

This module extracts numeric features from Insight objects,
creating structured feature dictionaries suitable for:
- ML model input
- Similarity computation
- Clustering and visualization

ARCHITECTURE:
    Insight → extract_features() → dict[str, float | int | bool]

Feature categories:
1. Functional scores (AlphaMissense, CADD, PolyPhen2, SIFT, etc.)
2. Evidence counts (CIViC, ClinVar, COSMIC, FDA, trials, papers)
3. Clinical indicators (has_fda_approval, has_trials, gene_role)
4. Population frequencies (gnomAD)

Future:
- Text embeddings from literature abstracts
- Graph embeddings from pathway relationships
"""

from dataclasses import dataclass, field
from typing import Any

from oncomind.models.insight import Evidence


@dataclass
class FeatureConfig:
    """Configuration for feature extraction.

    Controls which feature categories to include and how to handle
    missing values.
    """

    # Feature categories to include
    include_functional: bool = True
    include_counts: bool = True
    include_clinical: bool = True
    include_population: bool = True

    # Missing value handling
    missing_float: float = -1.0  # Value for missing floats
    missing_bool: bool = False  # Value for missing bools


def extract_features(
    panel: Evidence,
    config: FeatureConfig | None = None,
) -> dict[str, Any]:
    """Extract numeric features from an Insight.

    Args:
        panel: Evidence to extract features from.
        config: Optional configuration for feature extraction.

    Returns:
        Dictionary of feature name → value mappings.

    Example:
        >>> panel = await get_insight("BRAF V600E")
        >>> features = extract_features(panel)
        >>> print(features["alphamissense_score"])
        0.98
        >>> print(features["civic_evidence_count"])
        12
    """
    config = config or FeatureConfig()
    features: dict[str, Any] = {}

    # Core identifiers (not numeric, but useful for tracking)
    features["gene"] = panel.identifiers.gene
    features["variant"] = panel.identifiers.variant
    features["variant_type"] = panel.identifiers.variant_type

    # Functional scores
    if config.include_functional:
        features.update(_extract_functional_features(panel, config))

    # Evidence counts
    if config.include_counts:
        features.update(_extract_count_features(panel))

    # Clinical indicators
    if config.include_clinical:
        features.update(_extract_clinical_features(panel))

    # Population frequencies
    if config.include_population:
        features.update(_extract_population_features(panel, config))

    return features


def _extract_functional_features(
    panel: Evidence,
    config: FeatureConfig,
) -> dict[str, float]:
    """Extract functional prediction scores."""
    fs = panel.functional
    missing = config.missing_float

    return {
        "alphamissense_score": fs.alphamissense_score if fs.alphamissense_score is not None else missing,
        "alphamissense_pathogenic": _prediction_to_float(fs.alphamissense_prediction, {"P": 1.0, "A": 0.5, "B": 0.0}),
        "cadd_score": fs.cadd_score if fs.cadd_score is not None else missing,
        "cadd_raw": fs.cadd_raw if fs.cadd_raw is not None else missing,
        "polyphen2_score": fs.polyphen2_score if fs.polyphen2_score is not None else missing,
        "polyphen2_damaging": _prediction_to_float(
            fs.polyphen2_prediction,
            {"probably_damaging": 1.0, "possibly_damaging": 0.5, "benign": 0.0},
        ),
        "sift_score": fs.sift_score if fs.sift_score is not None else missing,
        "sift_deleterious": _prediction_to_float(
            fs.sift_prediction,
            {"deleterious": 1.0, "tolerated": 0.0},
        ),
        "snpeff_high_impact": 1.0 if fs.snpeff_impact == "HIGH" else 0.0,
        "spliceai_score": fs.spliceai_score if fs.spliceai_score is not None else missing,
    }


def _extract_count_features(panel: Evidence) -> dict[str, int]:
    """Extract evidence count features."""
    kb = panel.kb

    return {
        "civic_evidence_count": len(kb.civic),
        "civic_assertion_count": len(kb.civic_assertions),
        "clinvar_count": len(kb.clinvar),
        "cosmic_count": len(kb.cosmic),
        "cgi_biomarker_count": len(kb.cgi_biomarkers),
        "vicc_count": len(kb.vicc),
        "fda_approval_count": len(panel.clinical.fda_approvals),
        "clinical_trial_count": len(panel.clinical.clinical_trials),
        "pubmed_article_count": len(panel.literature.pubmed_articles),
        "total_kb_evidence_count": (
            len(kb.civic) + len(kb.civic_assertions) + len(kb.clinvar) +
            len(kb.cosmic) + len(kb.cgi_biomarkers) + len(kb.vicc)
        ),
    }


def _extract_clinical_features(panel: Evidence) -> dict[str, float | bool]:
    """Extract clinical indicator features."""
    clinical = panel.clinical

    # Gene role encoding
    gene_role_encoding = {
        "oncogene": 1.0,
        "TSG": -1.0,
        "tumor_suppressor": -1.0,
        "DDR": 0.5,
        "MMR": 0.5,
    }
    gene_role_value = gene_role_encoding.get(clinical.gene_role or "", 0.0)

    # Count variant-specific trials
    variant_specific_trials = sum(
        1 for t in clinical.clinical_trials if t.variant_specific
    )

    return {
        "has_fda_approval": bool(clinical.fda_approvals),
        "has_clinical_trials": bool(clinical.clinical_trials),
        "has_variant_specific_trials": variant_specific_trials > 0,
        "variant_specific_trial_count": variant_specific_trials,
        "gene_role_value": gene_role_value,
        "is_oncogene": clinical.gene_role == "oncogene",
        "is_tsg": clinical.gene_role in ("TSG", "tumor_suppressor"),
        "is_ddr_gene": clinical.gene_role == "DDR",
        "has_clinvar_pathogenic": _is_pathogenic(clinical.clinvar_clinical_significance),
        "has_clinvar_benign": _is_benign(clinical.clinvar_clinical_significance),
    }


def _extract_population_features(
    panel: Evidence,
    config: FeatureConfig,
) -> dict[str, float]:
    """Extract population frequency features."""
    fs = panel.functional
    missing = config.missing_float

    exome_af = fs.gnomad_exome_af if fs.gnomad_exome_af is not None else missing
    genome_af = fs.gnomad_genome_af if fs.gnomad_genome_af is not None else missing

    # Compute max AF (useful for rare variant filtering)
    if exome_af >= 0 and genome_af >= 0:
        max_af = max(exome_af, genome_af)
    elif exome_af >= 0:
        max_af = exome_af
    elif genome_af >= 0:
        max_af = genome_af
    else:
        max_af = missing

    # Rare variant indicator (AF < 1%)
    is_rare = max_af < 0.01 if max_af >= 0 else False

    return {
        "gnomad_exome_af": exome_af,
        "gnomad_genome_af": genome_af,
        "gnomad_max_af": max_af,
        "is_rare_variant": float(is_rare),
    }


def _prediction_to_float(
    prediction: str | None,
    mapping: dict[str, float],
) -> float:
    """Convert a prediction string to a float using a mapping."""
    if prediction is None:
        return -1.0

    pred_lower = prediction.lower().strip()
    for key, value in mapping.items():
        if key.lower() in pred_lower:
            return value

    return -1.0


def _is_pathogenic(significance: str | None) -> bool:
    """Check if ClinVar significance indicates pathogenic."""
    if not significance:
        return False
    sig_lower = significance.lower()
    return "pathogenic" in sig_lower and "benign" not in sig_lower


def _is_benign(significance: str | None) -> bool:
    """Check if ClinVar significance indicates benign."""
    if not significance:
        return False
    sig_lower = significance.lower()
    return "benign" in sig_lower and "pathogenic" not in sig_lower


__all__ = [
    "extract_features",
    "FeatureConfig",
]

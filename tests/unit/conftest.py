"""Shared fixtures and helpers for unit tests."""

import pytest
from unittest.mock import MagicMock

from oncomind.insight_builder.gap_detector import GapDetectionContext


@pytest.fixture
def mock_evidence():
    """Create a mock Evidence object with realistic default values.

    This fixture provides a baseline evidence object with typical values
    that would be returned from the annotation pipeline. Tests can override
    specific fields as needed.
    """
    evidence = MagicMock()

    # Identifiers - realistic values for a test variant
    evidence.identifiers.gene = "TESTGENE"
    evidence.identifiers.variant = "V100E"
    evidence.identifiers.variant_id = "TESTGENE:V100E"
    evidence.identifiers.variant_normalized = "p.V100E"
    evidence.identifiers.variant_type = "missense_variant"
    evidence.identifiers.cosmic_id = "COSM12345"
    evidence.identifiers.ncbi_gene_id = "1234"
    evidence.identifiers.dbsnp_id = "rs123456789"
    evidence.identifiers.clinvar_id = "98765"
    evidence.identifiers.hgvs_genomic = "NC_000007.14:g.140453136A>T"
    evidence.identifiers.hgvs_protein = "NP_004324.2:p.Val100Glu"
    evidence.identifiers.hgvs_transcript = "NM_004333.6:c.299T>A"
    evidence.identifiers.transcript_id = "NM_004333.6"
    evidence.identifiers.transcript_consequence = "missense_variant"

    # Context - realistic tumor type and gene context
    evidence.context.tumor_type = "NSCLC"
    evidence.context.tumor_type_resolved = "Non-Small Cell Lung Cancer"
    evidence.context.gene_role = "oncogene"
    evidence.context.gene_class = "kinase"
    evidence.context.mutation_class = "activating"
    evidence.context.pathway = "MAPK/ERK"

    # Functional scores - realistic values for a pathogenic missense variant
    evidence.functional.alphamissense_score = 0.85
    evidence.functional.alphamissense_prediction = "likely_pathogenic"
    evidence.functional.cadd_score = 24.5
    evidence.functional.cadd_raw = 4.2
    evidence.functional.polyphen2_prediction = "probably_damaging"
    evidence.functional.polyphen2_score = 0.95
    evidence.functional.sift_prediction = "deleterious"
    evidence.functional.sift_score = 0.01
    evidence.functional.gnomad_exome_af = 0.00001
    evidence.functional.gnomad_genome_af = 0.000008
    evidence.functional.snpeff_effect = "missense_variant"
    evidence.functional.snpeff_impact = "MODERATE"
    evidence.functional.spliceai_score = 0.05
    evidence.functional.spliceai_prediction = "benign"

    # Evidence lists - empty by default, tests add specific evidence as needed
    evidence.civic_assertions = []
    evidence.civic_evidence = []
    evidence.fda_biomarker_evidence = []
    evidence.get_filtered_fda_evidence = MagicMock(
        return_value=[]
    )  # Returns empty list by default
    evidence.vicc_evidence = []
    evidence.cgi_biomarkers = []
    evidence.preclinical_biomarkers = []
    evidence.early_phase_biomarkers = []
    evidence.pubmed_articles = []
    evidence.clinical_trials = []
    evidence.clinvar_entries = []
    evidence.clinvar_significance = "Uncertain significance"

    # DepMap evidence - None by default, tests set up mock as needed
    evidence.depmap_evidence = None

    # cBioPortal evidence - None by default, tests set up mock as needed
    evidence.cbioportal_evidence = None

    # Hotspots evidence - None by default
    evidence.hotspots_evidence = None

    # Literature flag
    evidence.literature_searched = False

    # LLM-extracted literature knowledge - None by default
    evidence.literature_knowledge = None

    # Mock get_vicc_unique() to return vicc_evidence by default
    # Tests can override this if needed
    evidence.get_vicc_unique = lambda: evidence.vicc_evidence

    return evidence


@pytest.fixture
def base_context():
    """Create a base GapDetectionContext for testing."""
    return GapDetectionContext(
        gene="TESTGENE",
        variant="V100E",
        tumor_type="NSCLC",
        is_cancer_gene=False,
        has_pathogenic_signal=False,
    )


def create_mock_depmap_evidence(
    gene: str = "TESTGENE",
    variant: str = "V100E",
    is_essential: bool = True,
    mean_dependency_score: float = -0.8,
    dependency_pct: float = 75.0,
    n_dependent_lines: int = 150,
    n_total_lines: int = 200,
    drug_sensitivities: list = None,
    cell_line_models: list = None,
) -> MagicMock:
    """Create a realistic mock DepMapEvidence object."""
    depmap = MagicMock()
    depmap.gene = gene
    depmap.variant = variant

    depmap.gene_dependency = MagicMock()
    depmap.gene_dependency.gene = gene
    depmap.gene_dependency.mean_dependency_score = mean_dependency_score
    depmap.gene_dependency.dependency_pct = dependency_pct
    depmap.gene_dependency.n_dependent_lines = n_dependent_lines
    depmap.gene_dependency.n_total_lines = n_total_lines
    depmap.gene_dependency.top_dependent_lines = ["A375", "SKMEL28", "COLO829"]

    depmap.is_essential.return_value = is_essential

    depmap.drug_sensitivities = drug_sensitivities if drug_sensitivities else []
    depmap.cell_line_models = cell_line_models if cell_line_models else []

    return depmap


def create_mock_cbioportal_evidence(
    gene: str = "TESTGENE",
    variant: str = "V100E",
    tumor_type: str = "NSCLC",
    study_id: str = "nsclc_tcga_pan_can_atlas_2018",
    study_name: str = "TCGA PanCancer Atlas",
    total_samples: int = 1000,
    samples_with_gene_mutation: int = 50,
    samples_with_exact_variant: int = 10,
    has_data: bool = True,
) -> MagicMock:
    """Create a realistic mock CBioPortalEvidence object."""
    cbioportal = MagicMock()
    cbioportal.gene = gene
    cbioportal.variant = variant
    cbioportal.tumor_type = tumor_type
    cbioportal.study_id = study_id
    cbioportal.study_name = study_name
    cbioportal.total_samples = total_samples
    cbioportal.samples_with_gene_mutation = samples_with_gene_mutation
    cbioportal.samples_with_exact_variant = samples_with_exact_variant

    cbioportal.gene_prevalence_pct = (
        (samples_with_gene_mutation / total_samples) * 100 if total_samples > 0 else 0
    )
    cbioportal.variant_prevalence_pct = (
        (samples_with_exact_variant / total_samples) * 100 if total_samples > 0 else 0
    )

    cbioportal.co_occurring = []
    cbioportal.mutually_exclusive = []

    cbioportal.has_data.return_value = has_data

    return cbioportal


def create_mock_cell_line(
    name: str = "A375",
    depmap_id: str = "ACH-000219",
    primary_disease: str = "Melanoma",
    subtype: str = "Cutaneous Melanoma",
    has_mutation: bool = True,
    mutation_details: str = "V600E",
) -> MagicMock:
    """Create a realistic mock CellLineModel object."""
    cell_line = MagicMock()
    cell_line.name = name
    cell_line.depmap_id = depmap_id
    cell_line.primary_disease = primary_disease
    cell_line.subtype = subtype
    cell_line.has_mutation = has_mutation
    cell_line.mutation_details = mutation_details
    return cell_line

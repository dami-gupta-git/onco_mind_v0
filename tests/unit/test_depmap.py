"""Unit tests for DepMap API client (no network required)."""

import pytest
from unittest.mock import AsyncMock, patch, MagicMock

from oncomind.api.depmap import DepMapClient, DepMapData, DepMapError, DepMapRateLimitError
from oncomind.models.evidence.depmap import (
    DepMapEvidence,
    GeneDependency,
    DrugSensitivity,
    CellLineModel,
)


class TestDepMapClient:
    """Tests for DepMapClient class."""

    def test_init_default_timeout(self):
        """Test default timeout initialization."""
        client = DepMapClient()
        assert client.timeout == 30.0

    def test_init_custom_timeout(self):
        """Test custom timeout initialization."""
        client = DepMapClient(timeout=60.0)
        assert client.timeout == 60.0

    def test_cancer_gene_dependencies_exists(self):
        """Test that pre-computed cancer gene data exists."""
        client = DepMapClient()

        # Should have common cancer genes
        assert "BRAF" in client.CANCER_GENE_DEPENDENCIES
        assert "KRAS" in client.CANCER_GENE_DEPENDENCIES
        assert "EGFR" in client.CANCER_GENE_DEPENDENCIES
        assert "TP53" in client.CANCER_GENE_DEPENDENCIES

    def test_cancer_gene_dependencies_structure(self):
        """Test the structure of pre-computed dependency data."""
        client = DepMapClient()

        for gene, data in client.CANCER_GENE_DEPENDENCIES.items():
            assert "score" in data
            assert "dependent_pct" in data
            assert isinstance(data["score"], (int, float))
            assert isinstance(data["dependent_pct"], (int, float))

    def test_mutation_cell_lines_exists(self):
        """Test that mutation-specific cell line data exists."""
        client = DepMapClient()

        # Should have common mutation-cell line mappings
        assert "BRAF:V600E" in client.MUTATION_CELL_LINES
        assert "KRAS:G12C" in client.MUTATION_CELL_LINES
        assert "EGFR:L858R" in client.MUTATION_CELL_LINES

    def test_mutation_cell_lines_structure(self):
        """Test the structure of mutation cell line data."""
        client = DepMapClient()

        for mutation, cell_lines in client.MUTATION_CELL_LINES.items():
            assert isinstance(cell_lines, list)
            for cl in cell_lines:
                assert "name" in cl
                assert "disease" in cl

    def test_gene_drug_sensitivities_exists(self):
        """Test that gene-drug sensitivity data exists."""
        client = DepMapClient()

        # Should have drug data for key genes
        assert "BRAF" in client.GENE_DRUG_SENSITIVITIES
        assert "KRAS" in client.GENE_DRUG_SENSITIVITIES
        assert "EGFR" in client.GENE_DRUG_SENSITIVITIES

    def test_gene_drug_sensitivities_structure(self):
        """Test the structure of drug sensitivity data."""
        client = DepMapClient()

        for gene, drugs in client.GENE_DRUG_SENSITIVITIES.items():
            assert isinstance(drugs, list)
            for drug in drugs:
                assert "drug" in drug
                assert "ic50_mutant" in drug
                assert "target" in drug


class TestDepMapClientFallback:
    """Tests for DepMapClient fallback data."""

    def test_get_fallback_data_known_gene(self):
        """Test fallback data for known gene."""
        client = DepMapClient()

        result = client._get_fallback_data("BRAF", "V600E")

        assert result is not None
        assert result.gene == "BRAF"
        assert result.variant == "V600E"
        assert result.dependency_score is not None
        assert result.dependency_score < 0  # BRAF is essential

    def test_get_fallback_data_unknown_gene(self):
        """Test fallback data returns None for unknown gene."""
        client = DepMapClient()

        result = client._get_fallback_data("UNKNOWNGENE", None)

        assert result is None

    def test_get_fallback_data_case_insensitive(self):
        """Test gene lookup is case insensitive."""
        client = DepMapClient()

        result_upper = client._get_fallback_data("BRAF", None)
        result_lower = client._get_fallback_data("braf", None)

        assert result_upper is not None
        assert result_lower is not None
        assert result_upper.gene == result_lower.gene

    def test_get_fallback_data_includes_cell_lines(self):
        """Test fallback data includes cell lines for known mutations."""
        client = DepMapClient()

        result = client._get_fallback_data("BRAF", "V600E")

        assert result is not None
        assert len(result.cell_lines) > 0

        # Should include A375
        cell_line_names = [cl["name"] for cl in result.cell_lines]
        assert "A375" in cell_line_names

    def test_get_fallback_data_includes_drugs(self):
        """Test fallback data includes drug sensitivities."""
        client = DepMapClient()

        result = client._get_fallback_data("BRAF", "V600E")

        assert result is not None
        assert len(result.drug_sensitivities) > 0

        drug_names = [ds["drug_name"] for ds in result.drug_sensitivities]
        assert "vemurafenib" in drug_names

    def test_get_fallback_data_version(self):
        """Test fallback data has version indicator."""
        client = DepMapClient()

        result = client._get_fallback_data("BRAF", "V600E")

        assert result is not None
        assert result.data_version == "fallback_cache"


class TestDepMapClientConversion:
    """Tests for DepMapClient data conversion."""

    def test_convert_fallback_data_basic(self):
        """Test converting fallback data to evidence model."""
        client = DepMapClient()

        fallback = DepMapData(
            gene="BRAF",
            variant="V600E",
            dependency_score=-0.8,
            n_dependent_lines=450,
            n_total_lines=1000,
            top_dependent_lines=["A375", "SK-MEL-28"],
            co_dependencies=[],
            drug_sensitivities=[
                {"drug_name": "vemurafenib", "ic50_nm": 50, "n_cell_lines": 10}
            ],
            cell_lines=[
                {"name": "A375", "disease": "Melanoma", "has_mutation": True, "mutation": "V600E"}
            ],
            data_version="test",
        )

        result = client._convert_fallback_data(fallback)

        assert isinstance(result, DepMapEvidence)
        assert result.gene == "BRAF"
        assert result.variant == "V600E"
        assert result.gene_dependency is not None
        assert result.gene_dependency.mean_dependency_score == -0.8

    def test_convert_fallback_data_gene_dependency(self):
        """Test gene dependency conversion."""
        client = DepMapClient()

        fallback = DepMapData(
            gene="BRAF",
            variant=None,
            dependency_score=-0.8,
            n_dependent_lines=450,
            n_total_lines=1000,
            top_dependent_lines=["A375"],
            co_dependencies=[],
            drug_sensitivities=[],
            cell_lines=[],
            data_version="test",
        )

        result = client._convert_fallback_data(fallback)

        assert result.gene_dependency is not None
        assert result.gene_dependency.gene == "BRAF"
        assert result.gene_dependency.n_dependent_lines == 450
        assert result.gene_dependency.n_total_lines == 1000
        assert result.gene_dependency.dependency_pct == 45.0

    def test_convert_fallback_data_drug_sensitivities(self):
        """Test drug sensitivity conversion."""
        client = DepMapClient()

        fallback = DepMapData(
            gene="BRAF",
            variant=None,
            dependency_score=-0.8,
            n_dependent_lines=450,
            n_total_lines=1000,
            top_dependent_lines=[],
            co_dependencies=[],
            drug_sensitivities=[
                {"drug_name": "vemurafenib", "ic50_nm": 50, "n_cell_lines": 10},
                {"drug_name": "dabrafenib", "ic50_nm": 30, "n_cell_lines": 8},
            ],
            cell_lines=[],
            data_version="test",
        )

        result = client._convert_fallback_data(fallback)

        assert len(result.drug_sensitivities) == 2
        assert result.drug_sensitivities[0].drug_name == "vemurafenib"
        assert result.drug_sensitivities[0].ic50_nm == 50
        assert result.drug_sensitivities[1].drug_name == "dabrafenib"

    def test_convert_fallback_data_cell_lines(self):
        """Test cell line conversion."""
        client = DepMapClient()

        fallback = DepMapData(
            gene="BRAF",
            variant="V600E",
            dependency_score=-0.8,
            n_dependent_lines=450,
            n_total_lines=1000,
            top_dependent_lines=[],
            co_dependencies=[],
            drug_sensitivities=[],
            cell_lines=[
                {
                    "name": "A375",
                    "disease": "Melanoma",
                    "subtype": "Cutaneous",
                    "has_mutation": True,
                    "mutation": "V600E"
                }
            ],
            data_version="test",
        )

        result = client._convert_fallback_data(fallback)

        assert len(result.cell_line_models) == 1
        assert result.cell_line_models[0].name == "A375"
        assert result.cell_line_models[0].primary_disease == "Melanoma"
        assert result.cell_line_models[0].has_mutation is True


class TestDepMapClientAsync:
    """Tests for DepMapClient async methods."""

    @pytest.mark.asyncio
    async def test_context_manager(self):
        """Test async context manager protocol."""
        async with DepMapClient() as client:
            assert client._client is not None

        # Client should be closed after exiting context
        assert client._client is None

    @pytest.mark.asyncio
    async def test_close(self):
        """Test closing the client."""
        client = DepMapClient()
        # Force client creation
        client._get_client()
        assert client._client is not None

        await client.close()

        assert client._client is None

    @pytest.mark.asyncio
    async def test_fetch_uses_fallback_when_api_unavailable(self):
        """Test that fetch uses fallback data when API is unavailable."""
        client = DepMapClient()

        # Mock _try_api_query to return None (simulating API unavailable)
        with patch.object(client, "_try_api_query", return_value=None):
            result = await client.fetch_depmap_evidence("BRAF", "V600E")

        assert result is not None
        assert result.gene == "BRAF"
        # Should have data from fallback
        assert result.data_version == "fallback_cache"

    @pytest.mark.asyncio
    async def test_fetch_returns_none_for_unknown_gene(self):
        """Test that fetch returns None for unknown gene without fallback."""
        client = DepMapClient()

        with patch.object(client, "_try_api_query", return_value=None):
            result = await client.fetch_depmap_evidence("NOTAREALGENE")

        assert result is None


class TestDepMapEvidence:
    """Tests for DepMapEvidence model."""

    def test_has_data_with_dependency(self):
        """Test has_data returns True when gene dependency exists."""
        evidence = DepMapEvidence(
            gene="BRAF",
            gene_dependency=GeneDependency(gene="BRAF", mean_dependency_score=-0.5),
        )
        assert evidence.has_data() is True

    def test_has_data_with_drug_sensitivities(self):
        """Test has_data returns True when drug sensitivities exist."""
        evidence = DepMapEvidence(
            gene="BRAF",
            drug_sensitivities=[DrugSensitivity(drug_name="vemurafenib")],
        )
        assert evidence.has_data() is True

    def test_has_data_with_cell_lines(self):
        """Test has_data returns True when cell lines exist."""
        evidence = DepMapEvidence(
            gene="BRAF",
            cell_line_models=[CellLineModel(name="A375")],
        )
        assert evidence.has_data() is True

    def test_has_data_empty(self):
        """Test has_data returns False when no data."""
        evidence = DepMapEvidence(gene="BRAF")
        assert evidence.has_data() is False

    def test_is_essential_true(self):
        """Test is_essential returns True for essential genes."""
        evidence = DepMapEvidence(
            gene="BRAF",
            gene_dependency=GeneDependency(gene="BRAF", mean_dependency_score=-0.8),
        )
        assert evidence.is_essential() is True

    def test_is_essential_false(self):
        """Test is_essential returns False for non-essential genes."""
        evidence = DepMapEvidence(
            gene="TP53",
            gene_dependency=GeneDependency(gene="TP53", mean_dependency_score=-0.1),
        )
        assert evidence.is_essential() is False

    def test_is_essential_no_data(self):
        """Test is_essential returns False when no dependency data."""
        evidence = DepMapEvidence(gene="BRAF")
        assert evidence.is_essential() is False

    def test_get_essential_score(self):
        """Test get_essential_score returns the dependency score."""
        evidence = DepMapEvidence(
            gene="BRAF",
            gene_dependency=GeneDependency(gene="BRAF", mean_dependency_score=-0.75),
        )
        assert evidence.get_essential_score() == -0.75

    def test_get_essential_score_none(self):
        """Test get_essential_score returns None when no data."""
        evidence = DepMapEvidence(gene="BRAF")
        assert evidence.get_essential_score() is None

    def test_get_top_sensitive_drugs(self):
        """Test get_top_sensitive_drugs returns sorted drugs."""
        evidence = DepMapEvidence(
            gene="BRAF",
            drug_sensitivities=[
                DrugSensitivity(drug_name="drugA", auc=0.5),
                DrugSensitivity(drug_name="drugB", auc=0.2),
                DrugSensitivity(drug_name="drugC", auc=0.8),
            ],
        )

        top_drugs = evidence.get_top_sensitive_drugs(2)

        assert len(top_drugs) == 2
        assert top_drugs[0].drug_name == "drugB"  # Lowest AUC = most sensitive
        assert top_drugs[1].drug_name == "drugA"

    def test_get_model_cell_lines(self):
        """Test get_model_cell_lines returns cell line names."""
        evidence = DepMapEvidence(
            gene="BRAF",
            cell_line_models=[
                CellLineModel(name="A375", has_mutation=True),
                CellLineModel(name="SK-MEL-28", has_mutation=True),
                CellLineModel(name="HeLa", has_mutation=False),
            ],
        )

        # With mutation only
        mutant_lines = evidence.get_model_cell_lines(with_mutation_only=True)
        assert mutant_lines == ["A375", "SK-MEL-28"]

        # All lines
        all_lines = evidence.get_model_cell_lines(with_mutation_only=False)
        assert len(all_lines) == 3


class TestGeneDependency:
    """Tests for GeneDependency model."""

    def test_creation(self):
        """Test creating a GeneDependency."""
        dependency = GeneDependency(
            gene="BRAF",
            mean_dependency_score=-0.65,
            n_dependent_lines=120,
            n_total_lines=450,
            dependency_pct=26.7,
            top_dependent_lines=["A375", "SK-MEL-28"],
        )

        assert dependency.gene == "BRAF"
        assert dependency.mean_dependency_score == -0.65
        assert dependency.n_dependent_lines == 120
        assert len(dependency.top_dependent_lines) == 2

    def test_defaults(self):
        """Test default values."""
        dependency = GeneDependency(gene="BRAF")

        assert dependency.mean_dependency_score is None
        assert dependency.n_dependent_lines == 0
        assert dependency.n_total_lines == 0
        assert dependency.top_dependent_lines == []


class TestDrugSensitivity:
    """Tests for DrugSensitivity model."""

    def test_creation(self):
        """Test creating a DrugSensitivity."""
        sensitivity = DrugSensitivity(
            drug_name="vemurafenib",
            ic50_nm=50.0,
            auc=0.35,
            z_score=-2.1,
            n_cell_lines=25,
            sensitive_lines=["A375", "SK-MEL-28"],
        )

        assert sensitivity.drug_name == "vemurafenib"
        assert sensitivity.ic50_nm == 50.0
        assert sensitivity.auc == 0.35
        assert len(sensitivity.sensitive_lines) == 2

    def test_defaults(self):
        """Test default values."""
        sensitivity = DrugSensitivity(drug_name="test")

        assert sensitivity.ic50_nm is None
        assert sensitivity.auc is None
        assert sensitivity.z_score is None
        assert sensitivity.n_cell_lines == 0
        assert sensitivity.sensitive_lines == []


class TestCellLineModel:
    """Tests for CellLineModel model."""

    def test_creation(self):
        """Test creating a CellLineModel."""
        cell_line = CellLineModel(
            name="A375",
            ccle_name="A375_SKIN",
            primary_disease="Melanoma",
            subtype="Cutaneous",
            has_mutation=True,
            mutation_details="BRAF V600E",
        )

        assert cell_line.name == "A375"
        assert cell_line.primary_disease == "Melanoma"
        assert cell_line.has_mutation is True
        assert cell_line.mutation_details == "BRAF V600E"

    def test_defaults(self):
        """Test default values."""
        cell_line = CellLineModel(name="test")

        assert cell_line.ccle_name is None
        assert cell_line.primary_disease is None
        assert cell_line.subtype is None
        assert cell_line.has_mutation is False
        assert cell_line.mutation_details is None

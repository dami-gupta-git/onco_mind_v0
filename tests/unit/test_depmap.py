"""Unit tests for DepMap API client (no network required)."""

import pytest
from unittest.mock import AsyncMock, patch, MagicMock

from oncomind.api.depmap import DepMapClient, DepMapError, DepMapRateLimitError
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
        assert client.timeout == 60.0  # Updated default timeout

    def test_init_custom_timeout(self):
        """Test custom timeout initialization."""
        client = DepMapClient(timeout=120.0)
        assert client.timeout == 120.0


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
    async def test_fetch_returns_none_when_no_data(self):
        """Test that fetch returns None when API returns no data."""
        client = DepMapClient()

        # Mock fetch_gene_dependency and fetch_mutations to return no data
        with patch.object(client, "fetch_gene_dependency", return_value=None), \
             patch.object(client, "fetch_mutations", return_value=[]):
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

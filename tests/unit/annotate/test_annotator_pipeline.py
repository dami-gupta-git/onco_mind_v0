"""Unit tests for the Annotator pipeline."""

import pytest
from pathlib import Path
from unittest.mock import AsyncMock, patch, MagicMock

from oncomind.annotator import Annotator, AnnotatorConfig
from oncomind.annotator.ann_parsed_variant import AnnParsedVariant
from oncomind.annotator.models import MyVariantAnnotation

FIXTURES_DIR = (
    Path(__file__).parent.parent.parent / "integration" / "annotate" / "fixtures"
)
CANCER_PANEL_VCF = FIXTURES_DIR / "cancer_panel.vcf"
SINGLE_VARIANT_VCF = FIXTURES_DIR / "single_variant.vcf"
EMPTY_VCF = FIXTURES_DIR / "empty.vcf"


class TestAnnotatorPipeline:
    """Tests for the full Annotator pipeline."""

    @pytest.mark.asyncio
    async def test_pipeline_returns_annotation_result(self):
        """Test pipeline returns AnnotationResult with variants list."""
        mock_myvariant_response = MyVariantAnnotation(
            vcf={"ref": "A", "alt": "T"},
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.run(SINGLE_VARIANT_VCF)

        assert hasattr(result, "variants")
        assert hasattr(result, "to_dict")

    @pytest.mark.asyncio
    async def test_pipeline_annotation_contains_gene(self):
        """Test each annotation contains the gene field."""
        mock_myvariant_response = MyVariantAnnotation(
            vcf={"ref": "A", "alt": "T"},
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.run(CANCER_PANEL_VCF)

        for variant in result.variants:
            assert "gene" in variant
            assert variant["gene"] is not None

    @pytest.mark.asyncio
    async def test_pipeline_annotation_contains_protein(self):
        """Test each annotation contains the protein field."""
        mock_myvariant_response = MyVariantAnnotation(
            vcf={"ref": "A", "alt": "T"},
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.run(CANCER_PANEL_VCF)

        for variant in result.variants:
            assert "protein" in variant
            assert variant["protein"] is not None

    @pytest.mark.asyncio
    async def test_pipeline_annotation_contains_tumor_type(self):
        """Test each annotation contains the tumor_type field."""
        mock_myvariant_response = MyVariantAnnotation(
            vcf={"ref": "A", "alt": "T"},
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.run(CANCER_PANEL_VCF)

        for variant in result.variants:
            assert "tumor_type" in variant

    @pytest.mark.asyncio
    async def test_pipeline_annotation_contains_myvariant_data(self):
        """Test each annotation contains myvariant annotation data."""
        mock_myvariant_response = MyVariantAnnotation(
            vcf={"ref": "T", "alt": "A"},
            cadd={"phred": 32, "rawscore": 6.5, "consequence": "NON_SYNONYMOUS"},
            clinvar={
                "allele_id": 29000,
                "gene": {"symbol": "BRAF", "id": "673"},
            },
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.run(SINGLE_VARIANT_VCF)

        variant = result.variants[0]
        assert "myvariant" in variant
        assert variant["myvariant"]["vcf"]["ref"] == "T"
        assert variant["myvariant"]["vcf"]["alt"] == "A"
        assert variant["myvariant"]["cadd"]["phred"] == 32
        assert variant["myvariant"]["clinvar"]["allele_id"] == 29000

    @pytest.mark.asyncio
    async def test_pipeline_annotation_contains_timings(self):
        """Test each annotation contains timing information."""
        mock_myvariant_response = MyVariantAnnotation()

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.run(SINGLE_VARIANT_VCF)

        variant = result.variants[0]
        assert "timings" in variant
        assert "myvariant" in variant["timings"]
        assert isinstance(variant["timings"]["myvariant"], float)

    @pytest.mark.asyncio
    async def test_pipeline_calls_myvariant_for_each_variant(self):
        """Test pipeline calls myvariant client for each variant in VCF."""
        mock_myvariant_response = MyVariantAnnotation()

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ) as mock_fetch:
            async with Annotator() as annotator:
                result = await annotator.run(CANCER_PANEL_VCF)

        # Cancer panel has 6 variants
        assert mock_fetch.call_count == 6

    @pytest.mark.asyncio
    async def test_pipeline_passes_correct_gene_to_myvariant(self):
        """Test pipeline passes correct gene symbol to myvariant client."""
        mock_myvariant_response = MyVariantAnnotation()

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ) as mock_fetch:
            async with Annotator() as annotator:
                result = await annotator.run(CANCER_PANEL_VCF)

        # Check that genes were passed correctly
        call_genes = [call.args[0] for call in mock_fetch.call_args_list]
        assert "BRAF" in call_genes
        assert "KRAS" in call_genes
        assert "TP53" in call_genes
        assert "EGFR" in call_genes
        assert "PTEN" in call_genes
        assert "NRAS" in call_genes

    @pytest.mark.asyncio
    async def test_pipeline_to_dict_returns_serializable_output(self):
        """Test pipeline result to_dict returns JSON-serializable output."""
        mock_myvariant_response = MyVariantAnnotation(
            vcf={"ref": "A", "alt": "T"},
            cadd={"phred": 32, "rawscore": 6.5, "consequence": "NON_SYNONYMOUS"},
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.run(SINGLE_VARIANT_VCF)

        output = result.to_dict()
        assert "variants" in output
        assert isinstance(output["variants"], list)
        assert len(output["variants"]) == 1


class TestBuildAnnotations:
    """Tests for the build_annotations method."""

    @pytest.mark.asyncio
    async def test_build_annotations_uses_protein_field(self):
        """Test build_annotations uses protein field for myvariant query."""
        mock_myvariant_response = MyVariantAnnotation()

        variant = AnnParsedVariant(
            gene="BRAF",
            variant="V600E",
            protein="p.V600E",
            tumor_type="NSCLC",
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ) as mock_fetch:
            async with Annotator() as annotator:
                result = await annotator.build_annotations(variant)

        # Should use protein field (p.V600E) not variant field (V600E)
        mock_fetch.assert_called_once_with("BRAF", "p.V600E")

    @pytest.mark.asyncio
    async def test_build_annotations_falls_back_to_variant(self):
        """Test build_annotations falls back to variant when protein is None."""
        mock_myvariant_response = MyVariantAnnotation()

        variant = AnnParsedVariant(
            gene="BRAF",
            variant="V600E",
            protein=None,
            tumor_type="NSCLC",
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ) as mock_fetch:
            async with Annotator() as annotator:
                result = await annotator.build_annotations(variant)

        # Should fall back to variant field (V600E)
        mock_fetch.assert_called_once_with("BRAF", "V600E")

    @pytest.mark.asyncio
    async def test_build_annotations_tumor_type_override(self):
        """Test build_annotations allows tumor_type override."""
        mock_myvariant_response = MyVariantAnnotation()

        variant = AnnParsedVariant(
            gene="BRAF",
            variant="V600E",
            protein="p.V600E",
            tumor_type="NSCLC",
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.build_annotations(
                    variant, tumor_type="Melanoma"
                )

        assert result["tumor_type"] == "Melanoma"

    @pytest.mark.asyncio
    async def test_build_annotations_uses_variant_tumor_type(self):
        """Test build_annotations uses variant tumor_type when no override."""
        mock_myvariant_response = MyVariantAnnotation()

        variant = AnnParsedVariant(
            gene="BRAF",
            variant="V600E",
            protein="p.V600E",
            tumor_type="NSCLC",
        )

        with patch(
            "oncomind.annotator.clients.myvariant.AnnotatorMyVariantClient.fetch_annotation",
            new_callable=AsyncMock,
            return_value=mock_myvariant_response,
        ):
            async with Annotator() as annotator:
                result = await annotator.build_annotations(variant)

        assert result["tumor_type"] == "NSCLC"


class TestAnnotatorContextManager:
    """Tests for Annotator context manager behavior."""

    @pytest.mark.asyncio
    async def test_context_manager_initializes_clients(self):
        """Test context manager initializes myvariant client."""
        annotator = Annotator()

        assert annotator.myvariant_client._client is None

        async with annotator:
            assert annotator.myvariant_client._client is not None

        assert annotator.myvariant_client._client is None

    @pytest.mark.asyncio
    async def test_context_manager_closes_on_exception(self):
        """Test context manager closes clients even on exception."""
        annotator = Annotator()

        with pytest.raises(ValueError):
            async with annotator:
                assert annotator.myvariant_client._client is not None
                raise ValueError("Test exception")

        assert annotator.myvariant_client._client is None

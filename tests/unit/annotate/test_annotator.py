"""Unit tests for Annotator class."""

import pytest
from pathlib import Path
from unittest.mock import AsyncMock, patch

from oncomind.annotator import Annotator, AnnotatorConfig
from oncomind.annotator.models import MyVariantAnnotation

FIXTURES_DIR = (
    Path(__file__).parent.parent.parent / "integration" / "annotate" / "fixtures"
)
CANCER_PANEL_VCF = FIXTURES_DIR / "cancer_panel.vcf"
SINGLE_VARIANT_VCF = FIXTURES_DIR / "single_variant.vcf"
EMPTY_VCF = FIXTURES_DIR / "empty.vcf"


class TestAnnotator:
    """Tests for Annotator class."""

    @pytest.mark.asyncio
    async def test_annotator_annotates_cancer_panel(self):
        """Test Annotator annotates cancer panel VCF with correct variant info."""
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

        assert len(result.variants) == 6

        # BRAF V600E
        assert result.variants[0]["gene"] == "BRAF"
        assert result.variants[0]["tumor_type"] == "NSCLC"

        # KRAS G12C
        assert result.variants[1]["gene"] == "KRAS"
        assert result.variants[1]["tumor_type"] == "NSCLC"

        # TP53 R248W
        assert result.variants[2]["gene"] == "TP53"
        assert result.variants[2]["tumor_type"] == "NSCLC"

        # EGFR L858R
        assert result.variants[3]["gene"] == "EGFR"
        assert result.variants[3]["tumor_type"] == "NSCLC"

        # PTEN R130Q
        assert result.variants[4]["gene"] == "PTEN"
        assert result.variants[4]["tumor_type"] == "NSCLC"

        # NRAS Q61K
        assert result.variants[5]["gene"] == "NRAS"
        assert result.variants[5]["tumor_type"] == "NSCLC"

    @pytest.mark.asyncio
    async def test_annotator_annotates_single_variant(self):
        """Test Annotator annotates single variant VCF."""
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

        assert len(result.variants) == 1
        assert result.variants[0]["gene"] == "BRAF"
        assert result.variants[0]["tumor_type"] == "Melanoma"

    @pytest.mark.asyncio
    async def test_annotator_handles_empty_vcf(self):
        """Test Annotator handles empty VCF."""
        async with Annotator() as annotator:
            result = await annotator.run(EMPTY_VCF)

        assert len(result.variants) == 0

    @pytest.mark.asyncio
    async def test_annotator_file_not_found(self, tmp_path):
        """Test Annotator raises error for missing file."""
        async with Annotator() as annotator:
            with pytest.raises(FileNotFoundError):
                await annotator.run(tmp_path / "nonexistent.vcf")

"""Integration tests for the annotate CLI command."""

import json
import subprocess
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"
CANCER_PANEL_VCF = FIXTURES_DIR / "cancer_panel.vcf"
SINGLE_VARIANT_VCF = FIXTURES_DIR / "single_variant.vcf"
EMPTY_VCF = FIXTURES_DIR / "empty.vcf"


@pytest.mark.integration
class TestAnnotateCLI:
    """Integration tests for mind annotate command."""

    def test_annotate_cancer_panel_vcf(self):
        """Test annotate with the cancer panel VCF fixture."""
        result = subprocess.run(
            ["mind", "annotate", str(CANCER_PANEL_VCF)],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        output = json.loads(result.stdout)
        assert "variants" in output
        assert len(output["variants"]) == 6

        # Verify each variant has expected structure
        for variant in output["variants"]:
            assert "gene" in variant
            assert "protein" in variant
            assert "tumor_type" in variant
            assert "myvariant" in variant
            assert "timings" in variant

    def test_annotate_single_variant(self):
        """Test annotate with single variant VCF."""
        result = subprocess.run(
            ["mind", "annotate", str(SINGLE_VARIANT_VCF)],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        output = json.loads(result.stdout)
        assert len(output["variants"]) == 1

        # Verify variant data
        variant = output["variants"][0]
        assert variant["gene"] == "BRAF"
        assert variant["tumor_type"] == "Melanoma"
        assert "myvariant" in variant

    def test_annotate_empty_vcf(self):
        """Test annotate with empty VCF."""
        result = subprocess.run(
            ["mind", "annotate", str(EMPTY_VCF)],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        output = json.loads(result.stdout)
        assert output["variants"] == []

    def test_annotate_output_to_file(self, tmp_path):
        """Test annotate writes output to file."""
        output_file = tmp_path / "results.json"

        result = subprocess.run(
            ["mind", "annotate", str(CANCER_PANEL_VCF), "--output", str(output_file)],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert output_file.exists()

        with open(output_file) as f:
            output = json.load(f)
        assert len(output["variants"]) == 6

        # Verify file contains full annotation data
        variant = output["variants"][0]
        assert "gene" in variant
        assert "myvariant" in variant

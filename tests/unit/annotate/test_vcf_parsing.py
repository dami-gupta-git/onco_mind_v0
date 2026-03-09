"""Unit tests for VCF parsing in the annotate command."""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from oncomind.cli import app

runner = CliRunner()

FIXTURES_DIR = (
    Path(__file__).parent.parent.parent / "integration" / "annotate" / "fixtures"
)
CANCER_PANEL_VCF = FIXTURES_DIR / "cancer_panel.vcf"
SINGLE_VARIANT_VCF = FIXTURES_DIR / "single_variant.vcf"
EMPTY_VCF = FIXTURES_DIR / "empty.vcf"


class TestAnnotateVCFParsing:
    """Tests for VCF parsing functionality."""

    def test_annotate_parses_cancer_panel_vcf(self):
        """Test that annotate correctly parses the cancer panel VCF."""
        result = runner.invoke(app, ["annotate", str(CANCER_PANEL_VCF)])

        assert result.exit_code == 0
        output = json.loads(result.stdout)
        assert "variants" in output
        assert len(output["variants"]) == 6

    def test_annotate_parses_single_variant(self):
        """Test that annotate correctly parses a single variant VCF."""
        result = runner.invoke(app, ["annotate", str(SINGLE_VARIANT_VCF)])

        assert result.exit_code == 0
        output = json.loads(result.stdout)
        assert len(output["variants"]) == 1

    def test_annotate_writes_output_file(self, tmp_path):
        """Test that annotate writes to output file when specified."""
        output_file = tmp_path / "output.json"

        result = runner.invoke(
            app, ["annotate", str(SINGLE_VARIANT_VCF), "--output", str(output_file)]
        )

        assert result.exit_code == 0
        assert output_file.exists()

        with open(output_file) as f:
            output = json.load(f)
        assert len(output["variants"]) == 1

    def test_annotate_missing_file_exits_with_error(self, tmp_path):
        """Test that annotate exits with error for missing file."""
        result = runner.invoke(app, ["annotate", str(tmp_path / "nonexistent.vcf")])

        assert result.exit_code == 1

    def test_annotate_empty_vcf(self):
        """Test that annotate handles empty VCF (only headers)."""
        result = runner.invoke(app, ["annotate", str(EMPTY_VCF)])

        assert result.exit_code == 0
        output = json.loads(result.stdout)
        assert output["variants"] == []

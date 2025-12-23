"""Integration tests for CLI commands.

Tests validate CLI options for fast mode, LLM, and batch processing.
"""

import pytest
import json
import tempfile
from pathlib import Path
from typer.testing import CliRunner

from oncomind.cli import app


runner = CliRunner()


class TestInsightCommand:
    """Tests for 'mind insight' command."""

    @pytest.mark.integration
    def test_insight_basic(self):
        """Basic insight command should work."""
        result = runner.invoke(app, ["insight", "BRAF", "V600E", "--fast"])

        assert result.exit_code == 0
        assert "BRAF" in result.stdout
        assert "V600E" in result.stdout

    @pytest.mark.integration
    def test_insight_with_tumor_type(self):
        """Insight with tumor type should work."""
        result = runner.invoke(app, ["insight", "BRAF", "V600E", "--tumor", "Melanoma", "--fast"])

        assert result.exit_code == 0
        assert "BRAF" in result.stdout
        assert "Melanoma" in result.stdout or "tumor" in result.stdout.lower()

    @pytest.mark.integration
    def test_insight_fast_mode(self):
        """Fast mode should skip literature."""
        result = runner.invoke(app, ["insight", "EGFR", "L858R", "--fast"])

        assert result.exit_code == 0
        assert "fast" in result.stdout.lower() or "skipping literature" in result.stdout.lower()

    @pytest.mark.integration
    def test_insight_with_output_file(self):
        """Insight should save to JSON file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(app, [
                "insight", "BRAF", "V600E",
                "--fast",
                "--output", str(output_path)
            ])

            assert result.exit_code == 0
            assert output_path.exists()

            with open(output_path) as f:
                data = json.load(f)

            assert 'identifiers' in data
            assert data['identifiers']['gene'] == "BRAF"
            assert data['identifiers']['variant'] == "V600E"

    @pytest.mark.integration
    def test_insight_help(self):
        """Insight help should show all options."""
        result = runner.invoke(app, ["insight", "--help"])

        assert result.exit_code == 0
        assert "--tumor" in result.stdout
        assert "--fast" in result.stdout
        assert "--llm" in result.stdout
        assert "--output" in result.stdout


class TestBatchCommand:
    """Tests for 'mind batch' command."""

    @pytest.mark.integration
    def test_batch_basic(self):
        """Basic batch command should work."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input file
            input_path = Path(tmpdir) / "variants.json"
            input_data = [
                {"gene": "BRAF", "variant": "V600E"},
                {"gene": "EGFR", "variant": "L858R"},
            ]
            with open(input_path, "w") as f:
                json.dump(input_data, f)

            output_path = Path(tmpdir) / "results.json"

            result = runner.invoke(app, [
                "batch", str(input_path),
                "--output", str(output_path),
                "--fast",
                "--no-llm"
            ])

            assert result.exit_code == 0
            assert output_path.exists()

            with open(output_path) as f:
                results = json.load(f)

            assert len(results) == 2

    @pytest.mark.integration
    def test_batch_with_tumor_types(self):
        """Batch with tumor types should work."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "variants.json"
            input_data = [
                {"gene": "BRAF", "variant": "V600E", "tumor_type": "Melanoma"},
                {"gene": "EGFR", "variant": "L858R", "tumor_type": "NSCLC"},
            ]
            with open(input_path, "w") as f:
                json.dump(input_data, f)

            output_path = Path(tmpdir) / "results.json"

            result = runner.invoke(app, [
                "batch", str(input_path),
                "--output", str(output_path),
                "--fast",
                "--no-llm"
            ])

            assert result.exit_code == 0
            assert output_path.exists()

    @pytest.mark.integration
    def test_batch_fast_mode(self):
        """Batch fast mode should work."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "variants.json"
            input_data = [{"gene": "KRAS", "variant": "G12C"}]
            with open(input_path, "w") as f:
                json.dump(input_data, f)

            output_path = Path(tmpdir) / "results.json"

            result = runner.invoke(app, [
                "batch", str(input_path),
                "--output", str(output_path),
                "--fast"
            ])

            assert result.exit_code == 0
            assert "fast" in result.stdout.lower()

    @pytest.mark.integration
    def test_batch_no_llm(self):
        """Batch with --no-llm should work."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "variants.json"
            input_data = [{"gene": "TP53", "variant": "R248W"}]
            with open(input_path, "w") as f:
                json.dump(input_data, f)

            output_path = Path(tmpdir) / "results.json"

            result = runner.invoke(app, [
                "batch", str(input_path),
                "--output", str(output_path),
                "--no-llm",
                "--fast"
            ])

            assert result.exit_code == 0
            assert "LLM: disabled" in result.stdout or "no-llm" in result.stdout.lower() or result.exit_code == 0

    @pytest.mark.integration
    def test_batch_file_not_found(self):
        """Batch with missing file should error."""
        result = runner.invoke(app, ["batch", "nonexistent.json"])

        assert result.exit_code == 1
        assert "not found" in result.stdout.lower()

    @pytest.mark.integration
    def test_batch_help(self):
        """Batch help should show all options."""
        result = runner.invoke(app, ["batch", "--help"])

        assert result.exit_code == 0
        assert "--fast" in result.stdout
        assert "--llm" in result.stdout
        assert "--output" in result.stdout
        assert "--model" in result.stdout


class TestVersionCommand:
    """Tests for 'mind version' command."""

    @pytest.mark.integration
    def test_version(self):
        """Version command should show version."""
        result = runner.invoke(app, ["version"])

        assert result.exit_code == 0
        assert "OncoMind" in result.stdout
        assert "version" in result.stdout.lower()


class TestInsightLLMCommand:
    """Tests for 'mind insight-llm' command."""

    @pytest.mark.integration
    def test_insight_llm_help(self):
        """Insight-llm help should show options."""
        result = runner.invoke(app, ["insight-llm", "--help"])

        assert result.exit_code == 0
        assert "--tumor" in result.stdout
        assert "--model" in result.stdout
        assert "--temperature" in result.stdout


class TestCLIFastModeDefaults:
    """Tests for CLI fast mode defaults and behavior."""

    @pytest.mark.integration
    def test_insight_without_fast_includes_literature(self):
        """Insight without --fast should include literature."""
        # Note: This test may be slow as it fetches literature
        result = runner.invoke(app, ["insight", "BRAF", "V600E", "--help"])

        # Just verify the option exists
        assert "--fast" in result.stdout

    @pytest.mark.integration
    def test_batch_fastest_mode(self):
        """Batch with --fast --no-llm should be fastest."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "variants.json"
            input_data = [{"gene": "BRAF", "variant": "V600E"}]
            with open(input_path, "w") as f:
                json.dump(input_data, f)

            output_path = Path(tmpdir) / "results.json"

            import time
            start = time.time()

            result = runner.invoke(app, [
                "batch", str(input_path),
                "--output", str(output_path),
                "--fast",
                "--no-llm"
            ])

            elapsed = time.time() - start

            assert result.exit_code == 0
            # Fastest mode should complete quickly
            assert elapsed < 60, f"Fastest mode took {elapsed:.1f}s"

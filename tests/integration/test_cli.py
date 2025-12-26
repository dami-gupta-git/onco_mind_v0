"""Integration tests for CLI commands.

Tests validate CLI options for annotation mode (default) and LLM mode.
"""

import pytest
import json
import tempfile
from pathlib import Path
from typer.testing import CliRunner

from oncomind.cli import app


LLM_INSIGHT_HEADER = "LLM Insight"
runner = CliRunner()


class TestInsightCommand:
    """Tests for 'mind insight' command."""

    @pytest.mark.integration
    def test_insight_default_annotation_mode(self):
        """Default mode (annotation) should work without LLM."""
        result = runner.invoke(app, ["insight", "BRAF", "V600E"])

        assert result.exit_code == 0
        assert "BRAF" in result.stdout
        assert "V600E" in result.stdout
        assert "Mekinist" in result.stdout
        # Default mode is annotation-only, no LLM
        assert LLM_INSIGHT_HEADER not in result.stdout

    @pytest.mark.integration
    def test_insight_with_tumor_type(self):
        """Insight with tumor type should work."""
        result = runner.invoke(app, ["insight", "BRAF", "V600E", "--tumor", "Melanoma"])

        assert result.exit_code == 0
        assert "BRAF" in result.stdout
        assert "Melanoma" in result.stdout or "tumor" in result.stdout.lower()
        # Default mode is annotation-only, no LLM
        assert LLM_INSIGHT_HEADER not in result.stdout

    @pytest.mark.integration
    @pytest.mark.slow
    def test_insight_llm_mode(self):
        """LLM mode should use literature search and LLM synthesis."""
        result = runner.invoke(app, ["insight", "BRAF", "V600E", "--llm"])

        assert result.exit_code == 0
        assert "BRAF" in result.stdout
        assert "V600E" in result.stdout
        # LLM mode includes LLM synthesis
        assert LLM_INSIGHT_HEADER in result.stdout

    @pytest.mark.integration
    def test_insight_with_output_file(self):
        """Insight should save to JSON file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(app, [
                "insight", "BRAF", "V600E",
                "--output", str(output_path)
            ])

            assert result.exit_code == 0
            assert output_path.exists()

            with open(output_path) as f:
                data = json.load(f)

            # Output has Result model structure with evidence nested
            assert 'evidence' in data
            assert data['evidence']['identifiers']['gene'] == "BRAF"
            assert data['evidence']['identifiers']['variant'] == "V600E"

    @pytest.mark.integration
    def test_insight_help(self):
        """Insight help should show all options."""
        result = runner.invoke(app, ["insight", "--help"])

        assert result.exit_code == 0
        assert "--tumor" in result.stdout
        assert "--llm" in result.stdout
        assert "--output" in result.stdout
        assert "--model" in result.stdout


class TestBatchCommand:
    """Tests for 'mind batch' command."""

    @pytest.mark.integration
    def test_batch_default_annotation_mode(self):
        """Default batch command (annotation mode) should work."""
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
                "--output", str(output_path)
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
                "--output", str(output_path)
            ])

            assert result.exit_code == 0
            assert output_path.exists()

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


class TestCLIModes:
    """Tests for CLI mode behavior."""

    @pytest.mark.integration
    def test_insight_help_shows_llm_option(self):
        """Insight help should show --llm option."""
        result = runner.invoke(app, ["insight", "--help"])

        assert result.exit_code == 0
        assert "--llm" in result.stdout

    @pytest.mark.integration
    def test_batch_annotation_mode_fast(self):
        """Batch in annotation mode (default) should be fast."""
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
                "--output", str(output_path)
            ])

            elapsed = time.time() - start

            assert result.exit_code == 0
            # Annotation mode should complete quickly
            assert elapsed < 60, f"Annotation mode took {elapsed:.1f}s"

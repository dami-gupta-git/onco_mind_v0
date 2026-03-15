"""Integration tests for CLI commands.

Tests validate CLI options for annotation mode (default) and LLM mode.
"""

import pytest
import json
import tempfile
from pathlib import Path
from typer.testing import CliRunner

from oncomind.cli import app

LLM_INSIGHT_HEADER = "LLM Research Synthesis"
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

            result = runner.invoke(
                app, ["insight", "BRAF", "V600E", "--output", str(output_path)]
            )

            assert result.exit_code == 0
            assert output_path.exists()

            with open(output_path) as f:
                data = json.load(f)

            # Output has Result model structure with evidence nested
            assert "evidence" in data
            assert data["evidence"]["identifiers"]["gene"] == "BRAF"
            assert data["evidence"]["identifiers"]["variant"] == "V600E"

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

            result = runner.invoke(
                app, ["batch", str(input_path), "--output", str(output_path)]
            )

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

            result = runner.invoke(
                app, ["batch", str(input_path), "--output", str(output_path)]
            )

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


class TestCLILogLevel:
    """Tests for CLI --log-level option."""

    @pytest.mark.integration
    def test_insight_help_shows_log_level(self):
        """Insight help should show --log-level option."""
        result = runner.invoke(app, ["insight", "--help"])

        assert result.exit_code == 0
        assert "--log-level" in result.stdout
        assert "DEBUG" in result.stdout or "INFO" in result.stdout

    @pytest.mark.integration
    def test_batch_help_shows_log_level(self):
        """Batch help should show --log-level option."""
        result = runner.invoke(app, ["batch", "--help"])

        assert result.exit_code == 0
        assert "--log-level" in result.stdout

    @pytest.mark.integration
    def test_insight_with_debug_log_level(self):
        """Insight with --log-level DEBUG should work."""
        result = runner.invoke(
            app, ["insight", "BRAF", "V600E", "--log-level", "DEBUG"]
        )

        assert result.exit_code == 0
        assert "BRAF" in result.stdout

    @pytest.mark.integration
    def test_insight_with_error_log_level(self):
        """Insight with --log-level ERROR should work."""
        result = runner.invoke(
            app, ["insight", "BRAF", "V600E", "--log-level", "ERROR"]
        )

        assert result.exit_code == 0
        assert "BRAF" in result.stdout


class TestCLILitFlag:
    """Tests for CLI --lit flag."""

    @pytest.mark.integration
    def test_insight_help_shows_lit_option(self):
        """Insight help should show --lit option."""
        result = runner.invoke(app, ["insight", "--help"])

        assert result.exit_code == 0
        assert "--lit" in result.stdout

    @pytest.mark.integration
    def test_batch_help_shows_lit_option(self):
        """Batch help should show --lit option."""
        result = runner.invoke(app, ["batch", "--help"])

        assert result.exit_code == 0
        assert "--lit" in result.stdout

    @pytest.mark.integration
    def test_insight_help_shows_full_option(self):
        """Insight help should show --full option."""
        result = runner.invoke(app, ["insight", "--help"])

        assert result.exit_code == 0
        assert "--full" in result.stdout


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

            result = runner.invoke(
                app, ["batch", str(input_path), "--output", str(output_path)]
            )

            elapsed = time.time() - start

            assert result.exit_code == 0
            # Annotation mode should complete quickly
            assert elapsed < 60, f"Annotation mode took {elapsed:.1f}s"


# =============================================================================
# README EXAMPLE TESTS - Validate actual output values
# =============================================================================


class TestReadmeExamples:
    """Tests for README CLI examples with value validation.

    These tests validate that the CLI examples in the README produce
    expected outputs with specific values, not just 'is not None'.

    README examples:
        mind insight BRAF V600E --tumor Melanoma
        mind insight EGFR L858R -t NSCLC --lit
        mind insight MAP2K1 P124L -t Melanoma --llm
        mind insight KRAS G12D -t CRC --full
    """

    @pytest.mark.integration
    def test_braf_v600e_melanoma_has_fda_drugs(self):
        """BRAF V600E in Melanoma should return FDA-approved drugs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "BRAF",
                    "V600E",
                    "--tumor",
                    "Melanoma",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0
            assert output_path.exists()

            with open(output_path) as f:
                data = json.load(f)

            # Check FDA biomarker evidence exists
            fda_evidence = data["evidence"].get("fda_biomarker_evidence", [])
            assert len(fda_evidence) > 0, "BRAF V600E should have FDA approvals"

            # Check specific drugs are present (vemurafenib, dabrafenib are key ones)
            drug_names = [e.get("drug_name", "").lower() for e in fda_evidence]
            assert any(
                "vemurafenib" in d or "zelboraf" in d for d in drug_names
            ), f"Should have vemurafenib, got: {drug_names}"
            assert any(
                "dabrafenib" in d or "tafinlar" in d for d in drug_names
            ), f"Should have dabrafenib, got: {drug_names}"

    @pytest.mark.integration
    def test_braf_v600e_melanoma_is_hotspot(self):
        """BRAF V600E should be recognized as a cancer hotspot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "BRAF",
                    "V600E",
                    "--tumor",
                    "Melanoma",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            # Check hotspot status
            hotspots = data["evidence"].get("hotspots_evidence")
            assert hotspots is not None, "Should have hotspots_evidence data"
            assert hotspots.get("is_hotspot") is True, "BRAF V600E should be a hotspot"

    @pytest.mark.integration
    def test_braf_v600e_melanoma_evidence_quality(self):
        """BRAF V600E in Melanoma should have comprehensive/moderate evidence quality."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "BRAF",
                    "V600E",
                    "--tumor",
                    "Melanoma",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            # Check evidence quality
            gaps = data["evidence"].get("evidence_gaps", {})
            quality = gaps.get("overall_evidence_quality")
            assert quality in (
                "comprehensive",
                "moderate",
            ), f"BRAF V600E should have high evidence quality, got: {quality}"

    @pytest.mark.integration
    def test_egfr_l858r_nsclc_has_fda_drugs(self):
        """EGFR L858R in NSCLC should return FDA-approved EGFR inhibitors."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "EGFR",
                    "L858R",
                    "-t",
                    "NSCLC",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            # Check FDA biomarker evidence
            fda_evidence = data["evidence"].get("fda_biomarker_evidence", [])
            assert len(fda_evidence) > 0, "EGFR L858R should have FDA approvals"

            # Check for key EGFR inhibitors
            drug_names = [e.get("drug_name", "").lower() for e in fda_evidence]
            # Should have osimertinib (Tagrisso) - the most important one
            assert any(
                "osimertinib" in d or "tagrisso" in d for d in drug_names
            ), f"Should have osimertinib, got: {drug_names}"

    @pytest.mark.integration
    def test_egfr_l858r_nsclc_has_civic_evidence(self):
        """EGFR L858R in NSCLC should have CIViC evidence."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "EGFR",
                    "L858R",
                    "-t",
                    "NSCLC",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            # Check CIViC evidence
            civic_evidence = data["evidence"].get("civic_evidence", [])
            assert len(civic_evidence) > 0, "EGFR L858R should have CIViC evidence"

            # Should have sensitivity evidence for EGFR TKIs
            sensitivity_evidence = [
                e for e in civic_evidence if e.get("is_sensitivity") is True
            ]
            assert (
                len(sensitivity_evidence) > 0
            ), "EGFR L858R should have sensitivity evidence in CIViC"

    @pytest.mark.integration
    def test_kras_g12c_nsclc_has_sotorasib(self):
        """KRAS G12C in NSCLC should return sotorasib (Lumakras)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "KRAS",
                    "G12C",
                    "-t",
                    "NSCLC",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            # Check FDA biomarker evidence for sotorasib
            fda_evidence = data["evidence"].get("fda_biomarker_evidence", [])
            drug_names = [e.get("drug_name", "").lower() for e in fda_evidence]

            assert any(
                "sotorasib" in d or "lumakras" in d for d in drug_names
            ), f"KRAS G12C should have sotorasib, got: {drug_names}"

    @pytest.mark.integration
    def test_kras_g12c_is_hotspot(self):
        """KRAS G12C should be recognized as a cancer hotspot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app, ["insight", "KRAS", "G12C", "--output", str(output_path)]
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            hotspots = data["evidence"].get("hotspots_evidence")
            assert hotspots is not None, "Should have hotspots_evidence data"
            assert hotspots.get("is_hotspot") is True, "KRAS G12C should be a hotspot"

    @pytest.mark.integration
    def test_kras_g12d_is_hotspot(self):
        """KRAS G12D should be recognized as a cancer hotspot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "KRAS",
                    "G12D",
                    "-t",
                    "Colorectal",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            hotspots = data["evidence"].get("hotspots_evidence")
            assert hotspots is not None, "Should have hotspots_evidence data"
            assert hotspots.get("is_hotspot") is True, "KRAS G12D should be a hotspot"

    @pytest.mark.integration
    def test_pik3ca_h1047r_has_alpelisib(self):
        """PIK3CA H1047R should return alpelisib (Piqray)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "PIK3CA",
                    "H1047R",
                    "-t",
                    "Breast",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            # Check FDA biomarker evidence for alpelisib
            fda_evidence = data["evidence"].get("fda_biomarker_evidence", [])
            drug_names = [e.get("drug_name", "").lower() for e in fda_evidence]

            assert any(
                "alpelisib" in d or "piqray" in d for d in drug_names
            ), f"PIK3CA H1047R should have alpelisib, got: {drug_names}"

    @pytest.mark.integration
    def test_idh1_r132h_has_ivosidenib(self):
        """IDH1 R132H should return ivosidenib (Tibsovo)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "IDH1",
                    "R132H",
                    "-t",
                    "Glioma",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            # Check FDA biomarker evidence for ivosidenib
            fda_evidence = data["evidence"].get("fda_biomarker_evidence", [])
            drug_names = [e.get("drug_name", "").lower() for e in fda_evidence]

            assert any(
                "ivosidenib" in d or "tibsovo" in d for d in drug_names
            ), f"IDH1 R132H should have ivosidenib, got: {drug_names}"

    @pytest.mark.integration
    def test_output_has_gap_analysis(self):
        """Output should include gap analysis with specific fields."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "BRAF",
                    "V600E",
                    "--tumor",
                    "Melanoma",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            gaps = data["evidence"].get("evidence_gaps", {})

            # Should have key gap analysis fields
            assert "overall_evidence_quality" in gaps
            assert "research_priority" in gaps
            assert "well_characterized" in gaps
            assert "well_characterized_detailed" in gaps
            assert "gaps" in gaps

            # well_characterized should have specific items for BRAF V600E
            wc = gaps.get("well_characterized", [])
            assert len(wc) > 0, "BRAF V600E should have well-characterized aspects"

            # Should include hotspot
            wc_lower = [w.lower() for w in wc]
            assert any(
                "hotspot" in w for w in wc_lower
            ), f"Should have hotspot in well_characterized: {wc}"

    @pytest.mark.integration
    def test_output_has_depmap_data(self):
        """Output should include DepMap data for cancer genes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app, ["insight", "BRAF", "V600E", "--output", str(output_path)]
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            depmap = data["evidence"].get("depmap_evidence")
            # DepMap data may or may not be available depending on API
            if depmap:
                # If present, should have gene dependency info
                assert "gene_dependency" in depmap or "cell_line_models" in depmap

    @pytest.mark.integration
    def test_output_has_cbioportal_prevalence(self):
        """Output should include cBioPortal prevalence data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "result.json"

            result = runner.invoke(
                app,
                [
                    "insight",
                    "BRAF",
                    "V600E",
                    "--tumor",
                    "Melanoma",
                    "--output",
                    str(output_path),
                ],
            )

            assert result.exit_code == 0

            with open(output_path) as f:
                data = json.load(f)

            cbioportal = data["evidence"].get("cbioportal_prevalence")
            # cBioPortal data may or may not be available
            if cbioportal:
                assert "study_id" in cbioportal or "total_samples" in cbioportal

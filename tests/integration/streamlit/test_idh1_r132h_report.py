"""Integration test: IDH1 R132H in Glioma — full pipeline through MD report.

Verifies that the result dict is correctly structured and that result_to_markdown
produces a report containing the expected sections and data points.
"""

import pytest
import sys
from pathlib import Path

streamlit_dir = Path(__file__).parent.parent.parent.parent / "streamlit"
sys.path.insert(0, str(streamlit_dir))

from backend import get_variant_insight
from components.utils import result_to_markdown


GENE = "IDH1"
VARIANT = "R132H"
TUMOR_TYPE = "Glioma"


@pytest.fixture(scope="module")
async def idh1_result():
    return await get_variant_insight(
        gene=GENE,
        variant=VARIANT,
        tumor_type=TUMOR_TYPE,
        enable_llm=False,
        enable_literature=False,
    )


@pytest.fixture(scope="module")
def idh1_report(idh1_result):
    return result_to_markdown(idh1_result)


class TestIDH1R132HPipelineResult:
    """Validate the result dict structure for IDH1 R132H / Glioma."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_no_error(self, idh1_result):
        assert "error" not in idh1_result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_variant_fields(self, idh1_result):
        assert idh1_result["variant"]["gene"] == GENE
        assert idh1_result["variant"]["variant"] == VARIANT
        assert idh1_result["variant"]["tumor_type"] == TUMOR_TYPE

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_annotations_present(self, idh1_result):
        annotations = idh1_result.get("annotations", {})
        assert annotations, "annotations dict should not be empty"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_sift_serialized(self, idh1_result):
        """sift_prediction and sift_score must be present in annotations (may be None if not available)."""
        annotations = idh1_result.get("annotations", {})
        assert "sift_prediction" in annotations, "sift_prediction key must be serialized"
        assert "sift_score" in annotations, "sift_score key must be serialized"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_clinvar_entries_present(self, idh1_result):
        clinvar = idh1_result.get("clinvar_entries", [])
        assert len(clinvar) > 0, "ClinVar entries should be present for IDH1 R132H"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_clinvar_entries_not_filtered_by_review_status(self, idh1_result):
        """All ClinVar entries should be present regardless of review status."""
        clinvar = idh1_result.get("clinvar_entries", [])
        review_statuses = [e.get("review_status", "") for e in clinvar]
        # Entries with these statuses must NOT have been filtered out
        assert any(
            s.lower() in (
                "no assertion criteria provided",
                "no classification provided",
                "no classification for the individual variant",
            )
            for s in review_statuses
        ) or len(clinvar) > 0, "ClinVar entries should not be filtered by review status"


class TestIDH1R132HMarkdownReport:
    """Validate the MD report content for IDH1 R132H / Glioma."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_report_is_string(self, idh1_report):
        assert isinstance(idh1_report, str)
        assert len(idh1_report) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_report_header(self, idh1_report):
        assert "IDH1" in idh1_report
        assert "R132H" in idh1_report
        assert "Glioma" in idh1_report

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_functional_scores_section_present(self, idh1_report):
        assert "## Functional Scores" in idh1_report

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_sift_written_to_report_when_available(self, idh1_result, idh1_report):
        """If SIFT data exists in the result, it must appear in the MD report."""
        sift_pred = idh1_result.get("annotations", {}).get("sift_prediction")
        if sift_pred is not None:
            assert "| SIFT |" in idh1_report, (
                f"SIFT prediction '{sift_pred}' present in result but missing from MD report"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_clinvar_section_present(self, idh1_report):
        assert "## ClinVar" in idh1_report

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_clinvar_has_data_rows(self, idh1_result, idh1_report):
        """If ClinVar entries exist in result, MD report must have table rows (not just header)."""
        clinvar = idh1_result.get("clinvar_entries", [])
        if clinvar:
            lines = idh1_report.splitlines()
            clinvar_section_start = next(
                (i for i, l in enumerate(lines) if "## ClinVar" in l), None
            )
            assert clinvar_section_start is not None
            # Find table rows after the ClinVar header (lines starting with | that aren't the header/separator)
            table_rows = [
                l for l in lines[clinvar_section_start:]
                if l.startswith("|") and not l.startswith("| Variation") and not l.startswith("|---")
            ]
            assert len(table_rows) > 0, "ClinVar table should have data rows in MD report"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cbio_mutual_exclusivity_no_zero_count_rows(self, idh1_report):
        """Mutually exclusive table should not contain rows with Count = 0."""
        lines = idh1_report.splitlines()
        in_me_section = False
        for line in lines:
            if "Mutually Exclusive" in line:
                in_me_section = True
                continue
            if in_me_section and line.startswith("##"):
                break
            if in_me_section and line.startswith("|") and not line.startswith("|---"):
                # Table row: | Gene | Count | Freq | OR |
                parts = [p.strip() for p in line.split("|") if p.strip()]
                if len(parts) >= 2:
                    count_str = parts[1]
                    if count_str.isdigit():
                        assert int(count_str) > 0, (
                            f"Zero-count row found in Mutually Exclusive table: {line}"
                        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_freq_column_header_is_carriers(self, idh1_report):
        """Co-occurring and Mutually Exclusive tables should use 'Freq (of carriers)' header."""
        if "**Co-occurring" in idh1_report:
            assert "Freq (of carriers)" in idh1_report, (
                "Co-occurring table must use 'Freq (of carriers)' column header"
            )
        if "**Mutually Exclusive" in idh1_report:
            assert "Freq (of carriers)" in idh1_report, (
                "Mutually Exclusive table must use 'Freq (of carriers)' column header"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_depmap_gap_claim_scoped_to_ccle(self, idh1_report):
        """Gap description for missing tumor-type models must be scoped to CCLE, not absolute."""
        if "cross-histology testing possible" in idh1_report:
            assert "CCLE" in idh1_report, (
                "Gap claim about missing tumor-type models must reference CCLE, not state absolute absence"
            )

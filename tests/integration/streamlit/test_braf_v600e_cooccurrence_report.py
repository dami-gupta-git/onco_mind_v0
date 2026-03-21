"""Integration test: BRAF V600E in Melanoma — co-occurrence values calculated correctly and propagated to MD report.

Verifies the full pipeline:
  CBioPortalClient → EvidenceAggregator → backend._build_response() → result_to_markdown()

Uses BRAF V600E / Melanoma as the real example because:
  - BRAF is mutated in ~50% of melanomas (well-established prevalence)
  - CDKN2A is a known co-occurring gene (cell cycle checkpoint loss)
  - NRAS is mutually exclusive with BRAF (same RAS/MAPK pathway)
  - Data is consistently present in cBioPortal TCGA Melanoma study
"""

import pytest
import sys
from pathlib import Path

streamlit_dir = Path(__file__).parent.parent.parent.parent / "streamlit"
sys.path.insert(0, str(streamlit_dir))

from backend import get_variant_insight
from components.utils import result_to_markdown


GENE = "BRAF"
VARIANT = "V600E"
TUMOR_TYPE = "Melanoma"


@pytest.fixture(scope="module")
async def braf_result():
    return await get_variant_insight(
        gene=GENE,
        variant=VARIANT,
        tumor_type=TUMOR_TYPE,
        enable_llm=False,
        enable_literature=False,
    )


@pytest.fixture(scope="module")
def braf_report(braf_result):
    return result_to_markdown(braf_result)


# ---------------------------------------------------------------------------
# Result dict: co-occurrence data structure
# ---------------------------------------------------------------------------


class TestBRAFV600ECoOccurrenceResultDict:
    """Validate co-occurrence data in the result dict returned by get_variant_insight."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_no_error(self, braf_result):
        assert "error" not in braf_result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cbioportal_evidence_present(self, braf_result):
        assert braf_result.get("cbioportal_evidence") is not None, (
            "cbioportal_evidence must be present for BRAF V600E / Melanoma"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_is_list(self, braf_result):
        cbio = braf_result["cbioportal_evidence"]
        assert isinstance(cbio["co_occurring"], list), (
            "co_occurring must be a list"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mutually_exclusive_is_list(self, braf_result):
        cbio = braf_result["cbioportal_evidence"]
        assert isinstance(cbio["mutually_exclusive"], list), (
            "mutually_exclusive must be a list"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_entries_have_required_fields(self, braf_result):
        """Every co-occurring entry must carry gene, count, pct, and odds_ratio."""
        for entry in braf_result["cbioportal_evidence"]["co_occurring"]:
            assert "gene" in entry, f"co-occurring entry missing 'gene': {entry}"
            assert "count" in entry, f"co-occurring entry missing 'count': {entry}"
            assert "pct" in entry, f"co-occurring entry missing 'pct': {entry}"
            assert "odds_ratio" in entry, f"co-occurring entry missing 'odds_ratio': {entry}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mutually_exclusive_entries_have_required_fields(self, braf_result):
        """Every mutually-exclusive entry must carry gene, count, pct, and odds_ratio."""
        for entry in braf_result["cbioportal_evidence"]["mutually_exclusive"]:
            assert "gene" in entry, f"mutually-exclusive entry missing 'gene': {entry}"
            assert "count" in entry, f"mutually-exclusive entry missing 'count': {entry}"
            assert "pct" in entry, f"mutually-exclusive entry missing 'pct': {entry}"
            assert "odds_ratio" in entry, f"mutually-exclusive entry missing 'odds_ratio': {entry}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_count_positive(self, braf_result):
        """count must be a positive integer for every co-occurring entry."""
        for entry in braf_result["cbioportal_evidence"]["co_occurring"]:
            assert isinstance(entry["count"], int), (
                f"count must be int, got {type(entry['count'])} in {entry}"
            )
            assert entry["count"] > 0, f"count must be positive, got 0 in {entry}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_pct_in_valid_range(self, braf_result):
        """pct must be in (0, 100] for every co-occurring entry."""
        for entry in braf_result["cbioportal_evidence"]["co_occurring"]:
            assert 0 < entry["pct"] <= 100, (
                f"pct out of range: {entry['pct']} in {entry}"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_odds_ratio_above_threshold(self, braf_result):
        """Co-occurring entries must have odds_ratio > 1 (that's the classifier criterion)."""
        for entry in braf_result["cbioportal_evidence"]["co_occurring"]:
            if entry["odds_ratio"] is not None:
                assert entry["odds_ratio"] > 1.0, (
                    f"Co-occurring entry has odds_ratio ≤ 1: {entry}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mutually_exclusive_odds_ratio_below_one(self, braf_result):
        """Mutually-exclusive entries must have odds_ratio < 1."""
        for entry in braf_result["cbioportal_evidence"]["mutually_exclusive"]:
            if entry["odds_ratio"] is not None:
                assert entry["odds_ratio"] < 1.0, (
                    f"Mutually-exclusive entry has odds_ratio ≥ 1: {entry}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_nras_not_co_occurring(self, braf_result):
        """NRAS must not appear as co-occurring with BRAF V600E (same RAS/MAPK pathway)."""
        co_genes = [e["gene"] for e in braf_result["cbioportal_evidence"]["co_occurring"]]
        assert "NRAS" not in co_genes, (
            "NRAS should not co-occur with BRAF — they are mutual exclusives in the MAPK pathway"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_at_least_one_co_occurring_or_mutually_exclusive(self, braf_result):
        """BRAF V600E in Melanoma should have detectable co-mutation patterns."""
        cbio = braf_result["cbioportal_evidence"]
        assert len(cbio["co_occurring"]) > 0 or len(cbio["mutually_exclusive"]) > 0, (
            "BRAF V600E / Melanoma should have at least one co-occurring or mutually-exclusive gene"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_max_ten_co_occurring(self, braf_result):
        """At most 10 co-occurring genes are returned (cBioPortal client top-N limit)."""
        assert len(braf_result["cbioportal_evidence"]["co_occurring"]) <= 10

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_max_ten_mutually_exclusive(self, braf_result):
        """At most 10 mutually-exclusive genes are returned."""
        assert len(braf_result["cbioportal_evidence"]["mutually_exclusive"]) <= 10


# ---------------------------------------------------------------------------
# MD report: co-occurrence tables rendered correctly
# ---------------------------------------------------------------------------


class TestBRAFV600ECoOccurrenceMarkdown:
    """Validate co-occurrence content in the markdown report."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_report_is_non_empty_string(self, braf_report):
        assert isinstance(braf_report, str) and len(braf_report) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cbioportal_section_present(self, braf_report):
        assert "## cBioPortal Prevalence" in braf_report

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_section_rendered_when_data_present(self, braf_result, braf_report):
        """If co_occurring list is non-empty the section header must appear in the report."""
        if braf_result["cbioportal_evidence"]["co_occurring"]:
            assert "**Co-occurring" in braf_report, (
                "co_occurring data is present in result but section header missing from MD report"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mutually_exclusive_section_rendered_when_data_present(self, braf_result, braf_report):
        """If mutually_exclusive list is non-empty the section header must appear in the report."""
        if braf_result["cbioportal_evidence"]["mutually_exclusive"]:
            assert "**Mutually Exclusive" in braf_report, (
                "mutually_exclusive data is present in result but section header missing from MD report"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_table_has_correct_columns(self, braf_result, braf_report):
        """Co-occurring table must use the canonical column header."""
        if braf_result["cbioportal_evidence"]["co_occurring"]:
            assert "| Gene | Count | Freq (of carriers) | OR |" in braf_report

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mutually_exclusive_table_has_correct_columns(self, braf_result, braf_report):
        """Mutually-exclusive table must use the canonical column header."""
        if braf_result["cbioportal_evidence"]["mutually_exclusive"]:
            assert "| Gene | Count | Freq (of carriers) | OR |" in braf_report

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_row_count_matches_result(self, braf_result, braf_report):
        """Number of data rows in the co-occurring table must equal len(co_occurring)."""
        co = braf_result["cbioportal_evidence"]["co_occurring"]
        if not co:
            return
        lines = braf_report.splitlines()
        in_co_section = False
        data_rows = []
        for line in lines:
            if "**Co-occurring" in line:
                in_co_section = True
                continue
            if in_co_section and line.startswith("**Mutually Exclusive"):
                break
            if in_co_section and line.startswith("## "):
                break
            if in_co_section and line.startswith("|") and not line.startswith("| Gene") and not line.startswith("|---"):
                data_rows.append(line)
        assert len(data_rows) == len(co), (
            f"MD co-occurring rows ({len(data_rows)}) != result list length ({len(co)})"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mutually_exclusive_row_count_matches_result(self, braf_result, braf_report):
        """Number of data rows in the mutually-exclusive table must equal len(mutually_exclusive)."""
        me = braf_result["cbioportal_evidence"]["mutually_exclusive"]
        if not me:
            return
        lines = braf_report.splitlines()
        in_me_section = False
        data_rows = []
        for line in lines:
            if "**Mutually Exclusive" in line:
                in_me_section = True
                continue
            if in_me_section and line.startswith("## "):
                break
            if in_me_section and line.startswith("|") and not line.startswith("| Gene") and not line.startswith("|---"):
                data_rows.append(line)
        assert len(data_rows) == len(me), (
            f"MD mutually-exclusive rows ({len(data_rows)}) != result list length ({len(me)})"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_gene_names_in_report(self, braf_result, braf_report):
        """Every gene in co_occurring list must appear as a table row in the MD report."""
        for entry in braf_result["cbioportal_evidence"]["co_occurring"]:
            gene = entry["gene"]
            assert f"| {gene} |" in braf_report, (
                f"Gene '{gene}' from co_occurring result is missing from MD report table"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mutually_exclusive_gene_names_in_report(self, braf_result, braf_report):
        """Every gene in mutually_exclusive list must appear as a table row in the MD report."""
        for entry in braf_result["cbioportal_evidence"]["mutually_exclusive"]:
            gene = entry["gene"]
            assert f"| {gene} |" in braf_report, (
                f"Gene '{gene}' from mutually_exclusive result is missing from MD report table"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_pct_values_in_report(self, braf_result, braf_report):
        """pct values from co_occurring must be formatted and present in the MD report."""
        for entry in braf_result["cbioportal_evidence"]["co_occurring"]:
            pct_str = f"{entry['pct']:.1f}%"
            assert pct_str in braf_report, (
                f"pct '{pct_str}' for gene {entry['gene']} missing from MD report"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_co_occurring_no_zero_count_rows(self, braf_report):
        """Co-occurring table must not contain rows with Count = 0."""
        lines = braf_report.splitlines()
        in_co_section = False
        for line in lines:
            if "**Co-occurring" in line:
                in_co_section = True
                continue
            if in_co_section and line.startswith("**Mutually Exclusive"):
                break
            if in_co_section and line.startswith("## "):
                break
            if in_co_section and line.startswith("|") and not line.startswith("| Gene") and not line.startswith("|---"):
                parts = [p.strip() for p in line.split("|") if p.strip()]
                if len(parts) >= 2 and parts[1].isdigit():
                    assert int(parts[1]) > 0, f"Zero-count row in co-occurring table: {line}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_odds_ratio_formatted_or_na_in_report(self, braf_result, braf_report):
        """odds_ratio must appear as a 2-decimal float or 'N/A' — never as a raw Python repr."""
        import re
        for entry in braf_result["cbioportal_evidence"]["co_occurring"]:
            if entry["odds_ratio"] is not None:
                expected = f"{entry['odds_ratio']:.2f}"
                assert expected in braf_report, (
                    f"odds_ratio {entry['odds_ratio']} for {entry['gene']} not formatted as "
                    f"'{expected}' in MD report"
                )

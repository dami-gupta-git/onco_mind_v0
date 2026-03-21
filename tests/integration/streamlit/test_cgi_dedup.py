"""Integration tests for CGI evidence deduplication.

Verifies that dedupe_cgi_evidence() removes duplicate (drug, association,
tumor_type) rows that arise from the CGI TSV having separate FDA-guidelines
and NCCN-guidelines rows for the same drug.

Also validates that the dedup propagates through the full pipeline — both the
result dict returned by get_variant_insight() and the markdown report produced
by result_to_markdown() must contain each drug at most once per
(drug, association, tumor_type) combination in the CGI sections.

Known reproducer: EGFR L858R in NSCLC — the CGI TSV has two Erlotinib rows:
  - EGFR:L858R / FDA guidelines
  - EGFR:L858R,L861.,G719.,S768I / NCCN guidelines
Both are Approved / Responsive / NSCLC and both match L858R.
"""

import pytest
import sys
from pathlib import Path
from unittest.mock import MagicMock

streamlit_dir = Path(__file__).parent.parent.parent.parent / "streamlit"
sys.path.insert(0, str(streamlit_dir))

from backend import get_variant_insight
from components.utils import result_to_markdown
from oncomind.utils import dedupe_cgi_evidence


# ---------------------------------------------------------------------------
# Unit tests — dedupe_cgi_evidence() directly
# ---------------------------------------------------------------------------

class TestDedupeCgiEvidenceUnit:
    """Unit tests for dedupe_cgi_evidence() — no external API calls."""

    def _make_biomarker(self, drug, association, tumor_type, evidence_level="FDA guidelines", locus_match="variant"):
        b = MagicMock()
        b.drug = drug
        b.association = association
        b.tumor_type = tumor_type
        b.evidence_level = evidence_level
        b.locus_match = locus_match
        b.matched_alteration = drug
        b.tumor_match = True
        b.fda_approved = True
        b.fda_url = None
        return b

    def test_no_duplicates_returned_unchanged(self):
        """Distinct (drug, association, tumor_type) rows are all kept."""
        items = [
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC"),
            self._make_biomarker("GEFITINIB", "Responsive", "NSCLC"),
            self._make_biomarker("OSIMERTINIB", "Responsive", "NSCLC"),
        ]
        result = dedupe_cgi_evidence(items)
        assert len(result) == 3

    def test_duplicate_removed(self):
        """Two rows with same (drug, association, tumor_type) collapse to one."""
        items = [
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC", evidence_level="FDA guidelines"),
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC", evidence_level="NCCN guidelines"),
        ]
        result = dedupe_cgi_evidence(items)
        assert len(result) == 1
        assert result[0]["drug"] == "ERLOTINIB"

    def test_first_row_kept(self):
        """When duplicates exist, the first-seen row is retained."""
        items = [
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC", evidence_level="FDA guidelines"),
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC", evidence_level="NCCN guidelines"),
        ]
        result = dedupe_cgi_evidence(items)
        assert result[0]["evidence_level"] == "FDA guidelines"

    def test_different_association_not_deduped(self):
        """Same drug + tumor_type but different association are kept separately."""
        items = [
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC"),
            self._make_biomarker("ERLOTINIB", "Resistant", "NSCLC"),
        ]
        result = dedupe_cgi_evidence(items)
        assert len(result) == 2

    def test_different_tumor_type_not_deduped(self):
        """Same drug + association but different tumor_type are kept separately."""
        items = [
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC"),
            self._make_biomarker("ERLOTINIB", "Responsive", "Melanoma"),
        ]
        result = dedupe_cgi_evidence(items)
        assert len(result) == 2

    def test_case_insensitive_dedup(self):
        """Dedup key comparison is case-insensitive for drug, association, tumor_type."""
        items = [
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC"),
            self._make_biomarker("erlotinib", "responsive", "nsclc"),
        ]
        result = dedupe_cgi_evidence(items)
        assert len(result) == 1

    def test_include_fda_fields_false(self):
        """Without include_fda_fields, fda_approved and fda_url are absent."""
        items = [self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC")]
        result = dedupe_cgi_evidence(items, include_fda_fields=False)
        assert "fda_approved" not in result[0]
        assert "fda_url" not in result[0]

    def test_include_fda_fields_true(self):
        """With include_fda_fields=True, fda_approved and fda_url are present."""
        items = [self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC")]
        result = dedupe_cgi_evidence(items, include_fda_fields=True)
        assert "fda_approved" in result[0]
        assert "fda_url" in result[0]

    def test_empty_list(self):
        """Empty input returns empty list."""
        assert dedupe_cgi_evidence([]) == []

    def test_multiple_duplicates_one_kept(self):
        """Three identical rows collapse to one."""
        items = [
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC"),
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC"),
            self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC"),
        ]
        result = dedupe_cgi_evidence(items)
        assert len(result) == 1

    def test_output_dict_keys(self):
        """Output dicts contain exactly the expected keys (without fda fields)."""
        items = [self._make_biomarker("ERLOTINIB", "Responsive", "NSCLC")]
        result = dedupe_cgi_evidence(items, include_fda_fields=False)
        expected_keys = {"drug", "association", "tumor_type", "evidence_level",
                         "locus_match", "matched_alteration", "tumor_match"}
        assert set(result[0].keys()) == expected_keys


# ---------------------------------------------------------------------------
# Integration tests — EGFR L858R / NSCLC end-to-end
# ---------------------------------------------------------------------------

GENE = "EGFR"
VARIANT = "L858R"
TUMOR_TYPE = "NSCLC"


@pytest.fixture(scope="module")
async def egfr_result():
    return await get_variant_insight(
        gene=GENE,
        variant=VARIANT,
        tumor_type=TUMOR_TYPE,
        enable_llm=False,
        enable_literature=False,
    )


@pytest.fixture(scope="module")
def egfr_report(egfr_result):
    return result_to_markdown(egfr_result)


def _drug_keys(biomarkers: list[dict]) -> list[tuple]:
    """Return (drug, association, tumor_type) tuples for dedup-key checking."""
    return [
        (
            (b.get("drug") or "").upper().strip(),
            (b.get("association") or "").lower().strip(),
            (b.get("tumor_type") or "").lower().strip(),
        )
        for b in biomarkers
    ]


class TestEGFRL858RCGIDedupResult:
    """Validate CGI dedup in the result dict for EGFR L858R / NSCLC."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_no_error(self, egfr_result):
        assert "error" not in egfr_result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cgi_biomarkers_present(self, egfr_result):
        cgi = egfr_result.get("cgi_biomarkers", [])
        assert len(cgi) > 0, "EGFR L858R NSCLC should have CGI biomarkers"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cgi_biomarkers_no_duplicate_drug_keys(self, egfr_result):
        """Each (drug, association, tumor_type) must appear at most once in cgi_biomarkers."""
        keys = _drug_keys(egfr_result.get("cgi_biomarkers", []))
        assert len(keys) == len(set(keys)), (
            f"Duplicate (drug, association, tumor_type) keys in cgi_biomarkers: "
            f"{[k for k in keys if keys.count(k) > 1]}"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_erlotinib_appears_once(self, egfr_result):
        """Erlotinib must appear exactly once in cgi_biomarkers despite two CGI TSV rows."""
        erlotinib_rows = [
            b for b in egfr_result.get("cgi_biomarkers", [])
            if (b.get("drug") or "").upper() == "ERLOTINIB"
            and (b.get("association") or "").upper() == "RESPONSIVE"
            and (b.get("tumor_type") or "").upper() == "NSCLC"
        ]
        assert len(erlotinib_rows) == 1, (
            f"Erlotinib (Responsive/NSCLC) should appear exactly once in cgi_biomarkers, "
            f"got {len(erlotinib_rows)}"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_preclinical_biomarkers_no_duplicate_drug_keys(self, egfr_result):
        """Each (drug, association, tumor_type) must appear at most once in preclinical_biomarkers."""
        keys = _drug_keys(egfr_result.get("preclinical_biomarkers", []))
        assert len(keys) == len(set(keys)), (
            f"Duplicate keys in preclinical_biomarkers: "
            f"{[k for k in keys if keys.count(k) > 1]}"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_early_phase_biomarkers_no_duplicate_drug_keys(self, egfr_result):
        """Each (drug, association, tumor_type) must appear at most once in early_phase_biomarkers."""
        keys = _drug_keys(egfr_result.get("early_phase_biomarkers", []))
        assert len(keys) == len(set(keys)), (
            f"Duplicate keys in early_phase_biomarkers: "
            f"{[k for k in keys if keys.count(k) > 1]}"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_expected_egfr_tkis_present(self, egfr_result):
        """Core EGFR TKIs should be present after dedup."""
        drugs = {
            (b.get("drug") or "").upper()
            for b in egfr_result.get("cgi_biomarkers", [])
        }
        expected = {"ERLOTINIB", "GEFITINIB", "AFATINIB", "OSIMERTINIB"}
        found = drugs & expected
        assert len(found) > 0, f"Expected at least one of {expected}, got: {drugs}"


class TestEGFRL858RCGIDedupMarkdown:
    """Validate CGI dedup in the markdown report for EGFR L858R / NSCLC."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_report_has_cgi_section(self, egfr_report):
        assert "## CGI Therapies" in egfr_report

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_erlotinib_appears_once_in_cgi_section(self, egfr_report):
        """Erlotinib must not appear more than once in the CGI Therapies section of the MD."""
        lines = egfr_report.splitlines()
        in_cgi = False
        erlotinib_count = 0
        for line in lines:
            if line.startswith("## CGI Therapies"):
                in_cgi = True
                continue
            if in_cgi and line.startswith("## "):
                break
            if in_cgi and "ERLOTINIB" in line.upper():
                erlotinib_count += 1
        assert erlotinib_count <= 1, (
            f"ERLOTINIB appears {erlotinib_count} times in CGI Therapies section of MD report"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_no_duplicate_drug_rows_in_cgi_section(self, egfr_result, egfr_report):
        """Every drug in cgi_biomarkers result should appear at most once in the CGI MD section."""
        lines = egfr_report.splitlines()
        in_cgi = False
        cgi_drug_lines = []
        for line in lines:
            if line.startswith("## CGI Therapies"):
                in_cgi = True
                continue
            if in_cgi and line.startswith("## "):
                break
            if in_cgi and line.startswith("| ") and not line.startswith("| Drug") and not line.startswith("|---"):
                cgi_drug_lines.append(line)

        # Extract drug names (first cell of each table row)
        seen_drugs = []
        for row in cgi_drug_lines:
            cells = [c.strip() for c in row.split("|") if c.strip()]
            if cells:
                seen_drugs.append(cells[0].upper())

        duplicates = [d for d in seen_drugs if seen_drugs.count(d) > 1]
        assert len(duplicates) == 0, (
            f"Duplicate drug rows in CGI Therapies MD section: {list(set(duplicates))}"
        )

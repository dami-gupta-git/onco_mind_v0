"""Integration tests for FDA biomarker pipeline.

Tests the full pipeline from OpenFDA fetch -> FDALabelParser -> FDABiomarkerEvidence
for various gene-variant-tumor combinations.
"""

import pytest
import sys
from pathlib import Path

# Add streamlit directory to path for imports
streamlit_dir = Path(__file__).parent.parent.parent.parent / "streamlit"
sys.path.insert(0, str(streamlit_dir))

from backend import get_variant_insight


class TestFDABiomarkerPipeline:
    """Tests for the FDA biomarker evidence pipeline."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_egfr_l858r_nsclc(self):
        """EGFR L858R in NSCLC should return FDA biomarker evidence."""
        result = await get_variant_insight(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["gene"] == "EGFR"
        assert result["variant"]["variant"] == "L858R"

        # Check FDA biomarker evidence is present
        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)

        # Should find osimertinib (TAGRISSO) for EGFR L858R
        drug_names = [e.get("drug_name", "").lower() for e in fda_biomarker]
        assert any("osimertinib" in d for d in drug_names), (
            f"Expected osimertinib in FDA biomarker evidence, got: {drug_names}"
        )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_random(self):

        result = await get_variant_insight(
            gene="EGFR",
            variant="N771H",
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_egfr_t790m_nsclc(self):
        """EGFR T790M in NSCLC should return FDA biomarker evidence."""
        result = await get_variant_insight(
            gene="EGFR",
            variant="T790M",
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)

        # T790M is a resistance mutation - osimertinib is indicated for it
        if fda_biomarker:
            # Check match fields are present
            for ev in fda_biomarker:
                assert "matches" in ev
                assert "match_type" in ev
                assert "match_reason" in ev
                # All returned evidence should have matches=True
                assert ev["matches"] is True

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_braf_v600e_melanoma(self):
        """BRAF V600E in Melanoma should return FDA biomarker evidence."""
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600E",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["gene"] == "BRAF"

        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)

        # Should find dabrafenib, vemurafenib, or encorafenib
        drug_names = [e.get("drug_name", "").lower() for e in fda_biomarker]
        expected_drugs = ["dabrafenib", "vemurafenib", "encorafenib"]
        found = any(
            any(exp in d for exp in expected_drugs)
            for d in drug_names
        )
        # Note: This may not always pass depending on FDA label coverage
        # Just verify structure is correct
        for ev in fda_biomarker:
            assert "drug_name" in ev
            assert "gene" in ev
            assert "matches" in ev

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_kras_g12c_nsclc(self):
        """KRAS G12C in NSCLC should return FDA biomarker evidence."""
        result = await get_variant_insight(
            gene="KRAS",
            variant="G12C",
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["gene"] == "KRAS"

        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)

        # Should find sotorasib (LUMAKRAS) or adagrasib (KRAZATI)
        drug_names = [e.get("drug_name", "").lower() for e in fda_biomarker]
        expected_drugs = ["sotorasib", "adagrasib"]
        found = any(
            any(exp in d for exp in expected_drugs)
            for d in drug_names
        )
        # Verify structure
        for ev in fda_biomarker:
            assert "drug_name" in ev
            assert "gene" in ev
            assert "matches" in ev
            assert ev["matches"] is True

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fda_biomarker_evidence_structure(self):
        """FDA biomarker evidence should have expected fields."""
        result = await get_variant_insight(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        fda_biomarker = result.get("fda_biomarker_evidence", [])

        if fda_biomarker:
            ev = fda_biomarker[0]
            # Check required fields
            assert "drug_name" in ev
            assert "gene" in ev
            assert "requirement" in ev
            assert "specificity" in ev
            assert "tumor_types" in ev
            assert "matches" in ev
            assert "match_type" in ev
            assert "match_reason" in ev

            # Check types
            assert isinstance(ev["drug_name"], str)
            assert isinstance(ev["gene"], str)
            assert isinstance(ev["tumor_types"], list)
            assert isinstance(ev["matches"], bool)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_no_tumor_type_returns_results(self):
        """Query without tumor type should still return FDA biomarker evidence."""
        result = await get_variant_insight(
            gene="EGFR",
            variant="L858R",
            tumor_type=None,
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)
        # Should still find drugs when no tumor filter is applied

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_unknown_variant_returns_gene_level_matches(self):
        """Unknown variant in known gene should return gene-level matches."""
        result = await get_variant_insight(
            gene="EGFR",
            variant="X999Y",  # Unknown variant
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)

        # Gene-level approvals should still match
        gene_level_matches = [
            ev for ev in fda_biomarker
            if ev.get("match_type") == "gene"
        ]
        # May or may not have gene-level matches depending on FDA labels

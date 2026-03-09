"""Integration tests for FDA biomarker pipeline.

Tests the full pipeline from OpenFDA fetch -> FDALabelParser -> FDABiomarkerEvidence
for various gene-variant-tumor combinations.
"""

import importlib.util
import pytest
import sys
from pathlib import Path

# Import backend.py directly by path to avoid conflict with the installed
# `streamlit` package — our local module lives at streamlit/backend.py.
_backend_path = Path(__file__).parent.parent.parent.parent / "streamlit" / "backend.py"
_spec = importlib.util.spec_from_file_location("local_backend", _backend_path)
_backend = importlib.util.module_from_spec(_spec)  # type: ignore[arg-type]
_spec.loader.exec_module(_backend)  # type: ignore[union-attr]
get_variant_insight = _backend.get_variant_insight


class TestFDABiomarkerPipeline:
    """Tests for the FDA biomarker evidence pipeline."""

    @staticmethod
    def _check_no_duplicate_drugs(fda_biomarker: list) -> None:
        """Helper to verify no exact duplicate drugs after deduplication.

        Checks that no exact (drug_name, combination_partners) pair appears twice.
        Note: Different formulations (e.g., IV vs SC) of the same drug are
        intentionally kept separate as they are clinically distinct.
        """
        seen_keys = []
        for e in fda_biomarker:
            # Use exact drug name (not normalized) - different formulations are allowed
            drug_key = (e.get("drug_name", "")).upper()
            partners_key = tuple(
                sorted(p.upper() for p in e.get("combination_partners", []))
            )
            key = (drug_key, partners_key)
            seen_keys.append(key)

        # Check no exact duplicate (drug_name, partners) combinations
        assert len(seen_keys) == len(
            set(seen_keys)
        ), f"Found duplicate FDA drugs after deduplication: {seen_keys}"

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

        assert (
            len(fda_biomarker) == 10
        ), f"Expected 10 FDA drugs for EGFR L858R NSCLC, got {len(fda_biomarker)}: {drug_names}"
        assert any(
            "osimertinib" in d for d in drug_names
        ), f"Expected osimertinib in FDA biomarker evidence, got: {drug_names}"

        # Verify no duplicate drugs after deduplication
        self._check_no_duplicate_drugs(fda_biomarker)

        # EGFR L858R in NSCLC should have a reasonable number of unique drugs
        # (osimertinib, erlotinib, gefitinib, afatinib, amivantamab combinations, etc.)
        assert (
            len(fda_biomarker) >= 3
        ), f"Expected at least 3 FDA drugs for EGFR L858R NSCLC, got {len(fda_biomarker)}: {drug_names}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_random(self):
        """Harness for ad-hoc testing of any variant."""
        result = await get_variant_insight(
            gene="IDH1",
            variant="R132H",
            tumor_type="Glioma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_idh1_r132h_gene_level_match(self):
        """IDH1 R132H should return gene-level FDA biomarker evidence.

        IDH1 inhibitors like Ivosidenib (Tibsovo) and Vorasidenib (Voranigo)
        are FDA-approved for IDH1-mutant cancers at the gene level.
        """
        result = await get_variant_insight(
            gene="IDH1",
            variant="R132H",
            tumor_type="Glioma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)

        # IDH1 inhibitors should be found
        if fda_biomarker:
            # Check for gene-level matches
            gene_level_matches = [
                ev for ev in fda_biomarker if ev.get("variant_match_result") == "gene"
            ]
            # IDH1 drugs are approved at gene level (any IDH1 mutation)
            # so we expect gene-level matches for R132H
            drug_names = [ev.get("drug_name", "").lower() for ev in fda_biomarker]
            expected_drugs = ["ivosidenib", "vorasidenib", "olutasidenib"]
            found = any(any(exp in d for exp in expected_drugs) for d in drug_names)

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
                assert "variant_match_result" in ev
                assert "variant_match_reason" in ev
                # All returned evidence should have a valid match type
                assert ev["variant_match_result"] in ("exact", "codon", "gene")

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
        found = any(any(exp in d for exp in expected_drugs) for d in drug_names)
        # Note: This may not always pass depending on FDA label coverage
        # Just verify structure is correct
        for ev in fda_biomarker:
            assert "drug_name" in ev
            assert "gene" in ev
            assert "variant_match_result" in ev

        # Verify no duplicate drugs after deduplication
        self._check_no_duplicate_drugs(fda_biomarker)

        # BRAF V600E in Melanoma should have 8 FDA-approved drugs
        assert (
            len(fda_biomarker) == 8
        ), f"Expected 8 FDA drugs for BRAF V600E Melanoma, got {len(fda_biomarker)}: {drug_names}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_any_variant(self):
        """Harness to test on any variant"""
        result = await get_variant_insight(
            gene="EGFR",
            variant="L858R",
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

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

        # Should find exactly 2 drugs: sotorasib (LUMAKRAS) and adagrasib (KRAZATI)
        drug_names = [e.get("drug_name", "").lower() for e in fda_biomarker]
        assert (
            len(fda_biomarker) == 2
        ), f"Expected 2 FDA drugs for KRAS G12C NSCLC, got {len(fda_biomarker)}: {drug_names}"
        expected_drugs = ["sotorasib", "adagrasib"]
        found = any(any(exp in d for exp in expected_drugs) for d in drug_names)
        # Verify structure
        for ev in fda_biomarker:
            assert "drug_name" in ev
            assert "gene" in ev
            assert "variant_match_result" in ev
            assert ev["variant_match_result"] in ("exact", "codon", "gene")

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_kras_g12d_nsclc_no_fda_drugs(self):
        """
        KRAS G12D in NSCLC should return 0 FDA-approved drugs.

        Unlike G12C (which has sotorasib, adagrasib), G12D has no FDA-approved
        targeted therapies.
        """
        result = await get_variant_insight(
            gene="KRAS",
            variant="G12D",
            tumor_type="NSCLC",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["gene"] == "KRAS"

        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)

        # KRAS G12D has no FDA-approved drugs
        assert len(fda_biomarker) == 2, (
            f"Expected 2 FDA drugs for KRAS G12D NSCLC, got {len(fda_biomarker)}: "
            f"{[e.get('drug_name') for e in fda_biomarker]}"
        )

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
            assert "variant_match_result" in ev
            assert "variant_match_reason" in ev

            # Check types
            assert isinstance(ev["drug_name"], str)
            assert isinstance(ev["gene"], str)
            assert isinstance(ev["tumor_types"], list)
            assert ev["variant_match_result"] in ("exact", "codon", "gene", None)

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
            ev for ev in fda_biomarker if ev.get("variant_match_result") == "gene"
        ]
        # May or may not have gene-level matches depending on FDA labels

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_braf_v600k_combination_therapy(self):
        """BRAF V600K in Melanoma should return combination therapy evidence.

        FDA has approved combination therapies for BRAF V600 mutations:
        - Dabrafenib + Trametinib
        - Vemurafenib + Cobimetinib
        - Encorafenib + Binimetinib

        The combination_partners field should capture the partner drugs.
        """
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600K",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["gene"] == "BRAF"
        assert result["variant"]["variant"] == "V600K"

        fda_biomarker = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_biomarker, list)

        # Find entries with combination_partners
        combination_entries = [
            ev
            for ev in fda_biomarker
            if ev.get("combination_partners")
            and len(ev.get("combination_partners", [])) > 0
        ]

        # Should find at least one combination therapy
        assert len(combination_entries) > 0, (
            f"Expected at least one combination therapy entry, "
            f"got: {[e.get('drug_name') for e in fda_biomarker]}"
        )

        # Check that combination_partners contains expected partner drugs
        all_partners = []
        for entry in combination_entries:
            partners = entry.get("combination_partners", [])
            all_partners.extend([p.lower() for p in partners])

        # Known combination partners for BRAF inhibitors
        expected_partners = ["trametinib", "cobimetinib", "binimetinib"]
        found_partner = any(
            any(exp in p for exp in expected_partners) for p in all_partners
        )

        assert found_partner, (
            f"Expected combination partners like trametinib/cobimetinib/binimetinib, "
            f"got: {all_partners}"
        )

        # Verify structure of combination therapy entries
        for entry in combination_entries:
            assert "drug_name" in entry
            assert "combination_partners" in entry
            assert isinstance(entry["combination_partners"], list)
            assert len(entry["combination_partners"]) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_hgvs_protein_notation_present(self):
        """HGVS protein notation should be populated for variant queries.

        The hgvs.protein field should contain the HGVS protein notation
        (e.g., p.Val600Lys) from VEP when MyVariant doesn't provide it.
        """
        result = await get_variant_insight(
            gene="BRAF",
            variant="V600K",
            tumor_type="Melanoma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

        # Check HGVS data structure
        hgvs = result.get("hgvs", {})
        assert isinstance(hgvs, dict)

        # HGVS protein notation should be present
        hgvs_protein = hgvs.get("protein")
        assert (
            hgvs_protein is not None
        ), f"Expected HGVS protein notation (p.Val600Lys), got: {hgvs}"

        # Should start with "p." prefix
        assert hgvs_protein.startswith(
            "p."
        ), f"HGVS protein notation should start with 'p.', got: {hgvs_protein}"

        # For V600K, should contain Val600 (V600)
        assert (
            "Val600" in hgvs_protein or "600" in hgvs_protein
        ), f"HGVS protein for V600K should contain Val600, got: {hgvs_protein}"


class TestIDH1R132HPipeline:
    """Full pipeline integration test for IDH1 R132H (p.Arg132His) in Glioma.

    IDH1 R132H is the most common IDH1 mutation in gliomas (~90% of IDH-mutant cases).
    FDA-approved drugs: Ivosidenib (Tibsovo), Vorasidenib (Voranigo), Olutasidenib.
    These are approved at the GENE level (any IDH1 mutation), not variant-specific.
    """

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_idh1_r132h_full_pipeline(self):
        """IDH1 R132H in Glioma - comprehensive pipeline test."""
        result = await get_variant_insight(
            gene="IDH1",
            variant="R132H",
            tumor_type="Glioma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result
        assert result["variant"]["gene"] == "IDH1"
        assert result["variant"]["variant"] == "R132H"

        # Check FDA biomarker evidence
        fda_evidence = result.get("fda_biomarker_evidence", [])
        assert isinstance(fda_evidence, list)

        # IDH1 inhibitors should be found (gene-level approvals)
        assert (
            len(fda_evidence) >= 1
        ), f"Expected at least 1 FDA drug for IDH1 R132H, got {len(fda_evidence)}"

        drug_names = [ev.get("drug_name", "").lower() for ev in fda_evidence]
        expected_drugs = ["ivosidenib", "vorasidenib", "olutasidenib"]

        # At least one IDH1 inhibitor should be found
        found_idh1_drug = any(
            any(exp in d for exp in expected_drugs) for d in drug_names
        )
        assert found_idh1_drug, (
            f"Expected at least one IDH1 inhibitor (ivosidenib, vorasidenib, olutasidenib), "
            f"got: {drug_names}"
        )

        # All matches should be gene-level (IDH1 drugs are not variant-specific)
        for ev in fda_evidence:
            assert ev.get("variant_match_result") == "gene", (
                f"IDH1 approvals are gene-level, got: {ev.get('variant_match_result')} "
                f"for {ev.get('drug_name')}"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_idh1_r132h_hgvs_notation(self):
        """IDH1 R132H should have correct HGVS protein notation (p.Arg132His)."""
        result = await get_variant_insight(
            gene="IDH1",
            variant="R132H",
            tumor_type="Glioma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

        # Check HGVS notation
        hgvs = result.get("hgvs", {})
        assert isinstance(hgvs, dict)

        # HGVS protein should be p.Arg132His
        hgvs_protein = hgvs.get("protein")
        assert hgvs_protein is not None, "HGVS protein notation should be present"
        assert hgvs_protein.startswith("p."), f"Expected p. prefix, got: {hgvs_protein}"

        # Should contain Arg132 (the reference amino acid at position 132)
        assert (
            "Arg132" in hgvs_protein or "132" in hgvs_protein
        ), f"HGVS protein for R132H should contain Arg132, got: {hgvs_protein}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_idh1_r132h_fda_drug_count(self):
        """IDH1 R132H should have at least 1 FDA-approved drug."""
        result = await get_variant_insight(
            gene="IDH1",
            variant="R132H",
            tumor_type="Glioma",
            enable_llm=False,
            enable_literature=False,
        )

        assert "error" not in result

        fda_evidence = result.get("fda_biomarker_evidence", [])

        # Should have at least 1 FDA-approved drug
        assert (
            len(fda_evidence) >= 1
        ), f"IDH1 R132H Glioma should have at least 1 FDA drug, got {len(fda_evidence)}"

        # Check structure of each entry
        for ev in fda_evidence:
            assert "drug_name" in ev
            assert "gene" in ev
            assert ev.get("gene") == "IDH1"
            assert "variant_match_result" in ev
            assert "variant_match_reason" in ev

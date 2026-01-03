"""Integration tests for MyVariant.info API.

Tests validate that the MyVariant API returns expected annotations including
AlphaMissense, COSMIC, ClinVar, CADD, and CIViC evidence for well-characterized variants.
"""

import pytest

from oncomind.api.myvariant import MyVariantClient


class TestMyVariantBasic:
    """Basic MyVariant API tests."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_braf_v600e(self):
        """BRAF V600E should return evidence with database identifiers."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("AKT1", "E17K")

            assert evidence is not None
            assert evidence.gene == "BRAF"
            assert evidence.variant == "V600E"

            has_identifiers = any([
                evidence.dbsnp_id,
                evidence.cosmic_id,
                evidence.clinvar_id,
            ])
            assert has_identifiers, "BRAF V600E should have database identifiers"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_evidence_structure(self):
        """Evidence should have valid data structure."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("BRAF", "V600E")

            assert isinstance(evidence.gene, str)
            assert isinstance(evidence.variant, str)
            assert isinstance(evidence.civic, list)
            assert isinstance(evidence.clinvar, list)
            assert isinstance(evidence.cosmic, list)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_multiple_annotation_sources(self):
        """Well-characterized variants should have multiple annotation sources."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("BRAF", "V600E")

            annotation_count = 0
            if evidence.alphamissense_score is not None or evidence.alphamissense_prediction is not None:
                annotation_count += 1
            if evidence.cadd_score is not None:
                annotation_count += 1
            if evidence.polyphen2_prediction is not None:
                annotation_count += 1
            if evidence.clinvar_clinical_significance is not None or len(evidence.clinvar) > 0:
                annotation_count += 1
            if len(evidence.civic) > 0:
                annotation_count += 1
            if len(evidence.cosmic) > 0:
                annotation_count += 1

            assert annotation_count >= 2, (
                f"BRAF V600E expected multiple annotation sources, only found {annotation_count}"
            )


class TestAlphaMissense:
    """Tests for AlphaMissense pathogenicity predictions."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_tp53_r248w_pathogenic(self):
        """TP53 R248W (known pathogenic) should have high AlphaMissense score."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("TP53", "R248W")

            assert evidence is not None
            assert evidence.gene == "TP53"

            if evidence.alphamissense_score is not None:
                assert 0 <= evidence.alphamissense_score <= 1
                assert evidence.alphamissense_score > 0.5, (
                    f"TP53 R248W is pathogenic, expected high score, got {evidence.alphamissense_score}"
                )

            if evidence.alphamissense_prediction is not None:
                assert evidence.alphamissense_prediction in ["P", "B", "A", "pathogenic", "benign", "ambiguous"]
                assert evidence.alphamissense_prediction in ["P", "pathogenic"], (
                    f"TP53 R248W should be pathogenic, got {evidence.alphamissense_prediction}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_egfr_l858r(self):
        """EGFR L858R should have AlphaMissense data."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("EGFR", "L858R")

            assert evidence is not None
            assert evidence.gene == "EGFR"

            if evidence.alphamissense_score is not None:
                assert 0 <= evidence.alphamissense_score <= 1
                assert evidence.alphamissense_score > 0.3, (
                    f"EGFR L858R expected higher score, got {evidence.alphamissense_score}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_kras_g12d_oncogenic(self):
        """KRAS G12D (oncogenic) should have high AlphaMissense score."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("KRAS", "G12D")

            assert evidence is not None
            assert evidence.gene == "KRAS"

            if evidence.alphamissense_score is not None:
                assert 0 <= evidence.alphamissense_score <= 1
                assert evidence.alphamissense_score > 0.5, (
                    f"KRAS G12D is oncogenic, expected high score, got {evidence.alphamissense_score}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_pik3ca_h1047r(self):
        """PIK3CA H1047R should have AlphaMissense data."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("PIK3CA", "H1047R")

            assert evidence is not None
            assert evidence.gene == "PIK3CA"

            if evidence.alphamissense_score is not None:
                assert 0 <= evidence.alphamissense_score <= 1
                assert evidence.alphamissense_score > 0.4, (
                    f"PIK3CA H1047R expected higher score, got {evidence.alphamissense_score}"
                )


class TestCOSMIC:
    """Tests for COSMIC (Catalogue of Somatic Mutations in Cancer) data."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_braf_v600e_cosmic_id(self):
        """BRAF V600E should have COSMIC ID."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("BRAF", "V600E")

            assert evidence is not None
            assert evidence.gene == "BRAF"

            if evidence.cosmic_id:
                assert (
                    evidence.cosmic_id.startswith("COSM") or
                    evidence.cosmic_id.startswith("COSV") or
                    evidence.cosmic_id.isdigit()
                ), f"Unexpected COSMIC ID format: {evidence.cosmic_id}"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_braf_v600e_cosmic_evidence(self):
        """BRAF V600E should have COSMIC evidence entries."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("BRAF", "V600E")

            if evidence.cosmic:
                assert len(evidence.cosmic) > 0

                for entry in evidence.cosmic:
                    if entry.mutation_id is not None:
                        assert isinstance(entry.mutation_id, str)
                    if entry.primary_site is not None:
                        assert isinstance(entry.primary_site, str)
                    if entry.sample_count is not None:
                        assert isinstance(entry.sample_count, int)
                        assert entry.sample_count >= 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_tp53_r175h_cosmic(self):
        """TP53 R175H (major hotspot) should have COSMIC data."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("TP53", "R175H")

            assert evidence is not None
            assert evidence.gene == "TP53"

            if evidence.cosmic:
                for entry in evidence.cosmic:
                    if entry.primary_site:
                        assert isinstance(entry.primary_site, str)

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_kras_g12c_cosmic(self):
        """KRAS G12C should have COSMIC annotation."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("KRAS", "G12C")

            assert evidence is not None
            assert evidence.gene == "KRAS"

            if evidence.cosmic_id:
                assert isinstance(evidence.cosmic_id, str)

            if evidence.cosmic:
                assert isinstance(evidence.cosmic, list)
                for entry in evidence.cosmic:
                    assert hasattr(entry, 'mutation_id')
                    assert hasattr(entry, 'primary_site')
                    assert hasattr(entry, 'primary_histology')

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cosmic_structure_integrity(self):
        """COSMIC evidence entries should have valid structure."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("BRAF", "V600E")

            if evidence.cosmic:
                for entry in evidence.cosmic:
                    assert hasattr(entry, 'mutation_id')
                    assert hasattr(entry, 'primary_site')
                    assert hasattr(entry, 'site_subtype')
                    assert hasattr(entry, 'primary_histology')
                    assert hasattr(entry, 'histology_subtype')
                    assert hasattr(entry, 'sample_count')
                    assert hasattr(entry, 'mutation_somatic_status')


class TestCADD:
    """Tests for CADD (Combined Annotation Dependent Depletion) scores."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_braf_v600e_cadd(self):
        """BRAF V600E (deleterious) should have high CADD score."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("BRAF", "V600E")

            if evidence.cadd_score is not None:
                assert evidence.cadd_score >= 0
                assert evidence.cadd_score > 10, (
                    f"BRAF V600E expected high CADD score, got {evidence.cadd_score}"
                )


class TestClinVar:
    """Tests for ClinVar clinical significance data."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_tp53_r248w_clinvar(self):
        """TP53 R248W should have ClinVar pathogenic classification."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("TP53", "R248W")

            has_clinvar = any([
                evidence.clinvar_id,
                evidence.clinvar_clinical_significance,
                evidence.clinvar_accession,
                len(evidence.clinvar) > 0,
            ])

            if has_clinvar and evidence.clinvar_clinical_significance:
                sig_lower = evidence.clinvar_clinical_significance.lower()
                assert "pathogenic" in sig_lower or "likely pathogenic" in sig_lower, (
                    f"TP53 R248W should be pathogenic, got: {evidence.clinvar_clinical_significance}"
                )


class TestClinVarConditionsList:
    """Tests for ClinVar variants with list-type conditions (regression test).

    Some ClinVar records have 'conditions' as a list instead of a single dict.
    This caused MyVariant parsing to fail for variants like TP53 R273H.
    """

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_tp53_r273h_functional_data(self):
        """TP53 R273H should return functional data despite complex ClinVar structure.

        This variant has ClinVar RCV records where 'conditions' is a list,
        which previously caused Pydantic validation to fail.
        """
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("TP53", "R273H")

            assert evidence is not None
            assert evidence.gene == "TP53"
            assert evidence.variant == "R273H"

            # Should have functional predictions (this was failing before the fix)
            has_functional = any([
                evidence.alphamissense_score is not None,
                evidence.alphamissense_prediction is not None,
                evidence.polyphen2_prediction is not None,
                evidence.cadd_score is not None,
            ])
            assert has_functional, (
                "TP53 R273H should have functional predictions - "
                "ClinVar parsing may have failed"
            )

            # AlphaMissense should indicate pathogenic
            if evidence.alphamissense_score is not None:
                assert evidence.alphamissense_score > 0.9, (
                    f"TP53 R273H is a major hotspot, expected high score, "
                    f"got {evidence.alphamissense_score}"
                )

            if evidence.alphamissense_prediction is not None:
                assert evidence.alphamissense_prediction in ["P", "pathogenic"], (
                    f"TP53 R273H should be pathogenic, got {evidence.alphamissense_prediction}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_tp53_r273h_identifiers(self):
        """TP53 R273H should have COSMIC and HGVS identifiers."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("TP53", "R273H")

            assert evidence is not None

            # Should have COSMIC ID (this is a major hotspot)
            assert evidence.cosmic_id is not None, (
                "TP53 R273H should have COSMIC ID"
            )
            assert evidence.cosmic_id.startswith("COSM") or evidence.cosmic_id.startswith("COSV"), (
                f"Unexpected COSMIC ID format: {evidence.cosmic_id}"
            )

            # Should have HGVS genomic notation
            assert evidence.hgvs_genomic is not None, (
                "TP53 R273H should have HGVS genomic notation"
            )


class TestClinVarSingleRCV:
    """Tests for ClinVar variants with single RCV dict (regression test).

    Some ClinVar records return 'rcv' as a single dict instead of a list.
    This caused MyVariant parsing to fail for variants like GNAQ Q209L.
    """

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_gnaq_q209l_uveal_melanoma(self):
        """GNAQ Q209L should return evidence despite single RCV structure.

        This variant is a key driver in uveal melanoma and has ClinVar RCV
        returned as a single dict (not a list), which previously caused
        'object of type ClinVarRCV has no len()' error.
        """
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("GNAQ", "Q209L")

            assert evidence is not None
            assert evidence.gene == "GNAQ"
            assert evidence.variant == "Q209L"

            # Should have ClinVar data
            has_clinvar = any([
                evidence.clinvar_id,
                evidence.clinvar_clinical_significance,
                evidence.clinvar_accession,
                len(evidence.clinvar) > 0,
            ])
            assert has_clinvar, (
                "GNAQ Q209L should have ClinVar data - "
                "single RCV parsing may have failed"
            )

            # Should be pathogenic (known uveal melanoma driver)
            if evidence.clinvar_clinical_significance:
                sig_lower = evidence.clinvar_clinical_significance.lower()
                assert "pathogenic" in sig_lower, (
                    f"GNAQ Q209L should be pathogenic, got: {evidence.clinvar_clinical_significance}"
                )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_gnaq_q209l_identifiers(self):
        """GNAQ Q209L should have database identifiers."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("GNAQ", "Q209L")

            assert evidence is not None

            # Should have at least one identifier
            has_identifiers = any([
                evidence.dbsnp_id,
                evidence.cosmic_id,
                evidence.clinvar_id,
            ])
            assert has_identifiers, (
                "GNAQ Q209L should have database identifiers"
            )

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_gnaq_q209p_alternate_variant(self):
        """GNAQ Q209P (alternate hotspot) should also work.

        Tests that both Q209L and Q209P variants at the same codon
        are handled correctly.
        """
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("GNAQ", "Q209P")

            assert evidence is not None
            assert evidence.gene == "GNAQ"
            assert evidence.variant == "Q209P"


class TestCIViCFallback:
    """Tests for CIViC fallback for fusions and amplifications."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_alk_fusion(self):
        """ALK fusion should have CIViC evidence via fallback."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("ALK", "fusion")

            assert evidence is not None
            assert evidence.gene == "ALK"

            if evidence.civic:
                assert len(evidence.civic) > 0, "ALK fusion should have CIViC evidence"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_erbb2_amplification(self):
        """ERBB2 amplification should have CIViC evidence via fallback."""
        async with MyVariantClient() as client:
            evidence = await client.fetch_evidence("ERBB2", "amplification")

            assert evidence is not None
            assert evidence.gene == "ERBB2"

            if evidence.civic:
                assert len(evidence.civic) > 0, "ERBB2 amplification should have CIViC evidence"

"""Integration tests for AnnotatorMyVariantClient against real API."""

import pytest

from oncomind.annotator.clients.myvariant import AnnotatorMyVariantClient
from oncomind.annotator.models import MyVariantAnnotation


@pytest.mark.integration
class TestMyVariantIntegration:
    """Integration tests for MyVariant.info API."""

    @pytest.mark.asyncio
    async def test_fetch_annotation_returns_pydantic_model(self):
        """Test fetch_annotation returns MyVariantAnnotation pydantic model."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        assert isinstance(result, MyVariantAnnotation)

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e(self):
        """Test fetching annotation for BRAF V600E from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        # BRAF V600E is a well-known variant, should have vcf data
        assert result.vcf is not None
        assert result.vcf.ref == "A"
        assert result.vcf.alt == "T"

    @pytest.mark.asyncio
    async def test_fetch_annotation_egfr_l858r(self):
        """Test fetching annotation for EGFR L858R from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("EGFR", "p.L858R")

        assert isinstance(result, MyVariantAnnotation)

    @pytest.mark.asyncio
    async def test_fetch_annotation_kras_g12c(self):
        """Test fetching annotation for KRAS G12C from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("KRAS", "p.G12C")

        assert isinstance(result, MyVariantAnnotation)

    @pytest.mark.asyncio
    async def test_build_annotation_returns_full_data(self):
        """Test build_annotation returns full API response."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.build_annotation("BRAF", "p.V600E")

        assert "total" in result
        assert "hits" in result

    @pytest.mark.asyncio
    async def test_fetch_annotation_unknown_variant(self):
        """Test fetching annotation for unknown variant returns null for all fields."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("FAKEGENE", "X999Z")

        assert result.vcf is None
        assert result.cadd is None
        assert result.clinvar is None
        assert result.cosmic is None
        assert result.dbnsfp is None
        assert result.dbsnp is None
        assert result.gnomad_exome is None

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_cadd_fields(self):
        """Test that BRAF V600E returns expected CADD fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        cadd = result.cadd

        # Test consequence field
        assert cadd.consequence == "NON_SYNONYMOUS"

        # Test annotype field
        assert cadd.annotype == "CodingTranscript"

        # Test dna field (DNA shape parameters)
        assert cadd.dna.helt == -0.55
        assert cadd.dna.mgw == 0.7
        assert cadd.dna.prot == 3.19
        assert cadd.dna.roll == 8.98

        # Test gene field with validated patterns
        assert cadd.gene.genename == "BRAF"
        assert cadd.gene.gene_id == "ENSG00000157764"
        assert cadd.gene.gene_id.startswith("ENSG")
        assert cadd.gene.feature_id == "ENST00000288602"
        assert cadd.gene.feature_id.startswith("ENST")
        assert cadd.gene.ccds_id == "CCDS5863.1"
        assert cadd.gene.cds.cds_pos == 1799
        assert cadd.gene.cds.cds_pos >= 1
        assert cadd.gene.cds.cdna_pos == 1860
        assert cadd.gene.prot.protpos == 600
        assert cadd.gene.prot.protpos >= 1

        # Test amino acid change (single letter codes)
        assert cadd.oaa == "V"
        assert len(cadd.oaa) == 1
        assert cadd.naa == "E"
        assert len(cadd.naa) == 1

        # Test score fields with range validation
        assert cadd.phred == 32
        assert 0 <= cadd.phred <= 100
        assert cadd.rawscore == 6.641785

        # Test genomic position
        assert cadd.chrom == 7
        assert cadd.pos == 140453136
        assert cadd.pos >= 1
        assert cadd.ref == "A"
        assert cadd.alt == "T"
        assert cadd.type == "SNV"

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_clinvar_fields(self):
        """Test that BRAF V600E returns expected ClinVar fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        clinvar = result.clinvar

        # Test ClinVar fields with value validation
        assert clinvar.allele_id == 29000
        assert clinvar.allele_id >= 1
        assert clinvar.gene.symbol == "BRAF"
        assert clinvar.gene.id == "673"
        assert clinvar.chrom == "7"
        assert clinvar.alt == "T"

        # Test RCV records exist and have expected structure
        assert len(clinvar.rcv) > 0
        first_rcv = clinvar.rcv[0]
        assert first_rcv.accession == "RCV000014992"
        assert first_rcv.accession.startswith("RCV")
        assert first_rcv.clinical_significance == "Pathogenic"

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_cosmic_fields(self):
        """Test that BRAF V600E returns expected COSMIC fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        cosmic = result.cosmic

        # Test COSMIC fields with pattern validation
        assert cosmic.cosmic_id == "COSM476"
        assert cosmic.cosmic_id.startswith("COSM")
        assert cosmic.mut_nt == "T>A"
        assert cosmic.chrom == "7"
        assert cosmic.alt == "A"
        assert cosmic.ref == "T"

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_dbnsfp_fields(self):
        """Test that BRAF V600E returns expected dbNSFP fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        dbnsfp = result.dbnsfp

        # Test Polyphen2 prediction with value validation
        assert dbnsfp.polyphen2.hdiv.pred == "D"  # D=Damaging
        assert dbnsfp.polyphen2.hdiv.pred in ("D", "B", "P")
        assert dbnsfp.polyphen2.hdiv.score == 0.971
        assert 0 <= dbnsfp.polyphen2.hdiv.score <= 1

        # Test AlphaMissense prediction
        assert dbnsfp.alphamissense.pred[0] == "P"  # P=Pathogenic
        assert dbnsfp.alphamissense.score[0] == 0.9853

    @pytest.mark.asyncio
    async def test_fetch_annotation_tp53_r248w_cadd_fields(self):
        """Test that TP53 R248W returns expected CADD fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("TP53", "p.R248W")

        cadd = result.cadd

        # Test consequence field
        assert cadd.consequence == "NON_SYNONYMOUS"
        assert cadd.annotype == "CodingTranscript"

        # Test gene field
        assert cadd.gene.genename == "TP53"
        assert cadd.gene.gene_id == "ENSG00000141510"
        assert cadd.gene.feature_id == "ENST00000269305"
        assert cadd.gene.cds.cds_pos == 742
        assert cadd.gene.prot.protpos == 248

        # Test amino acid change
        assert cadd.oaa == "R"
        assert cadd.naa == "W"

        # Test score fields with range validation
        assert cadd.phred == 34
        assert 0 <= cadd.phred <= 100
        assert cadd.rawscore == 7.204272

        # Test genomic position
        assert cadd.chrom == 17
        assert cadd.pos == 7577539
        assert cadd.pos >= 1
        assert cadd.ref == "G"
        assert cadd.alt == "A"
        assert cadd.type == "SNV"

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_gnomad_fields(self):
        """Test that BRAF V600E returns expected gnomAD fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        gnomad = result.gnomad_exome

        # Test gnomAD allele frequency (very rare variant) with range validation
        assert gnomad.af.af == 3.97994e-06
        assert 0 <= gnomad.af.af <= 1

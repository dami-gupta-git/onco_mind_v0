"""Integration tests for AnnotatorMyVariantClient against real API."""

import pytest

from oncomind.annotator.clients.myvariant import AnnotatorMyVariantClient
from oncomind.annotator.models import MyVariantAnnotation


@pytest.mark.integration
class TestMyVariantIntegration:
    """Integration tests for MyVariant.info API."""

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_returns_valid_annotation(self):
        """Test fetch_annotation returns valid MyVariantAnnotation for BRAF V600E."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        assert isinstance(result, MyVariantAnnotation)
        # BRAF V600E is well-characterized, should have VCF data
        assert result.vcf.ref == "A"
        assert result.vcf.alt == "T"
        # Should have CADD data for this known pathogenic variant
        assert result.cadd.gene.genename == "BRAF"
        assert result.cadd.consequence == "NON_SYNONYMOUS"

    @pytest.mark.asyncio
    async def test_build_annotation_returns_api_response_structure(self):
        """Test build_annotation returns full API response with hits."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.build_annotation("BRAF", "p.V600E")

        assert "total" in result
        assert result["total"] >= 1
        assert "hits" in result
        assert len(result["hits"]) >= 1
        # First hit should have variant ID
        assert "_id" in result["hits"][0]

    @pytest.mark.asyncio
    async def test_fetch_annotation_unknown_variant_returns_empty(self):
        """Test fetching annotation for unknown variant returns empty annotation."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("FAKEGENE123", "p.X999Z")

        assert isinstance(result, MyVariantAnnotation)
        # All fields should be None for non-existent variant
        assert result.vcf is None
        assert result.cadd is None
        assert result.clinvar is None
        assert result.dbnsfp is None
        assert result.dbsnp is None
        assert result.gnomad_exome is None

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_cadd_fields(self):
        """Test that BRAF V600E returns expected CADD structure from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        cadd = result.cadd

        # Test consequence - V600E is a missense mutation
        assert cadd.consequence == "NON_SYNONYMOUS"
        assert cadd.annotype == "CodingTranscript"

        # Test DNA shape parameters are present and numeric
        assert isinstance(cadd.dna.helt, float)
        assert isinstance(cadd.dna.mgw, float)
        assert isinstance(cadd.dna.prot, float)
        assert isinstance(cadd.dna.roll, float)

        # Test gene field - stable identifiers
        assert cadd.gene.genename == "BRAF"
        assert cadd.gene.gene_id.startswith("ENSG")
        assert cadd.gene.feature_id.startswith("ENST")
        assert cadd.gene.cds.cds_pos >= 1
        assert cadd.gene.prot.protpos == 600  # Position 600 is stable

        # Test amino acid change - V600E means V->E
        assert cadd.oaa == "V"
        assert cadd.naa == "E"

        # Test score fields - BRAF V600E is highly pathogenic
        assert 25 <= cadd.phred <= 40  # Expected range for known pathogenic
        assert cadd.rawscore > 5.0

        # Test genomic position
        assert cadd.chrom == 7
        assert cadd.pos >= 1
        assert cadd.ref == "A"
        assert cadd.alt == "T"
        assert cadd.type == "SNV"

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_clinvar_fields(self):
        """Test that BRAF V600E returns expected ClinVar structure from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        clinvar = result.clinvar

        # Test ClinVar has valid allele ID
        assert clinvar.allele_id >= 1
        assert clinvar.gene.symbol == "BRAF"
        assert clinvar.gene.id == "673"  # NCBI gene ID is stable
        assert clinvar.chrom == "7"
        assert clinvar.alt == "T"

        # Test RCV records exist with valid structure
        assert len(clinvar.rcv) > 0
        first_rcv = clinvar.rcv[0]
        assert first_rcv.accession.startswith("RCV")
        # BRAF V600E is pathogenic in cancer
        assert first_rcv.clinical_significance in (
            "Pathogenic",
            "Likely pathogenic",
            "Pathogenic/Likely pathogenic",
        )

    @pytest.mark.asyncio
    async def test_fetch_annotation_braf_v600e_dbnsfp_fields(self):
        """Test that BRAF V600E returns expected dbNSFP structure from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("BRAF", "p.V600E")

        dbnsfp = result.dbnsfp

        # Test Polyphen2 prediction - V600E should be predicted damaging
        assert dbnsfp.polyphen2.hdiv.pred in (
            "D",
            "P",
        )  # D=Damaging, P=Possibly damaging
        assert 0 <= dbnsfp.polyphen2.hdiv.score <= 1
        assert dbnsfp.polyphen2.hdiv.score > 0.8  # Should be high for pathogenic

        # Test AlphaMissense prediction - should predict pathogenic
        assert dbnsfp.alphamissense.pred[0] in ("P", "A")  # P=Pathogenic, A=Ambiguous
        assert 0 <= dbnsfp.alphamissense.score[0] <= 1
        assert (
            dbnsfp.alphamissense.score[0] > 0.8
        )  # Should be high for known pathogenic

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

        # Test score fields
        assert cadd.phred == 34.0
        assert cadd.rawscore == 7.204272

        # Test genomic position
        assert cadd.chrom == 17
        assert cadd.pos == 7577539
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

    @pytest.mark.asyncio
    async def test_fetch_annotation_pik3ca_h1047r_cadd_fields(self):
        """Test that PIK3CA H1047R returns expected CADD fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("PIK3CA", "p.H1047R")

        cadd = result.cadd

        # Test consequence - PIK3CA H1047R affects multiple transcripts
        assert cadd.consequence == ["NON_SYNONYMOUS", "DOWNSTREAM"]
        assert cadd.annotype == ["CodingTranscript", "Intergenic"]

        # Test gene field - first entry is the primary transcript
        gene = cadd.gene[0] if isinstance(cadd.gene, list) else cadd.gene
        assert gene.genename == "PIK3CA"
        assert gene.gene_id == "ENSG00000121879"
        assert gene.feature_id == "ENST00000263967"
        assert gene.cds.cds_pos == 3140
        assert gene.prot.protpos == 1047

        # Test amino acid change - H1047R means H->R
        assert cadd.oaa == "H"
        assert cadd.naa == "R"

        # Test score fields
        assert cadd.phred == 20.1
        assert cadd.rawscore == 2.588728

        # Test genomic position
        assert cadd.chrom == 3
        assert cadd.pos == 178952085
        assert cadd.ref == "A"
        assert cadd.alt == "G"
        assert cadd.type == "SNV"

    @pytest.mark.asyncio
    async def test_fetch_annotation_pik3ca_h1047r_clinvar_fields(self):
        """Test that PIK3CA H1047R returns expected ClinVar fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("PIK3CA", "p.H1047R")

        clinvar = result.clinvar

        # Test ClinVar structure
        assert clinvar.allele_id == 28691
        assert clinvar.gene.symbol == "PIK3CA"
        assert clinvar.gene.id == "5290"
        assert clinvar.chrom == "3"
        assert clinvar.alt == "G"

        # Test RCV record
        assert clinvar.rcv[0].accession == "RCV000014622"
        assert clinvar.rcv[0].clinical_significance == "Pathogenic"

    @pytest.mark.asyncio
    async def test_fetch_annotation_pik3ca_h1047r_gnomad_fields(self):
        """Test that PIK3CA H1047R returns expected gnomAD fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("PIK3CA", "p.H1047R")

        gnomad = result.gnomad_exome

        # Test gnomAD allele frequency
        assert gnomad.af.af == 4.02781e-06
        assert 0 <= gnomad.af.af <= 1

    @pytest.mark.asyncio
    async def test_fetch_annotation_egfr_l858r_cadd_fields(self):
        """Test that EGFR L858R returns expected CADD fields from real API.

        Note: EGFR L858R requires three-letter amino acid code fallback query.
        """
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("EGFR", "p.L858R")

        cadd = result.cadd

        # Test consequence field
        assert cadd.consequence == "NON_SYNONYMOUS"
        assert cadd.annotype == "CodingTranscript"

        # Test gene field
        assert cadd.gene.genename == "EGFR"
        assert cadd.gene.gene_id == "ENSG00000146648"
        assert cadd.gene.feature_id == "ENST00000275493"
        assert cadd.gene.cds.cds_pos == 2573
        assert cadd.gene.prot.protpos == 858

        # Test amino acid change
        assert cadd.oaa == "L"
        assert cadd.naa == "R"

        # Test score fields
        assert cadd.phred == 33.0
        assert cadd.rawscore == 6.929752

        # Test genomic position
        assert cadd.chrom == 7
        assert cadd.pos == 55259515
        assert cadd.ref == "T"
        assert cadd.alt == "G"
        assert cadd.type == "SNV"

    @pytest.mark.asyncio
    async def test_fetch_annotation_egfr_l858r_clinvar_fields(self):
        """Test that EGFR L858R returns expected ClinVar fields from real API."""
        async with AnnotatorMyVariantClient() as client:
            result = await client.fetch_annotation("EGFR", "p.L858R")

        clinvar = result.clinvar

        # Test ClinVar structure
        assert clinvar.allele_id == 31648
        assert clinvar.gene.symbol == "EGFR"
        assert clinvar.gene.id == "1956"
        assert clinvar.chrom == "7"
        assert clinvar.alt == "G"

        # Test RCV record - EGFR L858R has drug response significance
        rcv = clinvar.rcv[0] if isinstance(clinvar.rcv, list) else clinvar.rcv
        assert rcv.accession == "RCV000018084"
        assert rcv.clinical_significance == "drug response"


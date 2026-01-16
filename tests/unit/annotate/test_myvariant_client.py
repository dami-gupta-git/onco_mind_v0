"""Unit tests for AnnotatorMyVariantClient."""

import pytest
from unittest.mock import AsyncMock, patch, MagicMock

from oncomind.annotator.clients.myvariant import AnnotatorMyVariantClient
from oncomind.annotator.models import MyVariantAnnotation


class TestAnnotatorMyVariantClient:
    """Tests for AnnotatorMyVariantClient class."""

    def test_client_initialization(self):
        """Test client initializes with default timeout."""
        client = AnnotatorMyVariantClient()

        assert client.timeout == 30.0
        assert client._client is None

    def test_client_custom_timeout(self):
        """Test client accepts custom timeout."""
        client = AnnotatorMyVariantClient(timeout=60.0)

        assert client.timeout == 60.0

    @pytest.mark.asyncio
    async def test_context_manager_creates_client(self):
        """Test async context manager creates HTTP client."""
        client = AnnotatorMyVariantClient()

        async with client:
            assert client._client is not None

        assert client._client is None

    @pytest.mark.asyncio
    async def test_fetch_annotation_returns_pydantic_model(self):
        """Test fetch_annotation returns MyVariantAnnotation pydantic model."""
        client = AnnotatorMyVariantClient()

        mock_response = {
            "total": 1,
            "hits": [
                {
                    "_id": "chr7:g.140453136A>T",
                    "vcf": {"ref": "A", "alt": "T"},
                }
            ],
        }

        with patch.object(client, "build_annotation", new_callable=AsyncMock) as mock_build:
            mock_build.return_value = mock_response
            result = await client.fetch_annotation("BRAF", "p.V600E")

        assert isinstance(result, MyVariantAnnotation)

    @pytest.mark.asyncio
    async def test_fetch_annotation_returns_all_fields(self):
        """Test fetch_annotation returns vcf, cadd, clinvar, cosmic, dbnsfp, dbsnp, gnomad_exome."""
        client = AnnotatorMyVariantClient()

        mock_cadd = {
            "phred": 32.0,
            "rawscore": 6.5,
            "consequence": "NON_SYNONYMOUS",
        }
        mock_clinvar = {
            "allele_id": 29000,
            "gene": {"symbol": "BRAF", "id": "673"},
            "rcv": [{"accession": "RCV000012345", "clinical_significance": "Pathogenic"}],
        }
        mock_cosmic = {"cosmic_id": "COSM476"}
        mock_dbnsfp = {"sift": {"pred": "D", "score": 0.0}}
        mock_dbsnp = {"rsid": "rs113488022"}
        mock_gnomad = {"af": {"af": 0.00001}}

        mock_response = {
            "total": 1,
            "hits": [
                {
                    "_id": "chr7:g.140453136A>T",
                    "vcf": {"ref": "A", "alt": "T"},
                    "cadd": mock_cadd,
                    "clinvar": mock_clinvar,
                    "cosmic": mock_cosmic,
                    "dbnsfp": mock_dbnsfp,
                    "dbsnp": mock_dbsnp,
                    "gnomad_exome": mock_gnomad,
                }
            ],
        }

        with patch.object(client, "build_annotation", new_callable=AsyncMock) as mock_build:
            mock_build.return_value = mock_response
            result = await client.fetch_annotation("BRAF", "p.V600E")

        # Test VCF fields
        assert result.vcf.ref == "A"
        assert result.vcf.alt == "T"

        # Test CADD fields with value validation
        assert result.cadd.phred == 32.0
        assert 0 <= result.cadd.phred <= 100  # phred score range
        assert result.cadd.rawscore == 6.5
        assert result.cadd.consequence == "NON_SYNONYMOUS"

        # Test ClinVar fields
        assert result.clinvar.allele_id == 29000
        assert result.clinvar.allele_id >= 1  # must be positive
        assert result.clinvar.gene.symbol == "BRAF"
        assert result.clinvar.rcv[0].accession == "RCV000012345"
        assert result.clinvar.rcv[0].clinical_significance == "Pathogenic"

        # Test COSMIC fields
        assert result.cosmic.cosmic_id == "COSM476"
        assert result.cosmic.cosmic_id.startswith("COSM")  # COSMIC ID format

        # Test dbSNP fields
        assert result.dbsnp.rsid == "rs113488022"
        assert result.dbsnp.rsid.startswith("rs")  # rsID format

        # Test gnomAD fields
        assert result.gnomad_exome.af.af == 0.00001
        assert 0 <= result.gnomad_exome.af.af <= 1  # AF must be between 0 and 1

    @pytest.mark.asyncio
    async def test_fetch_annotation_no_hits(self):
        """Test fetch_annotation returns None for all fields when no hits."""
        client = AnnotatorMyVariantClient()

        mock_response = {"total": 0, "hits": []}

        with patch.object(client, "build_annotation", new_callable=AsyncMock) as mock_build:
            mock_build.return_value = mock_response
            result = await client.fetch_annotation("UNKNOWN", "X123Y")

        assert isinstance(result, MyVariantAnnotation)
        assert result.vcf is None
        assert result.cadd is None
        assert result.clinvar is None
        assert result.cosmic is None
        assert result.dbnsfp is None
        assert result.dbsnp is None
        assert result.gnomad_exome is None

    @pytest.mark.asyncio
    async def test_fetch_annotation_missing_fields_in_hit(self):
        """Test fetch_annotation handles missing fields in hit - returns None for missing."""
        client = AnnotatorMyVariantClient()

        mock_response = {
            "total": 1,
            "hits": [
                {
                    "_id": "chr7:g.140453136A>T",
                    "clinvar": {
                        "allele_id": 12345,
                        "gene": {"symbol": "BRAF", "id": "673"},
                    },
                }
            ],
        }

        with patch.object(client, "build_annotation", new_callable=AsyncMock) as mock_build:
            mock_build.return_value = mock_response
            result = await client.fetch_annotation("BRAF", "p.V600E")

        assert result.vcf is None
        assert result.cadd is None
        assert result.clinvar.allele_id == 12345
        assert result.cosmic is None
        assert result.dbnsfp is None
        assert result.dbsnp is None
        assert result.gnomad_exome is None

    @pytest.mark.asyncio
    async def test_build_annotation_tries_multiple_queries(self):
        """Test build_annotation tries multiple query strategies."""
        client = AnnotatorMyVariantClient()

        # First query returns no hits, second returns hits
        mock_responses = [
            {"total": 0, "hits": []},
            {"total": 1, "hits": [{"_id": "test", "vcf": {"ref": "A", "alt": "T"}}]},
        ]

        with patch.object(client, "_query", new_callable=AsyncMock) as mock_query:
            mock_query.side_effect = mock_responses
            result = await client.build_annotation("BRAF", "p.V600E")

        assert mock_query.call_count == 2
        assert result["total"] == 1

    @pytest.mark.asyncio
    async def test_fetch_annotation_cadd_contains_expected_fields(self):
        """Test fetch_annotation returns cadd with expected structure and validated values."""
        client = AnnotatorMyVariantClient()

        mock_cadd = {
            "phred": 32,
            "rawscore": 6.641785,
            "consequence": "NON_SYNONYMOUS",
            "annotype": "CodingTranscript",
            "dna": {
                "helt": -0.55,
                "mgw": 0.7,
                "prot": 3.19,
                "roll": 8.98,
            },
            "gene": {
                "genename": "BRAF",
                "gene_id": "ENSG00000157764",
                "feature_id": "ENST00000288602",
                "ccds_id": "CCDS5863.1",
                "cds": {
                    "cdna_pos": 1860,
                    "cds_pos": 1799,
                    "rel_cdna_pos": 0.75,
                    "rel_cds_pos": 0.78,
                },
                "prot": {
                    "protpos": 600,
                    "rel_prot_pos": 0.78,
                },
            },
            "oaa": "V",
            "naa": "E",
            "type": "SNV",
            "chrom": 7,
            "pos": 140453136,
            "ref": "A",
            "alt": "T",
        }

        mock_response = {
            "total": 1,
            "hits": [
                {
                    "_id": "chr7:g.140453136A>T",
                    "vcf": {"ref": "A", "alt": "T"},
                    "cadd": mock_cadd,
                }
            ],
        }

        with patch.object(client, "build_annotation", new_callable=AsyncMock) as mock_build:
            mock_build.return_value = mock_response
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
        assert cadd.gene.gene_id.startswith("ENSG")  # Ensembl gene ID format
        assert cadd.gene.feature_id == "ENST00000288602"
        assert cadd.gene.feature_id.startswith("ENST")  # Ensembl transcript format
        assert cadd.gene.ccds_id == "CCDS5863.1"
        assert cadd.gene.cds.cds_pos == 1799
        assert cadd.gene.cds.cds_pos >= 1  # position must be positive
        assert cadd.gene.prot.protpos == 600
        assert cadd.gene.prot.protpos >= 1  # protein position must be positive

        # Test amino acid change (single letter codes)
        assert cadd.oaa == "V"
        assert len(cadd.oaa) == 1
        assert cadd.naa == "E"
        assert len(cadd.naa) == 1

        # Test score fields with range validation
        assert cadd.phred == 32
        assert 0 <= cadd.phred <= 100  # phred score range
        assert cadd.rawscore == 6.641785

        # Test genomic position
        assert cadd.chrom == 7
        assert cadd.pos == 140453136
        assert cadd.pos >= 1  # genomic position must be positive
        assert cadd.ref == "A"
        assert cadd.alt == "T"

    @pytest.mark.asyncio
    async def test_close_method(self):
        """Test close method closes client."""
        client = AnnotatorMyVariantClient()

        async with client:
            assert client._client is not None

        await client.close()
        assert client._client is None

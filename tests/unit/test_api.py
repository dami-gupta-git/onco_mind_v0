"""Tests for API client."""

import pytest
from unittest.mock import AsyncMock, MagicMock, patch

from oncomind.api.myvariant import MyVariantAPIError, MyVariantClient
from oncomind.models.evidence import ClinVarEvidence


class TestMyVariantClient:
    """Tests for MyVariantClient."""

    @pytest.mark.asyncio
    async def test_context_manager(self):
        """Test async context manager."""
        async with MyVariantClient() as client:
            assert client._client is not None

        # Client should be closed after exit
        assert client._client is None

    @pytest.mark.asyncio
    async def test_fetch_evidence_no_results(self):
        """Test fetching evidence with no results."""
        client = MyVariantClient()

        with patch.object(client, "_query", new_callable=AsyncMock) as mock_query:
            mock_query.return_value = {"hits": []}

            evidence = await client.fetch_evidence("UNKNOWN", "X123Y")

            assert evidence.gene == "UNKNOWN"
            assert evidence.variant == "X123Y"
            assert not evidence.has_evidence()

        await client.close()

    @pytest.mark.asyncio
    async def test_fetch_evidence_with_clinvar(self):
        """Test fetching evidence with ClinVar data.

        Note: CIViC evidence is now fetched separately via CIViCClient.
        MyVariant is used for ClinVar, COSMIC, AlphaMissense, CADD, gnomAD.
        """
        client = MyVariantClient()

        mock_response = {
            "hits": [
                {
                    "_id": "test123",
                    "clinvar": {
                        "variant_id": 13961,
                        "rcv": [
                            {
                                "clinical_significance": "Pathogenic",
                                "review_status": "reviewed by expert panel",
                                "accession": "RCV000013961",
                            }
                        ],
                    },
                }
            ]
        }

        with patch.object(client, "_query", new_callable=AsyncMock) as mock_query:
            mock_query.return_value = mock_response

            evidence = await client.fetch_evidence("BRAF", "V600E")

            assert evidence.has_evidence()
            assert evidence.clinvar_id == "13961"
            assert evidence.clinvar_clinical_significance == "Pathogenic"
            # CIViC should be empty (fetched separately via CIViCClient)
            assert evidence.civic == []

        await client.close()

    @pytest.mark.asyncio
    async def test_civic_not_parsed_from_myvariant(self):
        """Test that CIViC data is not parsed from MyVariant (fetched separately).

        CIViC evidence is now fetched directly from the CIViC GraphQL API
        via CIViCClient, not from MyVariant.info.
        """
        client = MyVariantClient()

        # Even if MyVariant returns CIViC data, we ignore it
        mock_response = {
            "hits": [
                {
                    "_id": "test123",
                    "civic": {
                        "evidence_items": [
                            {
                                "evidence_type": "Predictive",
                                "evidence_level": "A",
                            }
                        ]
                    },
                }
            ]
        }

        with patch.object(client, "_query", new_callable=AsyncMock) as mock_query:
            mock_query.return_value = mock_response

            evidence = await client.fetch_evidence("BRAF", "V600E")

            # CIViC should be empty even if MyVariant returned data
            assert evidence.civic == []

        await client.close()

    @pytest.mark.asyncio
    async def test_parse_clinvar_evidence(self):
        """Test parsing ClinVar evidence."""
        client = MyVariantClient()

        clinvar_data = {
            "clinical_significance": "Pathogenic",
            "review_status": "reviewed by expert panel",
            "conditions": [{"name": "Cancer"}],
            "variation_id": "12345",
        }

        parsed = client._parse_clinvar_evidence(clinvar_data)

        assert len(parsed) == 1
        assert "Pathogenic" in parsed[0].clinical_significance
        assert "Cancer" in parsed[0].conditions

    @pytest.mark.asyncio
    async def test_api_error_handling(self):
        """Test API error handling."""
        client = MyVariantClient()

        with patch.object(client, "_query", new_callable=AsyncMock) as mock_query:
            mock_query.side_effect = MyVariantAPIError("API error")

            with pytest.raises(MyVariantAPIError):
                await client.fetch_evidence("BRAF", "V600E")

        await client.close()

    @pytest.mark.asyncio
    async def test_fetch_evidence_with_identifiers(self):
        """Test fetching evidence with database identifiers."""
        client = MyVariantClient()

        mock_response = {
            "hits": [
                {
                    "_id": "chr7:g.140453136A>T",
                    "cosmic": {"cosmic_id": "COSM476"},
                    "dbsnp": {
                        "rsid": "rs113488022",
                        "gene": {"geneid": 673},
                    },
                    "clinvar": {"variant_id": "13961"},
                    "civic": {},
                }
            ]
        }

        with patch.object(client, "_query", new_callable=AsyncMock) as mock_query:
            mock_query.return_value = mock_response

            evidence = await client.fetch_evidence("BRAF", "V600E")

            # Verify identifiers were extracted
            assert evidence.cosmic_id == "COSM476"
            assert evidence.ncbi_gene_id == "673"
            assert evidence.dbsnp_id == "rs113488022"
            assert evidence.clinvar_id == "13961"
            assert evidence.hgvs_genomic == "chr7:g.140453136A>T"

        await client.close()

    @pytest.mark.asyncio
    async def test_query_strategy_with_protein_notation(self):
        """Test that the client tries protein notation query first."""
        client = MyVariantClient()

        with patch.object(client, "_query", new_callable=AsyncMock) as mock_query:
            # First call returns results (protein notation query succeeds)
            mock_query.return_value = {
                "hits": [{"_id": "test", "civic": {}, "clinvar": {}, "cosmic": {}}]
            }

            await client.fetch_evidence("BRAF", "V600E")

            # Verify the first query used protein notation
            first_call_args = mock_query.call_args_list[0]
            # Should be called with "BRAF p.V600E"
            assert (
                "p.V600E" in first_call_args[0][0]
                or "BRAF p.V600E" == first_call_args[0][0]
            )

        await client.close()

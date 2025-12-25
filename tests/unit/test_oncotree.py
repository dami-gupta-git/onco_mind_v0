"""Tests for OncoTree API client."""

import pytest
from unittest.mock import AsyncMock, patch

from oncomind.api.oncotree import OncoTreeAPIError, OncoTreeClient


class TestOncoTreeClient:
    """Tests for OncoTreeClient."""

    @pytest.mark.asyncio
    async def test_context_manager(self):
        """Test async context manager."""
        async with OncoTreeClient() as client:
            assert client._client is not None

        # Client should be closed after exit
        assert client._client is None

    @pytest.mark.asyncio
    async def test_fetch_all_tumor_types_caching(self):
        """Test that fetching all tumor types uses caching."""
        client = OncoTreeClient()

        mock_response = [
            {"code": "NSCLC", "name": "Non-Small Cell Lung Cancer", "tissue": "Lung"},
            {"code": "MEL", "name": "Melanoma", "tissue": "Skin"},
        ]

        with patch.object(client, "_get_client") as mock_get_client:
            mock_http_client = AsyncMock()
            mock_http_client.get = AsyncMock()
            mock_http_client.get.return_value.raise_for_status = lambda: None
            mock_http_client.get.return_value.json = lambda: mock_response
            mock_get_client.return_value = mock_http_client

            # First call should fetch from API
            result1 = await client._fetch_all_tumor_types()
            assert len(result1) == 2

            # Second call should use cache (no additional API call)
            result2 = await client._fetch_all_tumor_types()
            assert len(result2) == 2

            # Verify only one API call was made
            assert mock_http_client.get.call_count == 1

        await client.close()

    @pytest.mark.asyncio
    async def test_get_tumor_type_by_code(self):
        """Test getting tumor type by code."""
        client = OncoTreeClient()

        mock_response = [
            {"code": "NSCLC", "name": "Non-Small Cell Lung Cancer", "tissue": "Lung", "level": 2},
            {"code": "LUAD", "name": "Lung Adenocarcinoma", "tissue": "Lung", "level": 3},
        ]

        with patch.object(client, "_fetch_all_tumor_types", new_callable=AsyncMock) as mock_fetch:
            mock_fetch.return_value = mock_response

            # Test exact match
            nsclc = await client.get_tumor_type_by_code("NSCLC")
            assert nsclc is not None
            assert nsclc["code"] == "NSCLC"
            assert nsclc["name"] == "Non-Small Cell Lung Cancer"

            # Test case-insensitive
            nsclc_lower = await client.get_tumor_type_by_code("nsclc")
            assert nsclc_lower is not None
            assert nsclc_lower["code"] == "NSCLC"

            # Test not found
            unknown = await client.get_tumor_type_by_code("UNKNOWN")
            assert unknown is None

        await client.close()

    @pytest.mark.asyncio
    async def test_resolve_tumor_type(self):
        """Test resolving user input to full tumor type name."""
        client = OncoTreeClient()

        mock_response = [
            {"code": "NSCLC", "name": "Non-Small Cell Lung Cancer", "tissue": "Lung"},
            {"code": "MEL", "name": "Melanoma", "tissue": "Skin"},
        ]

        with patch.object(client, "_fetch_all_tumor_types", new_callable=AsyncMock) as mock_fetch:
            mock_fetch.return_value = mock_response

            # Test code resolution
            resolved_nsclc = await client.resolve_tumor_type("NSCLC")
            assert resolved_nsclc == "Non-Small Cell Lung Cancer"

            # Test case-insensitive code resolution
            resolved_mel = await client.resolve_tumor_type("mel")
            assert resolved_mel == "Melanoma"

            # Test "CODE - Name" format
            resolved_formatted = await client.resolve_tumor_type("NSCLC - Non-Small Cell Lung Cancer")
            assert resolved_formatted == "Non-Small Cell Lung Cancer"

            # Test full name passthrough
            resolved_full = await client.resolve_tumor_type("Non-Small Cell Lung Cancer")
            assert resolved_full == "Non-Small Cell Lung Cancer"

            # Test unknown code (returns original)
            resolved_unknown = await client.resolve_tumor_type("UNKNOWN")
            assert resolved_unknown == "UNKNOWN"

        await client.close()

    @pytest.mark.asyncio
    async def test_api_error_handling(self):
        """Test API error handling."""
        client = OncoTreeClient()

        with patch.object(client, "_get_client") as mock_get_client:
            mock_http_client = AsyncMock()
            mock_http_client.get = AsyncMock()
            mock_http_client.get.side_effect = Exception("Network error")
            mock_get_client.return_value = mock_http_client

            # Should raise OncoTreeAPIError
            with pytest.raises(OncoTreeAPIError):
                await client._fetch_all_tumor_types()

        await client.close()

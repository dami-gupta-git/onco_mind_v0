"""Integration tests for OncoTree API.

Tests validate that the OncoTree API returns expected tumor type
classifications and codes from MSK's standardized cancer ontology.
"""

import pytest

from oncomind.api.oncotree import OncoTreeClient


class TestOncoTreeBasic:
    """Basic OncoTree API tests."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_fetch_all_tumor_types(self):
        """Should fetch all tumor types from OncoTree."""
        async with OncoTreeClient() as client:
            tumor_types = await client._fetch_all_tumor_types()

        assert len(tumor_types) > 500, "OncoTree should have 500+ tumor types"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_tumor_type_structure(self):
        """Tumor types should have expected structure."""
        async with OncoTreeClient() as client:
            tumor_types = await client._fetch_all_tumor_types()

        # Check first tumor type has expected fields
        for tumor in tumor_types[:10]:
            assert "code" in tumor, "Tumor type should have code"
            assert "name" in tumor, "Tumor type should have name"
            # Optional but common fields
            assert "mainType" in tumor or "tissue" in tumor


class TestOncoTreeCodeLookup:
    """Tests for tumor type code lookup."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_nsclc_code(self):
        """NSCLC code should return Non-Small Cell Lung Cancer."""
        async with OncoTreeClient() as client:
            tumor_type = await client.get_tumor_type_by_code("NSCLC")

        assert tumor_type is not None
        assert tumor_type["code"].upper() == "NSCLC"
        assert "lung" in tumor_type["name"].lower()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_luad_code(self):
        """LUAD code should return Lung Adenocarcinoma."""
        async with OncoTreeClient() as client:
            tumor_type = await client.get_tumor_type_by_code("LUAD")

        assert tumor_type is not None
        assert tumor_type["code"].upper() == "LUAD"
        assert "adenocarcinoma" in tumor_type["name"].lower() or "lung" in tumor_type["name"].lower()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_mel_code(self):
        """MEL code should return Melanoma."""
        async with OncoTreeClient() as client:
            tumor_type = await client.get_tumor_type_by_code("MEL")

        assert tumor_type is not None
        assert tumor_type["code"].upper() == "MEL"
        assert "melanoma" in tumor_type["name"].lower()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_crc_code(self):
        """CRC code should return Colorectal Cancer."""
        async with OncoTreeClient() as client:
            tumor_type = await client.get_tumor_type_by_code("CRC")

        assert tumor_type is not None
        assert tumor_type["code"].upper() == "CRC"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_brca_code(self):
        """BRCA code should return Breast Cancer."""
        async with OncoTreeClient() as client:
            tumor_type = await client.get_tumor_type_by_code("BRCA")

        assert tumor_type is not None
        assert tumor_type["code"].upper() == "BRCA"
        assert "breast" in tumor_type["name"].lower()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_case_insensitive(self):
        """Code lookup should be case insensitive."""
        async with OncoTreeClient() as client:
            upper = await client.get_tumor_type_by_code("NSCLC")
            lower = await client.get_tumor_type_by_code("nsclc")
            mixed = await client.get_tumor_type_by_code("Nsclc")

        assert upper is not None
        assert lower is not None
        assert mixed is not None
        assert upper["code"] == lower["code"] == mixed["code"]

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_unknown_code(self):
        """Unknown code should return None."""
        async with OncoTreeClient() as client:
            tumor_type = await client.get_tumor_type_by_code("NOTAREALCODE123")

        assert tumor_type is None


class TestOncoTreeResolveTumorType:
    """Tests for resolve_tumor_type method."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_resolve_code_to_name(self):
        """Should resolve OncoTree code to full name."""
        async with OncoTreeClient() as client:
            name = await client.resolve_tumor_type("NSCLC")

        assert "lung" in name.lower() or "cell" in name.lower()

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_resolve_already_full_name(self):
        """Should return full name unchanged if already a name."""
        async with OncoTreeClient() as client:
            name = await client.resolve_tumor_type("Non-Small Cell Lung Cancer")

        # Since it's not a code, should return as-is
        assert name == "Non-Small Cell Lung Cancer"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_resolve_code_dash_name_format(self):
        """Should extract name from 'CODE - Name' format."""
        async with OncoTreeClient() as client:
            name = await client.resolve_tumor_type("NSCLC - Non-Small Cell Lung Cancer")

        assert name == "Non-Small Cell Lung Cancer"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_resolve_empty_string(self):
        """Should return empty string for empty input."""
        async with OncoTreeClient() as client:
            name = await client.resolve_tumor_type("")

        assert name == ""


class TestOncoTreeCaching:
    """Tests for caching behavior."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cache_populated(self):
        """Cache should be populated after first fetch."""
        async with OncoTreeClient() as client:
            assert "all_tumor_types" not in client._cache

            await client._fetch_all_tumor_types()

            assert "all_tumor_types" in client._cache
            assert len(client._cache["all_tumor_types"]) > 0

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_cache_reused(self):
        """Subsequent fetches should use cache."""
        async with OncoTreeClient() as client:
            # First fetch
            first_result = await client._fetch_all_tumor_types()
            # Second fetch (should use cache)
            second_result = await client._fetch_all_tumor_types()

            assert first_result is second_result  # Same object from cache


class TestOncoTreeContextManager:
    """Tests for async context manager protocol."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_context_manager_creates_client(self):
        """Context manager should create HTTP client."""
        async with OncoTreeClient() as client:
            assert client._client is not None

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_context_manager_closes_client(self):
        """Context manager should close HTTP client on exit."""
        client = OncoTreeClient()
        async with client:
            pass
        assert client._client is None


class TestOncoTreeCommonCancerTypes:
    """Tests for common cancer type codes."""

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_common_codes_exist(self):
        """Common cancer type codes should exist in OncoTree."""
        common_codes = [
            "NSCLC", "LUAD", "LUSC",  # Lung
            "BRCA", "IDC", "ILC",  # Breast
            "CRC", "COAD", "READ",  # Colorectal
            "MEL", "SKCM",  # Melanoma
            "PAAD",  # Pancreatic
            "GBM", "LGG",  # Brain
            "PRAD",  # Prostate
            "OV",  # Ovarian
            "AML", "ALL",  # Leukemia
        ]

        async with OncoTreeClient() as client:
            for code in common_codes:
                tumor_type = await client.get_tumor_type_by_code(code)
                assert tumor_type is not None, f"Common code {code} should exist"

    @pytest.mark.integration
    @pytest.mark.asyncio
    async def test_tumor_types_have_tissue(self):
        """Tumor types should have tissue/organ information."""
        async with OncoTreeClient() as client:
            tumor_types = await client._fetch_all_tumor_types()

        # Most tumor types should have tissue info
        tissue_count = sum(1 for t in tumor_types if t.get("tissue"))
        assert tissue_count > 100, "Most tumor types should have tissue info"

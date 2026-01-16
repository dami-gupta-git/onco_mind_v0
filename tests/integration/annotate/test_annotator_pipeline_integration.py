"""Integration tests for the Annotator pipeline against real APIs."""

import pytest
from pathlib import Path

from oncomind.annotator import Annotator

FIXTURES_DIR = Path(__file__).parent / "fixtures"
SINGLE_VARIANT_VCF = FIXTURES_DIR / "single_variant.vcf"
EGFR_L858R_VCF = FIXTURES_DIR / "egfr_l858r.vcf"


@pytest.mark.integration
class TestAnnotatorPipelineIntegration:
    """Integration tests for the full Annotator pipeline."""

    @pytest.mark.asyncio
    async def test_pipeline_annotates_braf_v600e_from_vcf(self):
        """Test full pipeline annotates BRAF V600E from VCF file."""
        async with Annotator() as annotator:
            result = await annotator.run(SINGLE_VARIANT_VCF)

        assert len(result.variants) == 1
        variant = result.variants[0]

        # Test basic variant info
        assert variant["gene"] == "BRAF"
        assert variant["tumor_type"] == "Melanoma"

        # Test myvariant data is populated
        myvariant = variant["myvariant"]
        assert myvariant is not None

        # Test VCF data
        assert myvariant["vcf"] is not None
        assert myvariant["vcf"]["ref"] == "A"
        assert myvariant["vcf"]["alt"] == "T"

        # Test CADD data with value validation
        assert myvariant["cadd"] is not None
        assert myvariant["cadd"]["phred"] == 32
        assert 0 <= myvariant["cadd"]["phred"] <= 100
        assert myvariant["cadd"]["consequence"] == "NON_SYNONYMOUS"
        assert myvariant["cadd"]["gene"]["genename"] == "BRAF"
        assert myvariant["cadd"]["gene"]["gene_id"].startswith("ENSG")
        assert myvariant["cadd"]["oaa"] == "V"
        assert myvariant["cadd"]["naa"] == "E"

        # Test ClinVar data with value validation
        assert myvariant["clinvar"] is not None
        assert myvariant["clinvar"]["allele_id"] == 29000
        assert myvariant["clinvar"]["allele_id"] >= 1
        assert myvariant["clinvar"]["gene"]["symbol"] == "BRAF"
        assert len(myvariant["clinvar"]["rcv"]) > 0
        assert myvariant["clinvar"]["rcv"][0]["accession"].startswith("RCV")
        assert myvariant["clinvar"]["rcv"][0]["clinical_significance"] == "Pathogenic"

        # Test COSMIC data with pattern validation
        assert myvariant["cosmic"] is not None
        assert myvariant["cosmic"]["cosmic_id"] == "COSM476"
        assert myvariant["cosmic"]["cosmic_id"].startswith("COSM")

        # Test gnomAD data with range validation
        assert myvariant["gnomad_exome"] is not None
        assert 0 <= myvariant["gnomad_exome"]["af"]["af"] <= 1

        # Test timing data
        assert "timings" in variant
        assert "myvariant" in variant["timings"]
        assert isinstance(variant["timings"]["myvariant"], float)
        assert variant["timings"]["myvariant"] >= 0

    @pytest.mark.asyncio
    async def test_pipeline_result_is_json_serializable(self):
        """Test pipeline result can be serialized to JSON."""
        import json

        async with Annotator() as annotator:
            result = await annotator.run(SINGLE_VARIANT_VCF)

        output = result.to_dict()

        # Should not raise
        json_str = json.dumps(output)
        assert json_str is not None

        # Verify structure
        parsed = json.loads(json_str)
        assert "variants" in parsed
        assert len(parsed["variants"]) == 1
        assert parsed["variants"][0]["gene"] == "BRAF"

    @pytest.mark.asyncio
    async def test_pipeline_annotates_egfr_l858r_from_vcf(self):
        """Test full pipeline annotates EGFR L858R from VCF file.

        This tests the three-letter amino acid fallback since EGFR L858R
        is indexed as p.Leu858Arg in MyVariant.info.
        """
        async with Annotator() as annotator:
            result = await annotator.run(EGFR_L858R_VCF)

        assert len(result.variants) == 1
        variant = result.variants[0]

        # Test basic variant info
        assert variant["gene"] == "EGFR"
        assert variant["tumor_type"] == "NSCLC"

        # Test myvariant data is populated
        myvariant = variant["myvariant"]
        assert myvariant is not None

        # Test VCF data
        assert myvariant["vcf"] is not None
        assert myvariant["vcf"]["ref"] == "T"
        assert myvariant["vcf"]["alt"] == "G"

        # Test CADD data with value validation
        assert myvariant["cadd"] is not None
        assert myvariant["cadd"]["phred"] == 33
        assert 0 <= myvariant["cadd"]["phred"] <= 100
        assert myvariant["cadd"]["gene"]["genename"] == "EGFR"
        assert myvariant["cadd"]["gene"]["gene_id"].startswith("ENSG")
        assert myvariant["cadd"]["oaa"] == "L"
        assert myvariant["cadd"]["naa"] == "R"

        # Test ClinVar data with value validation
        assert myvariant["clinvar"] is not None
        assert myvariant["clinvar"]["allele_id"] >= 1
        assert myvariant["clinvar"]["gene"]["symbol"] == "EGFR"
        assert len(myvariant["clinvar"]["rcv"]) > 0
        assert myvariant["clinvar"]["rcv"][0]["accession"].startswith("RCV")

        # Test COSMIC data with pattern validation
        assert myvariant["cosmic"] is not None
        assert myvariant["cosmic"]["cosmic_id"].startswith("COSM")

        # Test timing data
        assert "timings" in variant
        assert "myvariant" in variant["timings"]
        assert isinstance(variant["timings"]["myvariant"], float)
        assert variant["timings"]["myvariant"] >= 0

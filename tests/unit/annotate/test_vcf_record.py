"""Unit tests for VCFRecord."""

from oncomind.annotator.vcf_record import VCFRecord


class TestVCFRecordToParsedVariant:
    """Tests for VCFRecord.to_parsed_variant()."""

    def test_snv_conversion(self):
        """Test SNV is correctly converted to ParsedVariant."""
        record = VCFRecord(
            chrom="7", pos=140753336, id="COSV54736915", ref="T", alt="A"
        )

        parsed = record.to_parsed_variant()

        assert parsed.gene == "UNKNOWN"
        assert parsed.variant == "7:140753336:T>A"
        assert parsed.variant_type == "snv"
        assert parsed.parse_confidence == 0.5

    def test_indel_conversion(self):
        """Test indel is correctly converted to ParsedVariant."""
        record = VCFRecord(
            chrom="7", pos=55242465, id=None, ref="GAATTAAGAGAAGC", alt="G"
        )

        parsed = record.to_parsed_variant()

        assert parsed.variant == "7:55242465:GAATTAAGAGAAGC>G"
        assert parsed.variant_type == "indel"

    def test_null_id_handled(self):
        """Test that null ID is handled."""
        record = VCFRecord(chrom="1", pos=12345, id=None, ref="C", alt="T")

        parsed = record.to_parsed_variant()

        assert parsed.variant == "1:12345:C>T"

    def test_tumor_type_from_info(self):
        """Test tumor_type is extracted from INFO field."""
        record = VCFRecord(
            chrom="7",
            pos=140753336,
            id="COSV54736915",
            ref="T",
            alt="A",
            info={"GENE": "BRAF", "TUMOR_TYPE": "NSCLC"},
        )

        parsed = record.to_parsed_variant()

        assert parsed.gene == "BRAF"
        assert parsed.tumor_type == "NSCLC"
        assert parsed.parse_confidence == 0.8

    def test_gene_from_info(self):
        """Test gene is extracted from INFO field."""
        record = VCFRecord(
            chrom="7",
            pos=140753336,
            id=None,
            ref="T",
            alt="A",
            info={"GENE": "BRAF"},
        )

        parsed = record.to_parsed_variant()

        assert parsed.gene == "BRAF"
        assert parsed.parse_confidence == 0.8

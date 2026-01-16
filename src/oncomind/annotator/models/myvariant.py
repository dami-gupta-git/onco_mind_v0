"""Pydantic models for MyVariant.info API responses."""

from typing import Any

from pydantic import BaseModel, Field


class VcfData(BaseModel):
    """VCF reference and alternate allele data."""

    ref: str = Field(min_length=1)
    alt: str = Field(min_length=1)


class CaddDna(BaseModel):
    """CADD DNA shape parameters."""

    helt: float
    mgw: float
    prot: float
    roll: float


class CaddCds(BaseModel):
    """CADD CDS position data."""

    cds_pos: int = Field(ge=1)
    cdna_pos: int = Field(ge=1)
    rel_cdna_pos: float | None = Field(default=None, ge=0, le=1)
    rel_cds_pos: float | None = Field(default=None, ge=0, le=1)


class CaddProt(BaseModel):
    """CADD protein position data."""

    protpos: int = Field(ge=1)
    rel_prot_pos: float | None = Field(default=None, ge=0, le=1)


class CaddGene(BaseModel):
    """CADD gene information."""

    genename: str = Field(min_length=1)
    gene_id: str = Field(pattern=r"^ENSG\d+$")
    feature_id: str = Field(pattern=r"^ENST\d+$")
    ccds_id: str | None = None
    cds: CaddCds | None = None
    prot: CaddProt | None = None


class CaddData(BaseModel):
    """CADD (Combined Annotation Dependent Depletion) scores."""

    phred: float = Field(ge=0, le=100)
    rawscore: float
    # API returns single value or list depending on number of transcripts
    consequence: str | list[str]
    annotype: str | list[str] | None = None
    dna: CaddDna | None = None
    gene: CaddGene | list[CaddGene] | None = None
    oaa: str | None = Field(default=None, min_length=1, max_length=1)
    naa: str | None = Field(default=None, min_length=1, max_length=1)
    type: str | None = None
    chrom: int | str | None = None
    pos: int | None = Field(default=None, ge=1)
    ref: str | None = None
    alt: str | None = None


class ClinvarGene(BaseModel):
    """ClinVar gene information."""

    symbol: str = Field(min_length=1)
    id: str


class ClinvarRcv(BaseModel):
    """ClinVar RCV record."""

    accession: str = Field(pattern=r"^RCV\d+$")
    clinical_significance: str


class ClinvarData(BaseModel):
    """ClinVar clinical variant data."""

    allele_id: int = Field(ge=1)
    gene: ClinvarGene | None = None
    chrom: str | None = None
    alt: str | None = None
    # API returns single dict or list depending on number of RCV records
    rcv: ClinvarRcv | list[ClinvarRcv] | None = None


class CosmicData(BaseModel):
    """COSMIC (Catalogue Of Somatic Mutations In Cancer) data."""

    cosmic_id: str = Field(pattern=r"^COSM\d+$")
    mut_nt: str | None = None
    chrom: str | None = None
    ref: str | None = None
    alt: str | None = None


class Polyphen2Hdiv(BaseModel):
    """PolyPhen2 HDIV prediction."""

    # API returns single value or list depending on number of transcripts
    pred: str | list[str]  # D=Damaging, B=Benign, P=Possibly damaging
    score: float | list[float]


class Polyphen2(BaseModel):
    """PolyPhen2 predictions."""

    hdiv: Polyphen2Hdiv | None = None


class AlphaMissense(BaseModel):
    """AlphaMissense pathogenicity prediction."""

    # API returns single value or list depending on number of transcripts
    pred: str | list[str]  # P=Pathogenic, B=Benign, A=Ambiguous
    score: float | list[float]


class DbnsfpData(BaseModel):
    """dbNSFP functional prediction data."""

    polyphen2: Polyphen2 | None = None
    alphamissense: AlphaMissense | None = None
    sift: dict[str, Any] | None = None


class DbsnpData(BaseModel):
    """dbSNP reference SNP data."""

    rsid: str = Field(pattern=r"^rs\d+$")


class GnomadAf(BaseModel):
    """gnomAD allele frequency."""

    af: float = Field(ge=0, le=1)


class GnomadExomeData(BaseModel):
    """gnomAD exome allele frequency data."""

    af: GnomadAf


class MyVariantAnnotation(BaseModel):
    """Complete MyVariant.info annotation response."""

    vcf: VcfData | None = None
    cadd: CaddData | None = None
    clinvar: ClinvarData | None = None
    cosmic: CosmicData | None = None
    dbnsfp: DbnsfpData | None = None
    dbsnp: DbsnpData | None = None
    gnomad_exome: GnomadExomeData | None = None

    @classmethod
    def from_hit(cls, hit: dict[str, Any] | None) -> "MyVariantAnnotation":
        """Create from a MyVariant API hit dictionary.

        Args:
            hit: Dictionary from API response hits array, or None

        Returns:
            MyVariantAnnotation with validated fields
        """
        if hit is None:
            return cls()

        return cls(
            vcf=hit.get("vcf"),
            cadd=hit.get("cadd"),
            clinvar=hit.get("clinvar"),
            cosmic=hit.get("cosmic"),
            dbnsfp=hit.get("dbnsfp"),
            dbsnp=hit.get("dbsnp"),
            gnomad_exome=hit.get("gnomad_exome"),
        )


__all__ = [
    "MyVariantAnnotation",
    "VcfData",
    "CaddData",
    "CaddDna",
    "CaddGene",
    "CaddCds",
    "CaddProt",
    "ClinvarData",
    "ClinvarGene",
    "ClinvarRcv",
    "CosmicData",
    "DbnsfpData",
    "Polyphen2",
    "Polyphen2Hdiv",
    "AlphaMissense",
    "DbsnpData",
    "GnomadExomeData",
    "GnomadAf",
]

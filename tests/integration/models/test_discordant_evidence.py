"""Tests for discordant evidence detection logic.

Tests the detection of conflicting evidence between different sources:
- ClinVar internal conflicts (pathogenic vs benign at same locus level)
- ClinVar vs CIViC conflicts (benign vs actionable at variant level)
- Cross-source drug response conflicts (sensitive vs resistant at variant level)
- FDA vs VICC conflicts (significant severity due to VICC reliability concerns)
"""

import pytest

from oncomind.insight_builder.gap_detector import _detect_discordant_evidence_internal


class TestDiscordantEvidenceDetection:
    """Tests for discordant evidence detection logic."""

    def test_no_conflict_when_same_source(self):
        """No conflict should be flagged when only CIViC data exists."""
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="KRAS:G12D",
                gene="KRAS",
                variant="G12D"
            ),
            civic_evidence=[
                CIViCEvidence(drugs=["Erlotinib"], clinical_significance="Resistance"),
                CIViCEvidence(drugs=["Erlotinib"], clinical_significance="Sensitivity"),
            ]
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)
        all_conflicts = high_conflicts + significant_conflicts

        # Should NOT flag intra-source conflict
        assert len(all_conflicts) == 0, f"Intra-source conflict should not be flagged: {all_conflicts}"

    def test_combination_therapy_not_flagged(self):
        """Combination therapies should not be mixed with monotherapy."""
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="KRAS:G12D",
                gene="KRAS",
                variant="G12D"
            ),
            civic_evidence=[
                # Monotherapy resistance
                CIViCEvidence(drugs=["Erlotinib"], clinical_significance="Resistance"),
                # Combination therapy sensitivity (different context)
                CIViCEvidence(drugs=["Erlotinib", "Sotorasib"], clinical_significance="Sensitivity"),
            ]
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)
        all_conflicts = high_conflicts + significant_conflicts

        # Should NOT flag conflict - different contexts (mono vs combo)
        assert len(all_conflicts) == 0, f"Combo therapy should not create conflict: {all_conflicts}"

    def test_clinvar_benign_vs_civic_actionable_conflict_variant_level(self):
        """ClinVar benign + CIViC actionable at VARIANT level should flag conflict."""
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCEvidence, CIViCAssertionEvidence
        from oncomind.models.evidence.clinvar import ClinVarEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers
        from oncomind.models.evidence.base import EvidenceLevel

        # Both ClinVar and CIViC at variant level - should flag conflict
        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="TEST:V123M",
                gene="TEST",
                variant="V123M"
            ),
            clinvar_entries=[
                ClinVarEvidence(
                    clinical_significance="Benign",
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ],
            clinvar_significance="Benign",
            civic_assertions=[
                CIViCAssertionEvidence(
                    assertion_type="PREDICTIVE",
                    therapies=["Erlotinib"],
                    is_sensitivity=True,
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ]
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)

        # Should flag ClinVar benign vs CIViC actionable conflict at variant level (HIGH severity)
        assert len(high_conflicts) >= 1, "Should flag ClinVar benign vs CIViC actionable conflict at variant level"
        assert any("benign" in c.lower() and "civic" in c.lower() for c in high_conflicts), \
            f"Should mention ClinVar benign vs CIViC conflict: {high_conflicts}"
        # Should mention it's at variant level
        assert any("variant level" in c.lower() for c in high_conflicts), \
            f"Should mention 'variant level' in conflict: {high_conflicts}"

    def test_clinvar_benign_gene_level_no_civic_conflict(self):
        """ClinVar benign at GENE level + CIViC actionable at VARIANT level should NOT flag conflict.

        Gene-level benign doesn't mean THIS variant is benign - it could be a different variant.
        """
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCAssertionEvidence
        from oncomind.models.evidence.clinvar import ClinVarEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers
        from oncomind.models.evidence.base import EvidenceLevel

        # ClinVar at gene level, CIViC at variant level - should NOT flag conflict
        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="TEST:V123M",
                gene="TEST",
                variant="V123M"
            ),
            clinvar_entries=[
                ClinVarEvidence(
                    clinical_significance="Benign",
                    locus_variant_match=EvidenceLevel(level="gene"),  # Gene level - different variant
                ),
            ],
            clinvar_significance="Benign",
            civic_assertions=[
                CIViCAssertionEvidence(
                    assertion_type="PREDICTIVE",
                    therapies=["Erlotinib"],
                    is_sensitivity=True,
                    locus_variant_match=EvidenceLevel(level="variant"),  # Variant level
                ),
            ]
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)
        all_conflicts = high_conflicts + significant_conflicts

        # Should NOT flag ClinVar vs CIViC conflict - different locus levels
        clinvar_civic_conflicts = [c for c in all_conflicts if "benign" in c.lower() and "civic" in c.lower()]
        assert len(clinvar_civic_conflicts) == 0, \
            f"Gene-level ClinVar benign should not conflict with variant-level CIViC: {all_conflicts}"

    def test_clinvar_pathogenic_no_civic_conflict(self):
        """ClinVar pathogenic should NOT conflict with CIViC actionable."""
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCEvidence
        from oncomind.models.evidence.clinvar import ClinVarEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers
        from oncomind.models.evidence.base import EvidenceLevel

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="TEST:V123M",
                gene="TEST",
                variant="V123M"
            ),
            clinvar_entries=[
                ClinVarEvidence(
                    clinical_significance="Pathogenic",
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ],
            clinvar_significance="Pathogenic",
            civic_evidence=[
                CIViCEvidence(
                    evidence_type="PREDICTIVE",
                    drugs=["Erlotinib"],
                    clinical_significance="Sensitivity",
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ]
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)
        all_conflicts = high_conflicts + significant_conflicts

        # Should NOT flag conflict - both agree variant is significant
        clinvar_civic_conflicts = [c for c in all_conflicts if "benign" in c.lower() and "civic" in c.lower()]
        assert len(clinvar_civic_conflicts) == 0, \
            f"Pathogenic ClinVar should not conflict with CIViC: {all_conflicts}"

    def test_clinvar_internal_conflict_variant_level(self):
        """ClinVar pathogenic + benign at VARIANT level should flag internal conflict."""
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.clinvar import ClinVarEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers
        from oncomind.models.evidence.base import EvidenceLevel

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="TEST:V123M",
                gene="TEST",
                variant="V123M"
            ),
            clinvar_entries=[
                ClinVarEvidence(
                    clinical_significance="Pathogenic",
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
                ClinVarEvidence(
                    clinical_significance="Benign",
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ],
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)

        # Should flag ClinVar internal conflict at variant level (HIGH severity)
        clinvar_conflicts = [c for c in high_conflicts if "clinvar" in c.lower() and "conflicting" in c.lower()]
        assert len(clinvar_conflicts) >= 1, \
            f"Should flag ClinVar internal conflict: {high_conflicts}"
        # Should mention variant level
        assert any("variant level" in c.lower() for c in clinvar_conflicts), \
            f"Should mention 'variant level' in conflict: {high_conflicts}"

    def test_clinvar_internal_conflict_gene_level_not_flagged(self):
        """ClinVar pathogenic + benign at GENE level (different variants) should NOT flag conflict.

        Different variants in the same gene can legitimately have different pathogenicity.
        """
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.clinvar import ClinVarEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers
        from oncomind.models.evidence.base import EvidenceLevel

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="TEST:V123M",
                gene="TEST",
                variant="V123M"
            ),
            clinvar_entries=[
                ClinVarEvidence(
                    clinical_significance="Pathogenic",
                    locus_variant_match=EvidenceLevel(level="gene"),  # Different variant
                ),
                ClinVarEvidence(
                    clinical_significance="Benign",
                    locus_variant_match=EvidenceLevel(level="gene"),  # Different variant
                ),
            ],
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)
        all_conflicts = high_conflicts + significant_conflicts

        # Should NOT flag ClinVar internal conflict - gene-level means different variants
        clinvar_conflicts = [c for c in all_conflicts if "clinvar" in c.lower() and "conflicting" in c.lower()]
        assert len(clinvar_conflicts) == 0, \
            f"Gene-level ClinVar should not flag conflict (different variants): {all_conflicts}"

    def test_cross_source_drug_conflict_variant_level(self):
        """Cross-source drug conflict at VARIANT level should flag HIGH conflict."""
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCEvidence
        from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers
        from oncomind.models.evidence.base import EvidenceLevel

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="TEST:V123M",
                gene="TEST",
                variant="V123M"
            ),
            civic_evidence=[
                CIViCEvidence(
                    drugs=["Erlotinib"],
                    clinical_significance="Sensitivity",
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ],
            cgi_biomarkers=[
                CGIBiomarkerEvidence(
                    drug="Erlotinib",
                    association="Resistance",
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ],
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)

        # Should flag cross-source drug conflict at variant level (HIGH severity)
        drug_conflicts = [c for c in high_conflicts if "erlotinib" in c.lower()]
        assert len(drug_conflicts) >= 1, \
            f"Should flag cross-source drug conflict: {high_conflicts}"
        # Should mention it's at variant level
        assert any("variant level" in c.lower() for c in drug_conflicts), \
            f"Should mention 'variant level' in conflict: {drug_conflicts}"

    def test_cross_source_drug_conflict_gene_level_not_flagged(self):
        """Cross-source drug conflict at GENE level should NOT flag conflict.

        Gene-level evidence may be for different variants that behave differently.
        """
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.civic import CIViCEvidence
        from oncomind.models.evidence.cgi import CGIBiomarkerEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers
        from oncomind.models.evidence.base import EvidenceLevel

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="TEST:V123M",
                gene="TEST",
                variant="V123M"
            ),
            civic_evidence=[
                CIViCEvidence(
                    drugs=["Erlotinib"],
                    clinical_significance="Sensitivity",
                    locus_variant_match=EvidenceLevel(level="gene"),  # Gene level
                ),
            ],
            cgi_biomarkers=[
                CGIBiomarkerEvidence(
                    drug="Erlotinib",
                    association="Resistance",
                    locus_variant_match=EvidenceLevel(level="gene"),  # Gene level
                ),
            ],
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)
        all_conflicts = high_conflicts + significant_conflicts

        # Should NOT flag cross-source drug conflict at gene level
        drug_conflicts = [c for c in all_conflicts if "erlotinib" in c.lower() and "conflict" in c.lower()]
        assert len(drug_conflicts) == 0, \
            f"Gene-level drug conflict should not be flagged: {all_conflicts}"

    def test_fda_vs_vicc_conflict_significant_severity(self):
        """FDA sensitive vs VICC resistant at variant level should flag SIGNIFICANT conflict.

        VICC data is less reliable than FDA, so this is SIGNIFICANT not HIGH.
        """
        from oncomind.models.evidence import Evidence
        from oncomind.models.evidence.fda_biomarker import FDABiomarkerEvidence, SpecificityLevel, BiomarkerRequirement
        from oncomind.models.evidence.vicc import VICCEvidence
        from oncomind.models.evidence.evidence import VariantIdentifiers
        from oncomind.models.evidence.base import EvidenceLevel

        evidence = Evidence(
            identifiers=VariantIdentifiers(
                variant_id="TEST:V123M",
                gene="TEST",
                variant="V123M"
            ),
            fda_biomarker_evidence=[
                FDABiomarkerEvidence(
                    drug_name="Erlotinib",
                    gene="TEST",
                    specificity=SpecificityLevel.VARIANT,
                    requirement=BiomarkerRequirement.REQUIRED_POSITIVE,
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ],
            vicc_evidence=[
                VICCEvidence(
                    drugs=["Erlotinib"],
                    response_type="Resistance",
                    locus_variant_match=EvidenceLevel(level="variant"),
                ),
            ],
        )

        high_conflicts, significant_conflicts = _detect_discordant_evidence_internal(evidence)

        # Should flag FDA vs VICC conflict in SIGNIFICANT (not HIGH)
        fda_vicc_conflicts = [c for c in significant_conflicts if "fda" in c.lower() and "vicc" in c.lower()]
        assert len(fda_vicc_conflicts) >= 1, \
            f"Should flag FDA vs VICC conflict as SIGNIFICANT: {significant_conflicts}"
        # Should NOT be in HIGH conflicts
        fda_vicc_high = [c for c in high_conflicts if "fda" in c.lower() and "vicc" in c.lower()]
        assert len(fda_vicc_high) == 0, \
            f"FDA vs VICC conflict should NOT be HIGH severity: {high_conflicts}"

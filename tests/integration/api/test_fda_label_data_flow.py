"""Integration tests for FDA label data flow.

Verifies that all fields (especially initial_approval_date) are correctly
passed through the entire data pipeline:
1. OpenFDA API fetch -> FDALabelData
2. FDALabelData -> FDALabelEvidence (in evidence_aggregator)
3. FDALabelEvidence -> backend serialization (for Streamlit UI)

Run with: pytest tests/integration/api/test_fda_label_data_flow.py -v
"""

import pytest
from unittest.mock import patch, MagicMock

from oncomind.api.fda_label_service import (
    fetch_drug_label_from_openfda,
    get_fda_labels_for_drugs,
    FDALabelData,
)
from oncomind.models.evidence.fda import FDALabelEvidence, ClinicalStudyEvidence


class TestFDALabelDataFlow:
    """Test that all fields flow through the entire pipeline."""

    def test_fetch_returns_initial_approval_date(self):
        """Verify fetch_drug_label_from_openfda extracts initial_approval_date."""
        # This test hits the real API - skip in CI if needed
        result = fetch_drug_label_from_openfda("capivasertib")

        assert result is not None, "Should return data for capivasertib"
        assert "initial_approval_date" in result, "Should have initial_approval_date key"
        assert result["initial_approval_date"] is not None, "initial_approval_date should not be None"
        # Format should be YYYY-MM-DD
        assert len(result["initial_approval_date"]) == 10, "Date should be in YYYY-MM-DD format"
        assert result["initial_approval_date"][4] == "-", "Date should have dash at position 4"

    def test_fda_label_data_has_initial_approval_date(self):
        """Verify FDALabelData dataclass has initial_approval_date field."""
        data = FDALabelData(
            drug="test_drug",
            gene="TEST",
            initial_approval_date="2023-11-14",
        )

        assert data.initial_approval_date == "2023-11-14"

    def test_fda_label_data_to_dict_includes_initial_approval_date(self):
        """Verify FDALabelData.to_dict() includes initial_approval_date."""
        data = FDALabelData(
            drug="test_drug",
            gene="TEST",
            initial_approval_date="2023-11-14",
        )

        result = data.to_dict()

        assert "initial_approval_date" in result
        assert result["initial_approval_date"] == "2023-11-14"

    def test_fda_label_evidence_has_initial_approval_date(self):
        """Verify FDALabelEvidence model has initial_approval_date field."""
        evidence = FDALabelEvidence(
            drug="test_drug",
            gene="TEST",
            initial_approval_date="2023-11-14",
        )

        assert evidence.initial_approval_date == "2023-11-14"

    def test_fda_label_evidence_model_dump_includes_initial_approval_date(self):
        """Verify FDALabelEvidence.model_dump() includes initial_approval_date."""
        evidence = FDALabelEvidence(
            drug="test_drug",
            gene="TEST",
            initial_approval_date="2023-11-14",
        )

        result = evidence.model_dump()

        assert "initial_approval_date" in result
        assert result["initial_approval_date"] == "2023-11-14"

    def test_backend_serialization_includes_initial_approval_date(self):
        """Verify backend serialization includes initial_approval_date.

        This simulates what backend.py does when serializing fda_labels.
        """
        # Create FDALabelEvidence with all fields
        evidence = FDALabelEvidence(
            drug="capivasertib",
            gene="AKT1",
            brand_name="TRUQAP",
            initial_approval_date="2023-11-14",
            clinical_studies=ClinicalStudyEvidence(
                trial_name="CAPItello-291",
                nct_id="NCT04305496",
            ),
        )

        # Simulate backend serialization (from backend.py lines 252-288)
        serialized = {
            "drug": evidence.drug,
            "gene": evidence.gene,
            "brand_name": evidence.brand_name,
            "generic_name": evidence.generic_name,
            "manufacturer": evidence.manufacturer,
            "indications_and_usage": evidence.indications_and_usage,
            "initial_approval_date": evidence.initial_approval_date,
            "clinical_studies": {
                "trial_name": evidence.clinical_studies.trial_name,
                "nct_id": evidence.clinical_studies.nct_id,
            } if evidence.clinical_studies else None,
            "last_label_update": evidence.last_label_update,
            "update_reason": evidence.update_reason,
        }

        assert serialized["initial_approval_date"] == "2023-11-14"
        assert serialized["clinical_studies"]["trial_name"] == "CAPItello-291"


class TestFDALabelDataFieldParity:
    """Ensure FDALabelData and FDALabelEvidence have matching fields."""

    def test_both_models_have_initial_approval_date(self):
        """Both FDALabelData and FDALabelEvidence must have initial_approval_date."""
        # FDALabelData (dataclass in fda_label_service.py)
        fda_data = FDALabelData(drug="test", gene="TEST")
        assert hasattr(fda_data, "initial_approval_date"), \
            "FDALabelData must have initial_approval_date field"

        # FDALabelEvidence (Pydantic model in models/evidence/fda.py)
        fda_evidence = FDALabelEvidence(drug="test", gene="TEST")
        assert hasattr(fda_evidence, "initial_approval_date"), \
            "FDALabelEvidence must have initial_approval_date field"

    def test_required_fields_present_in_both(self):
        """Verify key fields exist in both models."""
        required_fields = [
            "drug",
            "gene",
            "brand_name",
            "generic_name",
            "manufacturer",
            "indications_and_usage",
            "initial_approval_date",
            "clinical_studies",
            "mechanism_of_action",
            "adverse_reactions",
        ]

        fda_data = FDALabelData(drug="test", gene="TEST")
        fda_evidence = FDALabelEvidence(drug="test", gene="TEST")

        for field in required_fields:
            assert hasattr(fda_data, field), \
                f"FDALabelData missing field: {field}"
            assert hasattr(fda_evidence, field), \
                f"FDALabelEvidence missing field: {field}"


class TestEvidenceAggregatorConversion:
    """Test the conversion from FDALabelData to FDALabelEvidence."""

    def test_conversion_preserves_initial_approval_date(self):
        """Verify evidence_aggregator passes initial_approval_date."""
        from oncomind.insight_builder.evidence_aggregator import EvidenceAggregator
        from oncomind.api.fda_label_service import (
            FDALabelData,
            ClinicalStudyData,
            MechanismOfActionData,
        )

        # Create FDALabelData with initial_approval_date
        label_data = FDALabelData(
            drug="capivasertib",
            gene="AKT1",
            brand_name="TRUQAP",
            initial_approval_date="2023-11-14",
            clinical_studies=ClinicalStudyData(
                trial_name="CAPItello-291",
                nct_id="NCT04305496",
            ),
            mechanism_of_action=MechanismOfActionData(
                targets=["AKT1", "AKT2", "AKT3"],
            ),
        )

        # Create aggregator and convert
        aggregator = EvidenceAggregator.__new__(EvidenceAggregator)
        aggregator.tumor_type = "Breast Cancer"

        # Call the conversion method
        result = aggregator._convert_fda_labels_to_evidence([label_data])

        assert len(result) == 1
        evidence = result[0]

        # Verify initial_approval_date was passed through
        assert evidence.initial_approval_date == "2023-11-14", \
            "initial_approval_date must be passed from FDALabelData to FDALabelEvidence"
        assert evidence.drug == "capivasertib"
        assert evidence.brand_name == "TRUQAP"
        assert evidence.clinical_studies is not None
        assert evidence.clinical_studies.trial_name == "CAPItello-291"

"""Unit tests for FDA label parsing functions.

Tests the parsing functions in fda_label_service.py that extract
structured data from FDA label text fields.
"""

import pytest

from oncomind.api.fda_label_service import (
    parse_clinical_studies,
    parse_mechanism_of_action,
    parse_adverse_reactions,
    parse_recent_major_changes,
    ClinicalStudyData,
    MechanismOfActionData,
    AdverseReactionsData,
    RecentMajorChangesData,
)


class TestParseClinicalStudies:
    """Tests for parse_clinical_studies function."""

    def test_parse_none_returns_none(self):
        """Empty or None input should return None."""
        assert parse_clinical_studies(None) is None
        assert parse_clinical_studies("") is None

    def test_parse_nct_id(self):
        """Should extract NCT ID from clinical studies text."""
        text = "The efficacy was evaluated in CAPItello-291 (NCT04305496), a randomized trial."
        result = parse_clinical_studies(text)

        assert result is not None
        assert result.nct_id == "NCT04305496"

    def test_parse_trial_name(self):
        """Should extract trial name like CAPItello-291."""
        text = "The efficacy was evaluated in CAPItello-291 (NCT04305496), a randomized trial."
        result = parse_clinical_studies(text)

        assert result is not None
        assert result.trial_name == "CAPItello-291"

    def test_parse_patient_count(self):
        """Should extract patient count from biomarker-altered subset."""
        text = "Of these, 289 patients had tumors with eligible PIK3CA/AKT1/PTEN alterations."
        result = parse_clinical_studies(text)

        assert result is not None
        assert result.patients_n == 289

    def test_parse_hazard_ratio(self):
        """Should extract hazard ratio with confidence interval."""
        text = "Hazard Ratio (stratified) was 0.50 (0.38, 0.65) with p-value <0.0001."
        result = parse_clinical_studies(text)

        assert result is not None
        assert result.hazard_ratio == 0.50
        assert result.hazard_ratio_ci == (0.38, 0.65)

    def test_parse_orr(self):
        """Should extract objective response rates."""
        text = "The confirmed ORR was 23% versus 8% in the treatment and control arms."
        result = parse_clinical_studies(text)

        assert result is not None
        assert result.orr_treatment == 23.0
        assert result.orr_control == 8.0

    def test_parse_biomarker_breakdown(self):
        """Should extract biomarker breakdown percentages."""
        text = "The alterations were: PIK3CA: 76%, AKT1: 13%, PTEN: 11%."
        result = parse_clinical_studies(text)

        assert result is not None
        assert result.biomarker_breakdown == {"PIK3CA": 76.0, "AKT1": 13.0, "PTEN": 11.0}

    def test_parse_full_capivasertib_text(self):
        """Test parsing realistic capivasertib clinical studies text."""
        text = """
        The efficacy of TRUQAP in combination with fulvestrant was evaluated in
        CAPItello-291 (NCT04305496), a randomized trial of patients with HR-positive,
        HER2-negative breast cancer. Of the 708 patients enrolled, 289 patients had
        tumors with eligible PIK3CA/AKT1/PTEN alterations. The primary endpoint was
        PFS. Hazard Ratio was 0.50 (0.38, 0.65). Confirmed ORR was 22% vs 12%.
        The biomarker breakdown was PIK3CA: 76%, AKT1: 13%, PTEN: 11%.
        """
        result = parse_clinical_studies(text)

        assert result is not None
        assert result.nct_id == "NCT04305496"
        assert result.trial_name == "CAPItello-291"
        assert result.patients_n == 289
        assert result.hazard_ratio == 0.50
        assert result.hazard_ratio_ci == (0.38, 0.65)
        assert result.orr_treatment == 22.0
        assert result.orr_control == 12.0
        assert result.biomarker_breakdown.get("PIK3CA") == 76.0
        assert result.biomarker_breakdown.get("AKT1") == 13.0

    def test_no_meaningful_data_returns_none(self):
        """Text without meaningful clinical data should return None."""
        text = "This is some general text about cancer treatment without any study data."
        result = parse_clinical_studies(text)

        assert result is None


class TestParseMechanismOfAction:
    """Tests for parse_mechanism_of_action function."""

    def test_parse_none_returns_none(self):
        """Empty or None input should return None."""
        assert parse_mechanism_of_action(None) is None
        assert parse_mechanism_of_action("") is None

    def test_parse_akt_targets(self):
        """Should extract AKT isoform targets."""
        text = "Capivasertib is an inhibitor of AKT (AKT1, AKT2 and AKT3) serine/threonine kinase."
        result = parse_mechanism_of_action(text)

        assert result is not None
        assert result.targets == ["AKT1", "AKT2", "AKT3"]

    def test_parse_mechanism_description(self):
        """Should extract mechanism description."""
        text = "It inhibits phosphorylation of downstream substrates in the PI3K pathway."
        result = parse_mechanism_of_action(text)

        assert result is not None
        assert "phosphorylation" in result.mechanism.lower()

    def test_parse_preclinical_findings(self):
        """Should extract preclinical findings."""
        text = "In vitro, capivasertib reduced cell proliferation. In vivo, it reduced tumor growth."
        result = parse_mechanism_of_action(text)

        assert result is not None
        assert result.preclinical is not None
        assert "In vitro" in result.preclinical

    def test_parse_full_mechanism_text(self):
        """Test parsing realistic mechanism of action text."""
        text = """
        Capivasertib is an inhibitor of all isoforms of AKT (AKT1, AKT2 and AKT3)
        serine/threonine kinase, and inhibits phosphorylation of downstream substrates.
        In vitro, capivasertib reduced cell proliferation in AKT1-mutant cell lines.
        """
        result = parse_mechanism_of_action(text)

        assert result is not None
        assert "AKT1" in result.targets
        assert "AKT2" in result.targets
        assert "AKT3" in result.targets
        assert result.preclinical is not None

    def test_no_meaningful_data_returns_none(self):
        """Text without meaningful mechanism data should return None."""
        text = "This is general text without mechanism information."
        result = parse_mechanism_of_action(text)

        assert result is None


class TestParseAdverseReactions:
    """Tests for parse_adverse_reactions function."""

    def test_parse_none_returns_none(self):
        """Empty or None input should return None."""
        assert parse_adverse_reactions(None) is None
        assert parse_adverse_reactions("") is None

    def test_parse_serious_rate(self):
        """Should extract serious adverse reactions rate."""
        text = "Serious adverse reactions occurred in 35% of patients receiving TRUQAP."
        result = parse_adverse_reactions(text)

        assert result is not None
        assert result.serious_rate == 35.0

    def test_parse_discontinuation_rate(self):
        """Should extract discontinuation rate."""
        text = "Permanent discontinuation due to adverse reactions occurred in 13% of patients."
        result = parse_adverse_reactions(text)

        assert result is not None
        assert result.discontinuation_rate == 13.0

    def test_parse_common_toxicities(self):
        """Should extract common toxicities with percentages."""
        text = """
        The most common adverse reactions were diarrhea 77%, nausea 50%,
        fatigue 40%, rash 35%, and vomiting 30%.
        """
        result = parse_adverse_reactions(text)

        assert result is not None
        assert len(result.common_toxicities) > 0

        # Find diarrhea in results
        tox_dict = {name.lower(): pct for name, pct in result.common_toxicities}
        assert tox_dict.get("diarrhea") == 77.0
        assert tox_dict.get("nausea") == 50.0

    def test_filters_low_percentage_toxicities(self):
        """Should filter out toxicities with less than 10%."""
        text = "Common toxicities were diarrhea 50%, headache 5%."
        result = parse_adverse_reactions(text)

        assert result is not None
        tox_names = [name.lower() for name, _ in result.common_toxicities]
        assert "diarrhea" in tox_names
        assert "headache" not in tox_names

    def test_limits_to_top_10(self):
        """Should limit to top 10 toxicities."""
        text = """
        Toxicities: diarrhea 80%, nausea 70%, fatigue 60%, rash 55%, vomiting 50%,
        stomatitis 45%, hyperglycemia 40%, anemia 35%, constipation 30%, headache 25%,
        arthralgia 22%, myalgia 20%.
        """
        result = parse_adverse_reactions(text)

        assert result is not None
        assert len(result.common_toxicities) <= 10

    def test_parse_full_adverse_reactions_text(self):
        """Test parsing realistic adverse reactions text."""
        text = """
        Serious adverse reactions occurred in 35% of patients receiving TRUQAP.
        Permanent discontinuation due to adverse reactions occurred in 13% of patients.
        The most common adverse reactions (≥10%) were diarrhea 77%, nausea 50%,
        fatigue 40%, rash 35%, and hyperglycemia 30%.
        """
        result = parse_adverse_reactions(text)

        assert result is not None
        assert result.serious_rate == 35.0
        assert result.discontinuation_rate == 13.0
        assert len(result.common_toxicities) >= 5

    def test_no_meaningful_data_returns_none(self):
        """Text without meaningful adverse reaction data should return None."""
        text = "This is general text without adverse reaction information."
        result = parse_adverse_reactions(text)

        assert result is None


class TestParseRecentMajorChanges:
    """Tests for parse_recent_major_changes function."""

    def test_parse_none_returns_none(self):
        """Empty or None input should return None."""
        assert parse_recent_major_changes(None) is None
        assert parse_recent_major_changes("") is None

    def test_parse_date(self):
        """Should extract label update date."""
        text = "Warnings and Precautions, Hyperglycemia (5.1) 02/2025"
        result = parse_recent_major_changes(text)

        assert result is not None
        assert result.last_label_update == "02/2025"

    def test_parse_reason(self):
        """Should extract update reason."""
        text = "Warnings and Precautions, Hyperglycemia (5.1) 02/2025"
        result = parse_recent_major_changes(text)

        assert result is not None
        assert "Warnings" in result.update_reason

    def test_parse_single_digit_month(self):
        """Should parse dates with single-digit month."""
        text = "Dosage and Administration (2.1) 1/2025"
        result = parse_recent_major_changes(text)

        assert result is not None
        assert result.last_label_update == "1/2025"

    def test_no_meaningful_data_returns_none(self):
        """Text without date or reason should return None."""
        text = "Some general text without date information."
        result = parse_recent_major_changes(text)

        assert result is None


class TestClinicalStudyDataModel:
    """Tests for ClinicalStudyData dataclass."""

    def test_default_values(self):
        """All fields should default to None."""
        data = ClinicalStudyData()

        assert data.trial_name is None
        assert data.nct_id is None
        assert data.patients_n is None
        assert data.pfs_months_treatment is None
        assert data.pfs_months_control is None
        assert data.hazard_ratio is None
        assert data.hazard_ratio_ci is None
        assert data.orr_treatment is None
        assert data.orr_control is None
        assert data.biomarker_breakdown is None

    def test_set_values(self):
        """Should correctly store all values."""
        data = ClinicalStudyData(
            trial_name="CAPItello-291",
            nct_id="NCT04305496",
            patients_n=289,
            hazard_ratio=0.50,
            hazard_ratio_ci=(0.38, 0.65),
            orr_treatment=22.0,
            orr_control=12.0,
            biomarker_breakdown={"PIK3CA": 76.0, "AKT1": 13.0},
        )

        assert data.trial_name == "CAPItello-291"
        assert data.nct_id == "NCT04305496"
        assert data.patients_n == 289
        assert data.hazard_ratio == 0.50
        assert data.hazard_ratio_ci == (0.38, 0.65)
        assert data.orr_treatment == 22.0
        assert data.orr_control == 12.0
        assert data.biomarker_breakdown["PIK3CA"] == 76.0


class TestMechanismOfActionDataModel:
    """Tests for MechanismOfActionData dataclass."""

    def test_default_values(self):
        """Default values should be correct."""
        data = MechanismOfActionData()

        assert data.targets == []
        assert data.mechanism is None
        assert data.preclinical is None

    def test_set_values(self):
        """Should correctly store all values."""
        data = MechanismOfActionData(
            targets=["AKT1", "AKT2", "AKT3"],
            mechanism="Inhibits phosphorylation of downstream substrates",
            preclinical="In vitro, reduced cell proliferation",
        )

        assert data.targets == ["AKT1", "AKT2", "AKT3"]
        assert "phosphorylation" in data.mechanism
        assert "In vitro" in data.preclinical


class TestAdverseReactionsDataModel:
    """Tests for AdverseReactionsData dataclass."""

    def test_default_values(self):
        """Default values should be correct."""
        data = AdverseReactionsData()

        assert data.common_toxicities == []
        assert data.serious_rate is None
        assert data.discontinuation_rate is None

    def test_set_values(self):
        """Should correctly store all values."""
        data = AdverseReactionsData(
            common_toxicities=[("Diarrhea", 77.0), ("Nausea", 50.0)],
            serious_rate=35.0,
            discontinuation_rate=13.0,
        )

        assert len(data.common_toxicities) == 2
        assert data.common_toxicities[0] == ("Diarrhea", 77.0)
        assert data.serious_rate == 35.0
        assert data.discontinuation_rate == 13.0


class TestRecentMajorChangesDataModel:
    """Tests for RecentMajorChangesData dataclass."""

    def test_default_values(self):
        """All fields should default to None."""
        data = RecentMajorChangesData()

        assert data.last_label_update is None
        assert data.update_reason is None

    def test_set_values(self):
        """Should correctly store all values."""
        data = RecentMajorChangesData(
            last_label_update="02/2025",
            update_reason="Warnings and Precautions: Hyperglycemia",
        )

        assert data.last_label_update == "02/2025"
        assert "Hyperglycemia" in data.update_reason

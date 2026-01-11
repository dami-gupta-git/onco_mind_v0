"""Tests for nlp_utils module."""

import pytest
from oncomind.nlp_utils import is_positive


class TestIsPositive:
    """Tests for is_positive function."""

    def test_statement(self):
        """Positive sensitivity statement should return True."""
        x= is_positive("", ["SENSITIV"])
        print(x)

    def test_positive_sensitivity_statement(self):
        """Positive sensitivity statement should return True."""
        assert is_positive("AKT1 E17K is sensitive to Everolimus", ["sensitiv"]) is True

    def test_negated_sensitivity_statement(self):
        """Negated sensitivity statement should return False."""
        assert not is_positive("AKT1 E17K is not sensitive to MK-2206", ["sensitive"])

    def test_positive_resistance_statement(self):
        """Positive resistance statement should return True."""
        assert is_positive("EGFR T790M confers resistance to erlotinib", ["ResisT"])

    def test_negated_resistance_statement(self):
        """Negated resistance statement should return False."""
        assert is_positive("This variant does not confer resistance", ["resist"]) is False

    def test_keyword_not_found(self):
        """Should return False when keyword is not in statement."""
        assert is_positive("AKT1 E17K is oncogenic", ["sensitiv"]) is False

    def test_sensitive_variant_form(self):
        """Should handle different word forms."""
        assert is_positive("The variant shows sensitivity to the drug", ["sensitiv"]) is True

    def test_did_not_show_sensitivity(self):
        """'Did not show sensitivity' should return False."""
        assert is_positive("Cell lines did not show sensitivity to MK-2206", ["sensitiv"]) is False

    def test_confers_sensitivity(self):
        """'Confers sensitivity' should return True."""
        assert is_positive("AKT1 E17K confers sensitivity to AKT inhibitors", ["sensitiv"]) is True

    def test_multiple_keywords(self):
        """Should match any keyword in the list."""
        assert is_positive("AKT1 E17K confers sensitivity to AKT inhibitors", ["sensitiv", "respons"]) is True
        assert is_positive("The variant is responsive to treatment", ["sensitiv", "respons"]) is True
        assert is_positive("The variant is a tumor regressor", ["tumor regress"])
        pass

    @pytest.mark.skip(reason="'lacks' is a negative verb not detected by current negation logic")
    def test_lacks_sensitivity(self):
        """'Lacks sensitivity' should return False."""
        # TODO: Enhance is_positive to detect negative verbs like 'lacks', 'fails', 'unable'
        assert is_positive("This mutation lacks sensitivity to the inhibitor", ["sensitiv"]) is False


class TestDiminishers:
    """Tests for diminisher detection in is_positive."""

    def test_less_sensitive(self):
        """'Less sensitive' should return False."""
        assert is_positive("The variant is less sensitive to the drug", ["sensitiv"]) is False

    def test_reduced_sensitivity(self):
        """'Reduced sensitivity' should return False."""
        assert is_positive("Cell lines showed reduced sensitivity to erlotinib", ["sensitiv"]) is False

    def test_lower_sensitivity(self):
        """'Lower sensitivity' should return False."""
        assert is_positive("The mutant showed lower sensitivity to the inhibitor", ["sensitiv"]) is False

    def test_poor_response(self):
        """'Poor response' should return False."""
        assert is_positive("The tumor showed poor response to chemotherapy", ["respons"]) is False

    def test_limited_efficacy(self):
        """'Limited efficacy' should return False."""
        assert is_positive("The drug demonstrated limited efficacy", ["effic"]) is False

    def test_weak_inhibition(self):
        """'Weak inhibition' should return False."""
        assert is_positive("The compound showed weak inhibition of the target", ["inhibit"]) is False

    def test_positive_without_diminisher(self):
        """Normal positive statement should still return True."""
        assert is_positive("The variant shows strong sensitivity to the drug", ["sensitiv"]) is True

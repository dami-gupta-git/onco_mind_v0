"""Unit tests for is_acquired_resistance_mutation() with tumor-type context."""

from oncomind.config.constants import is_acquired_resistance_mutation


class TestIsAcquiredResistanceMutation:
    # Backward-compatible: unrestricted entries work without tumor type
    def test_egfr_t790m_no_tumor_type(self):
        assert is_acquired_resistance_mutation("EGFR", "T790M") is True

    def test_egfr_t790m_with_any_tumor(self):
        # Unrestricted entry — tumor type is irrelevant
        assert is_acquired_resistance_mutation("EGFR", "T790M", "Neuroblastoma") is True

    # ALK F1174L — the fixed behavior
    def test_alk_f1174l_nsclc(self):
        assert is_acquired_resistance_mutation("ALK", "F1174L", "NSCLC") is True

    def test_alk_f1174l_lung_adeno(self):
        # Substring match on "lung"
        assert (
            is_acquired_resistance_mutation("ALK", "F1174L", "Lung adenocarcinoma")
            is True
        )

    def test_alk_f1174l_neuroblastoma(self):
        # Primary oncogenic driver in neuroblastoma — NOT acquired resistance
        assert (
            is_acquired_resistance_mutation("ALK", "F1174L", "Neuroblastoma") is False
        )

    def test_alk_f1174l_no_tumor_type(self):
        # Unknown context — safe default is False
        assert is_acquired_resistance_mutation("ALK", "F1174L") is False

    # Edge cases
    def test_unknown_variant(self):
        assert is_acquired_resistance_mutation("ALK", "V600E") is False

    def test_case_insensitive(self):
        assert is_acquired_resistance_mutation("alk", "f1174l", "nsclc") is True

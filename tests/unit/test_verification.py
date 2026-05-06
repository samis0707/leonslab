"""Tests for hard verification gates (D7.4)."""
from __future__ import annotations

from primer_design.verification import count_recognition_sites


class TestRecognitionCounting:
    """Direct unit tests for the public recognition-site counter. The full
    verification path is exercised in tests/integration/."""

    def test_palindromic_motif_not_double_counted(self):
        # HindIII motif AAGCTT is a palindrome (rc(AAGCTT) == AAGCTT). One
        # physical site must count as 1, not 2.
        assert count_recognition_sites("GGGAAGCTTGGG", "HindIII") == 1

        # Two non-overlapping sites → 2.
        assert count_recognition_sites("AAGCTT" + "AAAAAA" + "AAGCTT", "HindIII") == 2

        # No site → 0.
        assert count_recognition_sites("GGGGGGGG", "HindIII") == 0

        # Same property for other palindromes used by v1: ClaI (ATCGAT),
        # EcoRI (GAATTC), ApaI (GGGCCC), KpnI (GGTACC).
        for motif, enzyme in [("ATCGAT", "ClaI"), ("GAATTC", "EcoRI"),
                              ("GGGCCC", "ApaI"), ("GGTACC", "KpnI")]:
            assert count_recognition_sites("AAA" + motif + "AAA", enzyme) == 1

    def test_zero_count_when_destroyed(self):
        # The site_destroyed convention guarantees zero recognition sites
        # after assembly; check the counter on a representative motif-free
        # sequence on both strands.
        seq = "ACGTACGTACGT" * 12
        assert count_recognition_sites(seq, "HindIII") == 0
        assert count_recognition_sites(seq, "ClaI") == 0
        assert count_recognition_sites(seq, "EcoRI") == 0

    def test_one_count_when_partial(self):
        # The site_partial_AAGCT convention regenerates exactly one HindIII
        # site at one assembly junction.
        assert count_recognition_sites("TTT" + "AAGCTT" + "AAA", "HindIII") == 1

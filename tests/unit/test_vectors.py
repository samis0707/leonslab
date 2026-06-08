"""Vector tests: cut detection and tail-rule calibration.

The tail-derivation calibration tests are the FIRST checkpoint of Phase 2.
They will fail until vectors are added to data/primer_design/vectors/.

When .gb files are present, these MUST produce the empirically correct tails:
    pEXG2 + HindIII + deletion  → P1 tail = "AATGTAAAGCAAGCT"
    pBBR1MCS2 + HindIII + expression → P1 tail = "CGGTATCGATAAGCT"
"""
from __future__ import annotations

from pathlib import Path

import pytest

from primer_design.exceptions import NotSingleCutter, UnknownEnzyme
from primer_design.vectors import (
    derive_tails,
    find_cut_position,
    get_unique_cutters,
    get_vector,
    reverse_complement,
)


@pytest.fixture
def vectors_present() -> bool:
    """Skip vector tests if vector files haven't been added yet (Phase 2 step 2)."""
    vectors_dir = Path(__file__).resolve().parents[2] / "data" / "primer_design" / "vectors"
    has_pexg2 = (vectors_dir / "pEXG2.gb").exists() or (vectors_dir / "pEXG2.fasta").exists()
    has_pbbr = (vectors_dir / "pBBR1MCS2.gb").exists() or (vectors_dir / "pBBR1MCS2.fasta").exists()
    return has_pexg2 and has_pbbr


class TestReverseComplement:

    def test_acgt(self):
        assert reverse_complement("ACGT") == "ACGT"

    def test_atcg(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        assert reverse_complement("AAGCTT") == "AAGCTT"


class TestCutDetection:

    def test_unknown_enzyme(self, vectors_present):
        if not vectors_present:
            pytest.skip("vectors not yet added")
        vector = get_vector("pEXG2")
        with pytest.raises(UnknownEnzyme):
            find_cut_position(vector, "FakeEnzyme")

    def test_pEXG2_HindIII_is_single_cutter(self, vectors_present):
        if not vectors_present:
            pytest.skip("vectors not yet added")
        vector = get_vector("pEXG2")
        nick = find_cut_position(vector, "HindIII")
        assert nick > 0


class TestTailCalibration:
    """These four tests are the empirical anchor for the entire pipeline."""

    def test_pEXG2_HindIII_deletion_p1_tail(self, vectors_present):
        if not vectors_present:
            pytest.skip("vectors not yet added")
        vector = get_vector("pEXG2")
        p1_tail, _ = derive_tails(vector, "HindIII", "deletion")
        assert p1_tail == "AATGTAAAGCAAGCT", \
            f"P1 tail mismatch: got {p1_tail!r}, expected 'AATGTAAAGCAAGCT'"

    def test_pEXG2_HindIII_deletion_p4_tail(self, vectors_present):
        if not vectors_present:
            pytest.skip("vectors not yet added")
        vector = get_vector("pEXG2")
        _, p4_tail = derive_tails(vector, "HindIII", "deletion")
        # From v1 working primer P4: tail = "CGACCTGCAGAAGCT"
        assert p4_tail == "CGACCTGCAGAAGCT", \
            f"P4 tail mismatch: got {p4_tail!r}, expected 'CGACCTGCAGAAGCT'"

    def test_pBBR1MCS2_HindIII_expression_p1_tail(self, vectors_present):
        if not vectors_present:
            pytest.skip("vectors not yet added")
        vector = get_vector("pBBR1MCS2")
        p1_tail, _ = derive_tails(vector, "HindIII", "expression")
        assert p1_tail == "CGGTATCGATAAGCT", \
            f"P1 tail mismatch: got {p1_tail!r}, expected 'CGGTATCGATAAGCT'"

    def test_pBBR1MCS2_HindIII_expression_p4_tail(self, vectors_present):
        if not vectors_present:
            pytest.skip("vectors not yet added")
        vector = get_vector("pBBR1MCS2")
        _, p4_tail = derive_tails(vector, "HindIII", "expression")
        assert p4_tail == "ATTCGATATCAAGCT", \
            f"P4 tail mismatch: got {p4_tail!r}, expected 'ATTCGATATCAAGCT'"


class TestUniqueCutters:

    def test_pEXG2_returns_list(self, vectors_present):
        if not vectors_present:
            pytest.skip("vectors not yet added")
        vector = get_vector("pEXG2")
        cutters = get_unique_cutters(vector)
        assert isinstance(cutters, list)
        assert "HindIII" in cutters

    def test_pBBR1MCS2_returns_13_cutters(self, vectors_present):
        if not vectors_present:
            pytest.skip("vectors not yet added")
        vector = get_vector("pBBR1MCS2")
        cutters = get_unique_cutters(vector)
        # User verified 13 single-cutters earlier in Phase 1
        assert len(cutters) == 13, f"Expected 13 unique cutters, got {len(cutters)}: {cutters}"

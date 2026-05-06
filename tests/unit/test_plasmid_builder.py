"""Unit tests for ``plasmid_builder.assemble``."""
from __future__ import annotations

import pytest

from primer_design.exceptions import PrimerDesignError
from primer_design.plasmid_builder import assemble
from primer_design.types import TailConvention, VectorRecord


def _vec(seq: str) -> VectorRecord:
    return VectorRecord(name="synth", sequence=seq, length=len(seq))


SITE_DESTROYED = TailConvention(
    name="site_destroyed",
    p1_tail_rule="left_arm_15nt_excl_recognition",
    p4_tail_rule="rc_right_arm_15nt_starting_at_second_nt_of_recognition",
    expected_recognition_count_in_final_plasmid=0,
    enzyme="HindIII", vector="pEXG2", application="deletion",
)


class TestAssembleSiteDestroyed:

    def test_length_identity(self):
        # vector: 100 nt, nick at 50, insert: 60 nt → final = 100 + 60 - 30 = 130
        vec = _vec("A" * 35 + "GATCAGGGCTAGCAA" + "T" * 50)
        nick = 50
        insert = "GATCAGGGCTAGCAA" + "C" * 30 + "TTTTTTTTTTTTTTT"  # left+body+right
        assert insert[:15] == vec.sequence[nick - 15 : nick]
        assert insert[-15:] == vec.sequence[nick : nick + 15]
        result = assemble(vec, insert, nick, SITE_DESTROYED)
        assert result.length == len(vec.sequence) + len(insert) - 30
        # Body of insert (between collapsed overlaps) sits at the nick
        assert result.sequence[:nick] == vec.sequence[:nick]
        assert result.sequence[nick : nick + 30] == "C" * 30
        assert result.sequence[nick + 30 :] == vec.sequence[nick:]

    def test_junction_positions(self):
        vec = _vec("ACGT" * 25 + "X" * 15 + "Y" * 40)  # nick=100; left tail 15 = "ACGT...ACGT"[-15:]
        # ensure left arm end matches insert head
        nick = 100
        left_tail = vec.sequence[nick - 15 : nick]
        right_tail = vec.sequence[nick : nick + 15]
        insert = left_tail + "Z" * 20 + right_tail
        result = assemble(vec, insert, nick, SITE_DESTROYED)
        assert result.junctions == [(nick, nick + 20)]

    def test_too_short_insert_raises(self):
        vec = _vec("N" * 100)
        with pytest.raises(PrimerDesignError):
            assemble(vec, "ACGT", 50, SITE_DESTROYED)

    def test_nick_too_close_to_end_raises(self):
        vec = _vec("N" * 100)
        with pytest.raises(PrimerDesignError):
            assemble(vec, "X" * 40, 5, SITE_DESTROYED)

    def test_unknown_convention_raises(self):
        vec = _vec("N" * 100)
        bad = TailConvention(
            name="some_other", p1_tail_rule="x", p4_tail_rule="y",
            expected_recognition_count_in_final_plasmid=0,
            enzyme="HindIII", vector="pEXG2", application="deletion",
        )
        with pytest.raises(NotImplementedError):
            assemble(vec, "X" * 40, 50, bad)

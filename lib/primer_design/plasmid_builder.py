"""Assemble the final circular plasmid by replacing the cut region of the vector
with the assembled insert.

For ``site_destroyed`` conventions: vector arms join directly to insert via the
15-nt In-Fusion overlap (which collapses), so:
    final = vector_left_arm + insert_minus_overlap_left + insert_body + insert_minus_overlap_right + vector_right_arm

Concretely we place the insert at the cut nick, and rely on the upstream tail
derivation to have made the junction sequences identical. Verification (§D7.4)
confirms this with strict equality checks.
"""
from __future__ import annotations

from .config import VECTOR_TAIL_LEN
from .types import CircularPlasmid, TailConvention, VectorRecord


def assemble(
    vector: VectorRecord,
    insert: str,
    cut_nick: int,
    convention: TailConvention,
) -> CircularPlasmid:
    """Build the linear representation of the final circular plasmid.

    Note on circular topology: pEXG2 / pBBR1MCS2 are circular. We carry a linear
    string with the cut at the start/end and stitch by string concatenation.
    The final length should equal ``len(vector) + len(insert) - 2 * VECTOR_TAIL_LEN``
    when the In-Fusion overlap collapses on both sides.

    Args:
        vector: parsed vector.
        insert: full insert from In-Fusion assembly. For deletion/tagging this is
            UP+DN with junction collapsed (length = len(UP) + len(DN) - 30).
            For expression it is the single PCR amplicon.
        cut_nick: 0-based top-strand position of the cut.
        convention: needed to know how many nt of the insert overlap with each arm.

    Returns:
        CircularPlasmid.
    """
    if len(insert) < 2 * VECTOR_TAIL_LEN:
        raise ValueError(
            f"Insert too short ({len(insert)} bp) for two {VECTOR_TAIL_LEN}-nt arms"
        )

    vec_seq = vector.sequence
    n = vector.length
    left_arm = vec_seq[(cut_nick - VECTOR_TAIL_LEN) % n : cut_nick]
    right_arm = vec_seq[cut_nick : cut_nick + VECTOR_TAIL_LEN]

    if insert[:VECTOR_TAIL_LEN] != left_arm:
        raise ValueError(
            "Insert 5' arm does not match vector left arm at cut site "
            f"(insert: {insert[:VECTOR_TAIL_LEN]!r}, vector: {left_arm!r})"
        )
    if insert[-VECTOR_TAIL_LEN:] != right_arm:
        raise ValueError(
            "Insert 3' arm does not match vector right arm at cut site "
            f"(insert: {insert[-VECTOR_TAIL_LEN:]!r}, vector: {right_arm!r})"
        )

    # Collapse both 15-nt overlaps: keep vector through the nick, splice in the
    # interior of the insert (= insert minus its two arm copies), keep the rest
    # of the vector. Equivalent to vec[:nick-15] + insert + vec[nick+15:].
    final_seq = (
        vec_seq[: cut_nick - VECTOR_TAIL_LEN]
        + insert
        + vec_seq[cut_nick + VECTOR_TAIL_LEN :]
    )

    insert_start = cut_nick - VECTOR_TAIL_LEN
    insert_end = insert_start + len(insert)
    return CircularPlasmid(
        name=f"{vector.name}_assembly",
        sequence=final_seq,
        length=len(final_seq),
        junctions=[(insert_start, insert_start + VECTOR_TAIL_LEN),
                   (insert_end - VECTOR_TAIL_LEN, insert_end)],
    )

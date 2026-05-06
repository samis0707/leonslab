"""Assemble the final circular plasmid by replacing the cut region of the vector
with the assembled insert.

The In-Fusion / Gibson overlap is 15 nt at each junction. The exact position of
those 15-nt windows in the vector depends on the tail rule chosen for that
(vector, enzyme, application) triple. For example:

    site_destroyed   (pEXG2|HindIII|deletion)   left arm = vec[nick-15:nick]
                                                right arm = vec[nick:nick+15]
    site_partial_AAGCT (pBBR1MCS2|HindIII|expression) left arm = vec[nick-11:nick+4]
                                                       right arm = vec[nick:nick+15]

Rather than hard-code the offsets per convention, we locate each arm by literal
match against the vector sequence. The insert's first / last 15 nt are required
to occur exactly once in the vector — which holds because the tail derivation
copies from the vector deterministically. The two vector positions then define
the splice; any nucleotide that lies in BOTH the left-arm range and the
right-arm range (e.g. the 4 nt of AAGCT shared by site_partial_AAGCT) is
collapsed once via the insert.
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

    Algorithm:
        1. Locate insert[:15] in vector.sequence  → left_arm at [L, L+15).
        2. Locate insert[-15:] in vector.sequence → right_arm at [R, R+15).
        3. Splice: final = vec[:L] + insert + vec[R+15:], assuming the vector
           is linear over [0, vec_len) and the cut lies inside [L, R+15).
           Equivalent to keeping vector outside the [L .. R+15) span and
           letting the insert fill that span (the insert already contains
           identical copies of the two 15-nt overlap windows, so they collapse
           naturally).

    Length: len(vec) - (R + 15 - L) + len(insert).
    """
    if len(insert) < 2 * VECTOR_TAIL_LEN:
        raise ValueError(
            f"Insert too short ({len(insert)} bp) for two {VECTOR_TAIL_LEN}-nt arms"
        )

    vec_seq = vector.sequence
    left_arm_seq = insert[:VECTOR_TAIL_LEN]
    right_arm_seq = insert[-VECTOR_TAIL_LEN:]

    L = _unique_index(vec_seq, left_arm_seq)
    if L < 0:
        raise ValueError(
            f"Insert 5' arm {left_arm_seq!r} not found uniquely in vector "
            f"{vector.name!r} (count={vec_seq.count(left_arm_seq)})"
        )
    R = _unique_index(vec_seq, right_arm_seq)
    if R < 0:
        raise ValueError(
            f"Insert 3' arm {right_arm_seq!r} not found uniquely in vector "
            f"{vector.name!r} (count={vec_seq.count(right_arm_seq)})"
        )

    # Sanity: the cut nick should lie within the spliced range [L .. R+15).
    if not (L <= cut_nick <= R + VECTOR_TAIL_LEN):
        raise ValueError(
            f"Cut nick {cut_nick} outside the insert overlap span "
            f"[{L}, {R + VECTOR_TAIL_LEN}) in {vector.name}"
        )
    if R < L:
        raise ValueError(
            f"Right arm precedes left arm in {vector.name}: L={L}, R={R}"
        )

    final_seq = vec_seq[:L] + insert + vec_seq[R + VECTOR_TAIL_LEN :]

    insert_start = L
    insert_end = insert_start + len(insert)
    return CircularPlasmid(
        name=f"{vector.name}_assembly",
        sequence=final_seq,
        length=len(final_seq),
        junctions=[(insert_start, insert_start + VECTOR_TAIL_LEN),
                   (insert_end - VECTOR_TAIL_LEN, insert_end)],
    )


def _unique_index(haystack: str, needle: str) -> int:
    """Return the unique 0-based position of ``needle`` in ``haystack``, or -1
    if it occurs zero or multiple times."""
    first = haystack.find(needle)
    if first < 0:
        return -1
    if haystack.find(needle, first + 1) >= 0:
        return -1
    return first

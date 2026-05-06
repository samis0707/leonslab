"""Assemble the final circular plasmid by replacing the cut region of the vector
with the assembled insert.

For a ``site_destroyed`` convention (deletion / tagging on pEXG2):
    P1.tail = vector[nick - 15 : nick]                   (forward, left arm)
    P4.tail = rc(vector[nick : nick + 15])               (RC of right arm)

The insert produced by In-Fusion carries those same 15 nt at each end (the P1
tail at the 5' end and rc(P4.tail) = vector[nick:nick+15] at the 3' end). When
the overlap collapses on both sides:

    final = vector[:nick] + insert[15:-15] + vector[nick:]

Length identity: ``len(final) = len(vector) + len(insert) - 30``.

For the ``site_partial_AAGCT`` convention (pBBR1MCS2 expression) the geometry is
asymmetric — the P1 tail straddles the recognition site rather than ending at the
nick. That case is implemented in step 8.
"""
from __future__ import annotations

from .config import VECTOR_TAIL_LEN
from .exceptions import PrimerDesignError
from .types import CircularPlasmid, TailConvention, VectorRecord


def assemble(
    vector: VectorRecord,
    insert: str,
    cut_nick: int,
    convention: TailConvention,
) -> CircularPlasmid:
    """Build the linear representation of the final circular plasmid.

    Args:
        vector: parsed vector (linear representation of the circular plasmid).
        insert: full insert from In-Fusion assembly. For deletion/tagging this is
            UP+DN with the 30-nt junction collapsed; for expression it is the
            single PCR amplicon.
        cut_nick: 0-based top-strand position of the restriction nick.
        convention: tells us which tail-derivation rule was used for P1/P4, which
            determines how the insert ends overlap the vector arms.

    Returns:
        CircularPlasmid in linear representation, with the 0-based positions of
        the two insert/vector junctions in ``junctions``.
    """
    overlap = VECTOR_TAIL_LEN

    if convention.name == "site_destroyed":
        # Symmetric 15-nt collapse on both sides.
        if not (overlap <= cut_nick <= len(vector.sequence) - overlap):
            raise PrimerDesignError(
                message=f"cut_nick {cut_nick} too close to vector ends for {overlap}-nt overlap",
                details={"cut_nick": cut_nick, "vector_len": len(vector.sequence)},
            )
        if len(insert) < 2 * overlap:
            raise PrimerDesignError(
                message=f"insert length {len(insert)} < 2*overlap ({2 * overlap})",
                details={"insert_len": len(insert)},
            )
        left_arm = vector.sequence[:cut_nick]
        right_arm = vector.sequence[cut_nick:]
        body = insert[overlap:-overlap]
        seq = left_arm + body + right_arm
        return CircularPlasmid(
            name=f"{vector.name}_assembly",
            sequence=seq,
            length=len(seq),
            junctions=[(cut_nick, cut_nick + len(body))],
        )

    raise NotImplementedError(
        f"plasmid assembly for convention {convention.name!r} is not implemented yet "
        f"(only 'site_destroyed' is supported in step 7)."
    )

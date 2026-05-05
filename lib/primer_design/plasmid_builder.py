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
    raise NotImplementedError(
        "Phase 2 step 10: implement plasmid assembly. "
        "Hardest part: getting the offsets right at both junctions. "
        "Validate against tests/integration/test_deletion_LB001_lasB.py expected "
        "final_plasmid_length_bp = 6155 and SHA-256 against pEXG2-lasB_delta_LB001.fasta."
    )

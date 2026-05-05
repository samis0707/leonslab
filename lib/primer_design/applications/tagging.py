"""Orchestrate the pEXG2 in-locus C-terminal tagging pipeline.

Per skill_v2 §5 and project_plan §4.2.
"""
from __future__ import annotations

from .. import gene_finder, primer_designer, storage_adapter, vectors, verification
from ..plasmid_builder import assemble
from ..tags import build_in_locus_cassette, get_tag, validate_tag_position
from ..types import DesignRequest, DesignResult


def run(request: DesignRequest) -> DesignResult:
    """End-to-end in-locus C-terminal tagging design.

    Steps (additional vs deletion):
        - Validate tag is in_locus_ok (≤36 nt cassette); raise TagTooLongForInLocus otherwise.
        - Build cassette (LINKER + tag DNA + new stop).
        - Junction overlap = len(cassette); split between P2 and P3 tails.

    Note:
        UP fragment ends BEFORE native stop codon (cds[:-3]); DN fragment starts
        AT native stop position (dn_flank[0:]). The cassette occupies the junction.
    """
    if request.tag is None:
        raise ValueError("tagging application requires a tag")
    tag = get_tag(request.tag)
    validate_tag_position(tag, request.tag_position or "C")
    cassette = build_in_locus_cassette(tag)   # raises TagTooLongForInLocus if needed

    raise NotImplementedError(
        "Phase 2 step 9: orchestrate per skill_v2 §5. "
        "Test target: tests/integration/test_tagging_PA14_lasR_His6.py"
    )

"""Orchestrate the pBBR1MCS2 plasmid-expression pipeline.

Per skill_v2 §6 and project_plan §4.3.
"""
from __future__ import annotations

from .. import gene_finder, primer_designer, storage_adapter, vectors, verification
from ..plasmid_builder import assemble
from ..tags import build_plasmid_fusion_cds, get_tag, validate_tag_position
from ..types import DesignRequest, DesignResult


def run(request: DesignRequest) -> DesignResult:
    """End-to-end plasmid-expression design.

    Steps:
        - Optional tag → build_plasmid_fusion_cds to construct the full ATG-to-stop
          coding sequence (with N-term or C-term tag).
        - 2-primer search:
            P1 = vector_p1_tail + AGGAGG + spacer + body_at_ATG
            P2 = vector_p2_tail + body_at_stop_RC
        - Single PCR amplicon, 2-fragment In-Fusion.
    """
    tag = None
    if request.tag is not None:
        tag = get_tag(request.tag)
        if request.tag_position is None:
            raise ValueError("tag_position required when tag is set")
        validate_tag_position(tag, request.tag_position)

    raise NotImplementedError(
        "Phase 2 step 8: orchestrate per skill_v2 §6. "
        "Test target: tests/integration/test_expression_PA14_lasR.py"
    )

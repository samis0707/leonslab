"""Orchestrate the pBBR1MCS2 plasmid-expression pipeline.

Per skill_v2 §6 and project_plan §4.3.
"""
from __future__ import annotations

from typing import Callable

from .. import gene_finder, primer_designer, vectors, verification
from ..config import OFFTARGET_PRODUCT_SIZE_MAX_BP
from ..exceptions import OffTargetDetected
from ..plasmid_builder import assemble
from ..tags import get_tag, validate_tag_position
from ..types import DesignRequest, DesignResult
from ..vectors import reverse_complement
from .deletion import _default_genome_loader


GenomeLoader = Callable[[str], bytes]


def run(
    request: DesignRequest,
    *,
    genome_loader: GenomeLoader | None = None,
) -> DesignResult:
    """End-to-end plasmid-expression design (skill_v2 §6).

    Steps:
        1. Load gene record.
        2. Resolve optional tag (validate position).
        3. Load vector + expression tail convention; locate cut nick.
        4. 2-primer search: P1 anchored at ATG, P2 anchored at stop.
        5. Build single PCR amplicon = P1.tail + coding_seq + RC(P2.tail).
        6. Off-target scan.
        7. Assemble plasmid.
        8. Verify.
    """
    gene = gene_finder.get_gene_record(request.isolate_id, request.gene)
    vector = vectors.get_vector(request.vector)
    convention = vectors.get_tail_convention(
        request.vector, request.enzyme, "expression"
    )
    cut_nick = vectors.find_cut_position(vector, request.enzyme)

    tag = None
    if request.tag is not None:
        tag = get_tag(request.tag)
        if request.tag_position is None:
            raise ValueError("tag_position required when tag is set")
        validate_tag_position(tag, request.tag_position)

    best, _alts = primer_designer.search_expression_primers(
        gene, vector, convention, tag, request.tag_position
    )

    # Full PCR amplicon (used for plasmid assembly): includes both vector
    # homology arms.
    amplicon = (
        best.p1.tail + best.coding_seq + reverse_complement(best.p2.tail)
    )
    # Reported "insert": the novel sequence carried into the plasmid, i.e.
    # everything except the left vector homology arm (which duplicates the
    # vector's left arm and gets collapsed at assembly). Equivalent to
    # RBS+spacer + coding_seq + right_vector_arm_RC.
    from ..config import VECTOR_TAIL_LEN
    insert = amplicon[VECTOR_TAIL_LEN:]

    loader = genome_loader or _default_genome_loader
    genome_bytes = loader(request.isolate_id)
    expected_products = _expected_products_for_expression(amplicon)
    off_target = primer_designer.off_target_scan(
        [best.p1, best.p2], genome_bytes, expected_products,
    )
    if not off_target.passed:
        raise OffTargetDetected(
            message="Off-target scan reports unintended PCR products",
            details={"violations": [vars(v) for v in off_target.violations]},
        )

    final_plasmid = assemble(vector, amplicon, cut_nick, convention)
    final_plasmid.name = f"{vector.name}_{gene.gene}_{gene.isolate_id}_expr"

    result = DesignResult(
        request=request,
        gene_record=gene,
        vector=vector,
        convention=convention,
        primer_set=best,
        up_amplicon=None,
        dn_amplicon=None,
        insert=insert,
        final_plasmid=final_plasmid,
        off_target=off_target,
    )
    verification.verify(result)
    return result


def _expected_products_for_expression(insert: str) -> dict:
    size = len(insert)
    tol = 150
    window = (
        max(0, size - tol),
        min(OFFTARGET_PRODUCT_SIZE_MAX_BP, size + tol),
    )
    # off_target_scan iterates ordered (a, b) pairs and may report the same
    # genomic locus under both (P1, P2) and (P2, P1) when both primers happen
    # to have matches on both strands. Accept either ordering as expected.
    return {("P1", "P2"): window, ("P2", "P1"): window}

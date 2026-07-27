"""Orchestrate the pBBR1MCS2 plasmid-expression pipeline.

Per skill_v2 §6 and project_plan §4.3.
"""
from __future__ import annotations

from typing import Callable

from .. import cross_isolate_checker, fixed_primers, gene_finder, primer_designer, vectors, verification
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

    best, alts = primer_designer.search_expression_primers(
        gene, vector, convention, tag, request.tag_position
    )

    loader = genome_loader or _default_genome_loader
    genome_bytes = loader(request.isolate_id)

    # A pre-validated (wet-lab) primer set takes priority over a fresh design
    # whenever its genomic bodies match this isolate's sequence. Only applies
    # to plain (untagged) complementation — a requested tag has no room in a
    # fixed body/tail.
    fixed_primer_set = None
    if tag is None:
        fixed = fixed_primers.build_fixed_expression_set(gene, vector, convention, genome_bytes)
        if fixed is not None:
            gene = fixed.gene_record
            fixed_primer_set = fixed.primer_set

    candidates: list = []
    seen: set = set()
    for cand in [fixed_primer_set, best, *alts]:
        if cand is None:
            continue
        key = (cand.p1.body, cand.p2.body)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(cand)

    chosen = amplicon = chosen_off = None
    last_violations: list = []
    for cand in candidates:
        # Full PCR amplicon: P1.tail + native_anneal_segment + RC(P2.tail).
        # The template is native genomic DNA (no tag), so the tag-encoding
        # nucleotides ride along inside the primer overhangs; reconstituting
        # the amplicon from primers + native template gives the same final
        # insert structure the assembly is supposed to produce.
        amp = cand.p1.tail + cand.native_anneal_segment + reverse_complement(cand.p2.tail)
        expected_products = _expected_products_for_expression(amp)
        off_target = primer_designer.off_target_scan(
            [cand.p1, cand.p2], genome_bytes, expected_products,
        )
        if off_target.passed:
            chosen, amplicon, chosen_off = cand, amp, off_target
            break
        last_violations = off_target.violations

    if chosen is None:
        raise OffTargetDetected(
            message=(
                "Off-target scan reports unintended PCR products for the "
                f"{len(candidates)} top-scoring primer alternatives."
            ),
            details={
                "candidates_tried": len(candidates),
                "last_violations": [vars(v) for v in last_violations],
            },
        )

    best = chosen
    insert = amplicon
    off_target = chosen_off

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
    if fixed_primer_set is not None and chosen is fixed_primer_set:
        result.warnings.append(
            f"PRE-VALIDATED PRIMER SET: using the fixed, wet-lab-validated {gene.gene} "
            f"complementation primers (not freshly designed) — their genomic bodies "
            f"matched {gene.isolate_id} exactly."
        )
    result.compatible_isolates = cross_isolate_checker.find_compatible_isolates(
        best, request.gene, request.isolate_id
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

"""Orchestrate the pEXG2 in-locus C-terminal tagging pipeline.

Per skill_v2 §5 and project_plan §4.2.
"""
from __future__ import annotations

from typing import Callable

from .. import cross_isolate_checker, gene_finder, primer_designer, vectors, verification
from ..config import OFFTARGET_PRODUCT_SIZE_MAX_BP
from ..exceptions import OffTargetDetected
from ..plasmid_builder import assemble
from ..tags import build_in_locus_cassette, get_tag, validate_tag_position
from ..types import DesignRequest, DesignResult
from ..vectors import reverse_complement
from .deletion import _default_genome_loader


GenomeLoader = Callable[[str], bytes]


def run(
    request: DesignRequest,
    *,
    genome_loader: GenomeLoader | None = None,
) -> DesignResult:
    """End-to-end in-locus C-terminal tagging design (skill_v2 §5).

    Steps:
        1. Validate tag (in_locus_ok and C-term).
        2. Build cassette = GGS-linker + tag DNA + new stop.
        3. Load gene record + vector + tagging tail convention.
        4. Search primer set (full CDS minus native stop, no (N, C) loop).
        5. Build UP / DN amplicons; the cassette spans the junction with the
           central overlap_len nt shared between them.
        6. Off-target scan.
        7. Assemble plasmid.
        8. Verify.
    """
    if request.tag is None:
        raise ValueError("tagging application requires a tag")
    tag = get_tag(request.tag)
    validate_tag_position(tag, request.tag_position or "C")
    build_in_locus_cassette(tag)  # raises TagTooLongForInLocus before any work

    gene = gene_finder.get_gene_record(request.isolate_id, request.gene)
    vector = vectors.get_vector(request.vector)
    convention = vectors.get_tail_convention(
        request.vector, request.enzyme, "tagging"
    )
    cut_nick = vectors.find_cut_position(vector, request.enzyme)

    best, alts = primer_designer.search_tagging_primers(
        gene, vector, convention, tag
    )

    loader = genome_loader or _default_genome_loader
    genome_bytes = loader(request.isolate_id)

    candidates: list = [best]
    seen = {(best.p1.body, best.p2.body, best.p3.body, best.p4.body)}
    for alt in alts:
        key = (alt.p1.body, alt.p2.body, alt.p3.body, alt.p4.body)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(alt)

    chosen = chosen_amps = chosen_off = None
    last_violations: list = []
    for cand in candidates:
        up_amplicon, dn_amplicon, insert = _build_amplicons(gene, cand)
        expected = _expected_products_for_tagging(up_amplicon, dn_amplicon, gene, cand)
        off_target = primer_designer.off_target_scan(
            [cand.p1, cand.p2, cand.p3, cand.p4], genome_bytes, expected,
        )
        if off_target.passed:
            chosen, chosen_amps, chosen_off = cand, (up_amplicon, dn_amplicon, insert), off_target
            break
        last_violations = off_target.violations

    if chosen is None:
        raise OffTargetDetected(
            message=(
                "Off-target scan reports unintended PCR products for the "
                f"{len(candidates)} top-scoring primer alternatives. Manual "
                "primer design or anchor tuning required."
            ),
            details={
                "candidates_tried": len(candidates),
                "last_violations": [vars(v) for v in last_violations],
            },
        )

    best = chosen
    up_amplicon, dn_amplicon, insert = chosen_amps
    off_target = chosen_off

    final_plasmid = assemble(vector, insert, cut_nick, convention)
    final_plasmid.name = (
        f"{vector.name}_{gene.gene}_{gene.isolate_id}_{tag.name}_C"
    )

    result = DesignResult(
        request=request,
        gene_record=gene,
        vector=vector,
        convention=convention,
        primer_set=best,
        up_amplicon=up_amplicon,
        dn_amplicon=dn_amplicon,
        insert=insert,
        final_plasmid=final_plasmid,
        off_target=off_target,
    )
    result.compatible_isolates = cross_isolate_checker.find_compatible_isolates(
        best, request.gene, request.isolate_id
    )
    verification.verify(result)
    return result


def _build_amplicons(gene, primer_set):
    """Reconstruct UP / DN amplicons + the joined insert.

    UP top strand: P1.tail + up_segment[p1_anchor:] + cassette[:overlap_left]
    DN top strand: cassette[len(cassette)-overlap_right:] + dn_segment[:p4_end]
                    + RC(P4.tail)

    The shared region between UP's 3' end and DN's 5' end is
    cassette[overlap_start:overlap_end] of length ``overlap_len``. After
    collapsing that overlap once, the insert reads::

        P1.tail + up_segment[p1_anchor:] + cassette + dn_segment[:p4_end]
                + RC(P4.tail)
    """
    cds_no_stop = gene.cds_seq[:-3]
    up_segment = gene.up_flank + cds_no_stop
    dn_segment = gene.dn_flank

    p1_anchor = up_segment.find(primer_set.p1.body)
    if p1_anchor < 0:
        raise ValueError(
            f"P1 body {primer_set.p1.body!r} not found in up_segment "
            f"({gene.gene}/{gene.isolate_id})"
        )

    p4_body_rc = reverse_complement(primer_set.p4.body)
    p4_anchor = dn_segment.find(p4_body_rc)
    if p4_anchor < 0:
        raise ValueError(
            f"P4 body RC {p4_body_rc!r} not found in dn_segment "
            f"({gene.gene}/{gene.isolate_id})"
        )
    p4_end = p4_anchor + len(p4_body_rc)

    cassette = primer_set.cassette
    overlap_left = primer_set.overlap_left
    overlap_right = primer_set.overlap_right
    overlap_len = overlap_left + overlap_right - len(cassette)
    assert overlap_len > 0, "overlap_left + overlap_right must exceed cassette length"

    up_amplicon = (
        primer_set.p1.tail
        + up_segment[p1_anchor:]
        + cassette[:overlap_left]
    )
    dn_amplicon = (
        cassette[len(cassette) - overlap_right :]
        + dn_segment[:p4_end]
        + reverse_complement(primer_set.p4.tail)
    )
    insert = up_amplicon + dn_amplicon[overlap_len:]
    return up_amplicon, dn_amplicon, insert


def _expected_products_for_tagging(
    up_amplicon: str, dn_amplicon: str, gene, primer_set
) -> dict:
    up_size = len(up_amplicon)
    dn_size = len(dn_amplicon)
    overlap_len = (
        primer_set.overlap_left + primer_set.overlap_right
        - len(primer_set.cassette)
    )
    tol = 150
    # WT locus product (P1 + P4) on the same isolate genome: full gene
    # including the native stop, plus the same flanking spans the UP and DN
    # amplicons cover. = (UP without cassette portion) + native CDS-tail-stop +
    # (DN without cassette portion).
    wt_size = (
        up_size - primer_set.overlap_left
        + 3                                   # native stop codon
        + dn_size - primer_set.overlap_right
    )
    return {
        ("P1", "P2"): (max(0, up_size - tol), up_size + tol),
        ("P3", "P4"): (max(0, dn_size - tol), dn_size + tol),
        ("P1", "P4"): (max(0, wt_size - tol),
                       min(OFFTARGET_PRODUCT_SIZE_MAX_BP, wt_size + tol)),
    }

"""Orchestrate the pEXG2 in-frame deletion pipeline.

Per skill_v2 §4 and project_plan §4.1.
"""
from __future__ import annotations

from .. import gene_finder, primer_designer, storage_adapter, vectors, verification
from ..config import JUNCTION_LEN_PER_PRIMER_DEFAULT
from ..exceptions import OffTargetDetected
from ..plasmid_builder import assemble
from ..types import (
    DeletionPrimerSet,
    DesignRequest,
    DesignResult,
    GeneRecord,
    OffTargetReport,
)
from ..vectors import reverse_complement


def run(request: DesignRequest) -> DesignResult:
    """End-to-end deletion design (skill_v2 §4)."""
    result = design_without_off_target(request)
    genome_bytes = storage_adapter.fetch_genome(request.isolate_id)
    result.off_target = _scan(result, genome_bytes)
    if not result.off_target.passed:
        raise OffTargetDetected(
            message="off-target scan reports unintended products",
            details={"violations": [vars(v) for v in result.off_target.violations]},
        )
    verification.verify(result)
    return result


def design_without_off_target(request: DesignRequest) -> DesignResult:
    """Search + assemble + verify (excluding off-target gates).

    Useful for environments without a genome available; the resulting
    ``DesignResult.off_target`` is an empty placeholder until a caller fills it
    by invoking ``_scan`` separately.
    """
    if request.application != "deletion":
        raise ValueError(f"deletion.run called with application={request.application!r}")

    gene = gene_finder.get_gene_record(request.isolate_id, request.gene)
    vector = vectors.get_vector(request.vector)
    convention = vectors.get_tail_convention(request.vector, request.enzyme, "deletion")

    primer_set, _top5 = primer_designer.search_deletion_primers(gene, vector, convention)
    p1_anchor, p4_anchor = _anchor_offsets(gene, primer_set)
    up_amp, dn_amp, insert = _build_amplicons(gene, primer_set, p1_anchor, p4_anchor)
    cut_nick = vectors.find_cut_position(vector, request.enzyme)
    final_plasmid = assemble(vector, insert, cut_nick, convention)

    return DesignResult(
        request=request,
        gene_record=gene,
        vector=vector,
        convention=convention,
        primer_set=primer_set,
        up_amplicon=up_amp,
        dn_amplicon=dn_amp,
        insert=insert,
        final_plasmid=final_plasmid,
        off_target=OffTargetReport(sites_per_primer={}, products=[], passed=True),
    )


def _scan(result: DesignResult, genome_bytes: bytes) -> OffTargetReport:
    ps: DeletionPrimerSet = result.primer_set      # type: ignore[assignment]
    p1_anchor, p4_anchor = _anchor_offsets(result.gene_record, ps)
    expected = _expected_off_target_products(result.gene_record, ps, p1_anchor, p4_anchor)
    return primer_designer.off_target_scan(
        [ps.p1, ps.p2, ps.p3, ps.p4], genome_bytes, expected,
    )


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _anchor_offsets(gene: GeneRecord, ps: DeletionPrimerSet) -> tuple[int, int]:
    p1_anchor = gene.up_flank.find(ps.p1.body)
    if p1_anchor < 0:
        raise ValueError(f"P1 body not found in up_flank ({gene.gene}/{gene.isolate_id})")
    fwd_p4 = reverse_complement(ps.p4.body)
    fwd_start = gene.dn_flank.find(fwd_p4)
    if fwd_start < 0:
        raise ValueError(f"P4 body not found in dn_flank ({gene.gene}/{gene.isolate_id})")
    p4_anchor = len(gene.dn_flank) - fwd_start - len(ps.p4.body)
    return p1_anchor, p4_anchor


def _build_amplicons(
    gene: GeneRecord,
    ps: DeletionPrimerSet,
    p1_anchor: int,
    p4_anchor: int,
) -> tuple[str, str, str]:
    """Reconstruct UP/DN amplicons (top strand) and the assembled In-Fusion insert."""
    L = len(gene.cds_seq) // 3
    up_segment = gene.up_flank + gene.cds_seq[: 3 * ps.N]
    dn_segment = gene.cds_seq[3 * (L - ps.C - 1) :] + gene.dn_flank

    j = JUNCTION_LEN_PER_PRIMER_DEFAULT  # 15
    up_amp = ps.p1.tail + up_segment[p1_anchor:] + dn_segment[:j]
    dn_seg_end = len(dn_segment) - p4_anchor
    dn_amp = up_segment[-j:] + dn_segment[:dn_seg_end] + reverse_complement(ps.p4.tail)
    insert = up_amp + dn_amp[2 * j :]
    return up_amp, dn_amp, insert


def _expected_off_target_products(
    gene: GeneRecord,
    ps: DeletionPrimerSet,
    p1_anchor: int,
    p4_anchor: int,
) -> dict[tuple[str, str], tuple[int, int] | None]:
    """Per skill_v2 §10. Genome amplicon sizes (no vector tails — tails are
    non-templated)."""
    up_to_locus = len(gene.up_flank) - p1_anchor
    dn_from_locus = len(gene.dn_flank) - p4_anchor
    p1p2 = up_to_locus + 3 * ps.N
    p3p4 = 3 * (ps.C + 1) + dn_from_locus
    p1p4 = up_to_locus + len(gene.cds_seq) + dn_from_locus
    tol = 5
    return {
        ("P1", "P2"): (p1p2 - tol, p1p2 + tol),
        ("P3", "P4"): (p3p4 - tol, p3p4 + tol),
        ("P1", "P4"): (p1p4 - tol, p1p4 + tol),
        ("P3", "P2"): None,
    }

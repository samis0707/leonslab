"""Orchestrate the pEXG2 in-frame deletion pipeline.

Per skill_v2 §4 and project_plan §4.1.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Callable

from .. import colony_pcr as _colony_pcr
from .. import cross_isolate_checker, gene_finder, primer_designer, storage_adapter, vectors, verification
from ..config import (
    JUNCTION_LEN_DEFAULT,
    JUNCTION_LEN_PER_PRIMER_DEFAULT,
    OFFTARGET_PRODUCT_SIZE_MAX_BP,
    VECTOR_TAIL_LEN,
)
from ..exceptions import OffTargetDetected
from ..plasmid_builder import assemble
from ..types import DesignRequest, DesignResult
from ..vectors import reverse_complement


GenomeLoader = Callable[[str], bytes]


def run(
    request: DesignRequest,
    *,
    genome_loader: GenomeLoader | None = None,
) -> DesignResult:
    """End-to-end deletion design.

    Steps (skill_v2 §4):
        1. Load gene record (bundled, no R2 fetch).
        2. Load vector + tail convention; locate cut nick.
        3. Search primer set (exhaustive (N, C) × P1 × P2 × P3 × P4).
        4. Build UP/DN amplicons and the joined insert.
        5. Lazy-fetch genome from R2 for off-target scan.
        6. Off-target scan; raise on unintended product.
        7. Assemble final plasmid.
        8. Verify (hard checks).
        9. Return DesignResult.
    """
    gene = gene_finder.get_gene_record(request.isolate_id, request.gene)
    vector = vectors.get_vector(request.vector)
    convention = vectors.get_tail_convention(
        request.vector, request.enzyme, "deletion"
    )
    cut_nick = vectors.find_cut_position(vector, request.enzyme)

    best, alts = primer_designer.search_deletion_primers(
        gene, vector, convention
    )

    loader = genome_loader or _default_genome_loader
    genome_bytes = loader(request.isolate_id)

    # Try the top-scoring set first; if it produces an off-target product,
    # walk the top-N alternatives and pick the first that scans cleanly.
    # `alts` is already sorted ascending by score (best last). We dedupe by
    # primer-tuple identity since the search returns near-duplicates when
    # different (N, C) pairs yield the same P1/P4 bodies.
    candidates: list = [best]
    seen = {(best.p1.body, best.p2.body, best.p3.body, best.p4.body)}
    for alt in alts:
        key = (alt.p1.body, alt.p2.body, alt.p3.body, alt.p4.body)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(alt)

    chosen = None
    chosen_amps = None
    chosen_off = None
    last_violations: list = []
    for cand in candidates:
        up_amplicon, dn_amplicon, insert = _build_amplicons(gene, cand)
        expected_products = _expected_products_for_deletion(
            up_amplicon, dn_amplicon, gene
        )
        off_target = primer_designer.off_target_scan(
            [cand.p1, cand.p2, cand.p3, cand.p4],
            genome_bytes,
            expected_products,
        )
        if off_target.passed:
            chosen = cand
            chosen_amps = (up_amplicon, dn_amplicon, insert)
            chosen_off = off_target
            break
        last_violations = off_target.violations

    if chosen is None:
        raise OffTargetDetected(
            message=(
                "Off-target scan reports unintended PCR products for the "
                f"{len(candidates)} top-scoring primer alternatives. Most "
                "likely cause: a paralog in this isolate's genome shares "
                "high homology with the lasR/lasB flank used for primer "
                "anchoring. Manual primer design or gene-specific anchor "
                "tuning required."
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
    final_plasmid.name = f"{vector.name}_{gene.gene}_delta_{gene.isolate_id}"

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
    if not gene.functional:
        reason = gene.truncation_reason or "non-functional allele in this isolate"
        result.warnings.append(
            f"NON-FUNCTIONAL ALLELE: {gene.gene} in {gene.isolate_id} appears to be "
            f"truncated/inactive at the protein level ({reason}). The primers below "
            f"target the nominal locus position (reconstructed from the closest "
            f"functional reference) and will work for allelic exchange, but the "
            f"native protein product is already absent or non-functional in this "
            f"isolate. Verify your experimental design accounts for this."
        )
    result.colony_pcr_primers = _colony_pcr.get_colony_pcr_primers(request.gene)
    result.compatible_isolates = cross_isolate_checker.find_compatible_isolates(
        best, request.gene, request.isolate_id
    )
    verification.verify(result)
    return result


def _build_amplicons(gene, primer_set):
    """Reconstruct UP / DN amplicons from the chosen primer set + gene record.

    UP amplicon (top strand):
        P1.tail + up_segment[p1_anchor:] + dn_segment[:15]
    DN amplicon (top strand):
        up_segment[-15:] + dn_segment[:p4_anchor_end] + RC(P4.tail)

    The 30-nt junction (p3_tail + p2_tail's RC = up_segment[-15:] + dn_segment[:15])
    is shared between the two amplicons, so the assembled insert collapses it once.
    """
    cds = gene.cds_seq
    L_codons = len(cds) // 3
    N = primer_set.N
    C = primer_set.C

    up_segment = gene.up_flank + cds[: 3 * N]
    dn_segment = cds[3 * (L_codons - C - 1) :] + gene.dn_flank

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

    overlap = JUNCTION_LEN_PER_PRIMER_DEFAULT  # 15

    up_amplicon = (
        primer_set.p1.tail
        + up_segment[p1_anchor:]
        + dn_segment[:overlap]
    )
    dn_amplicon = (
        up_segment[-overlap:]
        + dn_segment[:p4_end]
        + reverse_complement(primer_set.p4.tail)
    )
    insert = up_amplicon + dn_amplicon[JUNCTION_LEN_DEFAULT:]
    return up_amplicon, dn_amplicon, insert


def _expected_products_for_deletion(
    up_amplicon: str, dn_amplicon: str, gene
) -> dict:
    """Build the expected_products dict for off_target_scan.

    Note: _enumerate_products in primer_designer uses (sb_pos - sa_pos + 1) as
    the size estimate, which underestimates the true amplicon size by ~L_b - 1.
    Use a wide tolerance to absorb that.
    """
    up_size = len(up_amplicon)
    dn_size = len(dn_amplicon)
    tol = 150  # _enumerate_products underestimates by ~tail+body length
    # WT locus product (P1 + P4) on the same isolate genome: same genomic span as
    # UP + DN but with the full native CDS rather than the deleted scar.
    wt_size = up_size + dn_size - 30 + len(gene.cds_seq)
    return {
        ("P1", "P2"): (max(0, up_size - tol), up_size + tol),
        ("P3", "P4"): (max(0, dn_size - tol), dn_size + tol),
        ("P1", "P4"): (max(0, wt_size - tol), min(OFFTARGET_PRODUCT_SIZE_MAX_BP, wt_size + tol)),
    }


def _default_genome_loader(isolate_id: str) -> bytes:
    """Load genome bytes from R2, auto-loading credentials from api/design/.env
    if they aren't already set in the environment.
    """
    if "R2_ACCOUNT_ID" not in os.environ:
        _load_dotenv_if_present()
    return storage_adapter.fetch_genome(isolate_id)


def _load_dotenv_if_present() -> None:
    env_path = (
        Path(__file__).resolve().parents[3] / "api" / "design" / ".env"
    )
    if not env_path.exists():
        return
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        os.environ.setdefault(k.strip(), v.strip())

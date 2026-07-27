"""Wet-lab-validated primer sets that are preferred over a fresh computational
design whenever they fit the requested isolate.

Some (gene, application) combinations have primer sequences that were already
designed and used in the wet lab (see ``data/primer_design/fixed_primers.json``).
Only the isolate-annealing *body* of each primer is pinned by that file — the
vector-homology tail (and, for expression, the RBS) is always derived the same
way as for a computationally-designed primer (``vectors.derive_tails`` /
``config.RBS_TAIL_PART``), so a fixed set stays consistent with whichever
vector/enzyme/application convention the request specifies.

A fixed set is only used when every body matches the requested isolate's
curated sequence verbatim, at a genomically consistent position (in-frame
junction for deletion, correct ATG/stop anchoring for expression). If the
isolate's sequence differs at any of those positions, callers fall back to
``primer_designer.search_*_primers`` exactly as before.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path

from . import config as cfg
from . import primer_designer
from .primer_designer import _parse_fasta_bytes
from .types import DeletionPrimerSet, ExpressionPrimerSet, GeneRecord, Primer
from .vectors import derive_tails, reverse_complement


_DATA_PATH = Path(__file__).resolve().parents[2] / "data" / "primer_design" / "fixed_primers.json"

# How far past the curated (capped) flank to look, directly in the isolate's
# genome, when a fixed primer's anchor doesn't fall inside the curated window
# (e.g. lasB P4, which sits ~628 nt downstream of the stop — beyond the
# standard 600 nt dn_flank). Only used as a fallback; the fast path checks the
# curated gene-record flanks first with no genome access at all.
_WIDE_FLANK_LEN = 3000


@lru_cache(maxsize=1)
def _load_store() -> dict:
    if not _DATA_PATH.exists():
        return {}
    with open(_DATA_PATH) as f:
        return json.load(f)


def _get_bodies(gene: str, application: str) -> dict | None:
    entry = _load_store().get(f"{gene}|{application}")
    if entry is None:
        return None
    return {k: v.strip().upper() for k, v in entry.items() if k.endswith("_body")}


def _make_primer(name: str, role: str, tail: str, body: str, tail_kind: str) -> Primer:
    return Primer(
        name=name, role=role, tail=tail, body=body, tail_kind=tail_kind,
        tm_body_C=primer_designer.primer_body_tm(body),
        gc_body=primer_designer.primer_body_gc(body),
        length=len(tail) + len(body),
    )


def _wide_flanks(gene: GeneRecord, genome_bytes: bytes, span: int) -> tuple[str, str] | None:
    """Re-extract up/dn flanks for ``gene`` directly from the isolate's genome,
    wider than the curated flank stored in the gene record. Mirrors the
    strand-normalisation logic in scripts/build_gene_records.py::_strand_normalize.
    """
    contigs = dict(_parse_fasta_bytes(genome_bytes))
    seq = contigs.get(gene.contig_id)
    if seq is None:
        return None

    start_0based = gene.genome_start_1based - 1
    end_0based = gene.genome_end_1based

    if gene.original_strand == "+":
        wide_up = seq[max(0, start_0based - span):start_0based]
        wide_dn = seq[end_0based:min(len(seq), end_0based + span)]
    else:
        wide_up = reverse_complement(seq[end_0based:min(len(seq), end_0based + span)])
        wide_dn = reverse_complement(seq[max(0, start_0based - span):start_0based])
    return wide_up, wide_dn


# ===========================================================================
# Deletion (pEXG2)
# ===========================================================================

@dataclass
class FixedDeletionResult:
    primer_set: DeletionPrimerSet
    gene_record: GeneRecord   # == input gene, or a copy with widened flanks


def _locate_deletion_junction(
    gene: GeneRecord, up_flank: str, dn_flank: str,
    p1_body: str, p2_body: str, p3_body: str, p4_body: str,
) -> tuple[int, int] | None:
    """Return (N, C) if all four bodies anchor consistently in-frame, else None."""
    cds = gene.cds_seq
    L = len(cds) // 3
    full = up_flank + cds + dn_flank
    cds_start = len(up_flank)
    cds_end = cds_start + len(cds)

    if p1_body not in up_flank:
        return None
    if reverse_complement(p4_body) not in dn_flank:
        return None

    rc_p2 = reverse_complement(p2_body)
    p2_idx = full.find(rc_p2)
    if p2_idx < 0:
        return None
    offset_n = (p2_idx + len(rc_p2)) - cds_start
    if offset_n < 0 or offset_n % 3 != 0:
        return None
    N = offset_n // 3

    p3_idx = full.find(p3_body, cds_start)
    if p3_idx < 0 or p3_idx >= cds_end:
        return None
    offset_c = p3_idx - cds_start
    if offset_c % 3 != 0:
        return None
    C = L - 1 - offset_c // 3

    if N < 0 or C < 0 or N >= L or C + 1 >= L:
        return None

    scar = primer_designer.scar_orf(cds, N, C)
    if not primer_designer.assert_scar_valid(scar, N, C):
        return None
    return N, C


def build_fixed_deletion_set(
    gene: GeneRecord, vector, convention, genome_bytes: bytes | None = None,
) -> FixedDeletionResult | None:
    """Return a pre-validated DeletionPrimerSet for ``gene``, if one is on file
    and its bodies match this isolate's sequence; otherwise None."""
    bodies = _get_bodies(gene.gene, "deletion")
    if bodies is None:
        return None
    p1b, p2b, p3b, p4b = (bodies[k] for k in ("p1_body", "p2_body", "p3_body", "p4_body"))

    used_gene = gene
    junction = _locate_deletion_junction(gene, gene.up_flank, gene.dn_flank, p1b, p2b, p3b, p4b)
    if junction is None and genome_bytes is not None:
        wide = _wide_flanks(gene, genome_bytes, _WIDE_FLANK_LEN)
        if wide is not None:
            wide_up, wide_dn = wide
            junction = _locate_deletion_junction(gene, wide_up, wide_dn, p1b, p2b, p3b, p4b)
            if junction is not None:
                used_gene = replace(gene, up_flank=wide_up, dn_flank=wide_dn)
    if junction is None:
        return None
    N, C = junction

    p1_tail, p4_tail = derive_tails(vector, convention.enzyme, "deletion")
    cds = used_gene.cds_seq
    L = len(cds) // 3
    up_segment = used_gene.up_flank + cds[: 3 * N]
    dn_segment = cds[3 * (L - C - 1):] + used_gene.dn_flank
    overlap = cfg.JUNCTION_LEN_PER_PRIMER_DEFAULT
    p3_tail = up_segment[-overlap:]
    p2_tail = reverse_complement(dn_segment[:overlap])
    scar = primer_designer.scar_orf(cds, N, C)

    p1 = _make_primer("P1", "UP_Fwd", p1_tail, p1b, convention.name)
    p2 = _make_primer("P2", "UP_Rev", p2_tail, p2b, "junction_overlap")
    p3 = _make_primer("P3", "DN_Fwd", p3_tail, p3b, "junction_overlap")
    p4 = _make_primer("P4", "DN_Rev", p4_tail, p4b, convention.name)
    tms = [p1.tm_body_C, p2.tm_body_C, p3.tm_body_C, p4.tm_body_C]

    primer_set = DeletionPrimerSet(
        p1=p1, p2=p2, p3=p3, p4=p4, N=N, C=C, scar_dna=scar,
        score=-1.0, tm_spread=max(tms) - min(tms),
    )
    return FixedDeletionResult(primer_set=primer_set, gene_record=used_gene)


# ===========================================================================
# Expression (pBBR1MCS-2)
# ===========================================================================

@dataclass
class FixedExpressionResult:
    primer_set: ExpressionPrimerSet
    gene_record: GeneRecord   # == input gene, or a copy with widened flanks


def build_fixed_expression_set(
    gene: GeneRecord, vector, convention, genome_bytes: bytes | None = None,
) -> FixedExpressionResult | None:
    """Return a pre-validated ExpressionPrimerSet for ``gene``, if one is on
    file and its bodies match this isolate's sequence; otherwise None.

    Only applies to plain (untagged) complementation — a requested tag has no
    room in a fixed body/tail, so callers should not call this when a tag is
    requested.
    """
    bodies = _get_bodies(gene.gene, "expression")
    if bodies is None:
        return None
    p1b, p2b = bodies["p1_body"], bodies["p2_body"]
    rc_p2 = reverse_complement(p2b)

    def _locate(up_flank: str, dn_flank: str) -> tuple[int, int] | None:
        full = up_flank + gene.cds_seq + dn_flank
        p1_idx = full.find(p1b)
        if p1_idx < 0:
            return None
        p2_idx = full.find(rc_p2, p1_idx)
        if p2_idx < 0:
            return None
        return p1_idx, p2_idx + len(rc_p2)

    used_gene = gene
    span = _locate(gene.up_flank, gene.dn_flank)
    if span is None and genome_bytes is not None:
        wide = _wide_flanks(gene, genome_bytes, _WIDE_FLANK_LEN)
        if wide is not None:
            wide_up, wide_dn = wide
            span = _locate(wide_up, wide_dn)
            if span is not None:
                used_gene = replace(gene, up_flank=wide_up, dn_flank=wide_dn)
    if span is None:
        return None
    start, end = span

    p1_vector_tail, p2_vector_tail = derive_tails(vector, convention.enzyme, "expression")
    p1_tail = p1_vector_tail + cfg.CANONICAL_RBS
    p2_tail = p2_vector_tail

    full = used_gene.up_flank + used_gene.cds_seq + used_gene.dn_flank
    native_anneal_segment = full[start:end]

    p1 = _make_primer("P1", "INSERT_Fwd", p1_tail, p1b, convention.name)
    p2 = _make_primer("P2", "INSERT_Rev", p2_tail, p2b, convention.name)
    tms = [p1.tm_body_C, p2.tm_body_C]

    primer_set = ExpressionPrimerSet(
        p1=p1, p2=p2,
        coding_seq=used_gene.cds_seq,
        native_anneal_segment=native_anneal_segment,
        tag=None, tag_position=None,
        score=-1.0, tm_spread=max(tms) - min(tms),
    )
    return FixedExpressionResult(primer_set=primer_set, gene_record=used_gene)

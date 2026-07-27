"""Hand-curated, wet-lab-verified primer sets that pre-date the automated search.

Some constructs (e.g. ``pEXG2_lasB_PA14``) were built and sequence-verified in
Benchling before this tool existed. Where an isolate's flank sequence matches
one of these curated sets exactly, reuse it verbatim instead of running the
computed (N, C) x P1 x P2 x P3 x P4 search -- it is the primer set actually in
use in the lab. If the isolate's sequence diverges anywhere the fixed set
touches, fall back to the computed search unchanged.

Data lives in ``data/primer_design/fixed_primers/<gene>_<vector>.json``. Only
``lasB``/``pEXG2`` is populated so far; the same mechanism will be extended to
lasR/lasI (pEXG2) and lasB (pBBR1MCS2 complementation) next.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path

from . import primer_designer
from .primer_designer import assert_scar_valid, compute_score, primer_body_gc, primer_body_tm, scar_orf
from .types import DeletionPrimerSet, GeneRecord, Primer, VectorRecord
from .vectors import reverse_complement

# How far past the bundled dn_flank/up_flank we're willing to look on the
# isolate's own genome for a fixed primer anchored further out than the
# standard 600/800 bp window (see e.g. lasB P4, anchored 629 bp downstream).
EXTENDED_FLANK_SEARCH_BP = 3000


@dataclass
class FixedPrimerSpec:
    gene: str
    vector: str
    N: int
    C: int
    p1_tail: str
    p1_body: str
    p2_tail: str
    p2_body: str
    p3_tail: str
    p3_body: str
    p4_tail: str
    p4_body: str
    source: str
    note: str


_CACHE: dict[str, FixedPrimerSpec | None] = {}


def _fixed_primers_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "data" / "primer_design" / "fixed_primers"


def load_fixed_primer_spec(gene: str, vector: str) -> FixedPrimerSpec | None:
    """Return the curated primer set for (gene, vector), or None if none is defined."""
    key = f"{gene}|{vector}"
    if key in _CACHE:
        return _CACHE[key]

    path = _fixed_primers_dir() / f"{gene}_{vector}.json"
    if not path.exists():
        _CACHE[key] = None
        return None

    data = json.loads(path.read_text())
    spec = FixedPrimerSpec(
        gene=data["gene"],
        vector=data["vector"],
        N=data["N"],
        C=data["C"],
        p1_tail=data["p1"]["tail"], p1_body=data["p1"]["body"],
        p2_tail=data["p2"]["tail"], p2_body=data["p2"]["body"],
        p3_tail=data["p3"]["tail"], p3_body=data["p3"]["body"],
        p4_tail=data["p4"]["tail"], p4_body=data["p4"]["body"],
        source=data.get("source", ""),
        note=data.get("note", ""),
    )
    _CACHE[key] = spec
    return spec


def _find_contig(contigs: dict[str, str], contig_id: str) -> str | None:
    """Look up ``contig_id`` in a {contig_id: sequence} map, tolerating the
    isolate-id-prefixed naming some R2-uploaded whole genomes use (e.g. the
    bundled gene record's ``contig=NZ_CP104983.1`` vs. R2's
    ``PA14_NZ_CP104983.1``)."""
    if contig_id in contigs:
        return contigs[contig_id]
    for cid, seq in contigs.items():
        if cid.endswith(contig_id) or contig_id.endswith(cid):
            return seq
    return None


def _extend_dn_flank(gene: GeneRecord, genome_bytes: bytes, extra_bp: int) -> str | None:
    """Return a longer coding-strand-oriented dn_flank, fetched from the isolate's
    own genome, for primers anchored past the bundled ``gene.dn_flank`` window.

    Returns None if the gene's contig can't be found in ``genome_bytes``.
    """
    contigs = dict(primer_designer._parse_fasta_bytes(genome_bytes))
    seq = _find_contig(contigs, gene.contig_id)
    if seq is None:
        return None

    if gene.original_strand == "+":
        dn_end = min(len(seq), gene.genome_end_1based + extra_bp)
        return seq[gene.genome_end_1based : dn_end]
    else:
        up_start = max(0, gene.genome_start_1based - 1 - extra_bp)
        window = seq[up_start : gene.genome_start_1based - 1]
        return reverse_complement(window)


def try_build_fixed_primer_set(
    gene: GeneRecord,
    vector: VectorRecord,
    genome_bytes: bytes | None,
) -> tuple[DeletionPrimerSet | None, GeneRecord]:
    """Return ``(primer_set, gene_record_for_amplicon_building)`` if the curated
    fixed set for (gene.gene, vector.name) matches this isolate exactly.

    Returns ``(None, gene)`` if no fixed set is defined for this (gene, vector),
    or if the isolate's sequence diverges anywhere the fixed set anchors -- in
    both cases the caller should fall back to the computed search unchanged.
    """
    spec = load_fixed_primer_spec(gene.gene, vector.name)
    if spec is None:
        return None, gene

    if spec.p1_body not in gene.up_flank:
        return None, gene

    L = len(gene.cds_seq) // 3
    junction_pos = 3 * (L - spec.C - 1)
    if gene.cds_seq[junction_pos : junction_pos + len(spec.p3_body)] != spec.p3_body:
        return None, gene
    if junction_pos < 3 * spec.N:
        # N/C windows would overlap -- the fixed scar geometry doesn't fit this CDS.
        return None, gene

    gene_for_build = gene
    p4_body_rc = reverse_complement(spec.p4_body)
    if p4_body_rc not in gene.dn_flank:
        if genome_bytes is None:
            return None, gene
        extended = _extend_dn_flank(gene, genome_bytes, EXTENDED_FLANK_SEARCH_BP)
        if extended is None or p4_body_rc not in extended:
            return None, gene
        gene_for_build = replace(gene, dn_flank=extended)

    scar = scar_orf(gene.cds_seq, spec.N, spec.C)
    if not assert_scar_valid(scar, spec.N, spec.C):
        return None, gene

    def _primer(name: str, role: str, tail: str, body: str, tail_kind: str) -> Primer:
        return Primer(
            name=name, role=role, tail=tail, body=body, tail_kind=tail_kind,
            tm_body_C=primer_body_tm(body), gc_body=primer_body_gc(body),
            length=len(tail) + len(body),
        )

    p1 = _primer("P1", "UP_Fwd", spec.p1_tail, spec.p1_body, "fixed_verified")
    p2 = _primer("P2", "UP_Rev", spec.p2_tail, spec.p2_body, "fixed_verified")
    p3 = _primer("P3", "DN_Fwd", spec.p3_tail, spec.p3_body, "fixed_verified")
    p4 = _primer("P4", "DN_Rev", spec.p4_tail, spec.p4_body, "fixed_verified")

    tms = [p1.tm_body_C, p2.tm_body_C, p3.tm_body_C, p4.tm_body_C]
    gcs = [p1.gc_body, p2.gc_body, p3.gc_body, p4.gc_body]
    primer_set = DeletionPrimerSet(
        p1=p1, p2=p2, p3=p3, p4=p4,
        N=spec.N, C=spec.C,
        scar_dna=scar,
        score=compute_score(tms, gcs, body_lengths=[p.length for p in (p1, p2, p3, p4)]),
        tm_spread=max(tms) - min(tms),
    )
    return primer_set, gene_for_build

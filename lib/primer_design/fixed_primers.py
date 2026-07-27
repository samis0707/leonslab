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
from .config import RESTRICTION_SITES
from .primer_designer import _parse_fasta_bytes
from .types import DeletionPrimerSet, ExpressionPrimerSet, GeneRecord, Primer
from .vectors import derive_tails, forbidden_body_5prime_prefixes, reverse_complement


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


def _get_entry(gene: str, application: str) -> dict | None:
    return _load_store().get(f"{gene}|{application}")


def _get_bodies(gene: str, application: str) -> dict | None:
    entry = _get_entry(gene, application)
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
        # Some reference genomes (e.g. PA14) store contigs in R2 with the
        # isolate ID prefixed onto the RefSeq accession in the FASTA header
        # (e.g. "PA14_NZ_CP104983.1"), while the curated gene record kept the
        # bare accession ("NZ_CP104983.1") from the original source. Fall
        # back to a suffix match for that pattern.
        suffix = "_" + gene.contig_id
        matches = [s for cid, s in contigs.items() if cid.endswith(suffix)]
        if len(matches) == 1:
            seq = matches[0]
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
    expected_recognition_count: int   # verification target for this specific set;
                                       # usually == convention's, occasionally overridden
                                       # (see fixed_primers.json "expected_recognition_count")


def _deletion_bodies_fit(
    gene: GeneRecord, up_flank: str, dn_flank: str,
    p1_body: str, p2_body: str, p3_body: str, p4_body: str,
) -> bool:
    """True iff all four bodies anneal somewhere in the expected region of this
    isolate's sequence. P1/P4 must sit in the (possibly widened) flanks — the
    genuinely isolate-variable, non-coding regions. P2/P3 anneal within the
    CDS near the fixed N/C scar boundary (see ``build_fixed_deletion_set``);
    their exact annealing offset is a primer-design choice, not something the
    scar boundary is derived from, so this only checks they anneal somewhere
    in the CDS-adjacent region, not at an exact position.
    """
    if p1_body not in up_flank:
        return False
    if reverse_complement(p4_body) not in dn_flank:
        return False
    cds = gene.cds_seq
    if reverse_complement(p2_body) not in (up_flank + cds):
        return False
    if p3_body not in (cds + dn_flank):
        return False
    return True


def build_fixed_deletion_set(
    gene: GeneRecord, vector, convention, genome_bytes: bytes | None = None,
) -> FixedDeletionResult | None:
    """Return a pre-validated DeletionPrimerSet for ``gene``, if one is on file
    and its bodies match this isolate's sequence; otherwise None.

    The retained-codon scar (N, C) is a fixed design choice (which codons
    stay, i.e. the resulting protein scar) — it is stored directly in
    fixed_primers.json rather than inferred from where the P2/P3 bodies
    happen to anneal, since a hand-designed primer's body need not terminate
    exactly at the scar boundary (only its *tail*, rebuilt below from N/C,
    encodes the actual junction).
    """
    entry = _get_entry(gene.gene, "deletion")
    if entry is None or "N" not in entry or "C" not in entry:
        return None
    N, C = int(entry["N"]), int(entry["C"])
    bodies = {k: v.strip().upper() for k, v in entry.items() if k.endswith("_body")}
    p1b, p2b, p3b, p4b = (bodies[k] for k in ("p1_body", "p2_body", "p3_body", "p4_body"))

    L = len(gene.cds_seq) // 3
    if N < 0 or C < 0 or N >= L or C + 1 >= L:
        return None
    scar = primer_designer.scar_orf(gene.cds_seq, N, C)
    if not primer_designer.assert_scar_valid(scar, N, C):
        return None

    used_gene = gene
    fits = _deletion_bodies_fit(gene, gene.up_flank, gene.dn_flank, p1b, p2b, p3b, p4b)
    if not fits and genome_bytes is not None:
        wide = _wide_flanks(gene, genome_bytes, _WIDE_FLANK_LEN)
        if wide is not None:
            wide_up, wide_dn = wide
            if _deletion_bodies_fit(gene, wide_up, wide_dn, p1b, p2b, p3b, p4b):
                fits = True
                used_gene = replace(gene, up_flank=wide_up, dn_flank=wide_dn)
    if not fits:
        return None

    p1_tail, p4_tail = derive_tails(vector, convention.enzyme, "deletion")

    # A body's 5' end can accidentally regenerate the enzyme's recognition
    # site right at the vector junction once fused to the tail (e.g. a P4
    # body starting with "T" completes "...A" + tail's leading "AGCTT" into
    # "AAGCTT"). The dynamic search screens candidates for this; a fixed body
    # needs the same screen, since here it wasn't chosen to avoid it — unless
    # the JSON entry explicitly declares (and overrides the expected count
    # for) an accepted extra site, e.g. lasB's P4.
    expected_recognition_count = int(
        entry.get("expected_recognition_count", convention.expected_recognition_count_in_final_plasmid)
    )
    motif, _ = RESTRICTION_SITES[convention.enzyme]
    p1_forbidden, p4_forbidden = forbidden_body_5prime_prefixes(
        p1_tail, p4_tail, motif, expected_count=expected_recognition_count,
    )
    if any(p1b.startswith(p) for p in p1_forbidden):
        return None
    if any(p4b.startswith(p) for p in p4_forbidden):
        return None

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
    return FixedDeletionResult(
        primer_set=primer_set, gene_record=used_gene,
        expected_recognition_count=expected_recognition_count,
    )


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

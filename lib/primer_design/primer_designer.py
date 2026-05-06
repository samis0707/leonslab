"""Primer-design core: Tm calculation, hard filters, scoring, candidate generation,
search algorithms for the three applications.

Public APIs (see ``project_plan.md`` §3.4):
    search_deletion_primers   — exhaustive (N, C) × P1 × P2 × P3 × P4 search
    search_tagging_primers    — junction-fixed search with tag cassette in tails
    search_expression_primers — 2-primer optimization
    off_target_scan           — anchor-walk + body-extension scan against one genome

Implemented here:
    - Tm via Allawi & SantaLucia (Bio.SeqUtils.MeltingTemp.Tm_NN)
    - Hard body filters (D5.3)
    - Scoring function
    - Candidate enumeration (P1/P4 from flanks, P2/P3 from segment ends)
    - Tm-spread cutoff

Stubbed (Phase 2):
    - The orchestration loops in search_*_primers
    - off_target_scan walk and product enumeration
"""
from __future__ import annotations

from dataclasses import dataclass

from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_NN

from . import config as cfg
from .exceptions import NoCandidates, OffTargetDetected, PrimerDesignError
from .types import (
    DeletionPrimerSet,
    ExpressionPrimerSet,
    GeneRecord,
    OffTargetReport,
    OffTargetSite,
    Primer,
    Tag,
    TaggingPrimerSet,
    TailConvention,
    VectorRecord,
)
from .vectors import derive_tails, reverse_complement


# ===========================================================================
# Tm calculation
# ===========================================================================

def primer_body_tm(body: str) -> float:
    """Nearest-neighbor Tm (D5.1) for the primer body in In-Fusion buffer."""
    return float(Tm_NN(Seq(body), **cfg.TM_NN_PARAMS))


def primer_body_gc(body: str) -> float:
    return sum(b in "GC" for b in body) / len(body)


# ===========================================================================
# Hard filters (D5.3)
# ===========================================================================

def passes_hard_filters(body: str) -> bool:
    """True iff body meets all hard filters. See D5.3 / skill_v2 §7."""
    n = len(body)
    if not (cfg.BODY_LEN_MIN <= n <= cfg.BODY_LEN_MAX):
        return False

    # 3' G/C clamp
    if body[-1] not in cfg.CLAMP_LAST_BASE_OK:
        return False
    last5_gc = sum(b in "GC" for b in body[-5:])
    if not (cfg.CLAMP_LAST5_GC_MIN <= last5_gc <= cfg.CLAMP_LAST5_GC_MAX):
        return False
    if cfg.CLAMP_NO_4IDENT_LAST4 and len(set(body[-4:])) == 1:
        return False

    # 4-homopolymer anywhere
    for i in range(n - 3):
        if body[i] == body[i + 1] == body[i + 2] == body[i + 3]:
            return False

    # GC content
    gc = primer_body_gc(body)
    if not (cfg.GC_MIN <= gc <= cfg.GC_MAX):
        return False

    # Tm
    tm = primer_body_tm(body)
    if not (cfg.TM_HARD_MIN_C <= tm <= cfg.TM_HARD_MAX_C):
        return False

    # 3' self-dimer
    if max_3prime_self_dimer(body) > cfg.SELF_DIMER_MAX_3PRIME:
        return False

    return True


def max_3prime_self_dimer(body: str) -> int:
    """Length of the longest contiguous Watson-Crick match between the 3' end of
    ``body`` and its reverse complement, sliding the RC against ``body`` in 8-nt
    windows. Used to penalize primer-dimer formation.
    """
    rc = reverse_complement(body)
    best = 0
    for shift in range(0, len(body) - 4):
        # Compare body[-(8+shift):] against rc starting at offset
        a = body[max(0, len(body) - 8 - shift) :]
        b = rc[: len(a)]
        run = 0
        for x, y in zip(reversed(a), reversed(b)):
            if x == _complement(y):
                run += 1
            else:
                break
        best = max(best, run)
    return best


def _complement(base: str) -> str:
    return {"A": "T", "T": "A", "G": "C", "C": "G"}.get(base, "N")


# ===========================================================================
# Candidate enumeration
# ===========================================================================

@dataclass(frozen=True)
class Candidate:
    body: str
    tm: float
    gc: float
    anchor_offset: int    # 0-based position within the segment


def enumerate_p1_candidates(up_flank: str, p1_tail: str) -> list[Candidate]:
    """P1 candidates: forward primer anchored near 5' end of up_flank.

    For each (offset, length) within the configured offset/length windows, accept
    the body if it passes hard filters. Returns sorted-by-Tm list.
    """
    out: list[Candidate] = []
    for offset in range(cfg.P1_OFFSET_MIN, cfg.P1_OFFSET_MAX + 1):
        for L in range(cfg.BODY_LEN_MIN, cfg.BODY_LEN_MAX + 1):
            if offset + L > len(up_flank):
                break
            body = up_flank[offset : offset + L]
            if passes_hard_filters(body):
                out.append(
                    Candidate(body=body, tm=primer_body_tm(body),
                              gc=primer_body_gc(body), anchor_offset=offset)
                )
    out.sort(key=lambda c: c.tm)
    return out


def enumerate_p4_candidates(dn_flank: str, p4_tail: str) -> list[Candidate]:
    """P4 candidates: reverse primer anchored near 3' end of dn_flank.

    Body is RC of the forward sequence at the 3' anchor. Tm computed on the body
    sequence as it will appear in the primer (i.e., on the RC strand).
    """
    out: list[Candidate] = []
    for offset in range(cfg.P4_OFFSET_MIN, cfg.P4_OFFSET_MAX + 1):
        for L in range(cfg.BODY_LEN_MIN, cfg.BODY_LEN_MAX + 1):
            end = len(dn_flank) - offset
            start = end - L
            if start < 0:
                break
            fwd_segment = dn_flank[start:end]
            body = reverse_complement(fwd_segment)
            if passes_hard_filters(body):
                out.append(
                    Candidate(body=body, tm=primer_body_tm(body),
                              gc=primer_body_gc(body), anchor_offset=offset)
                )
    out.sort(key=lambda c: c.tm)
    return out


def enumerate_p2_candidates(up_segment: str) -> list[Candidate]:
    """P2 candidates: reverse primer anchored at 3' end of UP segment (= UP flank + retained N codons).

    Body is RC of the last L nt of up_segment. No offset window — the anchor is
    fixed at the 3' end of UP because the junction overlap defines the cut.
    """
    out: list[Candidate] = []
    for L in range(cfg.BODY_LEN_MIN, cfg.BODY_LEN_MAX + 1):
        if L > len(up_segment):
            break
        body = reverse_complement(up_segment[-L:])
        if passes_hard_filters(body):
            out.append(
                Candidate(body=body, tm=primer_body_tm(body),
                          gc=primer_body_gc(body), anchor_offset=0)
            )
    out.sort(key=lambda c: c.tm)
    return out


def enumerate_p3_candidates(dn_segment: str) -> list[Candidate]:
    """P3 candidates: forward primer anchored at 5' end of DN segment (= retained C codons + DN flank)."""
    out: list[Candidate] = []
    for L in range(cfg.BODY_LEN_MIN, cfg.BODY_LEN_MAX + 1):
        if L > len(dn_segment):
            break
        body = dn_segment[:L]
        if passes_hard_filters(body):
            out.append(
                Candidate(body=body, tm=primer_body_tm(body),
                          gc=primer_body_gc(body), anchor_offset=0)
            )
    out.sort(key=lambda c: c.tm)
    return out


# ===========================================================================
# Scoring
# ===========================================================================

def compute_score(tms: list[float], gcs: list[float]) -> float:
    """Lower is better. Hard reject elsewhere if tm_spread > TM_SPREAD_HARD_LIMIT_C."""
    tm_spread = max(tms) - min(tms)
    mean_tm = sum(tms) / len(tms)
    gc_spread = max(gcs) - min(gcs)
    gc_penalty = sum(
        0.1 * max(0.0, cfg.GC_PREFERRED_MIN - g) + 0.1 * max(0.0, g - cfg.GC_PREFERRED_MAX)
        for g in gcs
    )
    return tm_spread + 0.25 * abs(mean_tm - cfg.TM_TARGET_C) + 0.02 * gc_spread + gc_penalty


def closest_by_tm(pool: list[Candidate], target_tm: float) -> Candidate | None:
    """Binary-search-friendly lookup since pool is sorted by Tm."""
    if not pool:
        return None
    # Linear is fine; pools are small (<300 entries).
    return min(pool, key=lambda c: abs(c.tm - target_tm))


# ===========================================================================
# Search orchestration — STUBS
# ===========================================================================

# ===========================================================================
# Scar geometry (deletion application — D5.4 / skill_v2 §4.1)
# ===========================================================================

def search_scars(cds_len_codons: int, scar_max_total: int = cfg.SCAR_MAX_TOTAL_AA,
                 scar_min_per_side: int = cfg.SCAR_MIN_PER_SIDE) -> list[tuple[int, int]]:
    """Yield all (N, C) pairs with N, C ≥ 1 and N + C ≤ scar_max_total.

    N counts retained N-terminal amino-acid codons (start with ATG).
    C counts retained C-terminal amino-acid codons NOT including the stop codon
    — the stop is always carried along separately. Smallest scar = 2 aa total
    (one N + one C) + stop.
    """
    pairs: list[tuple[int, int]] = []
    for total in range(2 * scar_min_per_side, scar_max_total + 1):
        for N in range(scar_min_per_side, total - scar_min_per_side + 1):
            C = total - N
            pairs.append((N, C))
    return pairs


def scar_orf(cds: str, N: int, C: int) -> str:
    """Return retained N N-terminal codons + retained C C-terminal aa codons + stop.

    ``cds`` is assumed to include the stop codon at its 3' end. Resulting scar has
    length ``(N + C + 1) * 3`` nucleotides (the +1 is the stop).

    Note on indexing: ``L = len(cds) // 3`` is the codon count *including stop*,
    so the start of the retained C-terminal slice (which must include the stop)
    is at codon index ``L - C - 1``, not ``L - C``. The skill_v2 §4.1 worked
    example ("AL*" for LB001 lasB with C=2) confirms this off-by-one.
    """
    L = len(cds) // 3
    return cds[: 3 * N] + cds[3 * (L - C - 1) :]


def assert_scar_valid(scar: str, N: int, C: int) -> bool:
    """Return True iff ``scar`` is a valid in-frame ORF starting M, single trailing stop,
    and length matches the expected (N + C + 1) * 3 nucleotides."""
    if len(scar) != (N + C + 1) * 3:
        return False
    prot = str(Seq(scar).translate())
    if not prot.startswith("M"):
        return False
    if not prot.endswith("*"):
        return False
    if prot.count("*") != 1:
        return False
    return True


# ===========================================================================
# Top-K helper
# ===========================================================================

def _topk_insert(top: list, candidate, key=lambda x: x.score, k: int = 5) -> list:
    """Return updated top-k list (lowest score = best)."""
    merged = top + [candidate]
    merged.sort(key=key)
    return merged[:k]


# ===========================================================================
# Primer assembly
# ===========================================================================

def _make_primer(name: str, role: str, tail: str, candidate: Candidate, tail_kind: str) -> Primer:
    return Primer(
        name=name, role=role, tail=tail, body=candidate.body, tail_kind=tail_kind,
        tm_body_C=candidate.tm, gc_body=candidate.gc,
        length=len(tail) + len(candidate.body),
    )


# ===========================================================================
# Search orchestration
# ===========================================================================

def search_deletion_primers(
    gene: GeneRecord,
    vector: VectorRecord,
    convention: TailConvention,
) -> tuple[DeletionPrimerSet, list[DeletionPrimerSet]]:
    """Exhaustive (N, C) × P1 × P2 × P3 × P4 search per skill_v2 §4.3.

    P1 and P4 are anchored in the up-/down-flank and don't depend on the (N, C)
    choice, so their candidate pools are computed once. For each viable (N, C)
    pair we generate P2 / P3 pools (anchored at the UP/DN segment boundaries)
    and pick (P1, P4) closest in Tm to the (P2 + P3) average.

    Returns:
        (best_set, top5_alternatives) — alternatives sorted ascending by score.

    Raises:
        NoCandidates: if no tuple satisfies all hard filters and Tm-spread bound.
    """
    p1_tail, p4_tail = derive_tails(vector, convention.enzyme, "deletion")
    cds = gene.cds_seq
    L = len(cds) // 3                                          # total codons including stop

    p1_pool = enumerate_p1_candidates(gene.up_flank, p1_tail)
    p4_pool = enumerate_p4_candidates(gene.dn_flank, p4_tail)
    if not p1_pool:
        raise NoCandidates(message=f"P1: no body in up-flank passes hard filters "
                                   f"({gene.gene}/{gene.isolate_id}). Relax GC range or offset window.",
                           details={"primer": "P1", "gene": gene.gene, "isolate": gene.isolate_id})
    if not p4_pool:
        raise NoCandidates(message=f"P4: no body in dn-flank passes hard filters "
                                   f"({gene.gene}/{gene.isolate_id}). Relax GC range or offset window.",
                           details={"primer": "P4", "gene": gene.gene, "isolate": gene.isolate_id})

    junction_tail = cfg.JUNCTION_LEN_PER_PRIMER_DEFAULT       # 15 nt per primer

    best: DeletionPrimerSet | None = None
    top5: list[DeletionPrimerSet] = []

    for N, C in search_scars(L, cfg.SCAR_MAX_TOTAL_AA, cfg.SCAR_MIN_PER_SIDE):
        # The C-terminal slice must include the stop codon as its last codon.
        if N >= L or C + 1 >= L:
            continue
        scar = scar_orf(cds, N, C)
        if not assert_scar_valid(scar, N, C):
            continue

        up_segment = gene.up_flank + cds[: 3 * N]              # part the UP amplicon spans
        dn_segment = cds[3 * (L - C - 1) :] + gene.dn_flank    # part the DN amplicon spans
        if len(up_segment) < junction_tail or len(dn_segment) < junction_tail:
            continue

        # Junction overlap: 30 nt total (15 from P3 tail at end of UP, 15 from P2 tail at start of DN).
        p3_tail = up_segment[-junction_tail:]                  # last 15 nt of UP segment, fwd strand
        p2_tail = reverse_complement(dn_segment[:junction_tail])

        p2_pool = enumerate_p2_candidates(up_segment)
        p3_pool = enumerate_p3_candidates(dn_segment)
        if not p2_pool or not p3_pool:
            continue

        for p2c in p2_pool:
            for p3c in p3_pool:
                target_tm = (p2c.tm + p3c.tm) / 2.0
                p1c = closest_by_tm(p1_pool, target_tm)
                p4c = closest_by_tm(p4_pool, target_tm)
                if p1c is None or p4c is None:
                    continue
                tms = [p1c.tm, p2c.tm, p3c.tm, p4c.tm]
                spread = max(tms) - min(tms)
                if spread > cfg.TM_SPREAD_HARD_LIMIT_C:
                    continue
                gcs = [p1c.gc, p2c.gc, p3c.gc, p4c.gc]
                score = compute_score(tms, gcs)

                cand = DeletionPrimerSet(
                    p1=_make_primer("P1", "UP_Fwd", p1_tail, p1c, convention.name),
                    p2=_make_primer("P2", "UP_Rev", p2_tail, p2c, "junction_overlap"),
                    p3=_make_primer("P3", "DN_Fwd", p3_tail, p3c, "junction_overlap"),
                    p4=_make_primer("P4", "DN_Rev", p4_tail, p4c, convention.name),
                    N=N, C=C, scar_dna=scar,
                    score=score, tm_spread=spread,
                )
                top5 = _topk_insert(top5, cand)
                if best is None or score < best.score:
                    best = cand

    if best is None:
        raise NoCandidates(
            message=f"No (N, C, P1-P4) tuple satisfied all hard filters for "
                    f"{gene.gene}/{gene.isolate_id}. Consider widening tm_spread or "
                    f"GC bounds; common cause: overly GC-rich up/dn flank in this isolate.",
            details={"gene": gene.gene, "isolate": gene.isolate_id,
                     "p1_pool_size": len(p1_pool), "p4_pool_size": len(p4_pool)},
        )
    return best, top5


def search_tagging_primers(
    gene: GeneRecord,
    vector: VectorRecord,
    convention: TailConvention,
    tag: Tag,
) -> tuple[TaggingPrimerSet, list[TaggingPrimerSet]]:
    """Junction-fixed search: tag cassette occupies the P2/P3 junction overlap.

    See ``primer_design_skill_v2.md`` §5.2–5.3 for cassette splitting between P2
    and P3 tails (15+15 for His6, 18+18 for FLAG/His8).

    Raises:
        TagTooLongForInLocus: cassette > 36 nt (3xFLAG, HiBiT).
    """
    raise NotImplementedError(
        "Phase 2 step 9: implement per skill_v2 §5.3 pseudocode. "
        "Note that tail derivation differs from deletion: P2 and P3 tails are partly "
        "fixed by the cassette, so search degrees of freedom collapse."
    )


def search_expression_primers(
    gene: GeneRecord,
    vector: VectorRecord,
    convention: TailConvention,
    tag: Tag | None,
    tag_position: str | None,
) -> tuple[ExpressionPrimerSet, list[ExpressionPrimerSet]]:
    """2-primer search for plasmid expression.

    See ``primer_design_skill_v2.md`` §6.2.

    Note:
        - P1 tail = vector_p1_tail + RBS + spacer (29 nt total for HindIII)
        - P2 tail = vector_p2_tail (15 nt)
        - Tag handling: see ``tags.build_plasmid_fusion_cds``.
    """
    raise NotImplementedError(
        "Phase 2 step 8: implement per skill_v2 §6.2. "
        "Reuse passes_hard_filters / primer_body_tm / compute_score."
    )


# ===========================================================================
# Off-target scan
# ===========================================================================

def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def _parse_fasta_bytes(data: bytes) -> list[tuple[str, str]]:
    """Tiny FASTA parser. Returns list of (contig_id, sequence_uppercase)."""
    contigs: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for raw in data.splitlines():
        line = raw.decode("ascii", errors="ignore") if isinstance(raw, bytes) else raw
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                contigs.append((name, "".join(chunks).upper()))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        contigs.append((name, "".join(chunks).upper()))
    return contigs


def _scan_primer_sites(
    body: str, contigs: list[tuple[str, str]], primer_name: str,
) -> list[OffTargetSite]:
    """Return all (contig, position, strand) sites where the primer body binds
    within the configured anchor + body mismatch limits.

    Position semantics: ``position_0based`` is the 0-based start of the body's
    binding footprint on the top strand. For "+" strand the body matches the top
    strand directly at [pos, pos+L); for "-" strand the body matches the RC of
    the top strand at [pos, pos+L) — i.e., the primer would prime extension
    leftward from pos+L−1.
    """
    L = len(body)
    A = cfg.OFFTARGET_ANCHOR_LEN
    if L < A:
        return []
    anchor_fwd = body[-A:]
    body_rc = reverse_complement(body)
    anchor_rev = body_rc[-A:]

    sites: list[OffTargetSite] = []
    for contig_id, seq in contigs:
        n = len(seq)
        # Forward-strand binding: body matches seq directly. Anchor is at the 3' end.
        # For an extension running 5'→3' on the top strand, the anchor sits at
        # window_start + L − A. So slide window_start over [0, n − L].
        for ws in range(0, n - L + 1):
            window = seq[ws : ws + L]
            if _hamming(window[-A:], anchor_fwd) > cfg.OFFTARGET_ANCHOR_MAX_MM:
                continue
            mm = _hamming(window, body)
            if mm <= cfg.OFFTARGET_BODY_MAX_MM:
                sites.append(OffTargetSite(
                    primer_name=primer_name, contig_id=contig_id,
                    position_0based=ws, strand="+", body_mismatches=mm,
                ))
        # Reverse-strand binding: body matches RC of seq window.
        for ws in range(0, n - L + 1):
            window_rc = reverse_complement(seq[ws : ws + L])
            if _hamming(window_rc[-A:], anchor_fwd) > cfg.OFFTARGET_ANCHOR_MAX_MM:
                continue
            mm = _hamming(window_rc, body)
            if mm <= cfg.OFFTARGET_BODY_MAX_MM:
                sites.append(OffTargetSite(
                    primer_name=primer_name, contig_id=contig_id,
                    position_0based=ws, strand="-", body_mismatches=mm,
                ))
    return sites


def _enumerate_pcr_products(
    sites_per_primer: dict[str, list[OffTargetSite]],
    pairs: list[tuple[str, str]],
) -> list:
    """For each primer pair (A, B), enumerate amplicons where A primes on +
    strand at position p_a and B primes on − strand at position p_b on the same
    contig with size = (p_b + L_b) − p_a within [SIZE_MIN, SIZE_MAX].

    Returns a list of dicts (we don't bind to the OffTargetProduct dataclass
    here so we can also report the body lengths used; the DesignResult holds
    the structured products.)
    """
    from .types import OffTargetProduct
    out: list[OffTargetProduct] = []
    smin = cfg.OFFTARGET_PRODUCT_SIZE_MIN_BP
    smax = cfg.OFFTARGET_PRODUCT_SIZE_MAX_BP

    # Group sites by (primer_name, contig_id, strand) for quick lookup
    by_key: dict[tuple[str, str, str], list[OffTargetSite]] = {}
    for name, sites in sites_per_primer.items():
        for s in sites:
            by_key.setdefault((name, s.contig_id, s.strand), []).append(s)

    seen: set[tuple[str, str, str, int, int]] = set()
    for a, b in pairs:
        contigs = {k[1] for k in by_key if k[0] in (a, b)}
        for contig in contigs:
            fwd_a = by_key.get((a, contig, "+"), [])
            rev_b = by_key.get((b, contig, "-"), [])
            fwd_b = by_key.get((b, contig, "+"), [])
            rev_a = by_key.get((a, contig, "-"), [])
            for f in fwd_a:
                for r in rev_b:
                    end = r.position_0based + _site_len(r, sites_per_primer, b)
                    size = end - f.position_0based
                    if smin <= size <= smax and r.position_0based >= f.position_0based:
                        key = (a, b, contig, f.position_0based, end)
                        if key in seen:
                            continue
                        seen.add(key)
                        out.append(OffTargetProduct(
                            primer_pair=(a, b), contig_id=contig,
                            start_0based=f.position_0based, end_0based=end,
                            size_bp=size, is_expected=False,
                        ))
            # Same pair, swapped roles (B forward, A reverse) — same chemistry, count once.
            if a != b:
                for f in fwd_b:
                    for r in rev_a:
                        end = r.position_0based + _site_len(r, sites_per_primer, a)
                        size = end - f.position_0based
                        if smin <= size <= smax and r.position_0based >= f.position_0based:
                            key = (a, b, contig, f.position_0based, end)
                            if key in seen:
                                continue
                            seen.add(key)
                            out.append(OffTargetProduct(
                                primer_pair=(a, b), contig_id=contig,
                                start_0based=f.position_0based, end_0based=end,
                                size_bp=size, is_expected=False,
                            ))
    return out


def _site_len(site: OffTargetSite, sites_per_primer, primer_name: str) -> int:
    """Recover the body length the site refers to; the OffTargetSite object
    does not carry it directly, but every site for one primer uses the same body."""
    # Find any other site for the same primer to read length, fall back via lookup.
    return _PRIMER_BODY_LEN_CACHE.get(primer_name, 0)


_PRIMER_BODY_LEN_CACHE: dict[str, int] = {}


def off_target_scan(
    primers: list[Primer],
    genome_fasta_bytes: bytes,
    expected_products: dict[tuple[str, str], tuple[int, int] | None],
) -> OffTargetReport:
    """Walk-anchor + body-mismatch scan per skill_v2 §10.

    Args:
        primers: list of Primer objects to scan. Body sequences are matched;
            tails are non-templated and ignored.
        genome_fasta_bytes: raw FASTA bytes for the isolate's genome.
        expected_products: keys are primer-pair tuples like ``("P1", "P2")``;
            value is ``(size_min_bp, size_max_bp)`` for an expected product, or
            ``None`` if zero products are expected for that pair.

    Returns:
        OffTargetReport with ``passed = True`` iff every expected pair has
        exactly the expected outcome and no unexpected products exist.
    """
    contigs = _parse_fasta_bytes(genome_fasta_bytes)
    if not contigs:
        raise PrimerDesignError(
            error_code="empty_genome",
            message="genome FASTA produced no contigs",
            details={"bytes_len": len(genome_fasta_bytes)},
        )

    _PRIMER_BODY_LEN_CACHE.clear()
    sites_per_primer: dict[str, list[OffTargetSite]] = {}
    for p in primers:
        _PRIMER_BODY_LEN_CACHE[p.name] = len(p.body)
        sites_per_primer[p.name] = _scan_primer_sites(p.body, contigs, p.name)

    pairs_to_check = list(expected_products.keys())
    products = _enumerate_pcr_products(sites_per_primer, pairs_to_check)

    violations: list = []
    matched_expected: set[tuple[str, str]] = set()
    for prod in products:
        key = prod.primer_pair
        bounds = expected_products.get(key)
        if bounds is None:
            violations.append(prod)
            continue
        smin, smax = bounds
        if smin <= prod.size_bp <= smax:
            if key in matched_expected:
                # Multiple products in expected size range → unexpected duplicate
                violations.append(prod)
            else:
                matched_expected.add(key)
                # Mark as expected by replacing it
                from .types import OffTargetProduct
                idx = products.index(prod)
                products[idx] = OffTargetProduct(
                    primer_pair=prod.primer_pair, contig_id=prod.contig_id,
                    start_0based=prod.start_0based, end_0based=prod.end_0based,
                    size_bp=prod.size_bp, is_expected=True,
                )
        else:
            violations.append(prod)

    # Any expected pair we never matched is also a failure
    for key, bounds in expected_products.items():
        if bounds is not None and key not in matched_expected:
            from .types import OffTargetProduct
            violations.append(OffTargetProduct(
                primer_pair=key, contig_id="<missing>",
                start_0based=0, end_0based=0, size_bp=0, is_expected=False,
            ))

    return OffTargetReport(
        sites_per_primer=sites_per_primer,
        products=products,
        passed=len(violations) == 0,
        violations=violations,
    )

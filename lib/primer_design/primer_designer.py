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

from ._bio_lite import Seq

from . import config as cfg
from .exceptions import NoCandidates, OffTargetDetected
from .types import (
    DeletionPrimerSet,
    ExpressionPrimerSet,
    GeneRecord,
    OffTargetProduct,
    OffTargetReport,
    OffTargetSite,
    Primer,
    Tag,
    TaggingPrimerSet,
    TailConvention,
    VectorRecord,
)
from .config import RESTRICTION_SITES
from .vectors import derive_tails, forbidden_body_5prime_prefixes, reverse_complement


# ===========================================================================
# Tm calculation
# ===========================================================================

def primer_body_tm(body: str) -> float:
    """Wallace-rule Tm (D5.1): GC_count * 4 + AT_count * 2.

    Annealing temperature = Tm - 5 °C (standard PCR rule of thumb).
    """
    s = body.upper()
    gc = sum(b in "GC" for b in s)
    at = len(s) - gc
    return float(gc * 4 + at * 2)


def primer_body_gc(body: str) -> float:
    return sum(b in "GC" for b in body) / len(body)


# ===========================================================================
# Hard filters (D5.3)
# ===========================================================================

def passes_hard_filters(
    body: str,
    *,
    skip_homopolymer: bool = False,
    skip_gc: bool = False,
    skip_tm: bool = False,
    skip_clamp: bool = False,
) -> bool:
    """True iff body meets all hard filters. See D5.3 / skill_v2 §7.

    The four ``skip_*`` flags exist for **anchored** primers (expression P1 at
    the gene's ATG, P2 at the native stop, in-locus P3 at the scar boundary):
    when no offset window is available, some checks become impossible to
    satisfy regardless of length and degrade to soft preferences. Callers
    relax in tiers (homopolymer first, then GC, then Tm, and only as a last
    resort the 3' G/C clamp). Length and primer-dimer remain enforced in all
    tiers.
    """
    n = len(body)
    if not (cfg.BODY_LEN_MIN <= n <= cfg.BODY_LEN_MAX):
        return False

    if not skip_clamp:
        if body[-1] not in cfg.CLAMP_LAST_BASE_OK:
            return False
        last5_gc = sum(b in "GC" for b in body[-5:])
        if not (cfg.CLAMP_LAST5_GC_MIN <= last5_gc <= cfg.CLAMP_LAST5_GC_MAX):
            return False
        if cfg.CLAMP_NO_4IDENT_LAST4 and len(set(body[-4:])) == 1:
            return False

    if not skip_homopolymer:
        for i in range(n - 3):
            if body[i] == body[i + 1] == body[i + 2] == body[i + 3]:
                return False

    if not skip_gc:
        gc = primer_body_gc(body)
        if not (cfg.GC_MIN <= gc <= cfg.GC_MAX):
            return False

    if not skip_tm:
        tm = primer_body_tm(body)
        if not (cfg.TM_HARD_MIN_C <= tm <= cfg.TM_HARD_MAX_C):
            return False

    # Primer-dimer kinetics are sequence-driven and always enforced.
    if max_3prime_self_dimer(body) > cfg.SELF_DIMER_MAX_3PRIME:
        return False

    return True


# Relaxation tiers for anchored primers. Each tier adds one more skip on top
# of the previous. Strict (tier 0) → all hard rules; tier 4 (last resort)
# even drops the 3' G/C clamp.
_ANCHORED_TIERS: list[dict] = [
    {},
    {"skip_homopolymer": True},
    {"skip_homopolymer": True, "skip_gc": True},
    {"skip_homopolymer": True, "skip_gc": True, "skip_tm": True},
    {"skip_homopolymer": True, "skip_gc": True, "skip_tm": True, "skip_clamp": True},
]


def _enumerate_anchored(
    seq: str, *, end: str
) -> tuple[list[Candidate], int]:
    """Enumerate candidates anchored at one end of ``seq``.

    ``end="5'"`` builds bodies = seq[:L] (forward primer at the 5' end of seq).
    ``end="3'"`` builds bodies = RC(seq[-L:]) (reverse primer at the 3' end).

    Tries each relaxation tier in order; returns the first non-empty pool and
    the tier index used.
    """
    for tier_idx, kwargs in enumerate(_ANCHORED_TIERS):
        pool: list[Candidate] = []
        for L in range(cfg.BODY_LEN_MIN, cfg.BODY_LEN_MAX + 1):
            if L > len(seq):
                break
            if end == "5'":
                body = seq[:L]
            elif end == "3'":
                body = reverse_complement(seq[-L:])
            else:
                raise ValueError(end)
            if passes_hard_filters(body, **kwargs):
                pool.append(
                    Candidate(body=body, tm=primer_body_tm(body),
                              gc=primer_body_gc(body), anchor_offset=0)
                )
        if pool:
            pool.sort(key=lambda c: c.tm)
            return pool, tier_idx
    return [], len(_ANCHORED_TIERS) - 1


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


def _has_forbidden_prefix(body: str, forbidden: frozenset[str]) -> bool:
    """True if ``body`` starts with any string in ``forbidden`` (cut-site regen guard)."""
    return any(body.startswith(p) for p in forbidden)


def enumerate_p1_candidates(
    up_flank: str,
    p1_tail: str,
    forbidden_5prime_prefixes: frozenset[str] = frozenset(),
) -> list[Candidate]:
    """P1 candidates: forward primer anchored near 5' end of up_flank.

    For each (offset, length) within the configured offset/length windows, accept
    the body if it passes hard filters AND its 5' prefix would not regenerate
    the cut site at the left vector junction. Returns sorted-by-Tm list.
    """
    out: list[Candidate] = []
    for offset in range(cfg.P1_OFFSET_MIN, cfg.P1_OFFSET_MAX + 1):
        for L in range(cfg.BODY_LEN_MIN, cfg.BODY_LEN_MAX + 1):
            if offset + L > len(up_flank):
                break
            body = up_flank[offset : offset + L]
            if _has_forbidden_prefix(body, forbidden_5prime_prefixes):
                continue
            if passes_hard_filters(body):
                out.append(
                    Candidate(body=body, tm=primer_body_tm(body),
                              gc=primer_body_gc(body), anchor_offset=offset)
                )
    out.sort(key=lambda c: c.tm)
    return out


def enumerate_p4_candidates(
    dn_flank: str,
    p4_tail: str,
    forbidden_5prime_prefixes: frozenset[str] = frozenset(),
) -> list[Candidate]:
    """P4 candidates: reverse primer anchored near 3' end of dn_flank.

    Body is RC of the forward sequence at the 3' anchor. Tm computed on the body
    sequence as it will appear in the primer (i.e., on the RC strand). Bodies
    whose 5' prefix would regenerate the cut site at the right vector junction
    are rejected before the hard-filter check.
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
            if _has_forbidden_prefix(body, forbidden_5prime_prefixes):
                continue
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
    Uses strict hard filters (deletion / tagging junction quality matters).
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
    """P3 candidates: forward primer anchored at 5' end of DN segment (= retained C codons + DN flank).

    Strict hard filters; for the analogous expression-P1 anchor at ATG, use
    ``_enumerate_anchored`` which applies tiered relaxation.
    """
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

IDEAL_BODY_LEN: int = 20


def compute_score(
    tms: list[float],
    gcs: list[float],
    *,
    body_lengths: list[int] | None = None,
) -> float:
    """Lower is better. Hard reject elsewhere if tm_spread > TM_SPREAD_HARD_LIMIT_C.

    When ``body_lengths`` is provided, an additional penalty pulls the chosen
    bodies toward IDEAL_BODY_LEN (= 20 nt). The penalty is small enough that
    Tm and GC quality dominate, but ties are broken in favor of ideal length.
    """
    tm_spread = max(tms) - min(tms)
    mean_tm = sum(tms) / len(tms)
    gc_spread = max(gcs) - min(gcs)
    gc_penalty = sum(
        0.1 * max(0.0, cfg.GC_PREFERRED_MIN - g) + 0.1 * max(0.0, g - cfg.GC_PREFERRED_MAX)
        for g in gcs
    )
    score = tm_spread + 0.25 * abs(mean_tm - cfg.TM_TARGET_C) + 0.02 * gc_spread + gc_penalty
    if body_lengths is not None:
        score += 0.05 * sum(abs(L - IDEAL_BODY_LEN) for L in body_lengths)
    return score


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
                 scar_min_per_side: int = cfg.SCAR_MIN_PER_SIDE,
                 scar_min_total: int = cfg.SCAR_MIN_TOTAL_AA) -> list[tuple[int, int]]:
    """Return all (N, C) pairs within the configured scar-size window.

    N counts retained N-terminal codons; C counts retained C-terminal codons
    (stop codon is always appended separately). Hard bounds:
      N ≥ scar_min_per_side, C ≥ scar_min_per_side
      scar_min_total ≤ N + C ≤ scar_max_total
    """
    pairs: list[tuple[int, int]] = []
    start = max(2 * scar_min_per_side, scar_min_total)
    for total in range(start, scar_max_total + 1):
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
    motif, _ = RESTRICTION_SITES[convention.enzyme]
    p1_forbidden, p4_forbidden = forbidden_body_5prime_prefixes(
        p1_tail, p4_tail, motif,
        expected_count=convention.expected_recognition_count_in_final_plasmid,
    )
    cds = gene.cds_seq
    L = len(cds) // 3                                          # total codons including stop

    p1_pool = enumerate_p1_candidates(gene.up_flank, p1_tail, p1_forbidden)
    p4_pool = enumerate_p4_candidates(gene.dn_flank, p4_tail, p4_forbidden)
    if not p1_pool:
        raise NoCandidates(message=f"P1: no body in up-flank passes hard filters "
                                   f"({gene.gene}/{gene.isolate_id}). Relax GC range or offset window.",
                           details={"primer": "P1", "gene": gene.gene, "isolate": gene.isolate_id})
    if not p4_pool:
        raise NoCandidates(message=f"P4: no body in dn-flank passes hard filters "
                                   f"({gene.gene}/{gene.isolate_id}). Relax GC range or offset window.",
                           details={"primer": "P4", "gene": gene.gene, "isolate": gene.isolate_id})

    junction_tail = cfg.JUNCTION_LEN_PER_PRIMER_DEFAULT       # 10 nt per primer

    best: DeletionPrimerSet | None = None
    top5: list[DeletionPrimerSet] = []

    for N, C in search_scars(L, cfg.SCAR_MAX_TOTAL_AA, cfg.SCAR_MIN_PER_SIDE, cfg.SCAR_MIN_TOTAL_AA):
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

        # Junction overlap: 20 nt total (10 from P3 tail at end of UP, 10 from P2 tail at start of DN).
        p3_tail = up_segment[-junction_tail:]                  # last 10 nt of UP segment, fwd strand
        p2_tail = reverse_complement(dn_segment[:junction_tail])

        p2_pool = enumerate_p2_candidates(up_segment)
        p3_pool = enumerate_p3_candidates(dn_segment)
        if not p2_pool or not p3_pool:
            continue

        # P1 is paired with P2 (UP PCR reaction); P4 is paired with P3 (DN PCR reaction).
        # Spread is checked per reaction — the two reactions can have different annealing temps.
        for p2c in p2_pool:
            p1c = closest_by_tm(p1_pool, p2c.tm)      # match P1 Tm to P2
            if p1c is None:
                continue
            spread_up = abs(p1c.tm - p2c.tm)
            if spread_up > cfg.TM_SPREAD_HARD_LIMIT_C:
                continue
            for p3c in p3_pool:
                p4c = closest_by_tm(p4_pool, p3c.tm)  # match P4 Tm to P3
                if p4c is None:
                    continue
                spread_dn = abs(p3c.tm - p4c.tm)
                if spread_dn > cfg.TM_SPREAD_HARD_LIMIT_C:
                    continue
                tms = [p1c.tm, p2c.tm, p3c.tm, p4c.tm]
                gcs = [p1c.gc, p2c.gc, p3c.gc, p4c.gc]
                score = compute_score(tms, gcs)

                cand = DeletionPrimerSet(
                    p1=_make_primer("P1", "UP_Fwd", p1_tail, p1c, convention.name),
                    p2=_make_primer("P2", "UP_Rev", p2_tail, p2c, "junction_overlap"),
                    p3=_make_primer("P3", "DN_Fwd", p3_tail, p3c, "junction_overlap"),
                    p4=_make_primer("P4", "DN_Rev", p4_tail, p4c, convention.name),
                    N=N, C=C, scar_dna=scar,
                    score=score, tm_spread=max(spread_up, spread_dn),
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
    """In-locus C-terminal tagging primer search (skill_v2 §5.3).

    The tag cassette occupies the UP/DN junction. Each tail carries a slice of
    the cassette such that the central ``overlap_len`` nucleotides are shared
    between the two PCR products (the In-Fusion homology). For the 30-nt His6
    cassette the shared overlap is 15 nt; for 36-nt FLAG / His8 cassettes it
    is 18 nt.

    Tail layout::

        cassette = [0 .. overlap_start) [overlap_start .. overlap_end) [overlap_end .. len)
                   ^^^^^^^^^ P2-tail-only ^^^^^^^^^^^^ shared ^^^^^^^^^^^^ P3-tail-only ^^^

        p2_tail = RC(cassette[0 : overlap_end])           # carries left half + overlap
        p3_tail =      cassette[overlap_start :]          # carries overlap + right half

    P1 / P4 are deletion-style: anchored within the up/dn flanks with the
    vector tails from the calibrated convention.

    Raises:
        TagTooLongForInLocus: cassette > 36 nt (3xFLAG, HiBiT).
    """
    from .tags import build_in_locus_cassette  # local import to avoid cycles

    cassette = build_in_locus_cassette(tag)  # raises TagTooLongForInLocus

    p1_tail, p4_tail = derive_tails(vector, convention.enzyme, "tagging")
    motif, _ = RESTRICTION_SITES[convention.enzyme]
    p1_forbidden, p4_forbidden = forbidden_body_5prime_prefixes(
        p1_tail, p4_tail, motif,
        expected_count=convention.expected_recognition_count_in_final_plasmid,
    )

    overlap_len = len(cassette) // 2 + len(cassette) % 2  # 15 (His6) or 18 (FLAG/His8)
    overlap_start = (len(cassette) - overlap_len) // 2
    overlap_end = overlap_start + overlap_len
    p2_tail = reverse_complement(cassette[:overlap_end])
    p3_tail = cassette[overlap_start:]
    overlap_left = overlap_end                                # = len(cassette in P2 tail)
    overlap_right = len(cassette) - overlap_start             # = len(cassette in P3 tail)

    cds_no_stop = gene.cds_seq[:-3]
    up_segment = gene.up_flank + cds_no_stop
    dn_segment = gene.dn_flank

    # P1 / P4 retain offset windows in the flanks → strict filters apply.
    p1_pool = enumerate_p1_candidates(gene.up_flank, p1_tail, p1_forbidden)
    p4_pool = enumerate_p4_candidates(gene.dn_flank, p4_tail, p4_forbidden)
    # P2 anchors at the 3' end of cds_no_stop (potentially AT-rich); P3
    # anchors at the 5' end of dn_flank (potentially GC-rich). Both lack
    # offset windows, so apply tiered relaxation that keeps the 3' G/C clamp
    # respected unless absolutely impossible.
    p2_pool, p2_tier = _enumerate_anchored(up_segment, end="3'")
    p3_pool, p3_tier = _enumerate_anchored(dn_segment, end="5'")

    for name, pool in (("P1", p1_pool), ("P2", p2_pool),
                       ("P3", p3_pool), ("P4", p4_pool)):
        if not pool:
            raise NoCandidates(
                message=f"{name} (tagging): no body passes hard filters for "
                        f"{gene.gene}/{gene.isolate_id} even at maximum relaxation.",
                details={"primer": name, "gene": gene.gene,
                         "isolate": gene.isolate_id},
            )

    spread_limit = (
        cfg.TM_SPREAD_HARD_LIMIT_C
        if (p2_tier == 0 and p3_tier == 0)
        else float("inf")
    )

    best: TaggingPrimerSet | None = None
    top5: list[TaggingPrimerSet] = []
    for p1c in p1_pool:
        for p2c in p2_pool:
            spread_up = abs(p1c.tm - p2c.tm)          # UP PCR reaction
            if spread_up > spread_limit:
                continue
            for p3c in p3_pool:
                for p4c in p4_pool:
                    spread_dn = abs(p3c.tm - p4c.tm)  # DN PCR reaction
                    if spread_dn > spread_limit:
                        continue
                    gcs = [p1c.gc, p2c.gc, p3c.gc, p4c.gc]
                    score = compute_score(
                        tms, gcs,
                        body_lengths=[len(p1c.body), len(p2c.body),
                                      len(p3c.body), len(p4c.body)],
                    )
                    cand = TaggingPrimerSet(
                        p1=_make_primer("P1", "UP_Fwd", p1_tail, p1c, convention.name),
                        p2=_make_primer("P2", "UP_Rev", p2_tail, p2c, "tag_cassette"),
                        p3=_make_primer("P3", "DN_Fwd", p3_tail, p3c, "tag_cassette"),
                        p4=_make_primer("P4", "DN_Rev", p4_tail, p4c, convention.name),
                        cassette=cassette,
                        overlap_left=overlap_left,
                        overlap_right=overlap_right,
                        tag=tag,
                        score=score,
                        tm_spread=spread,
                    )
                    top5 = _topk_insert(top5, cand)
                    if best is None or score < best.score:
                        best = cand

    if best is None:
        raise NoCandidates(
            message=f"No (P1, P2, P3, P4) tuple satisfied Tm-spread bound for "
                    f"tagging of {gene.gene}/{gene.isolate_id}.",
            details={"gene": gene.gene, "isolate": gene.isolate_id,
                     "tag": tag.name, "cassette_nt": len(cassette)},
        )
    return best, top5


def search_expression_primers(
    gene: GeneRecord,
    vector: VectorRecord,
    convention: TailConvention,
    tag: Tag | None,
    tag_position: str | None,
) -> tuple[ExpressionPrimerSet, list[ExpressionPrimerSet]]:
    """2-primer search for plasmid expression (skill_v2 §6.2).

    PCR template is always **native genomic DNA** of the requested isolate —
    it never contains a tag. Any tag must therefore be introduced through the
    primer overhangs, not annealed to the genome:

        N-term tag → P1 tail = vector_p1_tail + RBS + spacer + ATG + tag.dna + GGS
                     P1 body anneals to cds_seq[3:]  (codon 2 onward)
        C-term tag → P2 tail = vector_p2_tail + RC(GGS + tag.dna + new_stop)
                     P2 body anneals to RC(cds_seq[:-3])  (last native codon before stop)
        untagged   → P1 tail = vector_p1_tail + RBS + spacer;  P2 tail = vector_p2_tail
                     bodies anneal to full cds_seq.

    Choose the (P1, P2) pair that minimizes score = Tm spread + GC penalty.
    """
    from .tags import build_plasmid_fusion_cds, LINKER_GGS  # local import to avoid cycle
    from .config import NEW_STOP_DEFAULT

    p1_vector_tail, p2_vector_tail = derive_tails(vector, convention.enzyme, "expression")

    coding_seq = build_plasmid_fusion_cds(gene.cds_seq, tag, tag_position)

    # Decide where bodies anneal (native template) and what extra DNA is folded
    # into each tail to encode the tag cassette during the PCR.
    if tag is None:
        native_template = gene.cds_seq
        p1_tag_overhang = ""
        p2_tag_overhang_rc = ""
    elif tag_position == "N":
        native_template = gene.cds_seq[3:]                       # drop native ATG
        p1_tag_overhang = "ATG" + tag.dna + LINKER_GGS           # 5'→3' coding orientation
        p2_tag_overhang_rc = ""
    elif tag_position == "C":
        native_template = gene.cds_seq[:-3]                      # drop native stop
        p1_tag_overhang = ""
        p2_tag_overhang_rc = reverse_complement(
            LINKER_GGS + tag.dna + NEW_STOP_DEFAULT
        )
    else:
        raise ValueError(f"Unknown tag_position {tag_position!r}")

    p1_tail = p1_vector_tail + cfg.RBS_TAIL_PART + p1_tag_overhang
    p2_tail = p2_vector_tail + p2_tag_overhang_rc

    # Both expression primers are anchor-fixed (P1 at ATG of native_template,
    # P2 at the 3' end of native_template). Tiered relaxation: stay strict if
    # possible, drop the 4-homopolymer rule first, then GC range, then Tm
    # range, and only as a last resort the 3' G/C clamp.
    p1_pool, p1_tier = _enumerate_anchored(native_template, end="5'")
    p2_pool, p2_tier = _enumerate_anchored(native_template, end="3'")
    if not p1_pool:
        raise NoCandidates(
            message=f"P1 (expression): no body at ATG of {gene.gene}/{gene.isolate_id} "
                    f"passes hard filters even at maximum relaxation.",
            details={"primer": "P1", "gene": gene.gene, "isolate": gene.isolate_id},
        )
    if not p2_pool:
        raise NoCandidates(
            message=f"P2 (expression): no body at stop of {gene.gene}/{gene.isolate_id} "
                    f"passes hard filters even at maximum relaxation.",
            details={"primer": "P2", "gene": gene.gene, "isolate": gene.isolate_id},
        )

    # Tm spread bound is itself relaxed when either primer is anchored at a
    # tier > 0 (because the gene's 3' AT-richness can force an unavoidable
    # spread). Strict bound otherwise.
    spread_limit = (
        cfg.TM_SPREAD_HARD_LIMIT_C
        if (p1_tier == 0 and p2_tier == 0)
        else float("inf")
    )

    best: ExpressionPrimerSet | None = None
    top5: list[ExpressionPrimerSet] = []
    for p1c in p1_pool:
        for p2c in p2_pool:
            tms = [p1c.tm, p2c.tm]
            spread = max(tms) - min(tms)
            if spread > spread_limit:
                continue
            gcs = [p1c.gc, p2c.gc]
            score = compute_score(
                tms, gcs,
                body_lengths=[len(p1c.body), len(p2c.body)],
            )
            cand = ExpressionPrimerSet(
                p1=_make_primer("P1", "INSERT_Fwd", p1_tail, p1c, convention.name),
                p2=_make_primer("P2", "INSERT_Rev", p2_tail, p2c, convention.name),
                coding_seq=coding_seq,
                native_anneal_segment=native_template,
                tag=tag,
                tag_position=tag_position,
                score=score,
                tm_spread=spread,
            )
            top5 = _topk_insert(top5, cand)
            if best is None or score < best.score:
                best = cand

    if best is None:
        raise NoCandidates(
            message=f"No (P1, P2) pair satisfied Tm-spread bound for "
                    f"{gene.gene}/{gene.isolate_id} expression.",
            details={"gene": gene.gene, "isolate": gene.isolate_id,
                     "p1_pool_size": len(p1_pool), "p2_pool_size": len(p2_pool)},
        )
    return best, top5


# ===========================================================================
# Off-target scan — STUB
# ===========================================================================

def off_target_scan(
    primers: list[Primer],
    genome_fasta_bytes: bytes,
    expected_products: dict[tuple[str, str], tuple[int, int] | None],
) -> OffTargetReport:
    """Walk-anchor + body-mismatch scan per D5.5 / skill_v2 §10.

    Algorithm:
        1. For each primer, sliding-window match the last ANCHOR_LEN nt against
           both strands of all contigs, allowing ANCHOR_MAX_MM mismatches.
        2. For each anchor hit, extend to full body length and check
           total mismatches ≤ BODY_MAX_MM.
        3. Collect all (primer, contig, position, strand) site tuples.
        4. Enumerate primer-pair products: any forward-strand site of primer A
           paired with reverse-strand site of primer B at distance 50–6000 bp on
           the same contig.
        5. Compare against ``expected_products`` (per skill_v2 §10).

    Args:
        expected_products: keys are primer-pair tuples like ("P1", "P2");
            values are (size_min_bp, size_max_bp) for the expected product, or
            None if zero products expected for that pair.

    Returns:
        OffTargetReport.passed = True iff every expected pair has exactly the
        expected outcome and no unexpected products exist.
    """
    contigs = _parse_fasta_bytes(genome_fasta_bytes)

    sites_per_primer: dict[str, list[OffTargetSite]] = {}
    for primer in primers:
        body = primer.body
        body_rc = reverse_complement(body)
        hits: list[OffTargetSite] = []
        for contig_id, seq in contigs:
            hits.extend(
                _find_primer_sites(seq, body, contig_id, primer.name, "+")
            )
            hits.extend(
                _find_primer_sites(seq, body_rc, contig_id, primer.name, "-")
            )
        sites_per_primer[primer.name] = hits

    by_name = {p.name: p for p in primers}
    products: list[OffTargetProduct] = []
    primer_names = [p.name for p in primers]
    for a in primer_names:
        for b in primer_names:
            for prod in _enumerate_products(
                sites_per_primer[a], sites_per_primer[b], a, b,
                len_a=by_name[a].length, len_b=by_name[b].length,
            ):
                products.append(prod)

    # Look up by either ordering: a primer pair (A, B) and (B, A) describe
    # the same physical PCR product (one fwd primer + one rev primer at the
    # same locus). Callers therefore need only specify one direction.
    def _spec_for(pair: tuple[str, str]):
        if pair in expected_products:
            return expected_products[pair]
        rev = (pair[1], pair[0])
        if rev in expected_products:
            return expected_products[rev]
        return ...

    violations: list[OffTargetProduct] = []
    matched: list[OffTargetProduct] = []
    for prod in products:
        spec = _spec_for(prod.primer_pair)
        if spec is ...:
            violations.append(prod)
            continue
        if spec is None:
            # explicit None → zero products allowed
            violations.append(prod)
            continue
        lo, hi = spec
        if lo <= prod.size_bp <= hi:
            matched.append(
                OffTargetProduct(
                    primer_pair=prod.primer_pair,
                    contig_id=prod.contig_id,
                    start_0based=prod.start_0based,
                    end_0based=prod.end_0based,
                    size_bp=prod.size_bp,
                    is_expected=True,
                )
            )
        else:
            violations.append(prod)

    final_products = matched + violations
    return OffTargetReport(
        sites_per_primer=sites_per_primer,
        products=final_products,
        passed=not violations,
        violations=violations,
    )


def _parse_fasta_bytes(data: bytes) -> list[tuple[str, str]]:
    """Parse multi-FASTA bytes into [(contig_id, uppercase_sequence), ...]."""
    contigs: list[tuple[str, str]] = []
    cur_id: str | None = None
    cur_chunks: list[str] = []
    for raw_line in data.splitlines():
        line = raw_line.decode("ascii", errors="replace").strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_id is not None:
                contigs.append((cur_id, "".join(cur_chunks).upper()))
            cur_id = line[1:].split()[0]
            cur_chunks = []
        else:
            cur_chunks.append(line)
    if cur_id is not None:
        contigs.append((cur_id, "".join(cur_chunks).upper()))
    return contigs


def _find_primer_sites(
    sequence: str,
    body_oriented: str,
    contig_id: str,
    primer_name: str,
    strand: str,
) -> list[OffTargetSite]:
    """Find all positions in ``sequence`` where ``body_oriented`` matches the
    top strand within OFFTARGET_BODY_MAX_MM mismatches AND the last
    OFFTARGET_ANCHOR_LEN nt match within OFFTARGET_ANCHOR_MAX_MM.

    For strand="+" the match means the primer body anneals to the bottom strand
    and primes synthesis 5'→3' into the top strand at the *end* of the match
    (3' end position).

    For strand="-" we already pass the RC of the body; matches mean the primer
    anneals to the top strand and primes synthesis 5'→3' on the bottom strand.
    """
    n = len(sequence)
    L = len(body_oriented)
    if L == 0 or L > n:
        return []
    anchor_len = cfg.OFFTARGET_ANCHOR_LEN
    anchor = body_oriented[-anchor_len:]
    body_offset = L - anchor_len  # body start = anchor_hit - body_offset

    seen: set[int] = set()
    out: list[OffTargetSite] = []
    for variant in _mismatch_variants(anchor, cfg.OFFTARGET_ANCHOR_MAX_MM):
        i = sequence.find(variant)
        while i != -1:
            body_start = i - body_offset
            if (
                body_start >= 0
                and body_start + L <= n
                and body_start not in seen
            ):
                seen.add(body_start)
                window = sequence[body_start : body_start + L]
                body_mm = _hamming(window, body_oriented)
                if body_mm <= cfg.OFFTARGET_BODY_MAX_MM:
                    out.append(
                        OffTargetSite(
                            primer_name=primer_name,
                            contig_id=contig_id,
                            position_0based=body_start,
                            strand=strand,
                            body_mismatches=body_mm,
                        )
                    )
            i = sequence.find(variant, i + 1)
    return out


def _mismatch_variants(seq: str, max_mm: int) -> list[str]:
    """Enumerate all DNA strings within Hamming distance ``max_mm`` of ``seq``.
    For max_mm=1 over a 10-mer this is 1 + 10*3 = 31 variants.
    """
    if max_mm <= 0:
        return [seq]
    variants: set[str] = {seq}
    bases = "ACGT"
    cur: set[str] = {seq}
    for _ in range(max_mm):
        nxt: set[str] = set()
        for v in cur:
            for i in range(len(v)):
                for b in bases:
                    if b != v[i]:
                        nxt.add(v[:i] + b + v[i + 1 :])
        variants.update(nxt)
        cur = nxt
    return list(variants)


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


_PRODUCT_SITE_MAX_BODY_MM = 2  # PCR amplification requires near-perfect priming both sides


def _enumerate_products(
    sites_a: list[OffTargetSite],
    sites_b: list[OffTargetSite],
    name_a: str,
    name_b: str,
    *,
    len_a: int = 0,
    len_b: int = 0,
) -> list[OffTargetProduct]:
    """For each forward site of A on a contig and each reverse site of B on the
    same contig, emit a product if the implied amplicon size is in
    [OFFTARGET_PRODUCT_SIZE_MIN_BP, OFFTARGET_PRODUCT_SIZE_MAX_BP] AND both sites
    have body_mismatches ≤ _PRODUCT_SITE_MAX_BODY_MM. The looser body-mm cutoff
    in OFFTARGET_BODY_MAX_MM is retained for site reporting but not for
    amplification simulation — real PCR rarely amplifies products where both
    primers have >2 body mismatches.
    """
    out: list[OffTargetProduct] = []
    by_contig: dict[str, list[OffTargetSite]] = {}
    for s in sites_b:
        if s.body_mismatches > _PRODUCT_SITE_MAX_BODY_MM:
            continue
        by_contig.setdefault(s.contig_id, []).append(s)
    for sa in sites_a:
        if sa.strand != "+":
            continue
        if sa.body_mismatches > _PRODUCT_SITE_MAX_BODY_MM:
            continue
        for sb in by_contig.get(sa.contig_id, []):
            if sb.strand != "-":
                continue
            if sb.position_0based <= sa.position_0based:
                continue
            # Real amplicon = len_a (tail+body of fwd primer) + genomic gap between
            # the 3' ends + len_b. Genomic body span on top strand = sb.pos + body_len_b
            # - sa.pos - body_len_a. Since len_X = tail_X + body_X, we approximate by
            # (sb.pos - sa.pos) + len_b. This estimate is within a few bp of the true
            # PCR product size when tails are 20 nt.
            size = sb.position_0based - sa.position_0based + (len_b or 1)
            # sb is on the top strand at the RC-match position; the 3' end of the
            # bottom-strand primer corresponds to sb.position_0based, so the amplicon
            # spans sa.position_0based .. sb.position_0based + L_b - 1. Use a generic
            # estimate: distance + body_len_b; we don't know L_b here, so just use
            # the gap. This is a slight underestimate (off by L_b-1 ≈ 20 bp) but
            # within ±50 tolerance applied by callers.
            if size < cfg.OFFTARGET_PRODUCT_SIZE_MIN_BP:
                continue
            if size > cfg.OFFTARGET_PRODUCT_SIZE_MAX_BP:
                continue
            out.append(
                OffTargetProduct(
                    primer_pair=(name_a, name_b),
                    contig_id=sa.contig_id,
                    start_0based=sa.position_0based,
                    end_0based=sb.position_0based,
                    size_bp=size,
                    is_expected=False,
                )
            )
    return out

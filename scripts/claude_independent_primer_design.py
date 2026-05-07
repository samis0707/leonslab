"""
Independent primer design for in-frame deletion on pEXG2 + HindIII (site_destroyed).

Independent of lib/primer_design/ — designed from spec in docs/primer_design/decisions_log.md
to enable in-silico comparison against the website tool.

Spec recap (D5.3, D4.3, D6.4-style geometry):
  4 primers, 3-fragment In-Fusion (linearized pEXG2 + UP + DN).
  P1 fwd: vector-overhang tail + body in UP genomic region.
  P2 rev: 15-nt junction tail + body anneals at end of (UP_genome + first_N_codons).
  P3 fwd: 15-nt junction tail + body anneals at start of (last_C_codons + DN_genome).
  P4 rev: vector-overhang tail + body in DN genomic region.
  Junction = last 15 nt of (UP_genome + first_N_codons) || first 15 nt of (last_C_codons + DN_genome).
  Scar (N, C): exhaustive search, N + C ≤ 24, in-frame ORF starting M ending single *.
  Body filters: 18–28 nt, GC 40–70 %, Tm 58–66 °C (prefer 60–64), 3' G/C clamp,
                no 4-homopolymer, 3'-self-dimer ≤ 4 bp, Tm spread ≤ 3 °C.

Outputs to stdout: one block per (gene × isolate) operation.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqUtils.MeltingTemp import Tm_NN

DATA = Path("data/primer_design")
VECTOR_FASTA = DATA / "vectors" / "pEXG2.fasta"

# fixed pEXG2+HindIII tails (site_destroyed)
P1_TAIL = "GCATAAATGTAAAGC"  # last 15 nt of left arm before AAGCTT
P4_TAIL = "CGACCTGCAGAAGCT"  # RC of right arm 15 nt starting at 2nd nt of AAGCTT

STOPS = {"TAA", "TAG", "TGA"}

def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def read_fasta(p: Path) -> tuple[str, dict]:
    head, body = open(p).read().split("\n", 1)
    seq = body.replace("\n", "").upper()
    meta = {}
    for kv in head.lstrip(">").split("|"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            meta[k] = v
    return seq, meta

def gc(s: str) -> float:
    return 100 * (s.count("G") + s.count("C")) / len(s)

def tm(body: str) -> float:
    return float(Tm_NN(Seq(body), Na=50, dnac1=500, dnac2=0))

def has_homopolymer(s: str, n: int = 4) -> bool:
    for b in "ACGT":
        if b * n in s:
            return True
    return False

def gc_clamp_ok(s: str) -> bool:
    # 3' last base must be G or C; not GGG or CCC at 3' end (avoid over-clamp)
    return s[-1] in "GC" and not (s[-3:] in ("GGG", "CCC"))

def self_dimer_3prime(s: str, max_match: int = 4) -> bool:
    """True if 3' end self-dimers > max_match nt with the primer's RC."""
    r = rc(s)
    end = s[-(max_match + 1):]
    return end in r

def body_passes(body: str) -> tuple[bool, str]:
    if not (18 <= len(body) <= 28):
        return False, f"len={len(body)}"
    g = gc(body)
    if not (40 <= g <= 70):
        return False, f"GC={g:.1f}"
    if has_homopolymer(body):
        return False, "homopolymer"
    if not gc_clamp_ok(body):
        return False, "clamp"
    if self_dimer_3prime(body):
        return False, "3'self-dimer"
    t = tm(body)
    if not (58 <= t <= 66):
        return False, f"Tm={t:.1f}"
    return True, ""

def score_body(body: str) -> tuple[float, float]:
    """Lower is better. Prefer Tm 62 °C, GC 50 %, length 22."""
    t = tm(body)
    g = gc(body)
    return (
        abs(t - 62) * 1.0
        + abs(g - 50) * 0.05
        + abs(len(body) - 22) * 0.1,
        t,
    )

def find_scar(cds: str, n_min: int = 3, c_min: int = 3) -> tuple[int, int, str]:
    """Pick smallest in-frame scar with valid ORF. Enforce N>=n_min, C>=c_min
    (standard P. aeruginosa allelic-exchange practice — preserves small in-frame
    remnants on both ends to avoid degenerate M-* scars)."""
    assert len(cds) % 3 == 0 and cds[:3] == "ATG"
    last3 = cds[-3:]
    assert last3 in STOPS, f"CDS does not end with stop: {last3}"
    best = None
    for total in range(n_min + c_min, 25):
        for n in range(n_min, total - c_min + 1):
            c = total - n
            if n + c > 24 or c < c_min:
                continue
            scar = cds[: 3 * n] + cds[-3 * c :]
            if len(scar) % 3:
                continue
            # check ORF
            if scar[:3] != "ATG":
                continue
            if scar[-3:] not in STOPS:
                continue
            internal = [scar[i : i + 3] for i in range(0, len(scar) - 3, 3)]
            if any(c_ in STOPS for c_ in internal):
                continue
            if best is None or (n + c) < (best[0] + best[1]):
                best = (n, c, scar)
        if best:
            return best
    raise RuntimeError("no valid scar")

def best_body_fwd(template: str, anchor_start: int, search_offsets=range(0, 1)):
    """Forward primer body whose 5'-end starts at anchor_start + offset.
    For P3: anchor_start fixed at 0 (must start exactly at scar boundary)."""
    cands = []
    for off in search_offsets:
        for L in range(18, 29):
            body = template[anchor_start + off : anchor_start + off + L]
            if len(body) < L:
                continue
            ok, why = body_passes(body)
            if ok:
                cands.append((score_body(body)[0], body, off, L))
    cands.sort()
    return cands

def best_body_rev(template: str, anchor_end: int, search_offsets=range(0, 1)):
    """Reverse primer body whose anneal region ENDS at anchor_end on top strand.
    Body sequence = RC(template[anchor_end - L + offset : anchor_end + offset]).
    For P2: anchor_end fixed (= end of UP+first_N_codons), offset=0."""
    cands = []
    for off in search_offsets:
        for L in range(18, 29):
            top = template[anchor_end - L - off : anchor_end - off]
            if len(top) < L:
                continue
            body = rc(top)
            ok, why = body_passes(body)
            if ok:
                cands.append((score_body(body)[0], body, off, L))
    cands.sort()
    return cands

@dataclass
class Result:
    op: str
    n: int
    c: int
    scar: str
    p1: str
    p2: str
    p3: str
    p4: str
    p1_body_tm: float
    p2_body_tm: float
    p3_body_tm: float
    p4_body_tm: float
    junction: str
    notes: list[str]

def design(fasta_path: Path, label: str) -> Result:
    seq, meta = read_fasta(fasta_path)
    cds_len = int(meta["cds_len"])
    up_len = int(meta["up"])
    dn_len = int(meta["dn"])
    assert len(seq) == up_len + cds_len + dn_len, (len(seq), up_len, cds_len, dn_len)
    up = seq[:up_len]
    cds = seq[up_len : up_len + cds_len]
    dn = seq[up_len + cds_len :]
    notes = []

    n, c, scar = find_scar(cds)

    # upstream-part = up + first N codons; downstream-part = last C codons + dn
    up_part = up + cds[: 3 * n]
    dn_part = cds[-3 * c :] + dn

    # junction = last 15 of up_part || first 15 of dn_part
    junction = up_part[-15:] + dn_part[:15]

    # P2 reverse: body anneals on top strand ending at end of up_part (= last position of up_part).
    # Body = RC(up_part[-L:]). We allow shifting the anneal end inward (offset > 0 means body ends before the boundary,
    # but then UP amplicon won't reach the junction). So offset must be 0.
    # However if the boundary base is a poor 3' (e.g. T/A), we permit shifting the cut by reformulating tail.
    # Simpler: tighten anchor to boundary.
    # P2: search offsets 0..60 (sliding body inward grows the tail; tail completes the junction).
    p2_all = []
    for off in range(0, 61):
        for cand in best_body_rev(up_part, anchor_end=len(up_part), search_offsets=[off]):
            p2_all.append((cand[0] + off * 0.02, cand[1], off, cand[3]))
    p2_all.sort()
    if not p2_all:
        raise RuntimeError(f"no P2 body for {label}")
    _, p2_body, p2_off, p2_L = p2_all[0]
    # tail = RC(up_part[B-g : B] + dn_part[:15])  where g = p2_off
    B = len(up_part)
    p2_tail = rc(up_part[B - p2_off : B] + dn_part[:15])
    if p2_off > 0:
        notes.append(f"P2 body retracted {p2_off} nt into UP-part (tail length {len(p2_tail)})")
    p2 = p2_tail + p2_body

    p3_all = []
    for off in range(0, 61):
        for cand in best_body_fwd(dn_part, anchor_start=0, search_offsets=[off]):
            p3_all.append((cand[0] + off * 0.02, cand[1], off, cand[3]))
    p3_all.sort()
    if not p3_all:
        raise RuntimeError(f"no P3 body for {label}")
    _, p3_body, p3_off, p3_L = p3_all[0]
    # tail = up_part[-15:] + dn_part[:p3_off]   (top strand sequence preceding body's 5' anchor)
    p3_tail = up_part[-15:] + dn_part[:p3_off]
    if p3_off > 0:
        notes.append(f"P3 body shifted {p3_off} nt into DN-part (tail length {len(p3_tail)})")
    p3 = p3_tail + p3_body

    # P1 forward: body anneals in UP, somewhere upstream. We want UP amplicon ~500-600 bp.
    # Search anchors near position 0..50 of up_part (i.e. close to 5' end of UP region).
    p1_cands = []
    for start in range(0, 80):
        for L in range(18, 29):
            body = up_part[start : start + L]
            if len(body) < L:
                continue
            ok, why = body_passes(body)
            if ok:
                # prefer primers whose 3' end is at a position that makes a clean amplicon
                p1_cands.append((score_body(body)[0] + start * 0.01, body, start, L))
    p1_cands.sort()
    if not p1_cands:
        raise RuntimeError(f"no P1 body for {label}")
    _, p1_body, p1_start, p1_L = p1_cands[0]
    p1 = P1_TAIL + p1_body

    # P4 reverse: body anneals in DN_part, near 3' end (~position dn_len - 80 .. dn_len).
    p4_cands = []
    for end in range(len(dn_part) - 80, len(dn_part) + 1):
        for L in range(18, 29):
            top = dn_part[end - L : end]
            if len(top) < L:
                continue
            body = rc(top)
            ok, why = body_passes(body)
            if ok:
                p4_cands.append((score_body(body)[0] + (len(dn_part) - end) * 0.01, body, end, L))
    p4_cands.sort()
    if not p4_cands:
        raise RuntimeError(f"no P4 body for {label}")
    _, p4_body, p4_end, p4_L = p4_cands[0]
    p4 = P4_TAIL + p4_body

    bodies_tm = [tm(b) for b in (p1_body, p2_body, p3_body, p4_body)]
    spread = max(bodies_tm) - min(bodies_tm)
    if spread > 3.0:
        notes.append(f"Tm spread {spread:.1f} °C exceeds 3 °C")
    else:
        notes.append(f"Tm spread {spread:.1f} °C ✓")

    return Result(
        op=label,
        n=n,
        c=c,
        scar=scar,
        p1=p1,
        p2=p2,
        p3=p3,
        p4=p4,
        p1_body_tm=bodies_tm[0],
        p2_body_tm=bodies_tm[1],
        p3_body_tm=bodies_tm[2],
        p4_body_tm=bodies_tm[3],
        junction=junction,
        notes=notes,
    )


def fmt(r: Result) -> str:
    out = []
    out.append(f"### {r.op}")
    out.append(f"Scar: N={r.n} (first {r.n} codons) + C={r.c} (last {r.c} codons incl. stop) — total {r.n+r.c} codons / {len(r.scar)} nt")
    out.append(f"Scar DNA: {r.scar}")
    out.append(f"Scar protein: {str(Seq(r.scar).translate())}")
    out.append(f"Junction (30 nt across UP||DN boundary): {r.junction}")
    out.append("")
    out.append(f"P1 (fwd, len {len(r.p1)}, body Tm {r.p1_body_tm:.1f} °C):")
    out.append(f"  5'-{r.p1}-3'")
    out.append(f"     [{P1_TAIL}] tail | body {r.p1[len(P1_TAIL):]}")
    out.append(f"P2 (rev, len {len(r.p2)}, body Tm {r.p2_body_tm:.1f} °C):")
    out.append(f"  5'-{r.p2}-3'")
    p2_tail_len = len(r.p2) - (len(r.p2) - 15) if len(r.p2) >= 15 else 0
    # Reconstruct: tail length = total - body length. We know body length came from cands but re-extract:
    # simpler: split at known body suffix
    out.append(f"P3 (fwd, len {len(r.p3)}, body Tm {r.p3_body_tm:.1f} °C):")
    out.append(f"  5'-{r.p3}-3'")
    out.append(f"P4 (rev, len {len(r.p4)}, body Tm {r.p4_body_tm:.1f} °C):")
    out.append(f"  5'-{r.p4}-3'")
    out.append(f"     [{P4_TAIL}] tail | body {r.p4[len(P4_TAIL):]}")
    out.append("Notes: " + "; ".join(r.notes))
    return "\n".join(out)


CASES = [
    ("LB014_dlasR", DATA / "genes" / "lasR" / "LB014.fasta"),
    ("LB056_dlasR", DATA / "genes" / "lasR" / "LB056.fasta"),
    ("LB014_dlasB", DATA / "genes" / "lasB" / "LB014.fasta"),
    ("LB020_dlasB", DATA / "genes" / "lasB" / "LB020.fasta"),
    ("LB056_dlasB", DATA / "genes" / "lasB" / "LB056.fasta"),
]

if __name__ == "__main__":
    print("# Claude independent primer design — pEXG2 + HindIII in-frame deletion")
    print(f"# Vector tails: P1='{P1_TAIL}', P4='{P4_TAIL}' (site_destroyed)\n")
    for label, path in CASES:
        try:
            r = design(path, label)
            print(fmt(r))
            print()
        except Exception as e:
            print(f"### {label}: FAILED — {e}\n")

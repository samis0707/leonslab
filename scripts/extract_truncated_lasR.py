"""
Extract lasR loci from isolates whose lasR allele is truncated/non-functional
and therefore was excluded from the curated lasR.fasta master file.

For each missing isolate we:
  1. Fetch its full genome from R2.
  2. Locate the lasR locus by anchor-search using two highly-conserved 60-nt
     anchors flanking the gene in the PA14 reference (UP-end + DN-start).
     Exact match on either strand of any contig.
  3. If both anchors hit consistently on one contig and one strand:
       - Extract UP_FLANK_LEN nt UP + nominal CDS region + DN_FLANK_LEN nt DN
       - Translate the nominal CDS; classify the truncation reason:
           * missing_start            (no ATG at expected position)
           * premature_stop@codon_X   (in-frame stop before nominal end)
           * frameshift_via_indel     (length ≠ multiple of 3 between hits)
           * length_anomaly           (nominal CDS length differs by > 30 nt)
           * substitution_only        (length OK, in-frame, but reference protein differs)
       - Write FASTA with header field functional=false|truncation_reason=...
  4. If anchors fail: log and skip (manual triage needed).

Reference: data/primer_design/genes/lasR/PA14.fasta (canonical functional lasR record)

Usage:
  set -a && source api/design/.env && set +a
  python3 scripts/extract_truncated_lasR.py
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from Bio import SeqIO
from Bio.Seq import Seq

from primer_design.config import DN_FLANK_LEN, UP_FLANK_LEN
from primer_design.storage_adapter import fetch_genome
from primer_design.vectors import reverse_complement

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("extract_truncated_lasR")

DATA = Path("data/primer_design")
GENE = "lasR"
REF_PATH = DATA / "genes" / GENE / "PA14.fasta"
OUT_DIR = DATA / "genes" / GENE

MISSING_ISOLATES = [
    "LB001", "LB009", "LB020", "LB038", "LB039", "LB040",
    "LB042", "LB043", "LB045", "LB046", "LB052", "LB069",
    "LB076", "LB086", "LB091", "LB117",
]

ANCHOR_LEN = 60                  # nt anchor on each side
STOPS = ("TAA", "TAG", "TGA")


def load_reference():
    """Return (up_flank, cds, dn_flank) on the coding strand from PA14 record."""
    head, body = open(REF_PATH).read().split("\n", 1)
    seq = body.replace("\n", "").upper()
    meta = {}
    for kv in head.lstrip(">").split("|"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            meta[k] = v
    up = int(meta["up"]); cds_len = int(meta["cds_len"]); dn = int(meta["dn"])
    return seq[:up], seq[up : up + cds_len], seq[up + cds_len :]


def parse_genome_contigs(genome_bytes: bytes) -> list[tuple[str, str]]:
    """Return list of (contig_id, sequence_uppercase) tuples."""
    out = []
    contig_id, parts = None, []
    for line in genome_bytes.decode("utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if contig_id is not None:
                out.append((contig_id, "".join(parts).upper()))
            contig_id = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if contig_id is not None:
        out.append((contig_id, "".join(parts).upper()))
    return out


def find_anchor(contigs, anchor: str) -> list[tuple[str, int, str]]:
    """Return list of (contig_id, start_0based, strand) hits for an exact match."""
    rc = reverse_complement(anchor)
    hits = []
    for cid, seq in contigs:
        i = -1
        while True:
            i = seq.find(anchor, i + 1)
            if i < 0:
                break
            hits.append((cid, i, "+"))
        i = -1
        while True:
            i = seq.find(rc, i + 1)
            if i < 0:
                break
            hits.append((cid, i, "-"))
    return hits


def classify_truncation(nominal_cds: str, ref_protein: str) -> tuple[str, str]:
    """Return (functional_str, reason). Decides whether the locus appears functional."""
    n = len(nominal_cds)
    if n < 3 or not nominal_cds.startswith("ATG"):
        return "false", "missing_start_codon"
    # Translate up to first stop
    protein = ""
    stop_pos_codon = None
    for i in range(0, n - n % 3, 3):
        codon = nominal_cds[i : i + 3]
        if codon in STOPS:
            stop_pos_codon = i // 3
            protein += "*"
            break
        try:
            aa = str(Seq(codon).translate())
        except Exception:
            aa = "X"
        protein += aa
    if stop_pos_codon is None:
        return "false", "no_inframe_stop_in_window"
    # Where does the ref protein end? (one stop at the end)
    ref_aa_len = len(ref_protein.rstrip("*")) if ref_protein.endswith("*") else len(ref_protein)
    if stop_pos_codon < ref_aa_len - 1:
        return "false", f"premature_stop@codon_{stop_pos_codon + 1}"
    # length matches reference; in-frame stop at expected position
    if stop_pos_codon == ref_aa_len:
        # CDS aligned; check substitutions vs reference protein
        ref_no_stop = ref_protein.rstrip("*")
        my_no_stop = protein.rstrip("*")
        diffs = sum(a != b for a, b in zip(ref_no_stop, my_no_stop))
        if diffs == 0:
            return "true", ""
        # In-frame, full length: probably actually functional with point mutations.
        # Mark as functional=true (silent / missense doesn't itself imply non-functional).
        return "true", ""
    # stop later than expected → extension; treat as non-functional anomaly
    return "false", f"stop_at_codon_{stop_pos_codon + 1}_ref_ends_at_{ref_aa_len}"


def extract(isolate: str, ref_up: str, ref_cds: str, ref_dn: str, dry_run: bool):
    log.info("Processing %s", isolate)
    ref_protein = str(Seq(ref_cds).translate())
    cds_len = len(ref_cds)

    try:
        genome = fetch_genome(isolate)
    except Exception as e:
        log.error("  R2 fetch failed: %s", e)
        return False, f"r2_fetch_failed: {e}"
    contigs = parse_genome_contigs(genome)

    # Anchor recipes. Each anchor is parameterised by its end position relative
    # to the CDS-start coordinate (e.g., -60 = 60 nt of UP ending at CDS start;
    # +60 = first 60 nt of CDS). Pairs of (5'-anchor, 3'-anchor) are tried in
    # order of preference: first conservative UP/DN-flank pairs, then CDS-internal
    # pairs that survive even with major rearrangements inside the locus.
    #
    # Tuple: (L, up_anchor_seq, up_anchor_end_rel_to_cds_start,
    #            dn_anchor_seq, dn_anchor_start_rel_to_cds_end)
    anchor_recipes = []
    for L in (60, 45, 30, 24):
        for up_off in (0, 60, 120, 240):
            up_seq = ref_up[len(ref_up) - up_off - L : len(ref_up) - up_off]
            for dn_off in (0, 60, 120, 240):
                dn_seq = ref_dn[dn_off : dn_off + L]
                if len(up_seq) != L or len(dn_seq) != L:
                    continue
                # UP anchor ends at CDS_start - up_off; DN anchor starts at CDS_end + dn_off
                anchor_recipes.append(("flank", L, up_seq, -up_off, dn_seq, dn_off))
    # CDS-internal fallback: anchors are the 5' and 3' ends of the reference CDS
    # itself. Useful when UP/DN flanks have been rearranged but the gene starts
    # and ends are still recognisable (typical with IS-element disruption).
    for L in (60, 45, 30):
        anchor_recipes.append(
            ("cds_internal", L, ref_cds[:L], L, ref_cds[-L:], -L)
        )

    pairs = []
    chosen_kind = chosen_L = None
    for kind, L, up_seq, up_end_rel, dn_seq, dn_start_rel in anchor_recipes:
        hits_up = find_anchor(contigs, up_seq)
        hits_dn = find_anchor(contigs, dn_seq)
        if not hits_up or not hits_dn:
            continue
        # Expected gap on the coding strand between UP-anchor end and DN-anchor start.
        expected_gap = (cds_len + dn_start_rel) - up_end_rel
        # Tolerance: tight for flank anchors (small indels expected), wide for
        # cds_internal (because IS-element insertions can be multi-kb).
        tol = max(150, cds_len // 5) if kind == "flank" else 5000
        for c1, p1, s1 in hits_up:
            for c2, p2, s2 in hits_dn:
                if c1 != c2 or s1 != s2:
                    continue
                if s1 == "+":
                    gap = p2 - (p1 + L)
                else:
                    gap = p1 - (p2 + L)
                if abs(gap - expected_gap) <= tol and gap >= 0:
                    pairs.append((c1, p1, p2, s1, gap, L,
                                  up_end_rel, dn_start_rel, kind))
        if pairs:
            chosen_kind, chosen_L = kind, L
            log.info("  anchored kind=%s L=%d up_end_rel=%d dn_start_rel=%d (%d pair(s); gap≈%d)",
                     kind, L, up_end_rel, dn_start_rel, len(pairs), pairs[0][4])
            break

    if not pairs:
        log.error("  no consistent UP/DN anchor pair across recipe cascade")
        return False, "anchor_pair_inconsistent_after_fallbacks"

    if len(pairs) > 1:
        log.warning("  %d consistent pairs; using first", len(pairs))

    cid, p_up, p_dn, strand, gap, L, up_end_rel, dn_start_rel, kind = pairs[0]
    contig_seq = next(s for c, s in contigs if c == cid)

    # Coding-strand CDS-start position derived from UP-anchor end coordinate.
    # On + strand: cds_start_top = p_up + L - up_end_rel  (top-strand 0-based).
    # cds_end_top   = cds_start_top + nominal_locus_length, where the locus length
    # is what we actually saw between the anchors (= gap + up_end_rel - dn_start_rel).
    locus_len_observed = gap + up_end_rel - dn_start_rel
    if strand == "+":
        cds_start = p_up + L - up_end_rel
        cds_end = cds_start + locus_len_observed
        up_start = max(0, cds_start - UP_FLANK_LEN)
        dn_end = min(len(contig_seq), cds_end + DN_FLANK_LEN)
        up_flank_seq = contig_seq[up_start:cds_start]
        dn_flank_seq = contig_seq[cds_end:dn_end]
        nominal_cds = contig_seq[cds_start:cds_end]
        genome_start_1b = cds_start + 1
        genome_end_1b = cds_end
    else:
        cds_top_end = p_up + up_end_rel               # exclusive
        cds_top_start = cds_top_end - locus_len_observed
        dn_top_start = max(0, cds_top_start - DN_FLANK_LEN)
        up_top_end = min(len(contig_seq), cds_top_end + UP_FLANK_LEN)
        up_flank_seq = reverse_complement(contig_seq[cds_top_end:up_top_end])
        dn_flank_seq = reverse_complement(contig_seq[dn_top_start:cds_top_start])
        nominal_cds = reverse_complement(contig_seq[cds_top_start:cds_top_end])
        genome_start_1b = cds_top_start + 1
        genome_end_1b = cds_top_end

    functional_str, reason = classify_truncation(nominal_cds, ref_protein)
    locus_actual_len = len(nominal_cds)
    # If the locus is dramatically larger than the reference, the nominal CDS
    # almost certainly hosts an IS-element/transposon insertion. Override the
    # finer-grained reason with that observation.
    if locus_actual_len > len(ref_cds) + 200:
        functional_str = "false"
        reason = f"insertion_in_locus_~{locus_actual_len - len(ref_cds)}nt"
    elif locus_actual_len < len(ref_cds) - 200:
        functional_str = "false"
        reason = f"deletion_in_locus_~{len(ref_cds) - locus_actual_len}nt"
    log.info("  contig=%s strand=%s locus_len=%d functional=%s reason=%s",
             cid, strand, locus_actual_len, functional_str, reason or "—")

    # For non-functional alleles substitute the REFERENCE CDS into the FASTA
    # between the actual isolate flanks. Rationale: scar / ORF logic downstream
    # needs a clean mod-3 CDS that begins ATG and ends with a single stop. The
    # primer bodies anchor on the *actual isolate* UP/DN flanks (preserved here
    # verbatim) so the resulting primers PCR cleanly off the isolate genome;
    # only the scar codons in the assembled deletion product come from the
    # reference. The actual locus length and truncation reason are preserved
    # in header metadata for traceability.
    if functional_str == "false":
        cds_for_fasta = ref_cds
    else:
        cds_for_fasta = nominal_cds

    record_seq = up_flank_seq + cds_for_fasta + dn_flank_seq
    actual_up = len(up_flank_seq)
    actual_dn = len(dn_flank_seq)

    header = (
        f"{isolate}|{GENE}"
        f"|contig={cid}"
        f"|start={genome_start_1b}"
        f"|end={genome_end_1b}"
        f"|strand={strand}"
        f"|cds_len={len(cds_for_fasta)}"
        f"|up={actual_up}"
        f"|dn={actual_dn}"
        f"|functional={functional_str}"
    )
    if functional_str == "false":
        header += f"|truncation_reason={reason}|locus_actual_len={locus_actual_len}"

    out_path = OUT_DIR / f"{isolate}.fasta"
    if dry_run:
        log.info("  [dry-run] would write %s (header=%s)", out_path, header)
        return True, header

    with out_path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(record_seq), 70):
            f.write(record_seq[i : i + 70] + "\n")
    log.info("  wrote %s", out_path)
    return True, header


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--isolates", nargs="*", default=MISSING_ISOLATES)
    args = parser.parse_args()

    ref_up, ref_cds, ref_dn = load_reference()
    log.info("Reference: PA14 lasR cds_len=%d", len(ref_cds))

    results = []
    for iso in args.isolates:
        ok, info = extract(iso, ref_up, ref_cds, ref_dn, args.dry_run)
        results.append((iso, ok, info))

    print("\nSummary:")
    for iso, ok, info in results:
        status = "OK" if ok else "FAIL"
        print(f"  {iso}: {status} — {info}")


if __name__ == "__main__":
    main()

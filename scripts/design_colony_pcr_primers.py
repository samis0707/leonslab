#!/usr/bin/env python3
"""Design conserved colony PCR verification primers for pEXG2 deletions.

Reads all per-isolate FASTA records for a gene, finds conserved primer sites
in the upstream buffer zone (outside) and downstream flank (inside), then
stores the best pair in data/primer_design/colony_pcr_primers.json.

Usage:
    python scripts/design_colony_pcr_primers.py --genes lasB lasR lasI
    python scripts/design_colony_pcr_primers.py --genes lasB --ref-isolate PAO1
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from primer_design._bio_lite import tm_nn, reverse_complement
from primer_design.config import (
    BODY_LEN_MIN, BODY_LEN_MAX,
    GC_MIN, GC_MAX,
    COLONY_PCR_TM_MIN_C, COLONY_PCR_TM_MAX_C, COLONY_PCR_TM_TARGET_C,
    COLONY_PCR_OUTSIDE_SEARCH_MAX, COLONY_PCR_INSIDE_SEARCH_MAX,
    TM_NN_PARAMS,
    POLYMERASE_OFFSETS_C,
    SUPPORTED_GENES,
)
from primer_design.gene_finder import get_gene_record, list_available_isolates

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger("design_colony_pcr_primers")

SCAR_ESTIMATE_BP = 20   # rough scar length for product size estimate
TAQMAN_OFFSET = POLYMERASE_OFFSETS_C["Taq"]   # -5.0 °C


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", nargs="+", default=list(SUPPORTED_GENES))
    parser.add_argument("--ref-isolate", default=None,
                        help="Preferred reference isolate (default: PAO1 > PA14 > first available)")
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "primer_design" / "colony_pcr_primers.json",
    )
    parser.add_argument("--max-mismatch", type=int, default=1,
                        help="Max allowed mismatches when checking conservation across isolates")
    args = parser.parse_args()

    store: dict = {}
    if args.output.exists():
        with open(args.output) as f:
            store = json.load(f)

    for gene in args.genes:
        log.info("Designing colony PCR primers for %s …", gene)
        isolates = list_available_isolates(gene)
        if not isolates:
            log.warning("No records found for %s — skipping", gene)
            continue

        ref = _pick_ref(isolates, args.ref_isolate)
        log.info("  Reference isolate: %s (%d total)", ref, len(isolates))

        ref_record = get_gene_record(ref, gene)
        all_records = [get_gene_record(iso, gene) for iso in isolates]

        outside = _design_outside(ref_record, all_records, args.max_mismatch, gene)
        inside = _design_inside(ref_record, all_records, args.max_mismatch, gene)

        if outside is None:
            log.error("  Could not find conserved outside primer for %s", gene)
            continue
        if inside is None:
            log.error("  Could not find conserved inside primer for %s", gene)
            continue

        # Estimate product sizes
        # Outside primer at position outside["flank_pos"] of up_flank
        # Inside primer at position inside["flank_pos"] of dn_flank (reverse primer, binds at pos:pos+len)
        out_pos = outside["flank_pos"]
        in_pos = inside["flank_pos"] + inside["length"]  # end of binding site
        up_len = len(ref_record.up_flank)
        del_product = (up_len - out_pos) + SCAR_ESTIMATE_BP + in_pos
        wt_product = del_product + len(ref_record.cds_seq)

        store[gene] = {
            "outside": {
                "body": outside["body"],
                "tm_body_C": round(outside["tm_C"], 2),
                "gc_body": round(outside["gc"], 4),
                "flank_pos": out_pos,
            },
            "inside": {
                "body": inside["body"],
                "tm_body_C": round(inside["tm_C"], 2),
                "gc_body": round(inside["gc"], 4),
                "flank_pos": in_pos,
            },
            "expected_deletion_product_bp_min": int(del_product * 0.9),
            "expected_deletion_product_bp_max": int(del_product * 1.1),
            "expected_wt_product_bp": wt_product,
            "note": (
                f"conserved across {outside['n_match']}/{len(all_records)} isolates (outside), "
                f"{inside['n_match']}/{len(all_records)} isolates (inside)"
            ),
        }
        log.info(
            "  OUT: %s  Tm=%.1f°C  pos=%d  conserved=%d/%d",
            outside["body"], outside["tm_C"], out_pos, outside["n_match"], len(all_records)
        )
        log.info(
            "  IN:  %s  Tm=%.1f°C  pos=%d  conserved=%d/%d",
            inside["body"], inside["tm_C"], in_pos, inside["n_match"], len(all_records)
        )
        log.info(
            "  Products: deletion ~%d–%d bp | WT ~%d bp",
            store[gene]["expected_deletion_product_bp_min"],
            store[gene]["expected_deletion_product_bp_max"],
            wt_product,
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(store, f, indent=2)
    log.info("Saved to %s", args.output)
    return 0


# ---------------------------------------------------------------------------
# Primer search helpers
# ---------------------------------------------------------------------------

def _design_outside(ref_record, all_records, max_mm: int, gene: str) -> dict | None:
    """Find best conserved forward primer in up_flank[0:COLONY_PCR_OUTSIDE_SEARCH_MAX]."""
    region = ref_record.up_flank[:COLONY_PCR_OUTSIDE_SEARCH_MAX]
    candidates = []
    for offset in range(len(region)):
        for length in range(BODY_LEN_MIN, BODY_LEN_MAX + 1):
            body = region[offset: offset + length]
            if len(body) < BODY_LEN_MIN:
                break
            if not _passes_filter(body):
                continue
            n_match = _check_conservation_fwd(body, offset, all_records, max_mm)
            tm = _tm(body)
            candidates.append({
                "body": body, "flank_pos": offset, "length": length,
                "tm_C": tm, "gc": _gc(body), "n_match": n_match,
                "score": abs(tm - COLONY_PCR_TM_TARGET_C) + 0.1 * (len(all_records) - n_match),
            })

    candidates = [c for c in candidates if c["n_match"] == len(all_records)]
    if not candidates:
        # Relax: allow up to max_mm mismatches total
        candidates = [c for c in _all_candidates_fwd(ref_record.up_flank[:COLONY_PCR_OUTSIDE_SEARCH_MAX])
                      if True]  # already filtered
        # Re-run with relaxed conservation
        candidates2 = []
        region = ref_record.up_flank[:COLONY_PCR_OUTSIDE_SEARCH_MAX]
        for offset in range(len(region)):
            for length in range(BODY_LEN_MIN, BODY_LEN_MAX + 1):
                body = region[offset: offset + length]
                if len(body) < BODY_LEN_MIN:
                    break
                if not _passes_filter(body):
                    continue
                n_match = _check_conservation_fwd(body, offset, all_records, max_mm)
                tm = _tm(body)
                candidates2.append({
                    "body": body, "flank_pos": offset, "length": length,
                    "tm_C": tm, "gc": _gc(body), "n_match": n_match,
                    "score": abs(tm - COLONY_PCR_TM_TARGET_C) + 0.1 * (len(all_records) - n_match),
                })
        candidates = sorted(candidates2, key=lambda c: (-c["n_match"], c["score"]))
        if not candidates:
            return None
        return candidates[0]

    candidates.sort(key=lambda c: c["score"])
    return candidates[0]


def _all_candidates_fwd(region: str) -> list[dict]:
    out = []
    for offset in range(len(region)):
        for length in range(BODY_LEN_MIN, BODY_LEN_MAX + 1):
            body = region[offset: offset + length]
            if len(body) < BODY_LEN_MIN:
                break
            if _passes_filter(body):
                out.append({"body": body, "offset": offset})
    return out


def _design_inside(ref_record, all_records, max_mm: int, gene: str) -> dict | None:
    """Find best conserved reverse primer in dn_flank[0:COLONY_PCR_INSIDE_SEARCH_MAX].

    The inside primer is a REVERSE primer: it binds the non-coding strand of dn_flank.
    body = RC(dn_flank[offset:offset+length]), so the primer sequence itself is the RC.
    """
    region = ref_record.dn_flank[:COLONY_PCR_INSIDE_SEARCH_MAX]
    candidates = []
    for offset in range(len(region)):
        for length in range(BODY_LEN_MIN, BODY_LEN_MAX + 1):
            template = region[offset: offset + length]
            if len(template) < BODY_LEN_MIN:
                break
            body = reverse_complement(template)   # the actual primer sequence
            if not _passes_filter(body):
                continue
            n_match = _check_conservation_rev(body, offset, length, all_records, max_mm)
            tm = _tm(body)
            candidates.append({
                "body": body, "flank_pos": offset, "length": length,
                "tm_C": tm, "gc": _gc(body), "n_match": n_match,
                "score": abs(tm - COLONY_PCR_TM_TARGET_C) + 0.1 * (len(all_records) - n_match),
            })

    best_match = max((c["n_match"] for c in candidates), default=0)
    candidates = [c for c in candidates if c["n_match"] == best_match]
    if not candidates:
        return None
    candidates.sort(key=lambda c: c["score"])
    return candidates[0]


# ---------------------------------------------------------------------------
# Conservation checks
# ---------------------------------------------------------------------------

def _check_conservation_fwd(body: str, offset: int, all_records, max_mm: int) -> int:
    """Count isolate records where body (fwd) matches up_flank at offset with ≤max_mm mismatches."""
    n = 0
    for rec in all_records:
        region = rec.up_flank[:COLONY_PCR_OUTSIDE_SEARCH_MAX + len(body)]
        # allow offset to shift by ±5 nt for small genomic rearrangements
        best = 0
        for delta in range(-5, 6):
            pos = offset + delta
            if pos < 0 or pos + len(body) > len(region):
                continue
            mm = sum(a != b for a, b in zip(body, region[pos: pos + len(body)]))
            best = max(best, len(body) - mm)
        if (len(body) - best) <= max_mm:
            n += 1
    return n


def _check_conservation_rev(body: str, offset: int, length: int, all_records, max_mm: int) -> int:
    """Count isolate records where RC(body) matches dn_flank at offset with ≤max_mm mismatches."""
    template = reverse_complement(body)
    n = 0
    for rec in all_records:
        region = rec.dn_flank[:COLONY_PCR_INSIDE_SEARCH_MAX + length]
        best = 0
        for delta in range(-5, 6):
            pos = offset + delta
            if pos < 0 or pos + length > len(region):
                continue
            mm = sum(a != b for a, b in zip(template, region[pos: pos + length]))
            best = max(best, length - mm)
        if (length - best) <= max_mm:
            n += 1
    return n


# ---------------------------------------------------------------------------
# Thermodynamic / filter helpers
# ---------------------------------------------------------------------------

def _tm(body: str) -> float:
    return tm_nn(body, **TM_NN_PARAMS) + TAQMAN_OFFSET


def _gc(body: str) -> float:
    return (body.count("G") + body.count("C")) / len(body)


def _passes_filter(body: str) -> bool:
    if len(body) < BODY_LEN_MIN or len(body) > BODY_LEN_MAX:
        return False
    gc = _gc(body)
    if gc < GC_MIN or gc > GC_MAX:
        return False
    tm = _tm(body)
    if tm < COLONY_PCR_TM_MIN_C or tm > COLONY_PCR_TM_MAX_C:
        return False
    if body[-1] not in ("G", "C"):
        return False
    for i in range(len(body) - 3):
        if len(set(body[i:i+4])) == 1:
            return False
    return True


def _pick_ref(isolates: list[str], preferred: str | None) -> str:
    if preferred and preferred in isolates:
        return preferred
    for ref in ("PAO1", "PA14"):
        if ref in isolates:
            return ref
    return isolates[0]


if __name__ == "__main__":
    sys.exit(main())

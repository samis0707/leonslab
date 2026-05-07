"""
Authoritative HindIII-regen sweep: invokes the TOOL's actual designer
(lib.primer_design.primer_designer.search_deletion_primers) on every
(gene × isolate) FASTA and checks whether the chosen primers regenerate
AAGCTT in the final plasmid.

Trap mechanism (pEXG2 + HindIII + site_destroyed, see comparison_round1.md §5):
  - P1 tail = CATAAATGTAAAGCA (ends in 'A'). If P1 body starts with AGCTT,
    the left vector junction reads ...AAAGCAAGCTT... → AAGCTT regenerated.
  - P4 tail (RC into top strand of insert) = AGCTTCTGCAGGTCG. If RC(P4_body)
    ends with 'A' (i.e. P4 body starts with 'T'), the right junction reads
    ...AAGCTTCTGCAG... → AAGCTT regenerated.

Output: docs/primer_design/in_silico_test/hindiii_trap_sweep.csv
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from primer_design.gene_finder import get_gene_record
from primer_design.vectors import get_tail_convention, get_vector
from primer_design.primer_designer import search_deletion_primers
from primer_design.exceptions import NoCandidates

DATA_GENES = Path("data/primer_design/genes")
OUT = Path("docs/primer_design/in_silico_test/hindiii_trap_sweep.csv")


def sweep_one(gene: str, isolate: str) -> dict:
    try:
        rec = get_gene_record(isolate, gene)
        vec = get_vector("pEXG2")
        conv = get_tail_convention("pEXG2", "HindIII", "deletion")
        best, _ = search_deletion_primers(rec, vec, conv)
    except NoCandidates as e:
        return {"isolate": isolate, "gene": gene, "result": "NO_CANDIDATES", "detail": str(e)[:120]}
    except Exception as e:
        return {"isolate": isolate, "gene": gene, "result": "ERROR", "detail": f"{type(e).__name__}: {e}"[:200]}

    p1_body = best.p1.body
    p4_body = best.p4.body
    p1_trap = p1_body.startswith("AGCTT")
    p4_trap = p4_body.startswith("T")
    insert_top_left  = "CATAAATGTAAAGCA" + p1_body[:6]
    from primer_design.vectors import reverse_complement
    insert_top_right = reverse_complement(p4_body)[-6:] + "AGCTTCTGCAGGTCG"
    regen_left  = "AAGCTT" in insert_top_left
    regen_right = "AAGCTT" in insert_top_right
    result = "TRAP" if (regen_left or regen_right) else "OK"
    return {
        "isolate": isolate, "gene": gene, "result": result,
        "N": best.N, "C": best.C,
        "p1_body": p1_body, "p4_body": p4_body,
        "p1_starts_AGCTT": p1_trap, "p4_starts_T": p4_trap,
        "regen_left": regen_left, "regen_right": regen_right,
        "score": round(best.score, 3), "tm_spread": round(best.tm_spread, 3),
        "detail": "",
    }


def main():
    rows = []
    for gene in ("lasR", "lasB"):
        for fasta in sorted((DATA_GENES / gene).glob("*.fasta")):
            isolate = fasta.stem
            rows.append(sweep_one(gene, isolate))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({k for r in rows for k in r.keys()})
    fieldnames = ["isolate", "gene", "result"] + [f for f in fieldnames if f not in {"isolate", "gene", "result"}]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    by_result: dict[str, list[str]] = {}
    for r in rows:
        by_result.setdefault(r["result"], []).append(f"{r['gene']}/{r['isolate']}")

    print(f"Total scanned: {len(rows)}")
    for k, lst in sorted(by_result.items()):
        print(f"  {k}: {len(lst)}")
        if k != "OK":
            for x in lst[:50]:
                print(f"    - {x}")
            if len(lst) > 50:
                print(f"    ... and {len(lst) - 50} more")


if __name__ == "__main__":
    main()

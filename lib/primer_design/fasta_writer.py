"""Combined annotated FASTA writer (skill_v2 §12.1).

Single file with: primers (tail/body markup), amplicons, insert, full circular
plasmid, scar/cassette ORF + translation. Produces a string for inline base64
encoding in the API response.
"""
from __future__ import annotations

from io import StringIO

from Bio.Seq import Seq

from .types import (
    DeletionPrimerSet,
    DesignResult,
    ExpressionPrimerSet,
    Primer,
    TaggingPrimerSet,
)


def write_fasta(result: DesignResult) -> str:
    """Return the combined FASTA file content as a single string."""
    buf = StringIO()
    app = result.request.application
    gene = result.gene_record.gene
    isolate = result.gene_record.isolate_id
    base = f"{gene}_{isolate}_{app}"

    # 1. Primers
    for p in _primers_in_order(result):
        _write_primer(buf, p, app, base)
        buf.write("\n")

    # 2. Amplicons (deletion / tagging only)
    if app in ("deletion", "tagging"):
        if result.up_amplicon is not None:
            _write_record(buf, f"{base}_UP_amplicon", result.up_amplicon,
                          extra=f"len={len(result.up_amplicon)} role=UP")
        if result.dn_amplicon is not None:
            _write_record(buf, f"{base}_DN_amplicon", result.dn_amplicon,
                          extra=f"len={len(result.dn_amplicon)} role=DN")

    # 3. Insert
    _write_record(buf, f"{base}_insert", result.insert,
                  extra=f"len={len(result.insert)} junctions_collapsed=true")

    # 4. Full plasmid
    plasmid = result.final_plasmid
    _write_record(buf, plasmid.name, plasmid.sequence,
                  extra=f"circular=true len={plasmid.length} validation=passed")

    # 5. Scar / cassette ORF
    if app == "deletion":
        assert isinstance(result.primer_set, DeletionPrimerSet)
        scar = result.primer_set.scar_dna
        protein = str(Seq(scar).translate())
        _write_record(buf, f"{base}_scar_ORF", scar,
                      extra=f"len={len(scar)} translation={protein}")
    elif app == "tagging":
        assert isinstance(result.primer_set, TaggingPrimerSet)
        cassette = result.primer_set.cassette
        protein = str(Seq(cassette).translate())
        _write_record(buf, f"{base}_tag_cassette", cassette,
                      extra=f"len={len(cassette)} translation={protein} tag={result.primer_set.tag.name}")
    elif app == "expression":
        assert isinstance(result.primer_set, ExpressionPrimerSet)
        cds = result.primer_set.coding_seq
        protein = str(Seq(cds).translate())
        _write_record(buf, f"{base}_expressed_CDS", cds,
                      extra=f"len={len(cds)} translation={protein}")

    return buf.getvalue()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _primers_in_order(result: DesignResult) -> list[Primer]:
    ps = result.primer_set
    if isinstance(ps, (DeletionPrimerSet, TaggingPrimerSet)):
        return [ps.p1, ps.p2, ps.p3, ps.p4]
    if isinstance(ps, ExpressionPrimerSet):
        return [ps.p1, ps.p2]
    raise TypeError(f"Unknown primer set type: {type(ps).__name__}")


def _write_primer(buf: StringIO, p: Primer, app: str, base: str) -> None:
    """Write a primer with annotated header.

    Header carries: app, role, len(tail), len(body), Tm, GC, tail_kind.
    Sequence is written with lowercase tail + UPPERCASE body for visual distinction.
    """
    buf.write(
        f">{p.name}_{base}_{p.role} app={app} tail_len={len(p.tail)} "
        f"body_len={len(p.body)} Tm_body={p.tm_body_C:.1f}C GC_body={int(p.gc_body * 100)}% "
        f"tail_kind={p.tail_kind}\n"
    )
    buf.write(p.annotated_sequence + "\n")


def _write_record(buf: StringIO, header: str, sequence: str, extra: str = "", line_width: int = 60) -> None:
    """Write one FASTA record with line-wrapping."""
    buf.write(f">{header}")
    if extra:
        buf.write(f"  {extra}")
    buf.write("\n")
    for i in range(0, len(sequence), line_width):
        buf.write(sequence[i : i + line_width] + "\n")

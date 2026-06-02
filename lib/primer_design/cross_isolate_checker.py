"""Check whether a primer set's binding sites exist identically in other isolates.

Only reads bundled gene FASTAs (no R2 / network access). Fast enough to run
inline during every design request (~300 string searches on ~2 400 nt strings).
"""
from __future__ import annotations

from pathlib import Path

from . import gene_finder
from .types import ExpressionPrimerSet, PrimerSet
from .vectors import reverse_complement


def find_compatible_isolates(
    primer_set: PrimerSet,
    gene: str,
    source_isolate_id: str,
    genes_dir: Path | None = None,
) -> list[str]:
    """Return sorted list of isolate IDs where all primer bodies match exactly.

    A primer body matches when it (or its reverse complement) appears verbatim
    in the other isolate's full gene sequence (up_flank + cds_seq + dn_flank).
    All four primers must match for deletion/tagging; both primers for expression.
    """
    if isinstance(primer_set, ExpressionPrimerSet):
        primers = [primer_set.p1, primer_set.p2]
    else:
        primers = [primer_set.p1, primer_set.p2, primer_set.p3, primer_set.p4]

    # Pre-compute both orientations once.
    search_targets = [
        (p.body.upper(), reverse_complement(p.body.upper())) for p in primers
    ]

    compatible: list[str] = []
    for isolate_id in gene_finder.list_available_isolates(gene, genes_dir):
        if isolate_id == source_isolate_id:
            continue
        try:
            record = gene_finder.get_gene_record(isolate_id, gene, genes_dir)
        except Exception:
            continue
        full_seq = (record.up_flank + record.cds_seq + record.dn_flank).upper()
        if all(fwd in full_seq or rev in full_seq for fwd, rev in search_targets):
            compatible.append(isolate_id)

    return sorted(compatible)

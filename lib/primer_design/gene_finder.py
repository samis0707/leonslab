"""Retrieve canonical gene records from bundled per-(gene, isolate) FASTAs.

The per-record FASTAs are produced once by ``scripts/build_gene_records.py`` and
shipped with the function deployment. At request time we just parse one small file
(no R2 fetch needed for the gene record itself).
"""
from __future__ import annotations

import re
from pathlib import Path

from ._bio_lite import Seq, parse_fasta as _parse_fasta

from .config import SUPPORTED_GENES
from .exceptions import CdsValidationError, GeneRecordNotFound
from .types import GeneRecord


def get_gene_record(isolate_id: str, gene: str, genes_dir: Path | None = None) -> GeneRecord:
    """Load the canonical strand-normalized record for one (gene, isolate) pair.

    Raises:
        GeneRecordNotFound: file missing (gene not present in this isolate).
        CdsValidationError: extracted CDS fails one of: length-mod-3, starts ATG,
            single internal stop, ends with stop.
    """
    if gene not in SUPPORTED_GENES:
        raise GeneRecordNotFound(
            message=f"Gene '{gene}' not in catalog. Supported: {SUPPORTED_GENES}"
        )

    genes_dir = genes_dir or _default_genes_dir()
    path = genes_dir / gene / f"{isolate_id}.fasta"
    if not path.exists():
        raise GeneRecordNotFound(
            message=f"No record for ({gene}, {isolate_id}) at {path}",
            details={"isolate_id": isolate_id, "gene": gene, "expected_path": str(path)},
        )

    record = next(_parse_fasta(path))
    meta = parse_header_metadata(record.description)
    sequence = str(record.seq).upper()
    up = int(meta["up"])
    dn = int(meta["dn"])

    if up + dn >= len(sequence):
        raise GeneRecordNotFound(
            message=f"Record header invalid for {path}: up={up}, dn={dn}, total={len(sequence)}"
        )

    cds = sequence[up : len(sequence) - dn]

    # functional defaults to True for legacy headers without the field
    functional = meta.get("functional", "true").lower() != "false"
    truncation_reason = meta.get("truncation_reason") or None

    if functional:
        _assert_cds_valid(cds, isolate_id, gene)
    # Non-functional records carry a *nominal* CDS (reference-derived) for
    # primer design; we deliberately skip strict CDS validation. The truncation
    # reason is propagated downstream so PDF/JSON output can warn the user.

    return GeneRecord(
        isolate_id=isolate_id,
        gene=gene,
        cds_seq=cds,
        up_flank=sequence[:up],
        dn_flank=sequence[len(sequence) - dn :],
        contig_id=meta["contig"],
        genome_start_1based=int(meta["start"]),
        genome_end_1based=int(meta["end"]),
        original_strand=meta["strand"],   # type: ignore[arg-type]
        functional=functional,
        truncation_reason=truncation_reason,
    )


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

# Header format produced by build_gene_records.py:
# >LB001|lasB|contig=LB001_00001|start=1848921|end=1850417|strand=+|cds_len=1497|up=600|dn=600
_HEADER_KV_RE = re.compile(r"(\w+)=([^\s|]+)")


def parse_header_metadata(header: str) -> dict[str, str]:
    """Parse the canonical header into a metadata dict.

    Required keys: contig, start, end, strand, cds_len, up, dn.
    The first two pipe-separated tokens (isolate_id and gene) are ignored here —
    they are recovered from the file path by ``get_gene_record``.
    """
    parts = header.split("|")
    meta = dict(_HEADER_KV_RE.findall(" ".join(parts)))
    required = {"contig", "start", "end", "strand", "cds_len", "up", "dn"}
    missing = required - meta.keys()
    if missing:
        raise GeneRecordNotFound(
            message=f"Canonical FASTA header missing keys: {missing}",
            details={"header": header},
        )
    return meta


# ---------------------------------------------------------------------------
# Validation invariants
# ---------------------------------------------------------------------------

def _assert_cds_valid(cds: str, isolate_id: str, gene: str) -> None:
    if len(cds) % 3 != 0:
        raise CdsValidationError(
            message=f"CDS length {len(cds)} not multiple of 3 ({gene}/{isolate_id})",
            details={"isolate_id": isolate_id, "gene": gene, "cds_len": len(cds)},
        )
    if not cds.startswith("ATG"):
        raise CdsValidationError(
            message=f"CDS does not start with ATG ({gene}/{isolate_id}): {cds[:10]}...",
            details={"isolate_id": isolate_id, "gene": gene},
        )
    protein = str(Seq(cds).translate())
    if protein.count("*") != 1 or not protein.endswith("*"):
        raise CdsValidationError(
            message=(
                f"CDS has {protein.count('*')} stop codons or not at end "
                f"({gene}/{isolate_id})"
            ),
            details={
                "isolate_id": isolate_id,
                "gene": gene,
                "stop_count": protein.count("*"),
                "protein_tail": protein[-15:],
            },
        )


# ---------------------------------------------------------------------------
# Discovery (used by frontend to populate isolate dropdown per gene)
# ---------------------------------------------------------------------------

def list_available_isolates(gene: str, genes_dir: Path | None = None) -> list[str]:
    """List isolate IDs for which a canonical record exists for ``gene``."""
    genes_dir = genes_dir or _default_genes_dir()
    gene_path = genes_dir / gene
    if not gene_path.exists():
        return []
    return sorted(p.stem for p in gene_path.glob("*.fasta"))


def list_supported_genes(genes_dir: Path | None = None) -> list[str]:
    """Return the configured gene catalog. Use ``list_available_isolates`` per gene
    to determine which (gene, isolate) pairs have curated records.
    """
    return list(SUPPORTED_GENES)


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

def _default_genes_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "data" / "primer_design" / "genes"

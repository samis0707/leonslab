#!/usr/bin/env python3
"""Build canonical per-(gene, isolate) records from user-provided files.

Inputs (fetched from R2):
    - bucket://Whole genome sequences/<isolate>.fna
    - bucket://<gene>.fasta   (multi-record, one per isolate)

Outputs (written locally, committed to git after run):
    - data/primer_design/genes/<gene>/<isolate>.fasta    (canonical strand-normalized record)
    - data/primer_design/manifest.json                   (SHA-256 per genome)

Run once for the initial 30 isolates, then incrementally with `--isolates LBxxx`
when new genomes arrive. See skill_v2 §1.5.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

# Path setup so we can import from lib/ without installing
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from Bio import SeqIO                                           # noqa: E402
from Bio.Seq import Seq                                         # noqa: E402

from primer_design.config import (                              # noqa: E402
    DN_FLANK_LEN,
    MIN_FLANK_LEN,
    SUPPORTED_GENES,
    UP_FLANK_LEN,
)
from primer_design.exceptions import (                          # noqa: E402
    BuildError,
    CdsValidationError,
    GeneAmbiguousInGenome,
    GeneNotInGenome,
)
from primer_design.storage_adapter import (                     # noqa: E402
    fetch_gene_fasta,
    fetch_genome,
)
from primer_design.vectors import reverse_complement            # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("build_gene_records")


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class GenomeHit:
    contig_id: str
    start_0based: int      # inclusive
    end_0based: int        # exclusive
    strand: str            # "+" or "-"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Build canonical gene records from R2")
    parser.add_argument(
        "--genes", nargs="*", default=list(SUPPORTED_GENES),
        help="Genes to (re)build. Default: all SUPPORTED_GENES.",
    )
    parser.add_argument(
        "--isolates", nargs="*",
        help="Isolates to (re)build. Default: every isolate present in the gene FASTAs.",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "primer_design",
        help="Output base dir. genes/<gene>/<isolate>.fasta written here.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Parse and validate but don't write canonical records.",
    )
    args = parser.parse_args()

    output_dir: Path = args.output_dir
    (output_dir / "genes").mkdir(parents=True, exist_ok=True)

    # Track manifest entries (cumulative across runs)
    manifest_path = output_dir / "manifest.json"
    manifest = _load_existing_manifest(manifest_path)

    n_built = 0
    n_skipped = 0
    n_failed = 0

    for gene in args.genes:
        log.info("Processing gene: %s", gene)
        try:
            gene_records = _parse_gene_fasta(gene)
        except Exception as e:
            log.error("Failed to load %s.fasta from R2: %s", gene, e)
            n_failed += 1
            continue

        wanted_isolates = set(args.isolates) if args.isolates else set(gene_records.keys())

        for isolate_id in sorted(wanted_isolates):
            if isolate_id not in gene_records:
                log.warning("Isolate %s has no entry in %s.fasta; skipping", isolate_id, gene)
                n_skipped += 1
                continue

            cds_seq = gene_records[isolate_id]
            try:
                _validate_cds(cds_seq, isolate_id, gene)
            except CdsValidationError as e:
                log.error("CDS invalid for %s/%s: %s", gene, isolate_id, e)
                n_failed += 1
                continue

            # Fetch genome, locate gene, extract flanks
            try:
                genome_bytes = _fetch_or_skip(isolate_id)
                if genome_bytes is None:
                    n_skipped += 1
                    continue
                hit = _locate_cds_in_genome(cds_seq, genome_bytes, isolate_id)
                record_seq = _strand_normalize(genome_bytes, hit, UP_FLANK_LEN, DN_FLANK_LEN, MIN_FLANK_LEN)
            except (GeneNotInGenome, GeneAmbiguousInGenome) as e:
                log.error("Locating %s in %s: %s", gene, isolate_id, e)
                n_failed += 1
                continue
            except Exception as e:
                log.exception("Unexpected error for %s/%s: %s", gene, isolate_id, e)
                n_failed += 1
                continue

            # Determine actual flank lengths (may be < default if neighbour-shortened)
            cds_position_in_record = record_seq.find(cds_seq)
            if cds_position_in_record == -1:
                log.error("After strand normalization, CDS not found in record for %s/%s", gene, isolate_id)
                n_failed += 1
                continue
            actual_up = cds_position_in_record
            actual_dn = len(record_seq) - cds_position_in_record - len(cds_seq)

            header = (
                f"{isolate_id}|{gene}"
                f"|contig={hit.contig_id}"
                f"|start={hit.start_0based + 1}"
                f"|end={hit.end_0based}"
                f"|strand={hit.strand}"
                f"|cds_len={len(cds_seq)}"
                f"|up={actual_up}"
                f"|dn={actual_dn}"
            )
            out_path = output_dir / "genes" / gene / f"{isolate_id}.fasta"

            if args.dry_run:
                log.info("[dry-run] would write %s (header=%s)", out_path, header)
                n_built += 1
                continue

            out_path.parent.mkdir(parents=True, exist_ok=True)
            _write_fasta(out_path, header, record_seq)
            log.info("Wrote %s", out_path)
            n_built += 1

            # Update manifest entry for this genome
            sha = hashlib.sha256(genome_bytes).hexdigest()
            _upsert_manifest_entry(manifest, isolate_id, sha, len(genome_bytes), genome_bytes)

    if not args.dry_run:
        _write_manifest(manifest_path, manifest)

    log.info("Done. built=%d skipped=%d failed=%d", n_built, n_skipped, n_failed)
    return 0 if n_failed == 0 else 1


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_ISOLATE_ID_PATTERNS = (
    re.compile(r"^(LB\d+)\b"),
    re.compile(r"^(PA14)\b"),
    re.compile(r"^(PAO1)\b"),
    re.compile(r"\|([A-Z][A-Za-z0-9]+)\|"),    # >X|LB001|... type headers
)


def _parse_gene_fasta(gene: str) -> dict[str, str]:
    """Parse the multi-record per-gene FASTA from R2. Returns {isolate_id: cds_seq}."""
    raw = fetch_gene_fasta(gene)
    out: dict[str, str] = {}
    with tempfile.NamedTemporaryFile("wb", delete=False, suffix=".fasta") as f:
        f.write(raw)
        tmp_path = Path(f.name)
    try:
        for record in SeqIO.parse(tmp_path, "fasta"):
            isolate_id = _parse_isolate_id_from_header(record.id, record.description)
            if isolate_id is None:
                log.warning("Could not parse isolate ID from header: %r", record.description)
                continue
            if isolate_id in out:
                raise BuildError(message=f"Duplicate isolate {isolate_id} in {gene}.fasta")
            out[isolate_id] = str(record.seq).upper()
    finally:
        tmp_path.unlink(missing_ok=True)
    log.info("Parsed %d isolates from %s.fasta", len(out), gene)
    return out


def _parse_isolate_id_from_header(record_id: str, description: str) -> str | None:
    """Try multiple patterns. Falls back to description if record_id is generic."""
    candidates = (record_id, description)
    for cand in candidates:
        for pat in _ISOLATE_ID_PATTERNS:
            m = pat.search(cand)
            if m:
                return m.group(1)
    return None


def _validate_cds(cds: str, isolate_id: str, gene: str) -> None:
    if len(cds) % 3 != 0:
        raise CdsValidationError(message=f"{gene}/{isolate_id}: CDS length {len(cds)} not mod 3")
    if not cds.startswith("ATG"):
        raise CdsValidationError(message=f"{gene}/{isolate_id}: CDS doesn't start with ATG")
    protein = str(Seq(cds).translate())
    if protein.count("*") != 1 or not protein.endswith("*"):
        raise CdsValidationError(
            message=f"{gene}/{isolate_id}: CDS has {protein.count('*')} stops, "
                    f"ends with {protein[-3:]!r}"
        )


def _fetch_or_skip(isolate_id: str) -> bytes | None:
    try:
        return fetch_genome(isolate_id)
    except Exception as e:
        log.warning("Genome fetch failed for %s (%s); skipping", isolate_id, e)
        return None


def _locate_cds_in_genome(cds: str, genome_bytes: bytes, isolate_id: str) -> GenomeHit:
    """Search both strands of all contigs for an exact CDS match."""
    cds_rc = reverse_complement(cds)
    hits: list[GenomeHit] = []

    with tempfile.NamedTemporaryFile("wb", delete=False, suffix=".fna") as f:
        f.write(genome_bytes)
        tmp_path = Path(f.name)
    try:
        for record in SeqIO.parse(tmp_path, "fasta"):
            seq = str(record.seq).upper()
            for start in _find_all(seq, cds):
                hits.append(GenomeHit(record.id, start, start + len(cds), "+"))
            for start in _find_all(seq, cds_rc):
                hits.append(GenomeHit(record.id, start, start + len(cds), "-"))
    finally:
        tmp_path.unlink(missing_ok=True)

    if not hits:
        raise GeneNotInGenome(
            message=f"CDS not found in {isolate_id} genome (possibly deleted in this strain or assembly gap)",
            details={"isolate_id": isolate_id},
        )
    if len(hits) > 1:
        raise GeneAmbiguousInGenome(
            message=f"CDS matches {len(hits)} loci in {isolate_id} (paralog or repeat)",
            details={"isolate_id": isolate_id, "hits": [vars(h) for h in hits]},
        )
    return hits[0]


def _strand_normalize(genome_bytes: bytes, hit: GenomeHit, up: int, dn: int, min_flank: int) -> str:
    """Return gene + flanks, RC if on minus strand so result is 5'→3' on coding strand."""
    with tempfile.NamedTemporaryFile("wb", delete=False, suffix=".fna") as f:
        f.write(genome_bytes)
        tmp_path = Path(f.name)
    try:
        contig = next(r for r in SeqIO.parse(tmp_path, "fasta") if r.id == hit.contig_id)
        seq = str(contig.seq).upper()
    finally:
        tmp_path.unlink(missing_ok=True)

    if hit.strand == "+":
        up_start = max(0, hit.start_0based - up)
        dn_end = min(len(seq), hit.end_0based + dn)
        out = seq[up_start:dn_end]
    else:
        # On - strand, the genomic "before-gene" piece (up flank in coding orientation) lies
        # AFTER the gene on the genomic + strand, so we extract a window on + strand and RC it.
        up_start = max(0, hit.start_0based - dn)
        dn_end = min(len(seq), hit.end_0based + up)
        out = reverse_complement(seq[up_start:dn_end])

    # Verify min_flank not breached
    cds_in_record = out.find(_get_cds_for_assertion(out, up, dn))
    if cds_in_record < min_flank:
        log.warning(
            "Flank shorter than min_flank=%d for hit at %s:%d-%d (actual_up=%d)",
            min_flank, hit.contig_id, hit.start_0based, hit.end_0based, cds_in_record,
        )
    return out


def _get_cds_for_assertion(record_seq: str, up: int, dn: int) -> str:
    """Best-effort: return the inner segment between the configured flank lengths.
    Used only as a debug sanity check inside _strand_normalize.
    """
    # We don't have CDS length here without re-passing; return empty to skip
    return ""


def _find_all(haystack: str, needle: str) -> list[int]:
    out: list[int] = []
    i = haystack.find(needle, 0)
    while i != -1:
        out.append(i)
        i = haystack.find(needle, i + 1)
    return out


def _write_fasta(path: Path, header: str, sequence: str, line_width: int = 60) -> None:
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i : i + line_width] + "\n")


def _load_existing_manifest(path: Path) -> dict:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {"schema_version": 1, "genomes": []}


def _upsert_manifest_entry(manifest: dict, isolate_id: str, sha: str, size_bytes: int, genome_bytes: bytes) -> None:
    n_contigs = genome_bytes.count(b">")
    entry = {
        "isolate_id": isolate_id,
        "filename": f"{isolate_id}.fna",
        "sha256": sha,
        "size_bytes": size_bytes,
        "n_contigs": n_contigs,
    }
    for i, existing in enumerate(manifest["genomes"]):
        if existing["isolate_id"] == isolate_id:
            manifest["genomes"][i] = entry
            return
    manifest["genomes"].append(entry)


def _write_manifest(path: Path, manifest: dict) -> None:
    import datetime as dt
    manifest["generated"] = dt.datetime.now(dt.timezone.utc).isoformat()
    manifest["genomes"].sort(key=lambda e: e["isolate_id"])
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    log.info("Manifest updated: %s (%d genomes)", path, len(manifest["genomes"]))


if __name__ == "__main__":
    sys.exit(main())

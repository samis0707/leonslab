#!/usr/bin/env python3
"""Build canonical gene records from a fully annotated GenBank reference genome.

Used for reference strains (PA14, PAO1) where we have a single annotated GBK
file rather than the (Prokka-annotated isolate fna + bucket-rooted gene fasta)
pair that ``build_gene_records.py`` consumes.

Pipeline per (gene, isolate) pair:
    1. Find the CDS feature whose /gene qualifier matches.
    2. Extract genomic CDS slice + ±UP_FLANK_LEN / DN_FLANK_LEN on the strand-
       normalized coding orientation.
    3. Validate (length-mod-3, starts ATG, single trailing stop).
    4. Write data/primer_design/genes/<gene>/<isolate>.fasta with the canonical
       header.
    5. Optionally derive the whole-genome FASTA (.fna) and upload it to R2 under
       ``Whole genome sequences/<isolate>.fna`` so storage_adapter.fetch_genome
       resolves the same way as for LB-prefixed isolates.
    6. Update data/primer_design/manifest.json with the genome SHA-256.

Inputs: a local annotated GenBank (.gbk/.gb) per isolate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from Bio import SeqIO                                            # noqa: E402

from primer_design.config import (                               # noqa: E402
    DN_FLANK_LEN,
    SUPPORTED_GENES,
    UP_FLANK_LEN,
)
from primer_design.exceptions import CdsValidationError, GeneNotInGenome  # noqa: E402
from primer_design.storage_adapter import _r2_client, _bucket_name       # noqa: E402
from primer_design.vectors import reverse_complement             # noqa: E402


def extract_gene(record, gene: str) -> tuple[int, int, str, str]:
    """Return (start_1based, end_1based, strand, gene_locus) for the CDS feature
    whose /gene qualifier matches ``gene``. Raises GeneNotInGenome if absent."""
    hits = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        if gene in feat.qualifiers.get("gene", []):
            hits.append(feat)
    if not hits:
        raise GeneNotInGenome(message=f"No CDS with /gene={gene} in {record.id}")
    if len(hits) > 1:
        raise GeneNotInGenome(
            message=f"Multiple /gene={gene} CDS in {record.id} ({len(hits)})"
        )
    feat = hits[0]
    start = int(feat.location.start) + 1   # GenBank locations are 0-based half-open
    end = int(feat.location.end)
    strand = "+" if feat.location.strand == 1 else "-"
    locus = feat.qualifiers.get("locus_tag", [None])[0]
    return start, end, strand, locus


def build_record(gbk_path: Path, isolate: str, gene: str, out_dir: Path) -> Path:
    record = next(SeqIO.parse(gbk_path, "genbank"))
    start_1, end_1, strand, _locus = extract_gene(record, gene)

    seq = str(record.seq).upper()
    cds_genomic = seq[start_1 - 1 : end_1]
    if strand == "+":
        coding_cds = cds_genomic
        up_genomic = seq[max(0, start_1 - 1 - UP_FLANK_LEN) : start_1 - 1]
        dn_genomic = seq[end_1 : end_1 + DN_FLANK_LEN]
        up_flank = up_genomic
        dn_flank = dn_genomic
    else:
        coding_cds = reverse_complement(cds_genomic)
        # On - strand, the coding 5' is at end_1 (genomic), 3' at start_1.
        up_genomic = seq[end_1 : end_1 + UP_FLANK_LEN]
        dn_genomic = seq[max(0, start_1 - 1 - DN_FLANK_LEN) : start_1 - 1]
        up_flank = reverse_complement(up_genomic)
        dn_flank = reverse_complement(dn_genomic)

    if len(coding_cds) % 3 != 0:
        raise CdsValidationError(
            message=f"CDS length {len(coding_cds)} not %3 ({gene}/{isolate})"
        )
    if not coding_cds.startswith("ATG"):
        raise CdsValidationError(
            message=f"CDS does not start with ATG ({gene}/{isolate})"
        )

    full_seq = up_flank + coding_cds + dn_flank
    header = (
        f"{isolate}|{gene}|contig={record.id}|start={start_1}|end={end_1}"
        f"|strand={strand}|cds_len={len(coding_cds)}"
        f"|up={len(up_flank)}|dn={len(dn_flank)}"
    )
    out_path = out_dir / gene / f"{isolate}.fasta"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(full_seq), 60):
            f.write(full_seq[i : i + 60] + "\n")
    return out_path


def derive_fna(gbk_path: Path, isolate: str, fna_path: Path) -> tuple[bytes, str]:
    """Write a single-record FASTA of the genome and return (bytes, sha256)."""
    record = next(SeqIO.parse(gbk_path, "genbank"))
    seq = str(record.seq).upper()
    fna_text = f">{isolate}_{record.id}\n"
    for i in range(0, len(seq), 60):
        fna_text += seq[i : i + 60] + "\n"
    data = fna_text.encode("ascii")
    fna_path.parent.mkdir(parents=True, exist_ok=True)
    fna_path.write_bytes(data)
    sha = hashlib.sha256(data).hexdigest()
    return data, sha


def upload_fna(isolate: str, data: bytes) -> None:
    key = f"Whole genome sequences/{isolate}.fna"
    _r2_client().put_object(
        Bucket=_bucket_name(), Key=key, Body=data, ContentType="text/x-fasta"
    )
    print(f"  uploaded R2://{key} ({len(data)} bytes)")


def update_manifest(isolate: str, sha: str, size: int, manifest_path: Path) -> None:
    manifest = json.loads(manifest_path.read_text())
    entries = [e for e in manifest.get("genomes", []) if e["isolate_id"] != isolate]
    entries.append({
        "isolate_id": isolate,
        "filename": f"{isolate}.fna",
        "sha256": sha,
        "size_bytes": size,
        "n_contigs": 1,
    })
    entries.sort(key=lambda e: e["isolate_id"])
    manifest["genomes"] = entries
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--isolate", required=True, help="e.g. PA14")
    parser.add_argument("--gbk", required=True, type=Path)
    parser.add_argument(
        "--genes", nargs="+", default=list(SUPPORTED_GENES),
        help="Gene names to extract; default: all SUPPORTED_GENES",
    )
    parser.add_argument(
        "--upload", action="store_true",
        help="Also derive the .fna, upload to R2, and update manifest.json",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    genes_dir = repo_root / "data" / "primer_design" / "genes"
    manifest_path = repo_root / "data" / "primer_design" / "manifest.json"

    for gene in args.genes:
        try:
            out = build_record(args.gbk, args.isolate, gene, genes_dir)
            print(f"  wrote {out.relative_to(repo_root)}")
        except (GeneNotInGenome, CdsValidationError) as e:
            print(f"  skip {gene}: {e}")

    if args.upload:
        fna_path = repo_root / "data" / "primer_design" / "_tmp_fna" / f"{args.isolate}.fna"
        data, sha = derive_fna(args.gbk, args.isolate, fna_path)
        upload_fna(args.isolate, data)
        update_manifest(args.isolate, sha, len(data), manifest_path)
        print(f"  manifest updated for {args.isolate} (sha={sha[:10]}…)")
        fna_path.unlink()


if __name__ == "__main__":
    main()

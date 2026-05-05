#!/usr/bin/env python3
"""One-time R2 upload of isolate genomes (.fna files).

Usage:
    python scripts/upload_genomes.py --source-dir /path/to/genomes
    python scripts/upload_genomes.py --source-dir ... --only LB200 LB201
    python scripts/upload_genomes.py --source-dir ... --overwrite

The user already uploaded the initial 30 genomes manually via the Cloudflare web UI
on 2026-05-04. This script is for incremental additions in Phase 2.

After upload, run ``scripts/build_gene_records.py`` to generate the canonical
per-(gene, isolate) records that ship with the function.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from primer_design.storage_adapter import upload_genome   # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("upload_genomes")


def main() -> int:
    parser = argparse.ArgumentParser(description="Upload .fna genomes to R2")
    parser.add_argument("--source-dir", type=Path, required=True,
                        help="Local directory containing LB001.fna, LB014.fna, ... files")
    parser.add_argument("--only", nargs="*",
                        help="Restrict to these isolate IDs. Default: all .fna in source-dir.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Re-upload even if key exists in bucket.")
    args = parser.parse_args()

    if not args.source_dir.exists():
        log.error("Source dir not found: %s", args.source_dir)
        return 1

    fna_files = sorted(args.source_dir.glob("*.fna"))
    if args.only:
        wanted = set(args.only)
        fna_files = [p for p in fna_files if p.stem in wanted]

    if not fna_files:
        log.error("No .fna files found in %s matching filter", args.source_dir)
        return 1

    log.info("Uploading %d genomes to R2", len(fna_files))
    n_ok = 0
    for fna in fna_files:
        try:
            sha = upload_genome(fna.stem, fna, overwrite=args.overwrite)
            log.info("  %s — sha256=%s", fna.name, sha[:16] + "...")
            n_ok += 1
        except Exception as e:
            log.error("  %s — FAILED: %s", fna.name, e)

    log.info("Done. Uploaded %d/%d.", n_ok, len(fna_files))
    log.info("Next step: python scripts/build_gene_records.py --isolates %s",
             " ".join(p.stem for p in fna_files))
    return 0 if n_ok == len(fna_files) else 1


if __name__ == "__main__":
    sys.exit(main())

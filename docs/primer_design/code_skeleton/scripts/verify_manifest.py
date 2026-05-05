#!/usr/bin/env python3
"""Verify that R2-bundled genomes match local manifest SHA-256.

Run after every genome upload or when investigating data integrity.
"""
from __future__ import annotations

import hashlib
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from primer_design.storage_adapter import fetch_genome   # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("verify_manifest")


def main() -> int:
    manifest_path = Path(__file__).resolve().parents[1] / "data" / "primer_design" / "manifest.json"
    if not manifest_path.exists():
        log.error("Manifest not found: %s", manifest_path)
        return 1

    with open(manifest_path) as f:
        manifest = json.load(f)

    n_ok = 0
    n_fail = 0
    for entry in manifest.get("genomes", []):
        isolate_id = entry["isolate_id"]
        try:
            data = fetch_genome.__wrapped__(isolate_id)  # bypass LRU cache
        except Exception:
            data = fetch_genome(isolate_id)
        actual = hashlib.sha256(data).hexdigest()
        if actual == entry["sha256"]:
            log.info("  OK    %s", isolate_id)
            n_ok += 1
        else:
            log.error("  FAIL  %s  expected=%s actual=%s", isolate_id, entry["sha256"][:16], actual[:16])
            n_fail += 1

    log.info("Verified %d/%d genomes", n_ok, n_ok + n_fail)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

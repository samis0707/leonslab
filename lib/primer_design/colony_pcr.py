"""Load pre-computed colony PCR verification primers.

Primers are designed once per gene (conserved across isolates) and stored in
data/primer_design/colony_pcr_primers.json. At request time we just look up the
entry — no computation needed.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from .types import ColonyPCRPrimerSet, Primer


_DATA_PATH = Path(__file__).resolve().parents[2] / "data" / "primer_design" / "colony_pcr_primers.json"


@lru_cache(maxsize=1)
def _load_store() -> dict:
    if not _DATA_PATH.exists():
        return {}
    with open(_DATA_PATH) as f:
        return json.load(f)


def get_colony_pcr_primers(gene: str) -> ColonyPCRPrimerSet | None:
    """Return the pre-computed colony PCR primer set for *gene*, or None if not stored."""
    store = _load_store()
    entry = store.get(gene)
    if entry is None:
        return None

    outside = Primer(
        name="OUT",
        role="Colony_Fwd",
        tail="",
        body=entry["outside"]["body"],
        tail_kind="colony_pcr",
        tm_body_C=entry["outside"]["tm_body_C"],
        gc_body=entry["outside"]["gc_body"],
        length=len(entry["outside"]["body"]),
    )
    inside = Primer(
        name="IN",
        role="Colony_Rev",
        tail="",
        body=entry["inside"]["body"],
        tail_kind="colony_pcr",
        tm_body_C=entry["inside"]["tm_body_C"],
        gc_body=entry["inside"]["gc_body"],
        length=len(entry["inside"]["body"]),
    )
    return ColonyPCRPrimerSet(
        gene=gene,
        outside=outside,
        inside=inside,
        expected_deletion_product_bp_min=entry["expected_deletion_product_bp_min"],
        expected_deletion_product_bp_max=entry["expected_deletion_product_bp_max"],
        note=entry.get("note", ""),
    )

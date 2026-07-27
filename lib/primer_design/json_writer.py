"""JSON parameter manifest for full reproducibility (D7.3 / project_plan §2.7)."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import uuid

from . import config as cfg
from ._bio_lite import translate_dna as _translate_dna
from .types import (
    ColonyPCRPrimerSet,
    DeletionPrimerSet,
    DesignResult,
    ExpressionPrimerSet,
    Primer,
    TaggingPrimerSet,
)


def write_manifest(result: DesignResult) -> dict:
    """Build the JSON manifest as a dict (caller serializes)."""
    ps = result.primer_set
    primers: list[Primer]
    if isinstance(ps, (DeletionPrimerSet, TaggingPrimerSet)):
        primers = [ps.p1, ps.p2, ps.p3, ps.p4]
    elif isinstance(ps, ExpressionPrimerSet):
        primers = [ps.p1, ps.p2]
    else:
        raise TypeError(type(ps).__name__)

    # Wallace-rule annealing temp: Ta = Tm_body - 5 (GC*4 + AT*2 - 5), lowest-Tm
    # primer governs so every primer in the set anneals at the reported Ta.
    annealing_temp = (
        min(p.tm_body_C for p in primers) - 5.0
        + cfg.POLYMERASE_OFFSETS_C[result.request.polymerase]
    )

    manifest: dict = {
        "tool_version": cfg.TOOL_VERSION,
        "job_id": str(uuid.uuid4()),
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "request": {
            "isolate_id": result.request.isolate_id,
            "gene": result.request.gene,
            "action": result.request.action,
            "tag": result.request.tag,
            "tag_position": result.request.tag_position,
            "use_plasmid_for_tag": result.request.use_plasmid_for_tag,
            "enzyme": result.request.enzyme,
            "polymerase": result.request.polymerase,
            "application": result.request.application,
            "vector": result.request.vector,
        },
        "data_hashes": {
            "vector_sha256": hashlib.sha256(result.vector.sequence.encode()).hexdigest(),
            "gene_record_sha256": hashlib.sha256(result.gene_record.cds_seq.encode()).hexdigest(),
        },
        "computed": {
            "primers": [_primer_to_dict(p) for p in primers],
            "annealing_temp_C_recommended": annealing_temp,
            "final_plasmid_length_bp": result.final_plasmid.length,
            "final_plasmid_sha256": hashlib.sha256(result.final_plasmid.sequence.encode()).hexdigest(),
            "off_target": {
                "passed": result.off_target.passed,
                "violations_count": len(result.off_target.violations),
            },
            "colony_pcr_primers": _colony_pcr_to_dict(result.colony_pcr_primers),
            "warnings": result.warnings,
        },
        "convention": {
            "name": result.convention.name,
            "p1_tail_rule": result.convention.p1_tail_rule,
            "p4_tail_rule": result.convention.p4_tail_rule,
            "expected_recognition_count_in_final_plasmid": result.convention.expected_recognition_count_in_final_plasmid,
        },
    }

    # Application-specific fields
    if isinstance(ps, DeletionPrimerSet):
        manifest["computed"]["scar"] = {
            "N": ps.N,
            "C": ps.C,
            "dna": ps.scar_dna,
            "translation": _translate_dna(ps.scar_dna),
        }
    elif isinstance(ps, TaggingPrimerSet):
        manifest["computed"]["tag_cassette"] = {
            "tag": ps.tag.name,
            "dna": ps.cassette,
            "translation": _translate_dna(ps.cassette),
            "overlap_left": ps.overlap_left,
            "overlap_right": ps.overlap_right,
        }
    elif isinstance(ps, ExpressionPrimerSet):
        manifest["computed"]["expression"] = {
            "tag": ps.tag.name if ps.tag else None,
            "tag_position": ps.tag_position,
            "coding_seq_len_nt": len(ps.coding_seq),
            "coding_seq_sha256": hashlib.sha256(ps.coding_seq.encode()).hexdigest(),
        }

    return manifest


def _primer_to_dict(p: Primer) -> dict:
    return {
        "name": p.name,
        "role": p.role,
        "sequence": p.sequence,
        "tail": p.tail,
        "body": p.body,
        "tail_len": len(p.tail),
        "body_len": len(p.body),
        "tm_body_C": round(p.tm_body_C, 2),
        "gc_body_pct": round(p.gc_body * 100, 1),
        "length": p.length,
        "tail_kind": p.tail_kind,
    }


def _colony_pcr_to_dict(cpcr) -> dict | None:
    if cpcr is None:
        return None
    return {
        "gene": cpcr.gene,
        "outside": {
            "sequence": cpcr.outside.body,
            "tm_body_C": round(cpcr.outside.tm_body_C, 2),
            "gc_body_pct": round(cpcr.outside.gc_body * 100, 1),
            "length": cpcr.outside.length,
        },
        "inside": {
            "sequence": cpcr.inside.body,
            "tm_body_C": round(cpcr.inside.tm_body_C, 2),
            "gc_body_pct": round(cpcr.inside.gc_body * 100, 1),
            "length": cpcr.inside.length,
        },
        "expected_deletion_product_bp_min": cpcr.expected_deletion_product_bp_min,
        "expected_deletion_product_bp_max": cpcr.expected_deletion_product_bp_max,
        "note": cpcr.note,
    }


def to_json_string(manifest: dict) -> str:
    return json.dumps(manifest, indent=2, ensure_ascii=False)

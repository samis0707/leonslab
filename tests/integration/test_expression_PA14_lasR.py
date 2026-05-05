"""Acceptance regression: PA14 lasR pBBR1MCS2 expression (skill_v2 §15.2).

Reproduces the empirical primer set used in the working PA14_lasR pBBR1MCS2
expression construct. Body length tolerated within ±2 nt; tails must be exact.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from primer_design.applications.expression import run as run_expression
from primer_design.types import DesignRequest


@pytest.fixture
def all_required_files_present() -> bool:
    root = Path(__file__).resolve().parents[2]
    required = [
        "data/primer_design/vectors/pBBR1MCS2.gb",
        "data/primer_design/genes/lasR/PA14.fasta",
    ]
    return all((root / r).exists() for r in required)


@pytest.fixture
def request_obj() -> DesignRequest:
    return DesignRequest(
        isolate_id="PA14",
        gene="lasR",
        action="express",
        enzyme="HindIII",
        polymerase="B7",
    )


# Expected from skill_v2 §15.2
EXPECTED = {
    "p1_tail": "CGGTATCGATAAGCT" + "AGGAGG" + "ACTTGTTC",   # 29 nt: vector + RBS + spacer
    "p2_tail": "ATTCGATATCAAGCT",                            # 15 nt: vector tail
    "p1_body_starts": "ATG",                                  # body anchored at native ATG
    "amplicon_length_bp_min": 740,                            # ±10 bp tolerance
    "amplicon_length_bp_max": 760,
    "AAGCTT_count_in_final_plasmid": 1,                       # site_partial_AAGCT regenerates one site
}


class TestPA14LasRExpression:

    def test_runs_without_error(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_expression(request_obj)
        assert result is not None

    def test_p1_tail_exact(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_expression(request_obj)
        assert result.primer_set.p1.tail == EXPECTED["p1_tail"]

    def test_p2_tail_exact(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_expression(request_obj)
        assert result.primer_set.p2.tail == EXPECTED["p2_tail"]

    def test_p1_body_starts_at_atg(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_expression(request_obj)
        assert result.primer_set.p1.body.startswith("ATG")

    def test_amplicon_size_in_range(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_expression(request_obj)
        size = len(result.insert)
        assert EXPECTED["amplicon_length_bp_min"] <= size <= EXPECTED["amplicon_length_bp_max"]

    def test_one_HindIII_site_in_final_plasmid(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_expression(request_obj)
        seq = result.final_plasmid.sequence
        assert seq.count("AAGCTT") == EXPECTED["AAGCTT_count_in_final_plasmid"]

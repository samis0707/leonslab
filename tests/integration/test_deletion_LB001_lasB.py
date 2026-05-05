"""Acceptance regression: LB001 ΔlasB pEXG2 HindIII (skill_v2 §15.1).

This is the canonical "must work" test for the deletion application. It encodes
the empirical primer set from v1, and Phase 2 is complete when this test passes.

Required state before running:
    - Vector files in data/primer_design/vectors/ (Phase 2 step 2)
    - Canonical record at data/primer_design/genes/lasB/LB001.fasta (Phase 2 step 1)
    - LB001.fna present in R2 (already done — user uploaded 2026-05-04)
"""
from __future__ import annotations

from pathlib import Path

import pytest

from primer_design.applications.deletion import run as run_deletion
from primer_design.types import DesignRequest


REQUIRED_FILES = [
    "data/primer_design/vectors/pEXG2.gb",
    "data/primer_design/genes/lasB/LB001.fasta",
]


@pytest.fixture
def all_required_files_present() -> bool:
    root = Path(__file__).resolve().parents[2]
    for rel in REQUIRED_FILES:
        if not (root / rel).exists():
            return False
    return True


@pytest.fixture
def request_obj() -> DesignRequest:
    return DesignRequest(
        isolate_id="LB001",
        gene="lasB",
        action="delete",
        enzyme="HindIII",
        polymerase="B7",
    )


# Expected values from skill_v2 §15.1 — empirically validated v1 working set
EXPECTED = {
    "N": 16,
    "C": 2,
    "scar_protein": "MKKVSTLDLLFVAIMGAL*",
    "p1_tail": "CATAAATGTAAAGCA",
    "p4_tail": "CGACCTGCAGAAGCT",
    "final_plasmid_length_bp": 6155,
    "up_amplicon_length_bp": 593,
    "dn_amplicon_length_bp": 538,
    "insert_length_bp": 1101,
    "AAGCTT_count_in_final_plasmid": 0,
}


class TestLB001LasBDeletion:

    def test_runs_without_error(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_deletion(request_obj)
        assert result is not None

    def test_scar_NC_match_v1(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_deletion(request_obj)
        assert result.primer_set.N == EXPECTED["N"]
        assert result.primer_set.C == EXPECTED["C"]

    def test_scar_translation(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        from Bio.Seq import Seq
        result = run_deletion(request_obj)
        protein = str(Seq(result.primer_set.scar_dna).translate())
        assert protein == EXPECTED["scar_protein"]

    def test_p1_tail_exact(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_deletion(request_obj)
        assert result.primer_set.p1.tail == EXPECTED["p1_tail"]

    def test_p4_tail_exact(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_deletion(request_obj)
        assert result.primer_set.p4.tail == EXPECTED["p4_tail"]

    def test_final_plasmid_length(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_deletion(request_obj)
        assert result.final_plasmid.length == EXPECTED["final_plasmid_length_bp"]

    def test_no_HindIII_in_final_plasmid(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_deletion(request_obj)
        seq = result.final_plasmid.sequence
        assert seq.count("AAGCTT") == EXPECTED["AAGCTT_count_in_final_plasmid"]

    def test_off_target_passes(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_deletion(request_obj)
        assert result.off_target.passed

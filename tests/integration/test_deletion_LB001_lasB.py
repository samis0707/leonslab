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

from primer_design.applications import deletion as deletion_app
from primer_design.applications.deletion import run as run_deletion
from primer_design.exceptions import GenomeNotFoundInR2
from primer_design.types import DesignRequest


REQUIRED_FILES = [
    "data/primer_design/genes/lasB/LB001.fasta",
]
# pEXG2 vector: either .gb (preferred) or .fasta (fallback) is fine.
VECTOR_FILES_ANY_OF = [
    "data/primer_design/vectors/pEXG2.gb",
    "data/primer_design/vectors/pEXG2.fasta",
]


@pytest.fixture
def all_required_files_present() -> bool:
    root = Path(__file__).resolve().parents[2]
    for rel in REQUIRED_FILES:
        if not (root / rel).exists():
            return False
    if not any((root / rel).exists() for rel in VECTOR_FILES_ANY_OF):
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


@pytest.fixture
def design_result(all_required_files_present, request_obj):
    """Run the full pipeline if a genome is available; otherwise fall back to
    the search + assembly path so the non-off-target assertions can still run."""
    if not all_required_files_present:
        pytest.skip("Phase 2 prerequisites not yet in place")
    try:
        return run_deletion(request_obj)
    except GenomeNotFoundInR2:
        return deletion_app.design_without_off_target(request_obj)


@pytest.fixture
def full_run_or_skip(all_required_files_present, request_obj):
    if not all_required_files_present:
        pytest.skip("Phase 2 prerequisites not yet in place")
    try:
        return run_deletion(request_obj)
    except GenomeNotFoundInR2:
        pytest.skip("Genome not available locally and R2 not configured")


class TestLB001LasBDeletion:

    def test_runs_without_error(self, design_result):
        assert design_result is not None

    def test_scar_NC_match_v1(self, design_result):
        assert design_result.primer_set.N == EXPECTED["N"]
        assert design_result.primer_set.C == EXPECTED["C"]

    def test_scar_translation(self, design_result):
        from Bio.Seq import Seq
        protein = str(Seq(design_result.primer_set.scar_dna).translate())
        assert protein == EXPECTED["scar_protein"]

    def test_p1_tail_exact(self, design_result):
        assert design_result.primer_set.p1.tail == EXPECTED["p1_tail"]

    def test_p4_tail_exact(self, design_result):
        assert design_result.primer_set.p4.tail == EXPECTED["p4_tail"]

    def test_up_amplicon_length(self, design_result):
        assert len(design_result.up_amplicon) == EXPECTED["up_amplicon_length_bp"]

    def test_dn_amplicon_length(self, design_result):
        assert len(design_result.dn_amplicon) == EXPECTED["dn_amplicon_length_bp"]

    def test_insert_length(self, design_result):
        assert len(design_result.insert) == EXPECTED["insert_length_bp"]

    def test_final_plasmid_length(self, design_result):
        assert design_result.final_plasmid.length == EXPECTED["final_plasmid_length_bp"]

    def test_no_HindIII_in_final_plasmid(self, design_result):
        seq = design_result.final_plasmid.sequence
        assert seq.count("AAGCTT") == EXPECTED["AAGCTT_count_in_final_plasmid"]

    def test_off_target_passes(self, full_run_or_skip):
        assert full_run_or_skip.off_target.passed

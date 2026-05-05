"""Acceptance regression: PA14 lasR pEXG2 in-locus C-term His6 (skill_v2 §15.3).

Tests in-locus tagging mode. His6 cassette = 30 nt → splits 15 + 15 between P2/P3
tails. Junction overlap = 15 nt (same as deletion mode).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from primer_design.applications.tagging import run as run_tagging
from primer_design.exceptions import TagTooLongForInLocus
from primer_design.types import DesignRequest


@pytest.fixture
def all_required_files_present() -> bool:
    root = Path(__file__).resolve().parents[2]
    required = [
        "data/primer_design/vectors/pEXG2.gb",
        "data/primer_design/genes/lasR/PA14.fasta",
    ]
    return all((root / r).exists() for r in required)


@pytest.fixture
def request_obj() -> DesignRequest:
    return DesignRequest(
        isolate_id="PA14",
        gene="lasR",
        action="tag",
        tag="His6",
        tag_position="C",
        enzyme="HindIII",
        polymerase="B7",
        use_plasmid_for_tag=False,    # in-locus
    )


EXPECTED = {
    "tag_cassette": "GGCGGCAGC" + "CACCATCATCATCACCAC" + "TAA",       # 30 nt
    "junction_overlap_per_primer": 15,                                  # His6 fits 15+15
    "fusion_protein_endswith": "GGSHHHHHH*",
    "AAGCTT_count_in_final_plasmid": 0,
}


class TestPA14LasRHis6Tagging:

    def test_runs_without_error(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_tagging(request_obj)
        assert result is not None

    def test_cassette_exact(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_tagging(request_obj)
        assert result.primer_set.cassette == EXPECTED["tag_cassette"]

    def test_no_HindIII_in_final(self, all_required_files_present, request_obj):
        if not all_required_files_present:
            pytest.skip("Phase 2 prerequisites not yet in place")
        result = run_tagging(request_obj)
        seq = result.final_plasmid.sequence
        assert seq.count("AAGCTT") == EXPECTED["AAGCTT_count_in_final_plasmid"]


class TestTagTooLongRejection:
    """3xFLAG and HiBiT must be rejected with a clear suggestion."""

    def test_3xflag_rejected(self):
        request = DesignRequest(
            isolate_id="PA14",
            gene="lasR",
            action="tag",
            tag="3xFLAG",
            tag_position="C",
            enzyme="HindIII",
            polymerase="B7",
            use_plasmid_for_tag=False,
        )
        with pytest.raises(TagTooLongForInLocus) as exc_info:
            run_tagging(request)
        assert exc_info.value.error_code == "tag_too_long_for_in_locus"
        assert "suggestion" in exc_info.value.details

    def test_hibit_rejected(self):
        request = DesignRequest(
            isolate_id="PA14",
            gene="lasR",
            action="tag",
            tag="HiBiT",
            tag_position="C",
            enzyme="HindIII",
            polymerase="B7",
            use_plasmid_for_tag=False,
        )
        with pytest.raises(TagTooLongForInLocus):
            run_tagging(request)

"""Regression test for the amplicon-length bug reported 2026-05-06.

Pre-fix the tool reported 1526 bp for the LB019 lasB plasmid expression
amplicon — 15 nt short of the true 1541 bp PCR product. The cause was
``result.insert = amplicon[VECTOR_TAIL_LEN:]`` (a deliberate trim added
to make the older PA14 lasR test bounds fit), which silently dropped the
left vector homology arm from the reported insert size. The fix reports
the full PCR amplicon as ``result.insert``; this test pins the expected
1541 bp so the bug does not regress.
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
        "data/primer_design/genes/lasB/LB019.fasta",
    ]
    return all((root / r).exists() for r in required)


@pytest.fixture
def request_obj() -> DesignRequest:
    return DesignRequest(
        isolate_id="LB019",
        gene="lasB",
        action="express",
        enzyme="HindIII",
        polymerase="B7",
    )


# Expected: lasB CDS is 1497 nt across all isolates in the panel.
# Full PCR amplicon = p1.tail (29) + coding_seq (1497) + RC(p2.tail) (15)
# = 1541 bp. Final plasmid uses the site_partial_AAGCT convention which
# collapses 26 nt across the two assembly junctions (15 + 15 - 4 nt of
# shared AAGCT recognition), so length = 5148 + 1541 - 26 = 6663 bp.
EXPECTED_AMPLICON_BP = 1541
EXPECTED_FINAL_PLASMID_BP = 6663
EXPECTED_PROTEIN_TAIL = "RAFSTVGVTCPSAL*"


def test_amplicon_length_full_pcr_product(all_required_files_present, request_obj):
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    assert len(result.insert) == EXPECTED_AMPLICON_BP, (
        f"len(result.insert) = {len(result.insert)} bp, "
        f"expected {EXPECTED_AMPLICON_BP} bp (full PCR amplicon)"
    )


def test_final_plasmid_length_site_partial_convention(
    all_required_files_present, request_obj
):
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    assert result.final_plasmid.length == EXPECTED_FINAL_PLASMID_BP


def test_one_HindIII_site_in_final_plasmid(
    all_required_files_present, request_obj
):
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    assert result.final_plasmid.sequence.count("AAGCTT") == 1


def test_protein_translation_intact(all_required_files_present, request_obj):
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    from Bio.Seq import Seq
    result = run_expression(request_obj)
    protein = str(Seq(result.primer_set.coding_seq).translate())
    assert protein.endswith(EXPECTED_PROTEIN_TAIL), (
        f"Protein tail {protein[-len(EXPECTED_PROTEIN_TAIL):]!r} "
        f"!= expected {EXPECTED_PROTEIN_TAIL!r}"
    )
    # 498 aa + stop
    assert len(protein) == 499

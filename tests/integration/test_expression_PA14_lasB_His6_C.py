"""Regression: PA14 lasB pBBR1MCS2 C-terminal His6 plasmid expression.

PCR template is **native genomic DNA** (no tag). The tag-encoding nucleotides
must therefore be folded into the P2 5' overhang, and P2's body must anneal to
the native lasB CDS (the codons immediately upstream of the native stop). Before
the fix, the algorithm anchored P2 inside the engineered tag cassette, leaving
the primer with ~0 nt of homology to the genomic template — it would not have
primed at all on the bench.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq

from primer_design.applications.expression import run as run_expression
from primer_design.config import LINKER_GGS, NEW_STOP_DEFAULT
from primer_design.tags import TAGS
from primer_design.types import DesignRequest
from primer_design.vectors import reverse_complement


@pytest.fixture
def all_required_files_present() -> bool:
    root = Path(__file__).resolve().parents[2]
    required = [
        "data/primer_design/vectors/pBBR1MCS2.gb",
        "data/primer_design/genes/lasB/PA14.fasta",
    ]
    return all((root / r).exists() for r in required)


@pytest.fixture
def request_obj() -> DesignRequest:
    return DesignRequest(
        isolate_id="PA14",
        gene="lasB",
        action="express",
        tag="His6",
        tag_position="C",
        enzyme="HindIII",
        polymerase="B7",
    )


# Vector tails for pBBR1MCS2 / HindIII / expression.
P1_VECTOR_TAIL = "CGGTATCGATAAGCT"
P2_VECTOR_TAIL = "ATTCGATATCAAGCT"
RBS_TAIL_PART = "AGGAGG" + "ACTTGTTC"  # 14 nt
# Tag cassette appended to native CDS (without native stop) on the coding strand.
TAG_CASSETTE = LINKER_GGS + TAGS["His6"].dna + NEW_STOP_DEFAULT  # 30 nt

# Native lasB CDS = 1497 nt; coding_seq = 1497 − 3 + 30 = 1524 nt.
EXPECTED_CODING_SEQ_NT = 1524
# PCR amplicon = P1.tail + native_template + RC(P2.tail)
#              = (15 + 14) + (1497 − 3) + (30 + 15)
#              = 29 + 1494 + 45 = 1568 nt
EXPECTED_AMPLICON_BP = 1568


def test_p2_body_anneals_to_native_cds(all_required_files_present, request_obj):
    """P2 body must be RC of a window of the NATIVE lasB CDS (without stop).

    Pre-fix: P2 body was RC of `cds[:-3] + GGS + His6 + TAA`, which has zero
    overlap with genomic DNA at the relevant locus.
    """
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    p2_body_target = reverse_complement(result.primer_set.p2.body)
    native_cds_no_stop = result.gene_record.cds_seq[:-3]
    # The RC of the body must be a suffix of the native CDS (minus stop).
    assert native_cds_no_stop.endswith(p2_body_target), (
        "P2 body does not anneal to native genomic DNA — it lies inside the "
        "engineered tag cassette, which is not present on the genome."
    )


def test_p2_tail_encodes_tag_cassette(all_required_files_present, request_obj):
    """P2 tail = vector tail + RC of full tag cassette (GGS + His6 + TAA)."""
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    expected_p2_tail = P2_VECTOR_TAIL + reverse_complement(TAG_CASSETTE)
    assert result.primer_set.p2.tail == expected_p2_tail


def test_p1_unchanged_for_c_term_tag(all_required_files_present, request_obj):
    """For a C-terminal tag, P1 stays identical to the untagged design."""
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    assert result.primer_set.p1.tail == P1_VECTOR_TAIL + RBS_TAIL_PART
    assert result.primer_set.p1.body.startswith("ATG")


def test_amplicon_length_matches_bench(all_required_files_present, request_obj):
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    assert len(result.insert) == EXPECTED_AMPLICON_BP


def test_coding_seq_includes_tag(all_required_files_present, request_obj):
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    coding_seq = result.primer_set.coding_seq
    assert len(coding_seq) == EXPECTED_CODING_SEQ_NT
    protein = str(Seq(coding_seq).translate())
    assert protein.startswith("M")
    assert protein.endswith("GGSHHHHHH*")
    assert protein.count("*") == 1


def test_amplicon_translates_with_tag(all_required_files_present, request_obj):
    """Amplicon contains a single in-frame ORF ending in GGS-His6-stop."""
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    insert = result.insert
    atg = insert.find("ATG", insert.find("AGGAGG"))
    # Translate from ATG to first stop on the same frame.
    orf_bp = (len(insert) - atg) // 3 * 3
    protein = str(Seq(insert[atg : atg + orf_bp]).translate())
    stop_idx = protein.find("*")
    assert stop_idx > 0
    assert protein[: stop_idx + 1].endswith("GGSHHHHHH*")


def test_one_HindIII_site_in_final_plasmid(all_required_files_present, request_obj):
    if not all_required_files_present:
        pytest.skip("Prerequisites not yet in place")
    result = run_expression(request_obj)
    assert result.final_plasmid.sequence.count("AAGCTT") == 1

"""Tests for the pre-validated ("fixed") primer-set lookup.

These only exercise the local, network-free matching logic (curated gene
records + vector GenBank files, both bundled in the repo) — no R2 genome
fetch, so they run fast and offline. The full lasB deletion path (which needs
a genome-wide search for P4, ~628 nt downstream of the stop — outside the
curated dn_flank window) is covered end-to-end by
tests/integration/test_deletion_LB001_lasB.py instead, since that requires
R2 access; here we only check the graceful, no-network-available fallback.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from primer_design import fixed_primers, gene_finder, vectors


REQUIRED = [
    "data/primer_design/vectors/pEXG2.gb",
    "data/primer_design/vectors/pBBR1MCS2.gb",
    "data/primer_design/genes/lasR/LB001.fasta",
    "data/primer_design/genes/lasI/LB001.fasta",
    "data/primer_design/genes/lasB/LB019.fasta",
]


@pytest.fixture
def all_required_files_present() -> bool:
    root = Path(__file__).resolve().parents[2]
    return all((root / rel).exists() for rel in REQUIRED)


class TestFixedDeletionSet:

    def test_lasR_matches_LB001_with_calibrated_scar(self, all_required_files_present):
        if not all_required_files_present:
            pytest.skip("required data files not present")
        gene = gene_finder.get_gene_record("LB001", "lasR")
        vector = vectors.get_vector("pEXG2")
        convention = vectors.get_tail_convention("pEXG2", "HindIII", "deletion")
        result = fixed_primers.build_fixed_deletion_set(gene, vector, convention)
        assert result is not None
        # Matches the N=6 / C=13 calibration noted in config.py.
        assert result.primer_set.N == 6
        assert result.primer_set.C == 13
        assert result.primer_set.p1.body == "CCGTTGCAGGCGCTGTTCGG"
        assert result.primer_set.p1.tail == "GCATAAATGTAAAGCAAGCT"
        assert result.expected_recognition_count == 0

    def test_lasI_matches_LB001_with_calibrated_scar(self, all_required_files_present):
        if not all_required_files_present:
            pytest.skip("required data files not present")
        gene = gene_finder.get_gene_record("LB001", "lasI")
        vector = vectors.get_vector("pEXG2")
        convention = vectors.get_tail_convention("pEXG2", "HindIII", "deletion")
        result = fixed_primers.build_fixed_deletion_set(gene, vector, convention)
        assert result is not None
        assert result.primer_set.N == 8
        assert result.primer_set.C == 6

    def test_gene_without_fixed_entry_returns_none(self, all_required_files_present):
        if not all_required_files_present:
            pytest.skip("required data files not present")
        gene = gene_finder.get_gene_record("LB001", "lasR")
        vector = vectors.get_vector("pEXG2")
        convention = vectors.get_tail_convention("pEXG2", "HindIII", "deletion")
        # rhlR has no fixed-primer entry on file.
        assert fixed_primers._get_bodies("rhlR", "deletion") is None

    def test_lasB_needs_genome_widening_and_returns_none_without_it(
        self, all_required_files_present,
    ):
        if not all_required_files_present:
            pytest.skip("required data files not present")
        gene = gene_finder.get_gene_record("LB001", "lasB")
        vector = vectors.get_vector("pEXG2")
        convention = vectors.get_tail_convention("pEXG2", "HindIII", "deletion")
        # P4 anneals ~628 nt downstream of the stop, outside the curated
        # 600 nt dn_flank, and no genome_bytes is supplied here — must not
        # raise, and must fall back cleanly (None) rather than guess.
        assert fixed_primers.build_fixed_deletion_set(gene, vector, convention) is None


class TestFixedExpressionSet:

    def test_lasB_matches_LB019(self, all_required_files_present):
        if not all_required_files_present:
            pytest.skip("required data files not present")
        gene = gene_finder.get_gene_record("LB019", "lasB")
        vector = vectors.get_vector("pBBR1MCS2")
        convention = vectors.get_tail_convention("pBBR1MCS2", "HindIII", "expression")
        result = fixed_primers.build_fixed_expression_set(gene, vector, convention)
        assert result is not None
        ps = result.primer_set
        assert ps.p1.body == "TGAACAAGATGAAGAAGGTTTCTA"
        assert ps.p2.body == "TTACAACGCGCTCGGGCAGG"
        # Vector arm auto-derived (20 nt) + bare RBS (no synthetic spacer —
        # the native upstream context between RBS and ATG lives in the body).
        assert ps.p1.tail == "GTCGACGGTATCGATAAGCT" + "AGGAGG"
        assert ps.native_anneal_segment.startswith("TGAACAAGATG")
        assert ps.coding_seq == gene.cds_seq

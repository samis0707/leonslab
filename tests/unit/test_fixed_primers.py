"""Unit tests for the curated fixed-primer override (lasB/pEXG2)."""
from __future__ import annotations

from dataclasses import replace

import pytest

from primer_design import fixed_primers, gene_finder, vectors

# 700 nt of PA14 genomic sequence starting right after the lasB stop codon
# (contig NZ_CP104983.1:1389012-1389711), i.e. the bundled 600 bp dn_flank plus
# 100 more nt. Contains the real P4 anchor (RC) at offset 629, just past the
# bundled window -- this is what exercises the genome-extension path.
_PA14_DN_EXTENSION = (
    "GCTCGGTGGCCCCGGCCGGCACTCCAGGAAGGAATGCCGGCCGGGGCCGCTCAAGCCGTCTTCCGCCAGGAGGACGGCTGCTTTATGTCGCTTAGGCCGTTGGCCTCCGCGAACCCGGTCTAAAGTTCAGGTGTGAGCATTTATTCCAAGACCGACCGGGAGTCCTGCCATGAGTCTGCTGTTCGAGCCTCTTAGCCTGCGTCAAATCACCTTGCCCAACCGCATCGCCGTATCGCCCATGTGCCAGTATTCGGCGCAGGAGGGCCTGGCCAACGACTGGCATCTCGTGCACCTGGGCAGTCGCGCGGTGGGCGGCGCCGGCCTGGTGATCGTCGAAGCCACCGCGGTGTTGCCCGAGGGGCGCATCACCGCCGACGACCTCGGCATCTGGAGCGACGCGCATGTCGAGCCGTTGCATCGCATCACCCGTTTCATCGAGTCCCAGGGCGCGGTCGCCGGGGTCCAGCTGGCCCACGCCGGGCGCAAGGCGAGTACCTGGCGGCCGTGGCTGGGCAAGCACGGCAGCGTGCCGATAGGCCAGGGCGGCTGGATACCGGTGGCGCCGTCGGCGATCCCGTTCGATCCCCAGCACACGACCCCCGAGGCTCTGAGCGAGGCGCAAATCGAGGCGCTGGTGCAGGCCTTCGTACGCGCCACCGAGCGCTCCCTGGCTGCCGGCTTCAAGGTCGCGGAGGTGCAT"
)


@pytest.fixture
def pa14_lasb():
    return gene_finder.get_gene_record("PA14", "lasB")


@pytest.fixture
def pexg2():
    return vectors.get_vector("pEXG2")


def _synthetic_genome_bytes(gene) -> bytes:
    """Build a minimal single-contig FASTA covering the gene's up_flank + cds +
    the extended dn region (700 nt), matching gene.contig_id, for the parts of
    try_build_fixed_primer_set that need to fetch beyond the bundled dn_flank."""
    seq = gene.up_flank + gene.cds_seq + _PA14_DN_EXTENSION
    return f">{gene.contig_id}\n{seq}\n".encode("ascii")


def _with_local_coords(gene):
    """Re-anchor genome_start/end_1based to this gene's own local up_flank/cds
    layout, so they index correctly into the small synthetic single-gene contig
    built by _synthetic_genome_bytes (rather than the real, much larger genome
    the bundled record's coordinates point into)."""
    up_len = len(gene.up_flank)
    cds_len = len(gene.cds_seq)
    return replace(
        gene,
        genome_start_1based=up_len + 1,
        genome_end_1based=up_len + cds_len,
    )


def test_spec_loads_for_lasb_pexg2():
    spec = fixed_primers.load_fixed_primer_spec("lasB", "pEXG2")
    assert spec is not None
    assert spec.N == 13
    assert spec.C == 12
    assert spec.p1_body == "CGCCTTCGATGCCGAAGTACG"


def test_spec_missing_for_unconfigured_pair():
    assert fixed_primers.load_fixed_primer_spec("lasR", "pEXG2") is None
    assert fixed_primers.load_fixed_primer_spec("lasB", "pBBR1MCS2") is None


def test_returns_none_without_genome_bytes(pa14_lasb, pexg2):
    """P4 anchors past the bundled 600 bp dn_flank, so without a genome to
    extend into, the fixed set must not apply (safe fallback to computed)."""
    primer_set, gene_used = fixed_primers.try_build_fixed_primer_set(
        pa14_lasb, pexg2, None
    )
    assert primer_set is None
    assert gene_used is pa14_lasb


def test_matches_pa14_with_extended_genome(pa14_lasb, pexg2):
    genome_bytes = _synthetic_genome_bytes(pa14_lasb)
    gene = _with_local_coords(pa14_lasb)
    primer_set, gene_used = fixed_primers.try_build_fixed_primer_set(
        gene, pexg2, genome_bytes
    )
    assert primer_set is not None
    assert primer_set.N == 13
    assert primer_set.C == 12
    assert primer_set.p1.body == "CGCCTTCGATGCCGAAGTACG"
    assert primer_set.p4.body == "TACGAAGGCCTGCACCAGCGC"
    assert len(gene_used.dn_flank) > len(pa14_lasb.dn_flank)

    from primer_design._bio_lite import translate_dna
    assert translate_dna(primer_set.scar_dna) == "MKKVSTLDLLFVAFSTVGVTCPSAL*"


def test_p1_mismatch_falls_back_to_none(pa14_lasb, pexg2):
    """If the isolate's up_flank doesn't carry the fixed P1 body verbatim
    (simulated here by truncating up_flank so it can't contain it), no match."""
    tampered = replace(pa14_lasb, up_flank=pa14_lasb.up_flank[-50:])
    genome_bytes = _synthetic_genome_bytes(pa14_lasb)
    primer_set, _ = fixed_primers.try_build_fixed_primer_set(
        tampered, pexg2, genome_bytes
    )
    assert primer_set is None


def test_find_contig_tolerates_isolate_prefix():
    contigs = {"PA14_NZ_CP104983.1": "ACGT"}
    assert fixed_primers._find_contig(contigs, "NZ_CP104983.1") == "ACGT"
    assert fixed_primers._find_contig(contigs, "NOPE") is None

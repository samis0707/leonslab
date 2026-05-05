"""Tests for gene record loading from canonical FASTA files."""
from __future__ import annotations

from pathlib import Path

import pytest

from primer_design.exceptions import CdsValidationError, GeneRecordNotFound
from primer_design.gene_finder import (
    list_available_isolates,
    list_supported_genes,
    parse_header_metadata,
)


class TestParseHeaderMetadata:

    def test_canonical_header(self):
        header = "LB001|lasB|contig=LB001_00001|start=1848921|end=1850417|strand=+|cds_len=1497|up=600|dn=600"
        meta = parse_header_metadata(header)
        assert meta["contig"] == "LB001_00001"
        assert meta["start"] == "1848921"
        assert meta["end"] == "1850417"
        assert meta["strand"] == "+"
        assert meta["cds_len"] == "1497"
        assert meta["up"] == "600"
        assert meta["dn"] == "600"

    def test_minus_strand(self):
        header = "PA14|lasR|contig=PA14_chrom|start=10000|end=11000|strand=-|cds_len=720|up=600|dn=600"
        meta = parse_header_metadata(header)
        assert meta["strand"] == "-"

    def test_missing_required_field_raises(self):
        header = "LB001|lasB|contig=LB001_00001"  # missing start/end/strand/etc.
        with pytest.raises(GeneRecordNotFound):
            parse_header_metadata(header)


class TestListSupportedGenes:

    def test_returns_configured_genes(self):
        genes = list_supported_genes()
        assert "lasB" in genes
        assert "lasR" in genes


class TestListAvailableIsolates:

    def test_empty_when_no_records_yet(self, tmp_path):
        # Phase 2 step 1 hasn't run yet → no records on disk
        result = list_available_isolates("lasB", genes_dir=tmp_path)
        assert result == []

    def test_finds_records(self, tmp_path):
        gene_dir = tmp_path / "lasB"
        gene_dir.mkdir()
        (gene_dir / "LB001.fasta").write_text(">LB001|lasB|...\nACGT\n")
        (gene_dir / "LB014.fasta").write_text(">LB014|lasB|...\nACGT\n")
        result = list_available_isolates("lasB", genes_dir=tmp_path)
        assert result == ["LB001", "LB014"]


# Note: Full integration tests of get_gene_record() against real canonical FASTA
# files live in tests/integration/ and run after Phase 2 step 1 builds them.

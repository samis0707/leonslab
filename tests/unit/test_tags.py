"""Tag-library tests. These are FULLY IMPLEMENTED and should pass immediately."""
from __future__ import annotations

import pytest
from Bio.Seq import Seq

from primer_design.config import LINKER_GGS, NEW_STOP_DEFAULT, TAGGING_MAX_CASSETTE_NT
from primer_design.exceptions import (
    InvalidTagPosition,
    TagTooLongForInLocus,
    UnknownTag,
)
from primer_design.tags import (
    TAGS,
    build_in_locus_cassette,
    build_plasmid_fusion_cds,
    get_tag,
    validate_tag_position,
)


class TestTagLibrary:

    def test_all_predefined_tags_translate_correctly(self):
        for tag in TAGS.values():
            translated = str(Seq(tag.dna).translate())
            assert translated == tag.protein, f"{tag.name}: {translated} != {tag.protein}"

    def test_all_predefined_tags_have_inframe_dna(self):
        for tag in TAGS.values():
            assert len(tag.dna) % 3 == 0, f"{tag.name}: DNA length {len(tag.dna)} not multiple of 3"

    def test_all_predefined_tag_dna_has_no_internal_stop(self):
        for tag in TAGS.values():
            translated = str(Seq(tag.dna).translate())
            assert "*" not in translated, f"{tag.name}: contains internal stop"

    def test_cassette_nt_is_consistent(self):
        for tag in TAGS.values():
            expected = len(LINKER_GGS) + len(tag.dna) + len(NEW_STOP_DEFAULT)
            assert tag.cassette_nt == expected, f"{tag.name}: {tag.cassette_nt} != {expected}"

    def test_in_locus_ok_threshold(self):
        for tag in TAGS.values():
            should_be_ok = tag.cassette_nt <= TAGGING_MAX_CASSETTE_NT
            assert tag.in_locus_ok is should_be_ok, f"{tag.name}: in_locus_ok={tag.in_locus_ok}"


class TestGetTag:

    def test_known(self):
        assert get_tag("His6").name == "His6"
        assert get_tag("FLAG").name == "FLAG"

    def test_unknown_raises(self):
        with pytest.raises(UnknownTag):
            get_tag("nonexistent_tag")


class TestBuildInLocusCassette:

    def test_his6_assembles_to_30nt(self):
        cassette = build_in_locus_cassette(get_tag("His6"))
        assert len(cassette) == 30
        assert cassette.startswith(LINKER_GGS)
        assert cassette.endswith(NEW_STOP_DEFAULT)
        protein = str(Seq(cassette).translate())
        assert protein == "GGSHHHHHH*"

    def test_flag_assembles_to_36nt(self):
        cassette = build_in_locus_cassette(get_tag("FLAG"))
        assert len(cassette) == 36
        protein = str(Seq(cassette).translate())
        assert protein == "GGSDYKDDDDK*"

    def test_his8_assembles_to_36nt(self):
        cassette = build_in_locus_cassette(get_tag("His8"))
        assert len(cassette) == 36
        protein = str(Seq(cassette).translate())
        assert protein == "GGSHHHHHHHH*"

    def test_3xflag_rejected_for_in_locus(self):
        with pytest.raises(TagTooLongForInLocus) as exc_info:
            build_in_locus_cassette(get_tag("3xFLAG"))
        assert "3xFLAG" in exc_info.value.message
        assert exc_info.value.details["limit"] == TAGGING_MAX_CASSETTE_NT

    def test_hibit_rejected_for_in_locus(self):
        with pytest.raises(TagTooLongForInLocus):
            build_in_locus_cassette(get_tag("HiBiT"))


class TestValidateTagPosition:

    def test_3xflag_rejected_n_term(self):
        with pytest.raises(InvalidTagPosition):
            validate_tag_position(get_tag("3xFLAG"), "N")

    def test_3xflag_ok_c_term(self):
        validate_tag_position(get_tag("3xFLAG"), "C")  # no error

    def test_hibit_rejected_n_term(self):
        with pytest.raises(InvalidTagPosition):
            validate_tag_position(get_tag("HiBiT"), "N")

    def test_his6_ok_both(self):
        validate_tag_position(get_tag("His6"), "N")
        validate_tag_position(get_tag("His6"), "C")


class TestBuildPlasmidFusionCds:

    def test_no_tag_returns_cds_unchanged(self):
        cds = "ATGAAACCGTTGTGA"
        assert build_plasmid_fusion_cds(cds, None, None) == cds

    def test_n_terminal_replaces_native_atg(self):
        cds = "ATGAAACCGTTGTGA"
        result = build_plasmid_fusion_cds(cds, get_tag("His6"), "N")
        # Result: ATG + His6_DNA + LINKER + (CDS without first ATG)
        expected = "ATG" + "CACCATCATCATCACCAC" + LINKER_GGS + "AAACCGTTGTGA"
        assert result == expected
        # Translates to: M + HHHHHH + GGS + native-protein-starting-after-M
        protein = str(Seq(result).translate())
        assert protein.startswith("MHHHHHHGGS")
        assert protein.endswith("*")

    def test_c_terminal_replaces_native_stop(self):
        cds = "ATGAAACCGTTGTGA"
        result = build_plasmid_fusion_cds(cds, get_tag("His6"), "C")
        # Result: CDS_without_stop + LINKER + His6_DNA + new_stop
        expected = "ATGAAACCGTTG" + LINKER_GGS + "CACCATCATCATCACCAC" + "TAA"
        assert result == expected
        protein = str(Seq(result).translate())
        assert protein.endswith("GGSHHHHHH*")

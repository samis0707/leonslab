"""Tag library (D6.2.v2).

Codon-optimized for *P. aeruginosa*. Linker is GGS (`GGCGGCAGC`).

Cassette length = ``len(LINKER_GGS) + len(tag.dna) + 3``  (3 = new stop codon).
``in_locus_ok`` is True iff cassette ≤ 36 nt — the maximum the 4-primer in-locus
design can encode in v1 (D6.4).
"""
from __future__ import annotations

from .config import LINKER_GGS, NEW_STOP_DEFAULT, TAGGING_MAX_CASSETTE_NT
from .exceptions import InvalidTagPosition, TagTooLongForInLocus, UnknownTag
from .types import Tag, TagPosition


def _cassette_nt(dna: str) -> int:
    return len(LINKER_GGS) + len(dna) + 3


# Pre-defined tag library.
# DNA sequences chosen for *P. aeruginosa* high-frequency codons.
TAGS: dict[str, Tag] = {
    "FLAG": Tag(
        name="FLAG",
        dna="GACTACAAGGACGACGATGACAAG",                                    # DYKDDDDK, 24 nt
        protein="DYKDDDDK",
        n_term_ok=True,
        c_term_ok=True,
        in_locus_ok=True,
        cassette_nt=_cassette_nt("GACTACAAGGACGACGATGACAAG"),              # 36
    ),
    "3xFLAG": Tag(
        name="3xFLAG",
        dna="GACTACAAGGACCACGACGGCGACTACAAGGATCATGATATCGATTACAAGGATGACGATGACAAG",  # 78 nt
        protein="DYKDHDGDYKDHDIDYKDDDDK",
        n_term_ok=False,
        c_term_ok=True,
        in_locus_ok=False,                                                 # cassette = 90 nt
        cassette_nt=_cassette_nt("GACTACAAGGACCACGACGGCGACTACAAGGATCATGATATCGATTACAAGGATGACGATGACAAG"),
    ),
    "His6": Tag(
        name="His6",
        dna="CACCATCATCATCACCAC",                                         # HHHHHH, 18 nt
        protein="HHHHHH",
        n_term_ok=True,
        c_term_ok=True,
        in_locus_ok=True,
        cassette_nt=_cassette_nt("CACCATCATCATCACCAC"),                   # 30
    ),
    "His8": Tag(
        name="His8",
        dna="CATCATCATCATCATCATCATCAT",                                   # HHHHHHHH, 24 nt
        protein="HHHHHHHH",
        n_term_ok=True,
        c_term_ok=True,
        in_locus_ok=True,
        cassette_nt=_cassette_nt("CATCATCATCATCATCATCATCAT"),             # 36
    ),
    "HiBiT": Tag(
        name="HiBiT",
        dna="GTGAGCGGCTGGCGGCTGTTCAAGAAGATCAGC",                          # VSGWRLFKKIS, 33 nt
        protein="VSGWRLFKKIS",
        n_term_ok=False,
        c_term_ok=True,
        in_locus_ok=False,                                                # cassette = 45 nt
        cassette_nt=_cassette_nt("GTGAGCGGCTGGCGGCTGTTCAAGAAGATCAGC"),
    ),
}


def get_tag(name: str) -> Tag:
    if name not in TAGS:
        raise UnknownTag(message=f"Unknown tag '{name}'. Available: {list(TAGS.keys())}")
    return TAGS[name]


def build_in_locus_cassette(tag: Tag) -> str:
    """Build the in-locus C-terminal tag cassette (LINKER + tag DNA + new stop).

    Returns the 30/36 nt cassette in coding orientation. Raises if tag too long.
    """
    if not tag.in_locus_ok:
        raise TagTooLongForInLocus(
            message=(
                f"Tag '{tag.name}' cassette is {tag.cassette_nt} nt, exceeds "
                f"{TAGGING_MAX_CASSETTE_NT} nt limit for 4-primer in-locus design."
            ),
            details={
                "tag": tag.name,
                "cassette_nt": tag.cassette_nt,
                "limit": TAGGING_MAX_CASSETTE_NT,
                "suggestion": "Use plasmid expression instead (action='express', tag_position='C').",
            },
        )
    cassette = LINKER_GGS + tag.dna + NEW_STOP_DEFAULT
    assert len(cassette) == tag.cassette_nt
    assert len(cassette) % 3 == 0
    return cassette


def validate_tag_position(tag: Tag, position: TagPosition) -> None:
    """Raise InvalidTagPosition if ``position`` is not allowed for ``tag``."""
    if position == "N" and not tag.n_term_ok:
        raise InvalidTagPosition(
            message=f"Tag '{tag.name}' is not valid for N-terminal fusion.",
            details={"tag": tag.name, "requested_position": position},
        )
    if position == "C" and not tag.c_term_ok:
        raise InvalidTagPosition(
            message=f"Tag '{tag.name}' is not valid for C-terminal fusion.",
            details={"tag": tag.name, "requested_position": position},
        )


def build_plasmid_fusion_cds(
    cds_seq: str,
    tag: Tag | None,
    position: TagPosition | None,
) -> str:
    """Construct the full ATG-to-stop CDS for plasmid expression with optional tag.

    Plasmid expression preserves stop codon at the very end. Native ATG is dropped
    when N-terminal tag is added; native stop is dropped when C-terminal tag is added.
    """
    if tag is None:
        return cds_seq

    validate_tag_position(tag, position)

    if position == "N":
        # M-tag-GGS-CDS_from_codon_2(stop preserved)
        return "ATG" + tag.dna + LINKER_GGS + cds_seq[3:]

    if position == "C":
        # CDS_to_last_codon-GGS-tag-newstop
        return cds_seq[:-3] + LINKER_GGS + tag.dna + NEW_STOP_DEFAULT

    raise InvalidTagPosition(message=f"Unknown tag_position {position!r}")

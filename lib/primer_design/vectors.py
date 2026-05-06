"""Vector library: load GenBank, find cut sites, derive primer tails.

Tail rules are named primitives (strings) that index into ``TAIL_RULES``. The
calibrated configuration lives in ``data/primer_design/tail_conventions.json``
and is loaded by ``get_tail_convention``.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

from Bio import SeqIO

from .config import RESTRICTION_SITES, SUPPORTED_VECTORS, VECTOR_TAIL_LEN
from .exceptions import (
    ConventionNotCalibrated,
    NotSingleCutter,
    UnknownEnzyme,
    UnknownVector,
)
from .types import Application, TailConvention, VectorRecord


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def reverse_complement(seq: str) -> str:
    """Return reverse complement of an uppercase DNA string. ACGT only."""
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def find_all(haystack: str, needle: str) -> list[int]:
    """Return all 0-based start positions of ``needle`` in ``haystack``, non-overlapping."""
    out: list[int] = []
    i = haystack.find(needle, 0)
    while i != -1:
        out.append(i)
        i = haystack.find(needle, i + 1)
    return out


# ---------------------------------------------------------------------------
# Vector loading
# ---------------------------------------------------------------------------

# Cache of parsed vectors. Cleared between deployments.
_VECTOR_CACHE: dict[str, VectorRecord] = {}


def get_vector(name: str, vectors_dir: Path | None = None) -> VectorRecord:
    """Load a vector from its bundled GenBank file.

    Args:
        name: ``"pEXG2"`` or ``"pBBR1MCS2"``.
        vectors_dir: override for tests; defaults to ``data/primer_design/vectors``.

    Raises:
        UnknownVector: if ``name`` is not in ``SUPPORTED_VECTORS``.
    """
    if name in _VECTOR_CACHE:
        return _VECTOR_CACHE[name]

    if name not in SUPPORTED_VECTORS:
        raise UnknownVector(message=f"Unknown vector '{name}'. Supported: {SUPPORTED_VECTORS}")

    vectors_dir = vectors_dir or _default_vectors_dir()
    gb_path = vectors_dir / f"{name}.gb"
    if not gb_path.exists():
        # Allow FASTA fallback during early Phase-2 work, before GenBank annotations exist
        fasta_path = vectors_dir / f"{name}.fasta"
        if fasta_path.exists():
            record = next(SeqIO.parse(fasta_path, "fasta"))
            features = []
        else:
            raise UnknownVector(
                message=f"Vector file not found: {gb_path} (or {fasta_path})",
                details={"vector": name, "search_dir": str(vectors_dir)},
            )
    else:
        record = next(SeqIO.parse(gb_path, "genbank"))
        features = list(record.features)

    vec = VectorRecord(
        name=name,
        sequence=str(record.seq).upper(),
        length=len(record.seq),
        features=features,
    )
    _VECTOR_CACHE[name] = vec
    return vec


def _default_vectors_dir() -> Path:
    """Resolve ``data/primer_design/vectors`` relative to the package."""
    return Path(__file__).resolve().parents[2] / "data" / "primer_design" / "vectors"


# ---------------------------------------------------------------------------
# Cut-site detection
# ---------------------------------------------------------------------------

def find_cut_position(vector: VectorRecord, enzyme: str) -> int:
    """Return the 0-based top-strand nick position for ``enzyme`` in ``vector``.

    Raises:
        UnknownEnzyme: if ``enzyme`` not in ``RESTRICTION_SITES``.
        NotSingleCutter: if zero or multiple sites in the vector.
    """
    if enzyme not in RESTRICTION_SITES:
        raise UnknownEnzyme(message=f"Unknown enzyme '{enzyme}'. Defined: {list(RESTRICTION_SITES)}")
    motif, nick_offset = RESTRICTION_SITES[enzyme]
    rc_motif = reverse_complement(motif)

    fwd_hits = find_all(vector.sequence, motif)
    # For palindromic motifs (most common 6-cutters), forward and rc are identical → don't double-count
    if motif == rc_motif:
        all_hits = fwd_hits
    else:
        rev_hits = find_all(vector.sequence, rc_motif)
        all_hits = sorted(set(fwd_hits) | set(rev_hits))

    if len(all_hits) != 1:
        raise NotSingleCutter(
            message=f"{enzyme} cuts {vector.name} {len(all_hits)} times; need exactly 1",
            details={"vector": vector.name, "enzyme": enzyme, "site_count": len(all_hits), "positions": all_hits},
        )
    return all_hits[0] + nick_offset


def get_unique_cutters(vector: VectorRecord) -> list[str]:
    """Enumerate enzymes from ``RESTRICTION_SITES`` that cut ``vector`` exactly once.

    Used by the frontend to populate the enzyme dropdown for a chosen vector.
    """
    out: list[str] = []
    for enzyme in RESTRICTION_SITES:
        try:
            find_cut_position(vector, enzyme)
            out.append(enzyme)
        except NotSingleCutter:
            continue
    return sorted(out)


# ---------------------------------------------------------------------------
# Tail derivation
# ---------------------------------------------------------------------------

# Each rule takes (vector, motif, nick_position_0based) and returns a 15-nt tail string.
# The four rules below cover the calibrated triples in tail_conventions.json.
# Adding a new rule = one new function here + a JSON entry referencing its name.

TailRuleFn = Callable[[VectorRecord, str, int], str]


def _rule_left_arm_15nt_excl_recognition(v: VectorRecord, motif: str, nick: int) -> str:
    """P1 tail for ``site_destroyed`` convention (pEXG2 deletion / tagging).

    Example (HindIII A^AGCTT, recognition starts at nick-1):
    Returns the 15 nt immediately before the nick (excludes recognition site entirely
    on the left arm side, since the nick is between A and AGCTT — "left arm" ends at nick).
    """
    assert nick >= VECTOR_TAIL_LEN, "vector too short before cut for 15 nt left tail"
    return v.sequence[nick - VECTOR_TAIL_LEN : nick]


def _rule_rc_right_arm_15nt_starting_at_second_nt_of_recognition(
    v: VectorRecord, motif: str, nick: int
) -> str:
    """P4 tail for ``site_destroyed`` (pEXG2 deletion / tagging).

    Returns RC of the 15 nt starting at the nick position. For HindIII, this means
    starting at the second A of AAGCTT (= AGCTT...). Captures 5 nt of recognition
    on the right arm + 10 nt past — but in RC, so the recognition is destroyed.
    """
    return reverse_complement(v.sequence[nick : nick + VECTOR_TAIL_LEN])


def _rule_rc_right_arm_15nt_incl_5nt_recognition(v: VectorRecord, motif: str, nick: int) -> str:
    """``site_partial_AAGCT`` tail derived from the right arm.

    Returns RC of the 15 nt starting at the nick position. Same physical span as
    the deletion P4 rule above (15 nt RC at nick). For Addgene's pBBR1MCS-2
    sequence (Addgene #85168, MCS oriented ...SalI-ClaI-HindIII-EcoRV-EcoRI...),
    this is the **P4** tail in the expression triple. The recognition site is
    regenerated at one junction in the final plasmid (the P1 vs P4 pairing
    contributes the matching half of `AAGCT` from each side).
    """
    return reverse_complement(v.sequence[nick : nick + VECTOR_TAIL_LEN])


def _rule_left_arm_15nt_incl_5nt_recognition(v: VectorRecord, motif: str, nick: int) -> str:
    """``site_partial_AAGCT`` tail derived from the left arm.

    Returns 15 nt of the forward strand: 10 nt of left-arm context immediately
    before the recognition site + the first 5 nt of the recognition site
    (for HindIII A^AGCTT this is `AAGCT`). For Addgene's pBBR1MCS-2 this is
    the **P1** tail in the expression triple (matches the empirical PA14
    lasB primer #1573 tail `CGGTATCGATAAGCT`).

    The slice is anchored on the recognition start (`s = nick - 1` for HindIII)
    rather than directly on `nick`, so the 5 nt of recognition included are
    the first 5 of the motif, not the last 5.
    """
    s = nick - 1                                              # recognition start
    return v.sequence[s - 10 : s + 5]                         # 10 nt context + AAGCT


TAIL_RULES: dict[str, TailRuleFn] = {
    "left_arm_15nt_excl_recognition": _rule_left_arm_15nt_excl_recognition,
    "rc_right_arm_15nt_starting_at_second_nt_of_recognition":
        _rule_rc_right_arm_15nt_starting_at_second_nt_of_recognition,
    "rc_right_arm_15nt_incl_5nt_recognition": _rule_rc_right_arm_15nt_incl_5nt_recognition,
    "left_arm_15nt_incl_5nt_recognition": _rule_left_arm_15nt_incl_5nt_recognition,
}


def derive_tails(vector: VectorRecord, enzyme: str, application: Application) -> tuple[str, str]:
    """Return (P1_tail, P4_tail) for the given vector / enzyme / application.

    Reads the calibrated convention from ``tail_conventions.json``, applies the named
    rules, and returns the literal 15-nt tail strings.
    """
    convention = get_tail_convention(vector.name, enzyme, application)
    motif, _ = RESTRICTION_SITES[enzyme]
    nick = find_cut_position(vector, enzyme)

    p1_rule = TAIL_RULES[convention.p1_tail_rule]
    p4_rule = TAIL_RULES[convention.p4_tail_rule]
    return p1_rule(vector, motif, nick), p4_rule(vector, motif, nick)


# ---------------------------------------------------------------------------
# Tail-convention loading
# ---------------------------------------------------------------------------

_CONVENTIONS_CACHE: dict[str, dict] | None = None


def get_tail_convention(vector: str, enzyme: str, application: Application) -> TailConvention:
    """Look up the calibrated convention for the (vector, enzyme, application) triple.

    If the triple is not calibrated, falls back to ``site_destroyed`` (deletion/tagging
    style) but emits no warning here — verification.py will catch misalignment.

    Raises:
        ConventionNotCalibrated: if neither the triple nor the fallback exists.
    """
    global _CONVENTIONS_CACHE
    if _CONVENTIONS_CACHE is None:
        path = _default_vectors_dir().parent / "tail_conventions.json"
        with open(path) as f:
            _CONVENTIONS_CACHE = json.load(f)

    key = f"{vector}|{enzyme}|{application}"
    fallback_key = "default|site_destroyed"

    if key in _CONVENTIONS_CACHE:
        entry = _CONVENTIONS_CACHE[key]
    elif fallback_key in _CONVENTIONS_CACHE:
        entry = _CONVENTIONS_CACHE[fallback_key]
    else:
        raise ConventionNotCalibrated(
            message=f"No tail convention for {key} and no fallback configured.",
            details={"key": key, "available": list(_CONVENTIONS_CACHE.keys())},
        )

    return TailConvention(
        name=entry["name"],
        p1_tail_rule=entry["p1_tail_rule"],
        p4_tail_rule=entry["p4_tail_rule"],
        expected_recognition_count_in_final_plasmid=entry["expected_recognition_count_in_final_plasmid"],
        enzyme=enzyme,
        vector=vector,
        application=application,
    )

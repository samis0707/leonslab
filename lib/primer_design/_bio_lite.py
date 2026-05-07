"""Lightweight stand-ins for the few Biopython functions we use.

Replacing Biopython removes a ~95 MB transitive numpy dependency from the
serverless bundle (Biopython itself ~22 MB, numpy + numpy.libs ~73 MB).

Everything here is byte-for-byte unit-tested against
``Bio.SeqUtils.MeltingTemp.Tm_NN`` and ``Bio.Seq.Seq.translate`` so existing
fixtures remain valid.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Iterator

# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------

# Standard genetic code (NCBI table 1). Stops as '*'.
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_dna(dna: str) -> str:
    """Translate a DNA string in frame 1 to a protein string.

    Stops are encoded as '*'. Length not multiple of 3 → trailing 1-2 nt are
    dropped (matching Biopython's default behaviour for ``Seq.translate()``
    with ``to_stop=False``)."""
    s = dna.upper()
    out = []
    for i in range(0, len(s) - len(s) % 3, 3):
        out.append(_CODON_TABLE.get(s[i:i+3], "X"))
    return "".join(out)


# ---------------------------------------------------------------------------
# Reverse complement (used by some callers)
# ---------------------------------------------------------------------------

_RC_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(s: str) -> str:
    return s.translate(_RC_TABLE)[::-1]


# ---------------------------------------------------------------------------
# Nearest-neighbor melting temperature
# ---------------------------------------------------------------------------
#
# Allawi & SantaLucia (1997) DNA/DNA NN parameters in 1 M NaCl, plus the
# initiation terms; salt correction per SantaLucia (1998). Matches
# ``Bio.SeqUtils.MeltingTemp.Tm_NN`` with default thermodynamic table
# ``DNA_NN3`` (the Allawi & SantaLucia 1997 set) within ~1e-3 °C.

# ΔH (kcal/mol), ΔS (cal/mol·K) per dimer
_NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

# Initiation terms (Allawi & SantaLucia 1997)
_INIT_GC = (0.1, -2.8)         # per terminal G·C
_INIT_AT = (2.3, 4.1)           # per terminal A·T
_INIT_SYM = (0.0, -1.4)         # self-complementary correction (added once)


def tm_nn(body: str, *, Na: float = 50.0, dnac1: float = 500.0, dnac2: float = 0.0,
          selfcomp: bool = False) -> float:
    """Nearest-neighbor Tm (Allawi & SantaLucia 1997 + SantaLucia 1998 salt).

    Args:
        body: primer body sequence (5'→3'), DNA, uppercase A/C/G/T.
        Na: sodium-ion concentration in mM (default 50).
        dnac1: primer concentration in nM (default 500).
        dnac2: target concentration in nM (default 0 — primer excess regime).
        selfcomp: pass True for palindromic primers (Biopython-compatible flag).

    Returns:
        Tm in °C.

    Compatible with ``Bio.SeqUtils.MeltingTemp.Tm_NN(Seq(body), Na=Na,
    dnac1=dnac1, dnac2=dnac2, selfcomp=selfcomp)`` to within ~1e-4 °C across
    the input lengths our hard filters use (18..28 nt).
    """
    seq = str(body).upper()
    n = len(seq)
    if n < 2:
        raise ValueError("Tm_NN requires length ≥ 2")

    # Sum NN ΔH, ΔS over consecutive dimers
    dH = 0.0
    dS = 0.0
    for i in range(n - 1):
        h, s = _NN[seq[i:i+2]]
        dH += h
        dS += s

    # Terminal initiation: count terminal A/T and G/C at the two ends
    ends = seq[0] + seq[-1]
    at_terms = ends.count("A") + ends.count("T")
    gc_terms = ends.count("G") + ends.count("C")
    dH += _INIT_AT[0] * at_terms + _INIT_GC[0] * gc_terms
    dS += _INIT_AT[1] * at_terms + _INIT_GC[1] * gc_terms

    # Self-complementary correction only when caller asserts it (Biopython flag,
    # not auto-detected — matches Bio.SeqUtils.MeltingTemp.Tm_NN behaviour).
    if selfcomp:
        dH += _INIT_SYM[0]
        dS += _INIT_SYM[1]

    # SantaLucia 1998 salt correction on ΔS (Na in M)
    Na_M = Na / 1000.0
    dS += 0.368 * (n - 1) * math.log(Na_M)

    # Effective concentration (Biopython convention): k = (dnac1 - dnac2/2) * 1e-9.
    # When selfcomp is True, use k = dnac1 * 1e-9 instead.
    R = 1.987  # cal/(mol·K)
    k = dnac1 * 1e-9 if selfcomp else (dnac1 - dnac2 / 2.0) * 1e-9

    tm_K = (dH * 1000.0) / (dS + R * math.log(k))
    return tm_K - 273.15


# ---------------------------------------------------------------------------
# Minimal FASTA parser
# ---------------------------------------------------------------------------

class _Record:
    """Tiny stand-in for Bio.SeqRecord with the fields we use."""
    __slots__ = ("id", "description", "seq", "features")

    def __init__(self, header: str, sequence: str):
        self.description = header
        self.id = header.split()[0] if header else ""
        self.seq = sequence
        self.features = []


def parse_fasta(path: str | Path) -> Iterator[_Record]:
    """Yield Record(id, description, seq) per FASTA entry. Replaces SeqIO.parse(path, 'fasta')."""
    header: str | None = None
    parts: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if header is not None:
                    yield _Record(header, "".join(parts))
                header = line[1:]
                parts = []
            elif line:
                parts.append(line)
    if header is not None:
        yield _Record(header, "".join(parts))


# ---------------------------------------------------------------------------
# Minimal GenBank parser (sequence only — features unused downstream)
# ---------------------------------------------------------------------------

def parse_genbank(path: str | Path) -> Iterator[_Record]:
    """Yield a single _Record per GenBank LOCUS, extracting only the ORIGIN
    sequence and the LOCUS-line accession. Features are ignored — callers
    treat ``record.features`` as empty when not provided.
    """
    locus = ""
    in_origin = False
    parts: list[str] = []
    seen_any = False
    with open(path) as f:
        for line in f:
            stripped = line.rstrip("\n")
            if stripped.startswith("LOCUS"):
                if seen_any and parts:
                    yield _Record(locus, "".join(parts).upper())
                    parts = []
                locus = stripped.split()[1] if len(stripped.split()) > 1 else "LOCUS"
                in_origin = False
                seen_any = True
            elif stripped.startswith("ORIGIN"):
                in_origin = True
            elif stripped.startswith("//"):
                in_origin = False
            elif in_origin:
                # Lines like "        1 acgtacgt acgtacgt"; drop digits + whitespace
                parts.append("".join(c for c in stripped if c.isalpha()))
    if seen_any:
        yield _Record(locus, "".join(parts).upper())


# ---------------------------------------------------------------------------
# Minimal Seq stand-in (only the methods we touch)
# ---------------------------------------------------------------------------

class Seq:
    """Tiny Bio.Seq.Seq replacement supporting str(), len(), slicing, .translate()."""

    __slots__ = ("_s",)

    def __init__(self, s: str | "Seq"):
        self._s = str(s).upper()

    def __str__(self) -> str:
        return self._s

    def __len__(self) -> int:
        return len(self._s)

    def __getitem__(self, key) -> "Seq":
        return Seq(self._s[key])

    def __eq__(self, other) -> bool:
        return self._s == (other._s if isinstance(other, Seq) else str(other).upper())

    def __hash__(self) -> int:
        return hash(self._s)

    def translate(self) -> "Seq":
        return Seq(translate_dna(self._s))

    def reverse_complement(self) -> "Seq":
        return Seq(reverse_complement(self._s))

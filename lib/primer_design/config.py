"""Constants for primer design.

All algorithmic parameters live here. Modules import what they need; nothing else
holds mutable global state. References ``D-IDs`` point to ``decisions_log.md``.
"""
from __future__ import annotations

from typing import Literal

# ============================================================================
# Tm calculation (D5.1 — Allawi & SantaLucia 1997, Na correction SantaLucia 1998)
# ============================================================================

TM_NN_PARAMS: dict[str, float] = {
    "Na": 50.0,         # mM
    "dnac1": 500.0,     # nM, primer concentration
    "dnac2": 0.0,       # nM, complementary strand concentration (0 for symmetric)
}

# ============================================================================
# Primer body filters (D5.3) — values tuned for *P. aeruginosa* GC ~66%
# ============================================================================

GC_MIN: float = 0.40
GC_MAX: float = 0.70
GC_PREFERRED_MIN: float = 0.45      # soft warning if outside
GC_PREFERRED_MAX: float = 0.55

BODY_LEN_MIN: int = 18
BODY_LEN_MAX: int = 28

TM_HARD_MIN_C: float = 58.0
TM_HARD_MAX_C: float = 66.0
TM_PREFERRED_MIN_C: float = 60.0
TM_PREFERRED_MAX_C: float = 64.0
TM_TARGET_C: float = 62.0           # used in scoring

TM_SPREAD_HARD_LIMIT_C: float = 3.0  # rejected above this
TM_SPREAD_SOFT_WARN_C: float = 2.0   # warning above this

SELF_DIMER_MAX_3PRIME: int = 4       # contiguous 3'-end self-pair length

# 3' GC clamp
CLAMP_LAST_BASE_OK: tuple[str, ...] = ("G", "C")
CLAMP_LAST5_GC_MIN: int = 1
CLAMP_LAST5_GC_MAX: int = 3
CLAMP_NO_4IDENT_LAST4: bool = True

# ============================================================================
# Scar (deletion application, D5.4)
# ============================================================================

SCAR_MAX_TOTAL_AA: int = 24          # N + C ≤ 24
SCAR_MIN_PER_SIDE: int = 1

# ============================================================================
# Junction overlap (deletion / tagging)
# ============================================================================

JUNCTION_LEN_DEFAULT: int = 30       # standard In-Fusion 3-fragment overlap (15 nt per primer tail)
JUNCTION_LEN_PER_PRIMER_DEFAULT: int = 15

# Tagging-specific (D6.4)
JUNCTION_LEN_HIS6: int = 30          # His6 cassette is exactly 30 nt
JUNCTION_LEN_FLAG_HIS8: int = 36     # FLAG and His8 cassettes = 36 nt
TAGGING_MAX_CASSETTE_NT: int = 36    # in-locus rejected above this in v1

# ============================================================================
# Vector arms (In-Fusion homology length per side)
# ============================================================================

VECTOR_TAIL_LEN: int = 15            # all four tail rules use 15 nt

# ============================================================================
# Flanks (gene record canonical extraction)
# ============================================================================

UP_FLANK_LEN: int = 800
DN_FLANK_LEN: int = 600
MIN_FLANK_LEN: int = 300             # shrink rather than overlap with neighbour gene

# Primer body anchor offset window (within flank, 5'-most start position)
P1_OFFSET_MIN: int = 200
P1_OFFSET_MAX: int = 350
P4_OFFSET_MIN: int = 0
P4_OFFSET_MAX: int = 150

# Colony PCR buffer in upstream flank (positions 0 – COLONY_PCR_BUFFER-1 reserved for outside primer)
COLONY_PCR_BUFFER: int = 200
COLONY_PCR_OUTSIDE_SEARCH_MAX: int = 150   # search the first 150 nt of the buffer for outside primer
COLONY_PCR_INSIDE_SEARCH_MAX: int = 200    # search the first 200 nt of dn_flank for inside primer
COLONY_PCR_TM_MIN_C: float = 53.0          # Taq Tm range (standard range – 5 °C)
COLONY_PCR_TM_MAX_C: float = 61.0
COLONY_PCR_TM_TARGET_C: float = 57.0

# ============================================================================
# RBS for plasmid expression (D6.1.v2 — fixed, no UI dropdown)
# ============================================================================

CANONICAL_RBS: str = "AGGAGG"
SYNTHETIC_SPACER_DEFAULT: str = "ACTTGTTC"   # 8 nt; matches PA14 lasR/lasI working primers
RBS_TAIL_PART: str = CANONICAL_RBS + SYNTHETIC_SPACER_DEFAULT   # 14 nt

# ============================================================================
# Tags (D6.2.v2)
# ============================================================================

LINKER_GGS: str = "GGCGGCAGC"        # GGS, 9 nt
NEW_STOP_DEFAULT: str = "TAA"        # preferred new stop after C-term tag

# ============================================================================
# Polymerases (D5.2)
# ============================================================================

PolymeraseName = Literal["B7", "Phusion", "Q5", "Taq"]

POLYMERASE_OFFSETS_C: dict[str, float] = {
    "B7": 0.0,        # Biozym B7 High Fidelity (Pfu-derived) — default
    "Phusion": 3.0,
    "Q5": 3.0,
    "Taq": -5.0,
}

POLYMERASE_DEFAULT: PolymeraseName = "B7"

# ============================================================================
# Off-target scan (D5.5)
# ============================================================================

OFFTARGET_ANCHOR_LEN: int = 10       # last N nt of body must match within mm cutoff
OFFTARGET_ANCHOR_MAX_MM: int = 1
OFFTARGET_BODY_MAX_MM: int = 4
OFFTARGET_PRODUCT_SIZE_MIN_BP: int = 50
OFFTARGET_PRODUCT_SIZE_MAX_BP: int = 6000

# ============================================================================
# Gene catalog (extend by adding to this set + adding gene FASTA in R2)
# ============================================================================

SUPPORTED_GENES: tuple[str, ...] = ("lasB", "lasR", "lasI", "rhlR", "aprA")

# ============================================================================
# Vector catalog
# ============================================================================

SUPPORTED_VECTORS: tuple[str, ...] = ("pEXG2", "pBBR1MCS2")

# ============================================================================
# Restriction enzymes (subset relevant to v1 vectors)
# Format: (recognition_site_5prime3prime, top-strand nick offset within motif)
# ============================================================================

RESTRICTION_SITES: dict[str, tuple[str, int]] = {
    "HindIII": ("AAGCTT", 1),     # A^AGCTT
    "EcoRI":   ("GAATTC", 1),     # G^AATTC
    "BamHI":   ("GGATCC", 1),
    "SalI":    ("GTCGAC", 1),
    "XhoI":    ("CTCGAG", 1),
    "KpnI":    ("GGTACC", 5),     # GGTAC^C
    "SacI":    ("GAGCTC", 5),
    "SpeI":    ("ACTAGT", 1),
    "XbaI":    ("TCTAGA", 1),
    "ApaI":    ("GGGCCC", 5),
    "ClaI":    ("ATCGAT", 2),
    "EcoRV":   ("GATATC", 3),     # blunt
    "SmaI":    ("CCCGGG", 3),     # blunt
}

# ============================================================================
# R2 storage layout (set by user 2026-05-04)
# ============================================================================

R2_ENDPOINT_PATTERN: str = "https://{account_id}.r2.cloudflarestorage.com"
R2_GENOMES_PREFIX: str = "Whole genome sequences/"   # note: spaces, case-sensitive
R2_GENE_FASTA_PREFIX: str = ""                       # gene FASTAs at bucket root

# Cache for warm Vercel function instances
GENOME_LRU_CACHE_SIZE: int = 4

# ============================================================================
# Tool version
# ============================================================================

TOOL_VERSION: str = "0.1.0"

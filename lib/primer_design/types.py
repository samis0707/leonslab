"""Dataclasses used across the primer-design pipeline.

All public APIs operate on these. Frozen where possible to prevent accidental
mutation in the orchestration layer.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal


Application = Literal["deletion", "tagging", "expression"]
TagPosition = Literal["N", "C"]
Strand = Literal["+", "-"]


# ============================================================================
# Request / response
# ============================================================================

@dataclass
class DesignRequest:
    """Validated input from POST /api/design."""
    isolate_id: str
    gene: str
    action: Literal["delete", "tag", "express"]
    enzyme: str
    polymerase: str = "B7"
    tag: str | None = None
    tag_position: TagPosition | None = None
    use_plasmid_for_tag: bool = False    # if True with action="tag", switches to expression with C-term tag

    @property
    def application(self) -> Application:
        if self.action == "delete":
            return "deletion"
        if self.action == "tag" and not self.use_plasmid_for_tag:
            return "tagging"
        return "expression"

    @property
    def vector(self) -> str:
        if self.application in ("deletion", "tagging"):
            return "pEXG2"
        return "pBBR1MCS2"


# ============================================================================
# Gene record (from canonical curated FASTA)
# ============================================================================

@dataclass(frozen=True)
class GeneRecord:
    """Strand-normalized record for one (gene, isolate) pair.

    All sequences are 5'→3' on the *coding strand*, regardless of original genomic strand.
    """
    isolate_id: str
    gene: str
    cds_seq: str           # starts with ATG, ends with stop, length multiple of 3
    up_flank: str          # may be < 600 nt if neighbour-shortened
    dn_flank: str          # may be < 600 nt if neighbour-shortened
    contig_id: str
    genome_start_1based: int
    genome_end_1based: int
    original_strand: Strand
    functional: bool = True
    """False if the allele is truncated / non-functional at the protein level
    (premature stop, frameshift, missing start, etc.). The CDS sequence is then
    a *reconstructed nominal* CDS based on the closest functional reference;
    primer design still works (user may want to delete or complement the
    truncated allele) but the output must surface this warning prominently."""

    truncation_reason: str | None = None
    """Free-text reason when ``functional`` is False
    (e.g. ``"premature_stop@codon_47"``, ``"frameshift_at_nt_213"``)."""


# ============================================================================
# Vector
# ============================================================================

@dataclass(frozen=True)
class VectorRecord:
    name: str              # "pEXG2" or "pBBR1MCS2"
    sequence: str          # uppercase, linear representation of circular vector
    length: int
    features: list[Any] = field(default_factory=list)   # SeqFeature instances from Biopython


@dataclass(frozen=True)
class TailConvention:
    """Per-(vector, enzyme, application) calibrated rule set."""
    name: str                                            # e.g. "site_destroyed"
    p1_tail_rule: str                                    # key into vectors.TAIL_RULES
    p4_tail_rule: str                                    # key into vectors.TAIL_RULES
    expected_recognition_count_in_final_plasmid: int     # 0, 1, ...
    enzyme: str
    vector: str
    application: Application


# ============================================================================
# Tag
# ============================================================================

@dataclass(frozen=True)
class Tag:
    name: str
    dna: str               # P. aeruginosa codon-optimized
    protein: str           # one-letter amino-acid string, no stop
    n_term_ok: bool
    c_term_ok: bool
    in_locus_ok: bool      # False if cassette > 36 nt
    cassette_nt: int       # full cassette: linker + dna + new stop


# ============================================================================
# Primer
# ============================================================================

@dataclass(frozen=True)
class Primer:
    name: str              # P1, P2, P3, P4
    role: str              # "UP_Fwd", "UP_Rev", "DN_Fwd", "DN_Rev", "INSERT_Fwd", "INSERT_Rev"
    tail: str              # 5' non-templated portion
    body: str              # 3' template-binding portion
    tail_kind: str         # convention name, for FASTA header annotation
    tm_body_C: float
    gc_body: float         # 0.0–1.0
    length: int            # = len(tail) + len(body)

    @property
    def sequence(self) -> str:
        return self.tail + self.body

    @property
    def annotated_sequence(self) -> str:
        """Tail in lowercase, body in UPPERCASE — for FASTA output."""
        return self.tail.lower() + self.body.upper()


# ============================================================================
# Primer sets per application
# ============================================================================

@dataclass
class DeletionPrimerSet:
    p1: Primer
    p2: Primer
    p3: Primer
    p4: Primer
    N: int
    C: int
    scar_dna: str          # = scar ORF including stop, multiple of 3
    score: float
    tm_spread: float


@dataclass
class TaggingPrimerSet:
    p1: Primer
    p2: Primer
    p3: Primer
    p4: Primer
    cassette: str          # GGS + tag DNA + new stop, full cassette in coding orientation
    overlap_left: int      # length of cassette in P2 tail (before junction-overlap region)
    overlap_right: int     # length of cassette in P3 tail (after junction-overlap region)
    tag: Tag
    score: float
    tm_spread: float


@dataclass
class ExpressionPrimerSet:
    p1: Primer
    p2: Primer
    coding_seq: str        # full insert CDS with optional tag, 5'→3', starts ATG, ends stop
    # Stretch of native CDS that is amplified from genomic template. Equals
    # ``coding_seq`` for the untagged case; for tagged designs the tag-encoding
    # nucleotides are carried in the primer overhangs (P1 tail for N-term,
    # P2 tail for C-term) rather than primed off the genome.
    native_anneal_segment: str
    tag: Tag | None
    tag_position: TagPosition | None
    score: float
    tm_spread: float


@dataclass(frozen=True)
class ColonyPCRPrimerSet:
    """Pre-computed conserved colony PCR primers for pEXG2 deletion verification."""
    gene: str
    outside: Primer          # forward primer upstream of P1 buffer zone
    inside: Primer           # reverse primer in downstream flank
    expected_deletion_product_bp_min: int
    expected_deletion_product_bp_max: int
    note: str = ""           # e.g. "conserved across N/M isolates"


PrimerSet = DeletionPrimerSet | TaggingPrimerSet | ExpressionPrimerSet


# ============================================================================
# Off-target scan
# ============================================================================

@dataclass(frozen=True)
class OffTargetSite:
    primer_name: str
    contig_id: str
    position_0based: int
    strand: Strand
    body_mismatches: int


@dataclass(frozen=True)
class OffTargetProduct:
    primer_pair: tuple[str, str]      # ("P1", "P2"), etc.
    contig_id: str
    start_0based: int
    end_0based: int
    size_bp: int
    is_expected: bool                  # True for the intended PCR product


@dataclass
class OffTargetReport:
    sites_per_primer: dict[str, list[OffTargetSite]]
    products: list[OffTargetProduct]
    passed: bool
    violations: list[OffTargetProduct] = field(default_factory=list)


# ============================================================================
# Final result
# ============================================================================

@dataclass
class CircularPlasmid:
    name: str
    sequence: str          # linear representation, 5'→3' on top strand
    length: int
    junctions: list[tuple[int, int]] = field(default_factory=list)   # 0-based positions of insert junctions


@dataclass
class DesignResult:
    """End-to-end result. Application-specific fields filled in by the orchestrator."""
    request: DesignRequest
    gene_record: GeneRecord
    vector: VectorRecord
    convention: TailConvention
    primer_set: PrimerSet
    up_amplicon: str | None              # None for expression
    dn_amplicon: str | None              # None for expression
    insert: str                           # full insert (UP+DN or single amplicon)
    final_plasmid: CircularPlasmid
    off_target: OffTargetReport
    warnings: list[str] = field(default_factory=list)
    colony_pcr_primers: ColonyPCRPrimerSet | None = None

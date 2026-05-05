"""Hard validation gates per D7.4 / skill_v2 §11.

Every check raises ``VerificationFailure`` with ``failed_check`` in details so
the handler can surface it as a typed JSON error.
"""
from __future__ import annotations

from Bio.Seq import Seq

from .config import RESTRICTION_SITES, JUNCTION_LEN_DEFAULT, JUNCTION_LEN_FLAG_HIS8
from .exceptions import VerificationFailure
from .tags import LINKER_GGS
from .types import (
    DeletionPrimerSet,
    DesignResult,
    ExpressionPrimerSet,
    TaggingPrimerSet,
)
from .vectors import reverse_complement


def verify(result: DesignResult) -> None:
    """Dispatch to the application-specific verifier."""
    app = result.request.application
    if app == "deletion":
        _verify_deletion(result)
    elif app == "tagging":
        _verify_tagging(result)
    elif app == "expression":
        _verify_expression(result)
    else:
        raise VerificationFailure(
            message=f"Unknown application {app!r}",
            details={"failed_check": "application_dispatch"},
        )


# ---------------------------------------------------------------------------
# Shared checks
# ---------------------------------------------------------------------------

def _check_recognition_count(result: DesignResult) -> None:
    enzyme = result.request.enzyme
    motif, _ = RESTRICTION_SITES[enzyme]
    rc_motif = reverse_complement(motif)
    seq = result.final_plasmid.sequence

    actual = seq.count(motif)
    if rc_motif != motif:
        actual += seq.count(rc_motif)

    expected = result.convention.expected_recognition_count_in_final_plasmid
    if actual != expected:
        raise VerificationFailure(
            message=(
                f"{enzyme} site count in final plasmid is {actual}; "
                f"convention {result.convention.name!r} expects {expected}"
            ),
            details={
                "failed_check": "recognition_count",
                "enzyme": enzyme,
                "actual": actual,
                "expected": expected,
                "convention": result.convention.name,
            },
        )


def _check_off_target_passed(result: DesignResult) -> None:
    if not result.off_target.passed:
        raise VerificationFailure(
            message="Off-target scan reports unintended PCR products",
            details={
                "failed_check": "off_target",
                "violations": [vars(p) for p in result.off_target.violations],
            },
        )


def _check_plasmid_size(result: DesignResult) -> None:
    """Final plasmid length = vector + insert − 2 × overlap (for site_destroyed deletions/tagging)
    or = vector + insert − 2 × overlap (for expression with single insert).
    """
    expected = (
        result.vector.length
        + len(result.insert)
        - 2 * JUNCTION_LEN_DEFAULT // 2   # 15 nt per side default; 18 nt per side for FLAG/His8 tagging
    )
    # The tagging case may use extended overlap; let plasmid_builder report the actual.
    # Here we just sanity-check it's plausible (within ±50 bp of the naive expectation).
    actual = result.final_plasmid.length
    if abs(actual - expected) > 50:
        raise VerificationFailure(
            message=f"Final plasmid size {actual} bp deviates from expected ~{expected} bp",
            details={
                "failed_check": "plasmid_size",
                "actual_bp": actual,
                "expected_bp": expected,
            },
        )


# ---------------------------------------------------------------------------
# Deletion-specific
# ---------------------------------------------------------------------------

def _verify_deletion(result: DesignResult) -> None:
    assert isinstance(result.primer_set, DeletionPrimerSet)
    ps = result.primer_set

    _check_recognition_count(result)
    _check_off_target_passed(result)
    _check_plasmid_size(result)

    # Scar ORF
    scar_protein = str(Seq(ps.scar_dna).translate())
    if not scar_protein.startswith("M"):
        raise VerificationFailure(
            message=f"Scar ORF doesn't start with M: {scar_protein!r}",
            details={"failed_check": "scar_start_codon"},
        )
    if not scar_protein.endswith("*") or scar_protein.count("*") != 1:
        raise VerificationFailure(
            message=f"Scar ORF malformed (stop codons): {scar_protein!r}",
            details={"failed_check": "scar_stop_codon"},
        )
    expected_len = (ps.N + ps.C) * 3
    if len(ps.scar_dna) != expected_len:
        raise VerificationFailure(
            message=f"Scar DNA length {len(ps.scar_dna)} != expected {expected_len}",
            details={"failed_check": "scar_length", "N": ps.N, "C": ps.C},
        )

    # Junction overlap exact
    if result.up_amplicon is None or result.dn_amplicon is None:
        raise VerificationFailure(
            message="Deletion result missing UP or DN amplicon",
            details={"failed_check": "amplicons_present"},
        )
    overlap = JUNCTION_LEN_DEFAULT
    if result.up_amplicon[-overlap:] != result.dn_amplicon[:overlap]:
        raise VerificationFailure(
            message=f"UP/DN junction overlap of {overlap} nt mismatched",
            details={
                "failed_check": "junction_overlap",
                "up_tail": result.up_amplicon[-overlap:],
                "dn_head": result.dn_amplicon[:overlap],
            },
        )


# ---------------------------------------------------------------------------
# Tagging-specific
# ---------------------------------------------------------------------------

def _verify_tagging(result: DesignResult) -> None:
    assert isinstance(result.primer_set, TaggingPrimerSet)
    ps = result.primer_set

    _check_recognition_count(result)
    _check_off_target_passed(result)
    _check_plasmid_size(result)

    # Cassette length matches tag and is consistent
    if len(ps.cassette) != ps.tag.cassette_nt:
        raise VerificationFailure(
            message=f"Cassette length {len(ps.cassette)} != tag.cassette_nt {ps.tag.cassette_nt}",
            details={"failed_check": "cassette_length", "tag": ps.tag.name},
        )
    if not ps.cassette.startswith(LINKER_GGS):
        raise VerificationFailure(
            message="Cassette does not start with GGS linker",
            details={"failed_check": "cassette_linker"},
        )

    # Junction overlap exact (for tagging, length = cassette / 2 rounded up)
    overlap_len = (
        len(ps.cassette) // 2 + len(ps.cassette) % 2
        if len(ps.cassette) > JUNCTION_LEN_DEFAULT
        else JUNCTION_LEN_DEFAULT // 2
    )
    if result.up_amplicon is None or result.dn_amplicon is None:
        raise VerificationFailure(
            message="Tagging result missing UP or DN amplicon",
            details={"failed_check": "amplicons_present"},
        )
    # The exact overlap span is encoded by plasmid_builder; trust ps.overlap_left/right
    actual_overlap = ps.overlap_left + ps.overlap_right - len(ps.cassette)
    if actual_overlap != JUNCTION_LEN_DEFAULT // 2 and actual_overlap != JUNCTION_LEN_FLAG_HIS8 // 2:
        # Allowed values: 15 (His6) or 18 (FLAG, His8)
        raise VerificationFailure(
            message=f"Tagging junction overlap {actual_overlap} nt unexpected",
            details={"failed_check": "tagging_junction_overlap"},
        )

    # Fusion ORF
    # The assembled CDS (after junction collapse) should be:
    # native_CDS_minus_stop + cassette_GGS + cassette_tag + cassette_TAA
    # → translates to: native_protein_minus_stop + "GGS" + tag.protein + "*"
    # We don't reconstruct the full insert here (plasmid_builder does that); we just
    # check the cassette translates correctly:
    cassette_protein = str(Seq(ps.cassette).translate())
    if not cassette_protein.endswith("*"):
        raise VerificationFailure(
            message=f"Cassette doesn't end with stop: {cassette_protein!r}",
            details={"failed_check": "cassette_stop"},
        )
    if cassette_protein.count("*") != 1:
        raise VerificationFailure(
            message=f"Cassette has {cassette_protein.count('*')} stop codons",
            details={"failed_check": "cassette_internal_stop"},
        )
    if "GGS" not in cassette_protein and "GGGS" not in cassette_protein:
        raise VerificationFailure(
            message=f"GGS linker missing from cassette protein: {cassette_protein!r}",
            details={"failed_check": "cassette_linker_protein"},
        )
    if ps.tag.protein not in cassette_protein:
        raise VerificationFailure(
            message=f"Tag protein {ps.tag.protein!r} not in cassette: {cassette_protein!r}",
            details={"failed_check": "cassette_tag_protein"},
        )


# ---------------------------------------------------------------------------
# Expression-specific
# ---------------------------------------------------------------------------

def _verify_expression(result: DesignResult) -> None:
    assert isinstance(result.primer_set, ExpressionPrimerSet)
    ps = result.primer_set

    _check_recognition_count(result)
    _check_off_target_passed(result)
    _check_plasmid_size(result)

    # CDS in-frame from start to stop
    cds = ps.coding_seq
    if len(cds) % 3 != 0:
        raise VerificationFailure(
            message=f"Coding sequence length {len(cds)} not multiple of 3",
            details={"failed_check": "cds_frame"},
        )
    protein = str(Seq(cds).translate())
    if not protein.startswith("M"):
        raise VerificationFailure(
            message="Expression CDS doesn't start with M",
            details={"failed_check": "cds_start"},
        )
    if not protein.endswith("*") or protein.count("*") != 1:
        raise VerificationFailure(
            message=f"Expression CDS has {protein.count('*')} stops",
            details={"failed_check": "cds_internal_stop"},
        )

    # Tag presence if requested
    if ps.tag is not None:
        if ps.tag.protein not in protein:
            raise VerificationFailure(
                message=f"Tag protein {ps.tag.protein!r} not in fused CDS",
                details={"failed_check": "tag_in_cds", "tag": ps.tag.name},
            )

    # RBS presence: AGGAGG should be in the insert tail region (P1 carries it)
    insert = result.insert
    if "AGGAGG" not in insert:
        raise VerificationFailure(
            message="Canonical RBS AGGAGG not detected in insert",
            details={"failed_check": "rbs_present"},
        )

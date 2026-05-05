"""Typed exceptions for the primer-design pipeline.

Each exception carries an ``error_code`` (snake_case string) used by ``handler.py``
to translate to the JSON error response. ``message`` is human-readable; ``details``
is structured data for debugging.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class PrimerDesignError(Exception):
    """Base for all tool errors."""
    error_code: str = "internal_error"
    message: str = ""
    details: dict = field(default_factory=dict)

    def __str__(self) -> str:
        return f"[{self.error_code}] {self.message}"


# ------------------------------------------------------------------
# Data layer
# ------------------------------------------------------------------

@dataclass
class GeneRecordNotFound(PrimerDesignError):
    error_code: str = "gene_record_not_found"


@dataclass
class GenomeNotFoundInR2(PrimerDesignError):
    error_code: str = "genome_not_found_in_r2"


@dataclass
class GenomeManifestMismatch(PrimerDesignError):
    """SHA-256 of fetched genome doesn't match manifest."""
    error_code: str = "genome_manifest_mismatch"


@dataclass
class CdsValidationError(PrimerDesignError):
    """CDS in user-provided FASTA fails one of: length-mod-3, starts-with-ATG, single-stop, ends-with-stop."""
    error_code: str = "cds_invalid"


@dataclass
class GeneNotInGenome(PrimerDesignError):
    """Build-time: gene CDS not findable in the isolate's genome (deletion in this strain, or assembly gap)."""
    error_code: str = "gene_not_in_genome"


@dataclass
class GeneAmbiguousInGenome(PrimerDesignError):
    """Build-time: gene CDS matches multiple distinct loci (paralog or repeat)."""
    error_code: str = "gene_ambiguous_in_genome"


# ------------------------------------------------------------------
# Vector / cut-site layer
# ------------------------------------------------------------------

@dataclass
class NotSingleCutter(PrimerDesignError):
    """Chosen restriction enzyme cuts the vector zero or multiple times."""
    error_code: str = "not_single_cutter"


@dataclass
class ConventionNotCalibrated(PrimerDesignError):
    """No tail convention defined for the (vector, enzyme, application) triple."""
    error_code: str = "convention_not_calibrated"


@dataclass
class UnknownEnzyme(PrimerDesignError):
    error_code: str = "unknown_enzyme"


@dataclass
class UnknownVector(PrimerDesignError):
    error_code: str = "unknown_vector"


# ------------------------------------------------------------------
# Tag layer
# ------------------------------------------------------------------

@dataclass
class TagTooLongForInLocus(PrimerDesignError):
    """Tag cassette exceeds 36 nt; cannot be encoded in 4-primer in-locus design."""
    error_code: str = "tag_too_long_for_in_locus"


@dataclass
class InvalidTagPosition(PrimerDesignError):
    """Tag does not support the requested N or C position (e.g., 3xFLAG-N)."""
    error_code: str = "invalid_tag_position"


@dataclass
class UnknownTag(PrimerDesignError):
    error_code: str = "unknown_tag"


# ------------------------------------------------------------------
# Algorithm layer
# ------------------------------------------------------------------

@dataclass
class NoCandidates(PrimerDesignError):
    """Primer pool empty after applying hard filters; relax_filters_suggestion in details."""
    error_code: str = "no_candidates"


@dataclass
class TmSpreadTooLarge(PrimerDesignError):
    """No (P1, P2, P3, P4) tuple has Tm spread ≤ 3 °C."""
    error_code: str = "tm_spread_too_large"


@dataclass
class OffTargetDetected(PrimerDesignError):
    """Off-target scan found unintended product. ``details["products"]`` lists them."""
    error_code: str = "off_target_detected"


# ------------------------------------------------------------------
# Verification layer (D7.4)
# ------------------------------------------------------------------

@dataclass
class VerificationFailure(PrimerDesignError):
    """One of the hard verification checks failed. ``details["failed_check"]`` names which."""
    error_code: str = "verification_failure"


# ------------------------------------------------------------------
# Request validation
# ------------------------------------------------------------------

@dataclass
class InvalidRequest(PrimerDesignError):
    """Request schema violation."""
    error_code: str = "invalid_request"


# ------------------------------------------------------------------
# Build-pipeline errors
# ------------------------------------------------------------------

@dataclass
class BuildError(PrimerDesignError):
    """Generic build-pipeline error (data ingestion)."""
    error_code: str = "build_error"

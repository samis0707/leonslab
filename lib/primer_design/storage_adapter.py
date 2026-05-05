"""Cloudflare R2 storage adapter.

Lazy-fetches isolate genomes (~7 MB each) from the private R2 bucket configured
via Vercel env vars. Caches by isolate ID on the warm function instance for
~minutes-scale lifetime; cold starts re-fetch.

Bucket layout (verified 2026-05-04):
    leonslab-pa-genomes/
    ├── Whole genome sequences/      ← .fna files (case-sensitive, contains spaces)
    │   ├── LB001.fna
    │   └── ...
    ├── lasB.fasta                   ← gene CDS lists at bucket root
    ├── lasR.fasta
    └── lasI.fasta                   (added later)
"""
from __future__ import annotations

import functools
import hashlib
import json
import logging
import os
from pathlib import Path

import boto3
from botocore.client import Config

from .config import (
    GENOME_LRU_CACHE_SIZE,
    R2_ENDPOINT_PATTERN,
    R2_GENE_FASTA_PREFIX,
    R2_GENOMES_PREFIX,
)
from .exceptions import GenomeManifestMismatch, GenomeNotFoundInR2

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# R2 client
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def _r2_client():
    """Return a cached boto3 S3-compatible client pointed at the R2 endpoint."""
    account_id = os.environ["R2_ACCOUNT_ID"]
    access_key = os.environ["R2_ACCESS_KEY_ID"]
    secret_key = os.environ["R2_SECRET_ACCESS_KEY"]
    endpoint = R2_ENDPOINT_PATTERN.format(account_id=account_id)

    return boto3.client(
        "s3",
        endpoint_url=endpoint,
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        config=Config(signature_version="s3v4", retries={"max_attempts": 3, "mode": "standard"}),
    )


def _bucket_name() -> str:
    return os.environ["R2_BUCKET_NAME"]


# ---------------------------------------------------------------------------
# Genome fetch (lazy + cached + manifest-verified)
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=GENOME_LRU_CACHE_SIZE)
def fetch_genome(isolate_id: str, manifest_path: Path | None = None) -> bytes:
    """Fetch a genome FASTA from R2 and verify its SHA-256 against the manifest.

    Args:
        isolate_id: e.g. ``"LB001"``. Maps to key ``Whole genome sequences/LB001.fna``.
        manifest_path: optional override; defaults to bundled manifest.

    Returns:
        Raw bytes of the .fna file.

    Raises:
        GenomeNotFoundInR2: key missing from bucket.
        GenomeManifestMismatch: SHA-256 disagrees with manifest entry.
    """
    key = f"{R2_GENOMES_PREFIX}{isolate_id}.fna"
    log.info("Fetching genome from R2: bucket=%s key=%s", _bucket_name(), key)

    try:
        resp = _r2_client().get_object(Bucket=_bucket_name(), Key=key)
    except _r2_client().exceptions.NoSuchKey as exc:
        raise GenomeNotFoundInR2(
            message=f"No genome '{isolate_id}' in R2 (key: {key})",
            details={"isolate_id": isolate_id, "key": key},
        ) from exc

    data = resp["Body"].read()
    expected_sha = _expected_sha256(isolate_id, manifest_path)
    if expected_sha is not None:
        actual_sha = hashlib.sha256(data).hexdigest()
        if actual_sha != expected_sha:
            raise GenomeManifestMismatch(
                message=f"Genome {isolate_id} SHA-256 mismatch",
                details={
                    "isolate_id": isolate_id,
                    "expected": expected_sha,
                    "actual": actual_sha,
                },
            )
    return data


def fetch_gene_fasta(gene: str) -> bytes:
    """Fetch a gene CDS list FASTA (e.g. ``lasB.fasta``) from R2 bucket root.

    Used by ``scripts/build_gene_records.py`` during build pipeline.
    """
    key = f"{R2_GENE_FASTA_PREFIX}{gene}.fasta"
    log.info("Fetching gene FASTA from R2: bucket=%s key=%s", _bucket_name(), key)
    try:
        resp = _r2_client().get_object(Bucket=_bucket_name(), Key=key)
    except _r2_client().exceptions.NoSuchKey as exc:
        raise GenomeNotFoundInR2(
            message=f"No gene FASTA '{gene}' in R2 (key: {key})",
            details={"gene": gene, "key": key},
        ) from exc
    return resp["Body"].read()


# ---------------------------------------------------------------------------
# Manifest helpers
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def _load_manifest(manifest_path: Path | None = None) -> dict:
    path = manifest_path or _default_manifest_path()
    if not path.exists():
        log.warning("Manifest not found at %s; SHA-256 verification disabled", path)
        return {"genomes": []}
    with open(path) as f:
        return json.load(f)


def _expected_sha256(isolate_id: str, manifest_path: Path | None = None) -> str | None:
    """Return expected SHA-256 from manifest, or None if no manifest entry."""
    manifest = _load_manifest(manifest_path)
    for entry in manifest.get("genomes", []):
        if entry["isolate_id"] == isolate_id:
            return entry["sha256"]
    return None


def _default_manifest_path() -> Path:
    return Path(__file__).resolve().parents[2] / "data" / "primer_design" / "manifest.json"


# ---------------------------------------------------------------------------
# Upload (used by scripts/upload_genomes.py only — not by request handler)
# ---------------------------------------------------------------------------

def upload_genome(isolate_id: str, file_path: Path, overwrite: bool = False) -> str:
    """Upload a genome to R2. Returns its SHA-256 hex digest.

    Used by the one-off ``upload_genomes.py`` script and by Phase-2 isolate-batch
    additions. Not used by the runtime handler.
    """
    key = f"{R2_GENOMES_PREFIX}{isolate_id}.fna"
    if not overwrite and _key_exists(key):
        log.info("Skipping %s — already present in R2 and overwrite=False", key)
        with open(file_path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()

    with open(file_path, "rb") as f:
        data = f.read()
    sha = hashlib.sha256(data).hexdigest()
    _r2_client().put_object(Bucket=_bucket_name(), Key=key, Body=data, ContentType="text/x-fasta")
    log.info("Uploaded %s (%d bytes, sha256=%s)", key, len(data), sha)
    return sha


def _key_exists(key: str) -> bool:
    try:
        _r2_client().head_object(Bucket=_bucket_name(), Key=key)
        return True
    except _r2_client().exceptions.ClientError:
        return False

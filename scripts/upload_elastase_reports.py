#!/usr/bin/env python3
"""One-time R2 upload of pre-generated elastase assay PDF reports.

The elastase assay calculator has a single integer input (number of samples),
so instead of running Python in a Vercel serverless function per request (which
hit an unrelated Vercel Python-builder packaging bug), we pre-render every
possible report locally and upload them to R2. The runtime path is then a tiny
Node.js function that streams the right object back — no Python at request time.

Usage:
    python scripts/upload_elastase_reports.py --min 1 --max 100
    python scripts/upload_elastase_reports.py --only 5 10 15 --overwrite
"""
from __future__ import annotations

import argparse
import functools
import logging
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from elastase_assay.pdf_writer import generate_assay_report  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("upload_elastase_reports")

R2_ENDPOINT_PATTERN = "https://{account_id}.r2.cloudflarestorage.com"
R2_PREFIX = "elastase-assay/reports/"


@functools.lru_cache(maxsize=1)
def _r2_client():
    import boto3
    from botocore.client import Config

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


def _key_exists(key: str) -> bool:
    try:
        _r2_client().head_object(Bucket=_bucket_name(), Key=key)
        return True
    except _r2_client().exceptions.ClientError:
        return False


def upload_report(num_samples: int, overwrite: bool = False) -> None:
    key = f"{R2_PREFIX}{num_samples}.pdf"
    if not overwrite and _key_exists(key):
        log.info("Skipping %s — already present in R2", key)
        return

    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as f:
        pdf_path = Path(f.name)
    try:
        generate_assay_report(num_samples, str(pdf_path))
        data = pdf_path.read_bytes()
    finally:
        pdf_path.unlink(missing_ok=True)

    _r2_client().put_object(Bucket=_bucket_name(), Key=key, Body=data, ContentType="application/pdf")
    log.info("Uploaded %s (%d bytes)", key, len(data))


def main() -> int:
    parser = argparse.ArgumentParser(description="Pre-generate and upload elastase assay PDF reports to R2")
    parser.add_argument("--min", type=int, default=1)
    parser.add_argument("--max", type=int, default=100)
    parser.add_argument("--only", type=int, nargs="*", help="Restrict to these sample counts")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    values = args.only if args.only else range(args.min, args.max + 1)
    n_ok = 0
    n_total = 0
    for n in values:
        n_total += 1
        try:
            upload_report(n, overwrite=args.overwrite)
            n_ok += 1
        except Exception:
            log.exception("Failed to upload report for num_samples=%d", n)

    log.info("Done: %d/%d uploaded", n_ok, n_total)
    return 0 if n_ok == n_total else 1


if __name__ == "__main__":
    raise SystemExit(main())

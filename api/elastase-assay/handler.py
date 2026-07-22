"""Vercel Python serverless function for POST /api/elastase-assay.

Reads {"num_samples": int}, computes the assay volumes, and returns the
generated report directly as a downloadable PDF.
"""
from __future__ import annotations

import json
import logging
import sys
import tempfile
from http.server import BaseHTTPRequestHandler
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from elastase_assay.pdf_writer import generate_assay_report  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
log = logging.getLogger("elastase_assay.handler")

MAX_SAMPLES = 1000


class handler(BaseHTTPRequestHandler):
    """Vercel Python runtime expects a class named ``handler`` exported from this module."""

    def do_POST(self):  # noqa: N802
        try:
            length = int(self.headers.get("content-length", 0))
            raw = self.rfile.read(length)
            body = json.loads(raw) if raw else {}
            num_samples = _validate_num_samples(body.get("num_samples"))

            with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as f:
                pdf_path = Path(f.name)
            try:
                generate_assay_report(num_samples, str(pdf_path))
                pdf_bytes = pdf_path.read_bytes()
            finally:
                pdf_path.unlink(missing_ok=True)

            self.send_response(200)
            self.send_header("Content-Type", "application/pdf")
            self.send_header("Content-Length", str(len(pdf_bytes)))
            self.send_header(
                "Content-Disposition",
                f'attachment; filename="Elastase_Assay_Report_{num_samples}_samples.pdf"',
            )
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            self.wfile.write(pdf_bytes)
        except ValueError as e:
            log.warning("Invalid request: %s", e)
            self._error(400, str(e))
        except Exception as e:  # noqa: BLE001
            log.exception("Unhandled error")
            self._error(500, str(e))

    def do_OPTIONS(self):  # noqa: N802 (CORS preflight)
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def _error(self, status: int, message: str) -> None:
        encoded = json.dumps({"status": "error", "message": message}, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(encoded)


def _validate_num_samples(value) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError("num_samples must be an integer")
    if value <= 0:
        raise ValueError("num_samples must be a positive integer")
    if value > MAX_SAMPLES:
        raise ValueError(f"num_samples must be at most {MAX_SAMPLES}")
    return value

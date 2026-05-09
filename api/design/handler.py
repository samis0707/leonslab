"""Vercel Python serverless function for POST /api/design.

Reads JSON body, validates against ``DesignRequest`` schema, dispatches to the
right application orchestrator, returns artefacts inline (base64 for binary,
JSON inline for the manifest).

Vercel auto-detects this file as a serverless function. The exported callable
is a WSGI-compatible handler.
"""
from __future__ import annotations

import base64
import io
import json
import logging
import sys
import tempfile
from http.server import BaseHTTPRequestHandler
from pathlib import Path

# Make /lib importable without installing as a package
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from primer_design.applications import dispatch     # noqa: E402
from primer_design.config import TOOL_VERSION       # noqa: E402
from primer_design.exceptions import (              # noqa: E402
    InvalidRequest,
    PrimerDesignError,
)
from primer_design.fasta_writer import write_fasta  # noqa: E402
from primer_design.genbank_writer import write_genbank  # noqa: E402
from primer_design.json_writer import to_json_string, write_manifest  # noqa: E402
from primer_design.pdf_writer import render_pdf     # noqa: E402
from primer_design.tags import TAGS as _TAG_LIBRARY  # noqa: E402
from primer_design.types import DesignRequest      # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
log = logging.getLogger("primer_design.handler")


# ---------------------------------------------------------------------------
# Vercel handler
# ---------------------------------------------------------------------------

class handler(BaseHTTPRequestHandler):
    """Vercel Python runtime expects a class named ``handler`` exported from this module.

    Each request creates an instance; ``do_POST`` handles the design request.
    """

    def do_POST(self):  # noqa: N802 (BaseHTTPRequestHandler convention)
        try:
            length = int(self.headers.get("content-length", 0))
            raw = self.rfile.read(length)
            body = json.loads(raw)
            request = _validate_request(body)
            response = run_design(request)
            self._respond(200, response)
        except PrimerDesignError as e:
            log.warning("PrimerDesignError: %s", e)
            self._respond(
                400,
                {
                    "status": "error",
                    "error_code": e.error_code,
                    "message": e.message,
                    "details": e.details,
                },
            )
        except Exception as e:                      # noqa: BLE001
            log.exception("Unhandled error")
            self._respond(
                500,
                {
                    "status": "error",
                    "error_code": "internal_error",
                    "message": str(e),
                },
            )

    def do_OPTIONS(self):  # noqa: N802 (CORS preflight)
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def _respond(self, status: int, body: dict) -> None:
        encoded = json.dumps(body, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(encoded)


# ---------------------------------------------------------------------------
# Pipeline (also callable from CLI / tests)
# ---------------------------------------------------------------------------

def run_design(request: DesignRequest) -> dict:
    """Run the full design pipeline; return the API response dict.

    Used by both the HTTP handler and unit/integration tests.
    """
    log.info("Designing: %s", request)
    result = dispatch(request)

    # Render artefacts
    fasta_str = write_fasta(result)
    genbank_str = write_genbank(result)
    manifest = write_manifest(result)

    # PDF goes to /tmp, then base64-encoded inline
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as f:
        pdf_path = Path(f.name)
    try:
        render_pdf(result, pdf_path)
        pdf_bytes = pdf_path.read_bytes()
    finally:
        pdf_path.unlink(missing_ok=True)

    return {
        "status": "ok",
        "tool_version": TOOL_VERSION,
        "job_id": manifest["job_id"],
        "artefacts": {
            "fasta": base64.b64encode(fasta_str.encode("utf-8")).decode("ascii"),
            "genbank": base64.b64encode(genbank_str.encode("utf-8")).decode("ascii"),
            "pdf": base64.b64encode(pdf_bytes).decode("ascii"),
            "json": manifest,
        },
        "summary": {
            "primers": [_short_primer(p) for p in manifest["computed"]["primers"]],
            "final_plasmid_length_bp": manifest["computed"]["final_plasmid_length_bp"],
            "annealing_temp_C": manifest["computed"]["annealing_temp_C_recommended"],
            "scar_or_cassette_translation": _summary_translation(manifest),
        },
        "warnings": result.warnings,
    }


# ---------------------------------------------------------------------------
# Request validation
# ---------------------------------------------------------------------------

_VALID_ACTIONS = ("delete", "tag", "express")
_VALID_TAGS = (None, *sorted(_TAG_LIBRARY.keys()))
_VALID_TAG_POSITIONS = (None, "N", "C")
_VALID_POLYMERASES = ("B7", "Phusion", "Q5", "Taq")


def _validate_request(body: dict) -> DesignRequest:
    """Validate and coerce the JSON body into a DesignRequest dataclass.

    Raises:
        InvalidRequest: with details about the specific violation.
    """
    required = ("isolate_id", "gene", "action", "enzyme")
    missing = [k for k in required if k not in body]
    if missing:
        raise InvalidRequest(message=f"Missing required fields: {missing}")

    action = body["action"]
    if action not in _VALID_ACTIONS:
        raise InvalidRequest(message=f"action must be one of {_VALID_ACTIONS}")

    tag = body.get("tag")
    if tag not in _VALID_TAGS:
        raise InvalidRequest(message=f"tag must be one of {_VALID_TAGS}")

    tag_position = body.get("tag_position")
    if tag_position not in _VALID_TAG_POSITIONS:
        raise InvalidRequest(message=f"tag_position must be one of {_VALID_TAG_POSITIONS}")

    polymerase = body.get("polymerase", "B7")
    if polymerase not in _VALID_POLYMERASES:
        raise InvalidRequest(message=f"polymerase must be one of {_VALID_POLYMERASES}")

    if action == "tag" and tag is None:
        raise InvalidRequest(message="action='tag' requires a tag")
    if action == "tag" and tag_position is None:
        # default to C for in-locus
        tag_position = "C"

    return DesignRequest(
        isolate_id=body["isolate_id"],
        gene=body["gene"],
        action=action,
        enzyme=body["enzyme"],
        polymerase=polymerase,
        tag=tag,
        tag_position=tag_position,
        use_plasmid_for_tag=bool(body.get("use_plasmid_for_tag", False)),
    )


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def _short_primer(p: dict) -> dict:
    return {
        "name": p["name"],
        "sequence": p["sequence"],
        "tail_len": p["tail_len"],
        "body_len": p["body_len"],
        "tm_body_C": p["tm_body_C"],
        "gc_body_pct": p["gc_body_pct"],
    }


def _summary_translation(manifest: dict) -> str:
    if "scar" in manifest["computed"]:
        return manifest["computed"]["scar"]["translation"]
    if "tag_cassette" in manifest["computed"]:
        return manifest["computed"]["tag_cassette"]["translation"]
    return ""

"""PDF report generator (3 templates: deletion, tagging, expression).

Implementation note for Phase 2:
    - Use ReportLab Platypus (Story / Paragraph / Table)
    - Common header and footer; section content varies per application
    - Match the visual style of v1 ``lasB_deletion_LB001_cloning_overview.pdf``:
        * Two-column header (Parameter | Value) for "Konstrukt-Überblick"
        * Bold tail in primer rows (use partial colored Run if ReportLab supports;
          otherwise lowercase tail / uppercase body convention)
        * Tables with borders, 9-10 pt body text
    - Use DejaVu fonts for °C, µ, Δ characters (register via pdfmetrics)

Section structure (from skill_v2 §12.2):

| Section | Deletion | Tagging | Expression |
|---|---|---|---|
| 1. Konstrukt-Überblick | Y (with scar) | Y (with cassette) | Y (with RBS+tag) |
| 2. Primer-Tabelle | 4 primers | 4 primers | 2 primers |
| 3. Off-Target-Check | 4 pairs | 4 pairs | 1 pair |
| 4. PCR-Bedingungen | per polymerase | per polymerase | per polymerase |
| 5. Vektor-Vorbereitung + In-Fusion | 3-fragment | 3-fragment, extended | 2-fragment |
| 6. Cloning workflow | conjugation + sucrose CS + Gm-loss screen | same | transformation + Km only |
| 7. Validation suggestions | Sanger junctions + colony PCR | + Western (anti-tag) | + Western + activity |
| 8. Ordering checklist | yes | yes | yes |
| 9. Hinweise (soft warnings) | yes | yes | yes |
"""
from __future__ import annotations

from pathlib import Path

from .types import DesignResult


def render_pdf(result: DesignResult, output_path: Path) -> None:
    """Dispatch to the application-specific renderer."""
    app = result.request.application
    if app == "deletion":
        render_deletion_pdf(result, output_path)
    elif app == "tagging":
        render_tagging_pdf(result, output_path)
    elif app == "expression":
        render_expression_pdf(result, output_path)
    else:
        raise ValueError(f"Unknown application {app!r}")


def render_deletion_pdf(result: DesignResult, output_path: Path) -> None:
    """Section 1–9 per skill_v2 §12.2 deletion column."""
    raise NotImplementedError("Phase 2 step 11: implement ReportLab template for deletion.")


def render_tagging_pdf(result: DesignResult, output_path: Path) -> None:
    """Section 1–9 per skill_v2 §12.2 tagging column."""
    raise NotImplementedError("Phase 2 step 11: implement ReportLab template for tagging.")


def render_expression_pdf(result: DesignResult, output_path: Path) -> None:
    """Section 1–9 per skill_v2 §12.2 expression column."""
    raise NotImplementedError("Phase 2 step 11: implement ReportLab template for expression.")


# ---------------------------------------------------------------------------
# Shared section helpers (TO IMPLEMENT)
# ---------------------------------------------------------------------------

def _render_header(story: list, result: DesignResult) -> None:
    """Title block + two-column construct overview table."""
    raise NotImplementedError


def _render_primer_table(story: list, result: DesignResult) -> None:
    """Table with columns: Primer name, Sequence (tail bold/lowercase + body),
    Length, Tm_body, GC_body."""
    raise NotImplementedError


def _render_off_target_table(story: list, result: DesignResult) -> None:
    """Table per primer-pair × outcome with PASS/FAIL coloring."""
    raise NotImplementedError


def _render_pcr_conditions(story: list, result: DesignResult) -> None:
    """Reaction setup + cycling table per polymerase. Annealing temp from
    json_writer.write_manifest()['computed']['annealing_temp_C_recommended']."""
    raise NotImplementedError


def _render_warnings(story: list, result: DesignResult) -> None:
    """Soft warnings list (D7.5). Empty if none."""
    raise NotImplementedError

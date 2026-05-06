"""PDF report generator (3 templates: deletion, tagging, expression).

A4 single-page layout. Visual style follows ``vorlage für pdf.pdf``:
    - Navy section heads with monospace primer table
    - Tail bold + body regular in the Sequenz column
    - Two-column "Parameter | Wert" overview
Compact: dropped ordering checklist, validation suggestions, soft warnings.

Polymerase data assumes Biozym B7 High Fidelity (only polymerase reported per
user request 2026-05-06): Ta = Tm_body, elongation 60 s/kb.
"""
from __future__ import annotations

from datetime import date
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)

from .types import DeletionPrimerSet, DesignResult, ExpressionPrimerSet, TaggingPrimerSet
from .verification import count_recognition_sites


# ---------------------------------------------------------------------------
# Fonts (DejaVu — has Δ, °, µ, →)
# ---------------------------------------------------------------------------

_FONTS_REGISTERED = False


def _register_fonts() -> None:
    global _FONTS_REGISTERED
    if _FONTS_REGISTERED:
        return
    base = Path("/usr/share/fonts/truetype/dejavu")
    pdfmetrics.registerFont(TTFont("DejaVu", str(base / "DejaVuSans.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVu-Bold", str(base / "DejaVuSans-Bold.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVuMono", str(base / "DejaVuSansMono.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVuMono-Bold", str(base / "DejaVuSansMono-Bold.ttf")))
    _FONTS_REGISTERED = True


# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------

NAVY = colors.HexColor("#1F3A5F")
ACCENT = colors.HexColor("#3A6EA5")
SOFT_GREY = colors.HexColor("#ECEFF2")
DARK_TEXT = colors.HexColor("#1A1A1A")
MUTED = colors.HexColor("#5A6470")


# ---------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------

def project_name(result: DesignResult) -> str:
    """Return the construct identifier per user spec.

    deletion:   ``{vector}-Δ{gene}_{isolate}``      e.g. pEXG2-ΔlasB_LB001
    tagging:    ``{vector}-{gene}-{tag}-C_{isolate}`` e.g. pEXG2-lasR-His6-C_PA14
    expression: ``{vector}-{gene}_{isolate}`` (with -{tag}-{N|C} suffix if tagged)
    """
    g = result.gene_record.gene
    iso = result.gene_record.isolate_id
    vec = result.vector.name
    app = result.request.application
    if app == "deletion":
        return f"{vec}-Δ{g}_{iso}"
    if app == "tagging":
        ps = result.primer_set
        tag_name = ps.tag.name if isinstance(ps, TaggingPrimerSet) else "tag"
        return f"{vec}-{g}-{tag_name}-C_{iso}"
    if isinstance(result.primer_set, ExpressionPrimerSet) and result.primer_set.tag:
        return (
            f"{vec}-{g}-{result.primer_set.tag.name}"
            f"-{result.primer_set.tag_position}_{iso}"
        )
    return f"{vec}-{g}_{iso}"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

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
    _render(result, output_path)


def render_tagging_pdf(result: DesignResult, output_path: Path) -> None:
    _render(result, output_path)


def render_expression_pdf(result: DesignResult, output_path: Path) -> None:
    _render(result, output_path)


# ---------------------------------------------------------------------------
# Styles
# ---------------------------------------------------------------------------

def _make_styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()["Normal"]
    return {
        "title": ParagraphStyle(
            "Title", parent=base, fontName="DejaVu-Bold",
            fontSize=15, leading=18, textColor=NAVY, spaceAfter=2,
        ),
        "subtitle": ParagraphStyle(
            "Sub", parent=base, fontName="DejaVu",
            fontSize=8.5, leading=11, textColor=MUTED, spaceAfter=8,
        ),
        "section": ParagraphStyle(
            "Section", parent=base, fontName="DejaVu-Bold",
            fontSize=10.5, leading=13, textColor=NAVY,
            spaceBefore=8, spaceAfter=3,
        ),
        "body": ParagraphStyle(
            "Body", parent=base, fontName="DejaVu",
            fontSize=8.5, leading=11, textColor=DARK_TEXT,
        ),
        "small": ParagraphStyle(
            "Small", parent=base, fontName="DejaVu",
            fontSize=7.5, leading=9.5, textColor=MUTED,
        ),
        "mono": ParagraphStyle(
            "Mono", parent=base, fontName="DejaVuMono",
            fontSize=7.8, leading=10, textColor=DARK_TEXT,
        ),
        "kv_key": ParagraphStyle(
            "KvKey", parent=base, fontName="DejaVu-Bold",
            fontSize=8.5, leading=10.5, textColor=DARK_TEXT,
        ),
        "kv_val": ParagraphStyle(
            "KvVal", parent=base, fontName="DejaVu",
            fontSize=8.5, leading=10.5, textColor=DARK_TEXT,
        ),
        "th": ParagraphStyle(
            "Th", parent=base, fontName="DejaVu-Bold",
            fontSize=8.5, leading=10.5, textColor=colors.whitesmoke,
        ),
    }


# ---------------------------------------------------------------------------
# Renderer
# ---------------------------------------------------------------------------

def _render(result: DesignResult, output_path: Path) -> None:
    _register_fonts()
    styles = _make_styles()

    doc = BaseDocTemplate(
        str(output_path), pagesize=A4,
        leftMargin=14 * mm, rightMargin=14 * mm,
        topMargin=18 * mm, bottomMargin=14 * mm,
        title=project_name(result), author="leonslab primer-design",
    )
    frame = Frame(
        doc.leftMargin, doc.bottomMargin,
        doc.width, doc.height, id="main",
        leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0,
    )

    pname = project_name(result)
    today = date.today().isoformat()

    def on_page(canvas, _doc):
        canvas.saveState()
        canvas.setFont("DejaVu", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(14 * mm, A4[1] - 10 * mm, pname)
        canvas.drawRightString(A4[0] - 14 * mm, A4[1] - 10 * mm, today)
        canvas.setStrokeColor(ACCENT)
        canvas.setLineWidth(0.6)
        canvas.line(14 * mm, A4[1] - 12 * mm, A4[0] - 14 * mm, A4[1] - 12 * mm)
        canvas.setFont("DejaVu", 7)
        canvas.setFillColor(MUTED)
        canvas.drawString(14 * mm, 8 * mm, pname)
        canvas.drawRightString(
            A4[0] - 14 * mm, 8 * mm,
            "leonslab primer-design • B7 High Fidelity",
        )
        canvas.restoreState()

    doc.addPageTemplates([PageTemplate(id="A4", frames=[frame], onPage=on_page)])

    story = []
    _section_title_block(story, result, styles)
    _section_konstrukt_ueberblick(story, result, styles)
    _section_primer_tabelle(story, result, styles)
    _section_pcr_conditions(story, result, styles)
    _section_construct_summary(story, result, styles)

    doc.build(story)


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

def _section_title_block(story, result, styles):
    g = result.gene_record.gene
    iso = result.gene_record.isolate_id
    vec = result.vector.name
    app = result.request.application
    if app == "deletion":
        title = (
            f"In-frame deletion of <i>{g}</i> in {iso} via {vec} allelic exchange"
        )
    elif app == "tagging":
        ps = result.primer_set
        tag_name = ps.tag.name if isinstance(ps, TaggingPrimerSet) else "tag"
        title = (
            f"In-locus C-terminal {tag_name} tagging of <i>{g}</i> in {iso} via {vec}"
        )
    else:
        title = f"Plasmid expression of <i>{g}</i> from {iso} in {vec}"
    story.append(Paragraph(title, styles["title"]))
    story.append(Paragraph(
        f"Primer design • PCR / Assembly protocol • {date.today().isoformat()}",
        styles["subtitle"],
    ))


def _section_konstrukt_ueberblick(story, result, styles):
    story.append(Paragraph("1   Konstrukt-Überblick", styles["section"]))

    g = result.gene_record
    vec = result.vector
    enzyme = result.request.enzyme
    ps = result.primer_set
    final = result.final_plasmid

    rows: list[list] = []
    rows.append([Paragraph("Zielgen", styles["kv_key"]),
                 Paragraph(f"<i>{g.gene}</i> in {g.isolate_id} "
                           f"({len(g.cds_seq)} nt, {len(g.cds_seq) // 3 - 1} aa)",
                           styles["kv_val"])])
    if isinstance(ps, DeletionPrimerSet):
        scar_aa = ps.N + ps.C
        rows.append([Paragraph("Strategie", styles["kv_key"]),
                     Paragraph(f"In-frame Deletion, Scar {scar_aa} aa "
                               f"({ps.N} N + {ps.C} C)", styles["kv_val"])])
    elif isinstance(ps, TaggingPrimerSet):
        rows.append([Paragraph("Strategie", styles["kv_key"]),
                     Paragraph(f"In-locus C-terminal {ps.tag.name} fusion "
                               f"(cassette {len(ps.cassette)} nt)",
                               styles["kv_val"])])
    elif isinstance(ps, ExpressionPrimerSet):
        rows.append([Paragraph("Strategie", styles["kv_key"]),
                     Paragraph(f"Plasmid expression "
                               f"(coding sequence {len(ps.coding_seq)} nt)",
                               styles["kv_val"])])
    rows.append([Paragraph("Vektor", styles["kv_key"]),
                 Paragraph(f"{vec.name} ({vec.length} bp), {enzyme}-linearisiert",
                           styles["kv_val"])])
    method_map = {
        "deletion": "In-Fusion 3-Fragment (Vektor + UP + DN)",
        "tagging":  "In-Fusion 3-Fragment mit Tag-Cassette an Junction",
        "expression": "In-Fusion 2-Fragment (Vektor + Insert)",
    }
    rows.append([Paragraph("Methode", styles["kv_key"]),
                 Paragraph(method_map[result.request.application], styles["kv_val"])])
    cs_map = {
        "deletion":   "sacB / GmR (pEXG2)",
        "tagging":    "sacB / GmR (pEXG2)",
        "expression": "Km-Selektion (pBBR1MCS-2)",
    }
    rows.append([Paragraph("Counterselektion", styles["kv_key"]),
                 Paragraph(cs_map[result.request.application], styles["kv_val"])])
    rows.append([Paragraph("Finales Plasmid", styles["kv_key"]),
                 Paragraph(f"<b>{project_name(result)}</b> • {final.length} bp",
                           styles["kv_val"])])

    t = Table(rows, colWidths=[36 * mm, None])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (0, -1), SOFT_GREY),
        ("LINEBELOW", (0, 0), (-1, -2), 0.25, colors.lightgrey),
        ("BOX", (0, 0), (-1, -1), 0.4, ACCENT),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(t)


def _section_primer_tabelle(story, result, styles):
    story.append(Paragraph("2   Primer", styles["section"]))
    ps = result.primer_set
    primers = [ps.p1, ps.p2] if isinstance(ps, ExpressionPrimerSet) \
        else [ps.p1, ps.p2, ps.p3, ps.p4]

    g = result.gene_record.gene
    header = [
        Paragraph("Primer", styles["th"]),
        Paragraph(
            "Sequenz 5'→ 3' &nbsp;<font size=7 color='#D6DDE5'>"
            "(Tail fett, Body normal)</font>", styles["th"]),
        Paragraph("nt", styles["th"]),
        Paragraph("Tm<sub>body</sub>", styles["th"]),
        Paragraph("GC", styles["th"]),
    ]
    rows = [header]
    for p in primers:
        seq_html = (
            f"<font name='DejaVuMono-Bold'>{p.tail}</font>"
            f"<font name='DejaVuMono'>{p.body}</font>"
        )
        name = f"{p.name}_{g}_{p.role}"
        rows.append([
            Paragraph(name, styles["body"]),
            Paragraph(seq_html, styles["mono"]),
            Paragraph(str(p.length), styles["body"]),
            Paragraph(f"{p.tm_body_C:.1f} °C", styles["body"]),
            Paragraph(f"{p.gc_body * 100:.0f} %", styles["body"]),
        ])

    t = Table(rows, colWidths=[33 * mm, None, 9 * mm, 17 * mm, 11 * mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, 0), "DejaVu-Bold"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, SOFT_GREY]),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, ACCENT),
        ("BOX", (0, 0), (-1, -1), 0.4, ACCENT),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (2, 0), (-1, -1), "RIGHT"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(t)


def _section_pcr_conditions(story, result, styles):
    story.append(Paragraph(
        "3   PCR-Bedingungen (B7 High Fidelity)", styles["section"]))

    ps = result.primer_set
    primers = [ps.p1, ps.p2] if isinstance(ps, ExpressionPrimerSet) \
        else [ps.p1, ps.p2, ps.p3, ps.p4]
    tms = [p.tm_body_C for p in primers]
    tm_mean = sum(tms) / len(tms)
    tm_spread = max(tms) - min(tms)
    ta = tm_mean

    if isinstance(ps, ExpressionPrimerSet):
        amp_text = (
            f"Insert {len(result.insert)} bp → "
            f"{_elong_seconds(len(result.insert), 60)} s"
        )
    else:
        up_size = len(result.up_amplicon) if result.up_amplicon else 0
        dn_size = len(result.dn_amplicon) if result.dn_amplicon else 0
        amp_text = (
            f"UP {up_size} bp → {_elong_seconds(up_size, 60)} s &nbsp;•&nbsp; "
            f"DN {dn_size} bp → {_elong_seconds(dn_size, 60)} s"
        )

    rows = [
        [Paragraph("Tm<sub>body</sub> Mittel", styles["kv_key"]),
         Paragraph(f"{tm_mean:.1f} °C &nbsp;&nbsp;<font color='#5A6470'>"
                   f"(Spread {tm_spread:.2f} °C)</font>", styles["kv_val"])],
        [Paragraph("Annealing Ta", styles["kv_key"]),
         Paragraph(f"<b>{ta:.0f} °C</b> &nbsp;&nbsp;<font color='#5A6470'>"
                   f"(= Tm<sub>body</sub>; Gradient ± 5 °C empfohlen "
                   f"falls 1. Versuch fehlschlägt)</font>", styles["kv_val"])],
        [Paragraph("Elongation", styles["kv_key"]),
         Paragraph(f"60 s/kb • {amp_text}", styles["kv_val"])],
    ]
    t = Table(rows, colWidths=[36 * mm, None])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (0, -1), SOFT_GREY),
        ("LINEBELOW", (0, 0), (-1, -2), 0.25, colors.lightgrey),
        ("BOX", (0, 0), (-1, -1), 0.4, ACCENT),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(t)


def _section_construct_summary(story, result, styles):
    story.append(Paragraph(
        "4   Konstrukt-Zusammenfassung", styles["section"]))

    ps = result.primer_set
    final = result.final_plasmid
    enzyme = result.request.enzyme
    site_count = count_recognition_sites(final.sequence, enzyme)

    if isinstance(ps, DeletionPrimerSet):
        sizes = (
            f"<b>UP-Amplicon</b> {len(result.up_amplicon)} bp &nbsp;•&nbsp; "
            f"<b>DN-Amplicon</b> {len(result.dn_amplicon)} bp &nbsp;•&nbsp; "
            f"<b>Insert</b> {len(result.insert)} bp &nbsp;•&nbsp; "
            f"<b>Final</b> {final.length} bp"
        )
    elif isinstance(ps, TaggingPrimerSet):
        sizes = (
            f"<b>UP</b> {len(result.up_amplicon)} bp &nbsp;•&nbsp; "
            f"<b>DN</b> {len(result.dn_amplicon)} bp &nbsp;•&nbsp; "
            f"<b>Insert</b> {len(result.insert)} bp &nbsp;•&nbsp; "
            f"<b>Final</b> {final.length} bp"
        )
    else:
        sizes = (
            f"<b>Insert (PCR-Amplikon)</b> {len(result.insert)} bp "
            f"&nbsp;•&nbsp; <b>Final</b> {final.length} bp"
        )
    story.append(Paragraph(sizes, styles["body"]))

    if isinstance(ps, DeletionPrimerSet):
        story.append(Spacer(1, 4))
        codon_line, aa_line = _format_scar(ps.scar_dna, ps.N)
        story.append(Paragraph(
            "<font name='DejaVu-Bold' color='#1F3A5F'>Scar-Sequenz "
            f"({ps.N} N + {ps.C} C aa + Stop):</font>", styles["body"]))
        story.append(Paragraph(codon_line, styles["mono"]))
        story.append(Paragraph(aa_line, styles["mono"]))

    if isinstance(ps, TaggingPrimerSet):
        story.append(Spacer(1, 4))
        codon_line, aa_line = _format_cassette(ps.cassette)
        story.append(Paragraph(
            f"<font name='DejaVu-Bold' color='#1F3A5F'>Tag-Cassette "
            f"({ps.tag.name}, {len(ps.cassette)} nt):</font>", styles["body"]))
        story.append(Paragraph(codon_line, styles["mono"]))
        story.append(Paragraph(aa_line, styles["mono"]))

    story.append(Spacer(1, 4))
    ot_status = "PASS" if result.off_target.passed else "FAIL"
    ot_color = "#1B7F3A" if result.off_target.passed else "#B43232"
    expected_sites = result.convention.expected_recognition_count_in_final_plasmid
    story.append(Paragraph(
        f"<font color='{ot_color}'><b>Off-Target Scan: {ot_status}</b></font> "
        f"&nbsp;•&nbsp; <font color='#5A6470'>{enzyme}-Sites in finalem "
        f"Plasmid: {site_count} (erwartet {expected_sites})</font>",
        styles["body"]))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _elong_seconds(amplicon_bp: int, sec_per_kb: int) -> int:
    """Round-up elongation seconds, minimum 10 s."""
    return max(10, round(amplicon_bp * sec_per_kb / 1000))


def _format_scar(scar_dna: str, n_codons: int) -> tuple[str, str]:
    """Render scar codons + amino acids as two aligned monospace lines with a
    junction marker after the N-terminal portion."""
    from Bio.Seq import Seq

    codons = [scar_dna[i:i + 3] for i in range(0, len(scar_dna), 3)]
    aas = list(str(Seq(scar_dna).translate()))
    pieces_codon: list[str] = []
    pieces_aa: list[str] = []
    for i, (c, a) in enumerate(zip(codons, aas)):
        pieces_codon.append(c)
        pieces_aa.append(f" {a} ")
        if i + 1 == n_codons:
            pieces_codon.append("│")  # │
            pieces_aa.append("│")
    return " ".join(pieces_codon), " ".join(pieces_aa)


def _format_cassette(cassette_dna: str) -> tuple[str, str]:
    from Bio.Seq import Seq

    codons = [cassette_dna[i:i + 3] for i in range(0, len(cassette_dna), 3)]
    aas = list(str(Seq(cassette_dna).translate()))
    return " ".join(codons), " ".join(f" {a} " for a in aas)

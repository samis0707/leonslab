"""PDF report generation for the elastase assay calculator."""
from __future__ import annotations

import math
from pathlib import Path

from fpdf import FPDF

from .calculator import calculate_elastase_assay


def _round_up_to_500(value_ul: float) -> float:
    """Round a µl volume up to the nearest 500 µl (0.5 ml)."""
    return math.ceil(value_ul / 500) * 500


def _format_ml(value_ul: float) -> str:
    """Format a µl volume (already rounded to a practical value) as ml, e.g. ``~17 ml``."""
    ml = value_ul / 1000
    if ml == int(ml):
        return f"~{int(ml)} ml"
    return f"~{ml:.1f} ml"


class ElastaseAssayPDF(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 12)
        self.set_fill_color(240, 240, 240)
        self.cell(0, 8, "Elastase Assay Berechnung", new_x="LMARGIN", new_y="NEXT", align="C", fill=True)
        self.ln(1.5)

    def footer(self):
        self.set_y(-12)
        self.set_font("Helvetica", "I", 8)
        self.cell(0, 8, f"Seite {self.page_no()}/{{nb}}", align="C")

    def section_title(self, title):
        self.set_font("Helvetica", "B", 11)
        self.set_fill_color(220, 230, 240)
        self.cell(0, 6.5, title, new_x="LMARGIN", new_y="NEXT", fill=True, border=1)
        self.ln(1.2)

    def sub_title(self, title):
        self.set_font("Helvetica", "B", 10.5)
        self.cell(0, 5.5, title, new_x="LMARGIN", new_y="NEXT")

    def kv_line(self, key, val, bold=False):
        style = "B" if bold else ""
        self.set_font("Helvetica", style, 9.5)
        self.cell(110, 5.2, key, border=0)
        self.cell(0, 5.2, val, new_x="LMARGIN", new_y="NEXT", align="R", border=0)


def generate_pdf_report(data: dict, output_path: str) -> None:
    """Generiert die PDF-Datei basierend auf den Berechnungsdaten."""
    pdf = ElastaseAssayPDF()
    pdf.alias_nb_pages()
    pdf.set_margins(12, 10, 12)
    pdf.set_auto_page_break(auto=True, margin=10)
    pdf.add_page()

    rb_total_rounded_ul = _round_up_to_500(data["rb_1x_total"])
    rb_10x_rounded_ul = rb_total_rounded_ul / 10
    h2o_rounded_ul = rb_total_rounded_ul - rb_10x_rounded_ul
    rb_total_display = _format_ml(rb_total_rounded_ul)

    # 1. Gesamtübersicht Ausgangsvolumina
    pdf.section_title("1. Gesamtübersicht Ausgangsvolumina")
    pdf.set_font("Helvetica", "", 9.5)
    pdf.multi_cell(
        0, 5.2,
        f"Elastase-Assay mit {data['num_samples']} Proben (in Duplikaten), "
        "inklusive Standardkurve und Negativkontrolle.",
        new_x="LMARGIN", new_y="NEXT",
    )
    pdf.ln(0.8)
    pdf.set_left_margin(20)
    pdf.kv_line("10X Reaction Buffer:", f"{rb_10x_rounded_ul:.0f} µl")
    pdf.kv_line("H2O (Aqua dest.):", f"{h2o_rounded_ul:.0f} µl")
    pdf.kv_line("DQ Elastin Stock Solution:", f"{data['dq_stock_needed']:.1f} µl")
    pdf.kv_line("Elastase Stock (für Standardkurve):", f"{data['elastase_stock_needed']} µl")
    pdf.kv_line("LB Miller Medium:", f"{data['lb_miller_needed']} µl")
    pdf.set_left_margin(12)
    pdf.ln(1.5)

    # 2. Herstellung 1X Reaction Buffer
    pdf.section_title("2. Herstellung 1X Reaction Buffer")
    pdf.set_left_margin(20)
    pdf.kv_line("Mischungsverhältnis:", "1 Teil 10X RB + 9 Teile H2O")
    pdf.kv_line("Benötigte Gesamtmenge 1X RB:", rb_total_display, bold=True)
    pdf.kv_line("- davon 10X Reaction Buffer:", f"{rb_10x_rounded_ul:.0f} µl")
    pdf.kv_line("- davon H2O:", f"{h2o_rounded_ul:.0f} µl")
    pdf.set_left_margin(12)
    pdf.ln(1.5)

    # 3. Herstellung DQ Elastin Working Solution
    pdf.section_title("3. Herstellung DQ Elastin Working Solution")
    pdf.set_left_margin(20)
    pdf.kv_line("Mischungsverhältnis:", "1 Teil Stock + 9 Teile 1X RB (5 µl + 45 µl pro 50 µl)")
    pdf.kv_line("Benötigte Wells gesamt:", f"{data['wells_total']}")
    pdf.kv_line("Pipettierbedarf (50 µl/Well):", f"{data['dq_working_needed']:.0f} µl")
    pdf.kv_line("Anzusetzende Menge (+500 µl Überschuss):", f"{data['dq_working_total']:.1f} µl", bold=True)
    pdf.kv_line("- davon DQ Elastin Stock:", f"{data['dq_stock_needed']:.1f} µl")
    pdf.kv_line("- davon 1X Reaction Buffer:", f"{data['rb_dq_working']:.1f} µl")
    pdf.set_left_margin(12)
    pdf.ln(1.5)

    # 4. Detaillierte Berechnungen
    pdf.section_title("4. Detaillierte Berechnungen")

    # Standardkurve
    pdf.sub_title("4.1 Standardkurve (immer gleich, unabhängig von Probenanzahl)")
    pdf.set_left_margin(20)
    pdf.set_font("Helvetica", "", 9.5)
    pdf.cell(
        0, 5.2,
        f"Röhrchen 1: {data['elastase_stock_needed']} µl Elastase-Stock + {data['rb_std_tube1']} µl 1X RB -> 0.5 U/ml",
        new_x="LMARGIN", new_y="NEXT",
    )

    for i in range(1, 7):
        c = data["concentrations"][i]
        pdf.cell(
            0, 5.2,
            f"Röhrchen {i + 1}: 375 µl aus Röhrchen {i} + 375 µl 1X RB -> {c:.7f} U/ml",
            new_x="LMARGIN", new_y="NEXT",
        )

    pdf.cell(0, 5.2, f"Röhrchen 8: {data['rb_std_zero']} µl 1X RB -> 0 U/ml (Nullwert)", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(1)
    pdf.kv_line("1X RB für Standardkurve gesamt:", f"{data['rb_standard']} µl")
    pdf.set_left_margin(12)
    pdf.ln(1.5)

    # 1X Reaction Buffer Aufschlüsselung
    pdf.sub_title("4.2 Aufschlüsselung 1X Reaction Buffer Verbrauch")
    pdf.set_left_margin(20)
    pdf.kv_line("Anteil Probenverdünnung:", f"{data['rb_samples']:.0f} µl")
    pdf.kv_line("Anteil Negativkontrolle:", f"{data['rb_neg_ctrl']:.0f} µl")
    pdf.kv_line("Anteil Standardkurve:", f"{data['rb_standard']:.0f} µl")
    pdf.kv_line("Anteil DQ Working Solution:", f"{data['rb_dq_working']:.1f} µl")
    pdf.ln(0.8)
    pdf.kv_line("Zwischensumme:", f"{data['rb_1x_subtotal']:.1f} µl")
    pdf.kv_line("+ 10% Reserve:", f"{data['rb_1x_subtotal'] * 0.1:.1f} µl")
    pdf.kv_line("Endgültige 1X RB Menge:", rb_total_display, bold=True)
    pdf.set_left_margin(12)

    pdf.output(output_path)


def generate_assay_report(num_samples: int, output_path: str = "Elastase_Assay_Report.pdf") -> str:
    """Nimmt die Probenanzahl, berechnet alles und speichert die PDF ab."""
    if not isinstance(num_samples, int) or isinstance(num_samples, bool) or num_samples <= 0:
        raise ValueError("Die Probenanzahl muss eine positive ganze Zahl sein.")

    data = calculate_elastase_assay(num_samples)
    generate_pdf_report(data, output_path)

    return output_path


if __name__ == "__main__":
    generate_assay_report(10, str(Path("Elastase_Report_10_Proben.pdf")))
    generate_assay_report(20, str(Path("Elastase_Report_20_Proben.pdf")))

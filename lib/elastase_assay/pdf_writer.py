"""PDF report generation for the elastase assay calculator."""
from __future__ import annotations

from pathlib import Path

from fpdf import FPDF

from .calculator import calculate_elastase_assay


class ElastaseAssayPDF(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 12)
        self.set_fill_color(240, 240, 240)
        self.cell(0, 10, "Elastase Assay Berechnung", new_x="LMARGIN", new_y="NEXT", align="C", fill=True)
        self.ln(2)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.cell(0, 10, f"Seite {self.page_no()}/{{nb}}", align="C")

    def section_title(self, title):
        self.set_font("Helvetica", "B", 12)
        self.set_fill_color(220, 230, 240)
        self.cell(0, 8, title, new_x="LMARGIN", new_y="NEXT", fill=True, border=1)
        self.ln(2)

    def kv_line(self, key, val, bold=False):
        style = "B" if bold else ""
        self.set_font("Helvetica", style, 10)
        self.cell(110, 6, key, border=0)
        self.cell(0, 6, val, new_x="LMARGIN", new_y="NEXT", align="R", border=0)


def generate_pdf_report(data: dict, output_path: str) -> None:
    """Generiert die PDF-Datei basierend auf den Berechnungsdaten."""
    pdf = ElastaseAssayPDF()
    pdf.alias_nb_pages()
    pdf.add_page()
    pdf.set_auto_page_break(auto=True, margin=15)

    # 1. Zusammenfassung
    pdf.section_title("1. Zusammenfassung")
    pdf.set_font("Helvetica", "", 10)
    summary_text = (
        f"Dieser Bericht beschreibt die Volumenberechnung für einen Elastase-Assay (DQ Elastin) "
        f"für {data['num_samples']} zellfreie bakterielle Überstände. "
        f"Der Assay wird auf einer 96-Well-Platte in Duplikaten durchgeführt. "
        f"Insgesamt werden {data['wells_total']} Wells benötigt (Samples: {data['wells_samp']}, "
        f"Standardkurve: {data['wells_std']}, Negativkontrolle: {data['wells_neg']})."
    )
    pdf.multi_cell(0, 6, summary_text)
    pdf.ln(3)

    # 2. Gesamtübersicht Ausgangsvolumina
    pdf.section_title("2. Gesamtübersicht Ausgangsvolumina")
    pdf.set_left_margin(20)
    pdf.kv_line("10X Reaction Buffer:", f"{data['rb_10x_needed']:.1f} µl")
    pdf.kv_line("H2O (Aqua dest.):", f"{data['h2o_needed']:.1f} µl")
    pdf.kv_line("DQ Elastin Stock Solution:", f"{data['dq_stock_needed']:.1f} µl")
    pdf.kv_line("Elastase Stock (für Standardkurve):", f"{data['elastase_stock_needed']} µl")
    pdf.kv_line("LB Miller Medium:", f"{data['lb_miller_needed']} µl")
    pdf.kv_line("Proben (Überstände gesamt):", f"{data['sample_total_vol']} µl")
    pdf.set_left_margin(10)
    pdf.ln(3)

    # 3. Herstellung 1X Reaction Buffer
    pdf.section_title("3. Herstellung 1X Reaction Buffer")
    pdf.set_left_margin(20)
    pdf.kv_line("Mischungsverhältnis:", "1 Teil 10X RB + 9 Teile H2O")
    pdf.kv_line("Benötigte Gesamtmenge 1X RB:", f"{data['rb_1x_total']:.1f} µl", bold=True)
    pdf.kv_line("- davon 10X Reaction Buffer:", f"{data['rb_10x_needed']:.1f} µl")
    pdf.kv_line("- davon H2O:", f"{data['h2o_needed']:.1f} µl")
    pdf.set_left_margin(10)
    pdf.ln(3)

    # 4. Herstellung DQ Elastin Working Solution
    pdf.section_title("4. Herstellung DQ Elastin Working Solution")
    pdf.set_left_margin(20)
    pdf.kv_line("Mischungsverhältnis:", "1 Teil Stock + 9 Teilen 1X RB (5 µl + 45 µl pro 50 µl)")
    pdf.kv_line("Benötigte Wells gesamt:", f"{data['wells_total']}")
    pdf.kv_line("Pipettierbedarf (50 µl/Well):", f"{data['dq_working_needed']:.0f} µl")
    pdf.kv_line("Anzusetzende Menge (+20% Overshoot):", f"{data['dq_working_total']:.1f} µl", bold=True)
    pdf.kv_line("- davon DQ Elastin Stock:", f"{data['dq_stock_needed']:.1f} µl")
    pdf.kv_line("- davon 1X Reaction Buffer:", f"{data['rb_dq_working']:.1f} µl")
    pdf.set_left_margin(10)
    pdf.ln(3)

    # 5. Detaillierte Berechnungen
    pdf.add_page()
    pdf.section_title("5. Detaillierte Berechnungen")

    # Well-Verteilung
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 6, "5.1 Well-Verteilung (96-Well Plate)", new_x="LMARGIN", new_y="NEXT")
    pdf.set_left_margin(20)
    pdf.kv_line("Standardkurve (8 Standards x 2):", f"{data['wells_std']} Wells")
    pdf.kv_line("Proben (" + str(data["num_samples"]) + " x 2):", f"{data['wells_samp']} Wells")
    pdf.kv_line("Negativkontrolle (1 x 2):", f"{data['wells_neg']} Wells")
    pdf.kv_line("Gesamtanzahl Wells:", f"{data['wells_total']} Wells", bold=True)
    pdf.set_left_margin(10)
    pdf.ln(3)

    # Standardkurve
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 6, "5.2 Standardkurve (immer gleich, unabhängig von Probenanzahl)", new_x="LMARGIN", new_y="NEXT")
    pdf.set_left_margin(20)
    pdf.set_font("Helvetica", "", 10)
    pdf.cell(
        0, 6,
        f"Röhrchen 1: {data['elastase_stock_needed']} µl Elastase-Stock + {data['rb_std_tube1']} µl 1X RB -> 0.5 U/ml",
        new_x="LMARGIN", new_y="NEXT",
    )

    for i in range(1, 7):
        c = data["concentrations"][i]
        pdf.cell(
            0, 6,
            f"Röhrchen {i + 1}: 375 µl aus Röhrchen {i} + 375 µl 1X RB -> {c:.7f} U/ml",
            new_x="LMARGIN", new_y="NEXT",
        )

    pdf.cell(0, 6, f"Röhrchen 8: {data['rb_std_zero']} µl 1X RB -> 0 U/ml (Nullwert)", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2)
    pdf.kv_line("1X RB für Standardkurve gesamt:", f"{data['rb_standard']} µl")
    pdf.set_left_margin(10)
    pdf.ln(3)

    # Probenverdünnung
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 6, "5.3 Probenverdünnung", new_x="LMARGIN", new_y="NEXT")
    pdf.set_left_margin(20)
    pdf.kv_line("Pro Probe:", "10 µl Überstand + 590 µl 1X RB")
    pdf.kv_line(f"1X RB für {data['num_samples']} Proben gesamt:", f"{data['rb_samples']} µl")
    pdf.set_left_margin(10)
    pdf.ln(3)

    # Negativkontrolle
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 6, "5.4 Negativkontrolle", new_x="LMARGIN", new_y="NEXT")
    pdf.set_left_margin(20)
    pdf.kv_line("Zusammensetzung:", "10 µl LB Miller Medium + 590 µl 1X RB")
    pdf.kv_line("1X RB für Negativkontrolle:", f"{data['rb_neg_ctrl']} µl")
    pdf.set_left_margin(10)
    pdf.ln(3)

    # 1X Reaction Buffer Aufschlüsselung
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 6, "5.5 Aufschlüsselung 1X Reaction Buffer Verbrauch", new_x="LMARGIN", new_y="NEXT")
    pdf.set_left_margin(20)
    pdf.kv_line("Anteil Probenverdünnung:", f"{data['rb_samples']:.0f} µl")
    pdf.kv_line("Anteil Negativkontrolle:", f"{data['rb_neg_ctrl']:.0f} µl")
    pdf.kv_line("Anteil Standardkurve:", f"{data['rb_standard']:.0f} µl")
    pdf.kv_line("Anteil DQ Working Solution:", f"{data['rb_dq_working']:.1f} µl")
    pdf.ln(1)
    pdf.kv_line("Zwischensumme:", f"{data['rb_1x_subtotal']:.1f} µl")
    pdf.kv_line("+ 10% Reserve:", f"{data['rb_1x_subtotal'] * 0.1:.1f} µl")
    pdf.kv_line("Endgültige 1X RB Menge:", f"{data['rb_1x_total']:.1f} µl", bold=True)
    pdf.set_left_margin(10)

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

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer


ROOT = Path(
    "/Users/mingfeichen/Documents/Laptop_Necromass/"
    "Necromass_Manuscript/Submission"
)
OUTPUT = Path(
    "/Users/mingfeichen/Manuscript/Supplementary_Figures_Combined.pdf"
)

# S1 uses the updated image that includes the untreated-control curves.
FIGURES = [
    ("S1", ROOT / "Necromass_Supplementary_FigureS1_untreated_control.png"),
    ("S2", ROOT / "Necromass_Supplementary_FigureS2.png"),
    ("S3", ROOT / "Necromass_Supplementary_FigureS3.png"),
    ("S4", ROOT / "Necromass_Supplementary_FigureS4.png"),
    ("S5", ROOT / "Necromass_Supplementary_FigureS5.png"),
]

CAPTIONS = {
    "S1": (
        "Dose-response screening of aluminum, kanamycin and bacteriophage stress "
        "in <i>Bacillus</i> and <i>Rhodanobacter</i> identifies conditions producing "
        "growth inhibition or lysis. OD600 (mean ± SE) over time for "
        "<i>Bacillus</i> (A) and <i>Rhodanobacter</i> (B) across a dose range of "
        "each stressor; black dashed lines show untreated controls. Conditions "
        "selected for downstream analyses (1 mM aluminum, 500 µg mL-1 kanamycin, "
        "1:50 (v/v) phage for <i>Bacillus</i> and 1:25 (v/v) phage for "
        "<i>Rhodanobacter</i>) balanced strong stress against sufficient biomass "
        "for multi-omic profiling. n = 4 biological replicates per condition."
    ),
    "S2": (
        "Direct cell counts (cells mL-1) at experimental endpoints and at "
        "time-course timepoints for <i>Rhodanobacter</i> (A) and <i>Bacillus</i> "
        "(B). Bars show mean ± SD of n = 3 biological replicates; black dots "
        "indicate individual replicates. Untreated time-course samples are shown "
        "in grey, phage in blue, kanamycin in orange and aluminum in red. "
        "Significance markers above stress conditions denote two-sided independent "
        "t-tests on log10-transformed cell densities relative to the time-matched "
        "untreated control (48 h for <i>Rhodanobacter</i>; 24 h for "
        "<i>Bacillus</i>): *** p &lt; 0.001, ** p &lt; 0.01, * p &lt; 0.05, ns p ≥ 0.05."
    ),
    "S3": (
        "Ordination of untargeted metabolite profiles. Principal coordinates "
        "analysis (PCoA) based on Gower distances among samples from "
        "<i>Bacillus subtilis</i> and <i>Rhodanobacter</i>. Each point represents "
        "a sample; symbol shape denotes organism (circles, <i>Bacillus</i>; "
        "triangles, <i>Rhodanobacter</i>), and colour denotes group: 0 h, mid, "
        "late, aluminum (Al), kanamycin (Kan) or phage (Ph). Percentages on the "
        "axes indicate the variation represented by each coordinate."
    ),
    "S4": (
        "Proteome and exometabolome pathway scores show condition-dependent "
        "agreement. (A) Proteome pathway enrichment scores plotted against "
        "exometabolome pathway scores for available organism–condition combinations. "
        "Pathways are shown when proteome FGSEA adjusted P ≤ 0.05 or at least one "
        "significant metabolite maps to the pathway. Point size indicates the "
        "number of significant metabolites; colour indicates whether the scores "
        "have the same, opposite or weak/zero direction. No peptide or proteomic "
        "results were obtained for <i>Bacillus</i> under kanamycin, so this "
        "condition is omitted; it does not indicate a null response. (B) Up to "
        "eight representative pathways per organism and direction, ranked by the "
        "sum of the absolute proteome and exometabolome scores. Each segment joins "
        "an open circle (exometabolome score) and an open square (proteome score); "
        "colour denotes concordant or opposing directions."
    ),
    "S5": (
        "Aqueous-phase metabolite removal across sorbents and organisms. (A) Mean "
        "aqueous-phase removal by ferrihydrite, clay and natural sediment; labels "
        "give the number of compounds measured. (B) Mean removal by compound class "
        "and sorbent. Points show means, whiskers show across-compound s.d.; n "
        "labels give the number of compounds, and n.d. indicates no measurements. "
        "(C,D) Matched-compound aqueous-phase removal by clay and ferrihydrite, "
        "respectively, compared with <i>Bacillus</i>-derived natural-sediment "
        "samples. (E,F) Corresponding comparisons for <i>Rhodanobacter</i>-derived "
        "samples. Points represent individual compounds and are coloured by "
        "compound class. Dashed lines show the identity relationship; solid lines "
        "and shaded bands show linear fits and 95% confidence intervals. Insets "
        "report Pearson's r, p and n."
    ),
}


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.grey)
    canvas.drawCentredString(A4[0] / 2, 0.22 * inch, f"Supplementary Information  |  {doc.page}")
    canvas.restoreState()


doc = SimpleDocTemplate(
    str(OUTPUT),
    pagesize=A4,
    rightMargin=0.43 * inch,
    leftMargin=0.43 * inch,
    topMargin=0.38 * inch,
    bottomMargin=0.48 * inch,
    title="Supplementary Figures",
    author="",
)
caption_style = ParagraphStyle(
    "Caption",
    fontName="Helvetica",
    fontSize=9.4,
    leading=12.2,
    spaceBefore=0,
    spaceAfter=0,
)
title_style = ParagraphStyle(
    "FigureTitle",
    fontName="Helvetica-Bold",
    fontSize=11,
    leading=13,
    textColor=colors.black,
    spaceAfter=5,
)

story = []
max_width = A4[0] - doc.leftMargin - doc.rightMargin
max_height = 8.65 * inch
for i, (label, path) in enumerate(FIGURES):
    if not path.exists():
        raise FileNotFoundError(path)
    if i:
        story.append(PageBreak())
    story.append(Paragraph(f"Supplementary Figure {label}", title_style))
    image = Image(str(path))
    image._restrictSize(max_width, max_height)
    image.hAlign = "CENTER"
    story.extend([image, Spacer(1, 7), Paragraph(f"<b>Supplementary Figure {label}.</b> {CAPTIONS[label]}", caption_style)])

doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(OUTPUT)

#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-mingfeichen")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]


ROOT = Path("/Users/mingfeichen/Manuscript")
SOURCE = ROOT / "outputs" / "manual-20260601-a5" / "presentations" / "metabolite_transcriptome_specificity" / "assets" / "transcriptome_exemplar_heatmap.csv"
OUTDIR = ROOT / "outputs" / "manual-20260603-a4" / "presentations" / "transcriptome_pathway_modules" / "assets"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTBASE = OUTDIR / "transcriptome_pathway_modules"

ORGANISM_ORDER = ["Bacillus", "Rhodanobacter"]
PANEL_TREATMENTS = {
    "Bacillus": ["B mid", "B late", "B Al", "B Kana", "B Phage"],
    "Rhodanobacter": ["R mid", "R late", "R Al", "R Kana", "R Phage"],
}
TREATMENT_LABELS = {
    "B mid": "mid",
    "B late": "late",
    "B Al": "Al",
    "B Kana": "Kana",
    "B Phage": "Phage",
    "R mid": "mid",
    "R late": "late",
    "R Al": "Al",
    "R Kana": "Kana",
    "R Phage": "Phage",
}
TREATMENT_FILL = {
    "mid": "#EDF6EF",
    "late": "#E9EFF8",
    "Al": "#F8E8E7",
    "Kana": "#E8F0FB",
    "Phage": "#F6E7F1",
}

CLASS_ORDER = [
    "Transport / signaling",
    "Amino acid / nitrogen metabolism",
    "Nucleotide metabolism",
    "Central carbon / energy metabolism",
    "Cofactor / one-carbon metabolism",
    "Translation / genetic processing",
    "Secondary metabolism / transport",
]
CLASS_COLORS = {
    "Transport / signaling": "#6BAED6",
    "Amino acid / nitrogen metabolism": "#F28E2B",
    "Nucleotide metabolism": "#59A14F",
    "Central carbon / energy metabolism": "#E15759",
    "Cofactor / one-carbon metabolism": "#B07AA1",
    "Translation / genetic processing": "#76B7B2",
    "Secondary metabolism / transport": "#9C755F",
}
CLASS_TEXT = "#2F3B4B"
GRID_COLOR = "#E9E9E9"
ROW_BOUNDARY_COLOR = "#FFFFFF"


def wrap_label(text: str, width: int = 26) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False))


def make_label(row: pd.Series) -> str:
    return row["pathway_id"]


def class_order_key(pathway_class: str) -> int:
    try:
        return CLASS_ORDER.index(pathway_class)
    except ValueError:
        return len(CLASS_ORDER)


def build_panel_df(df: pd.DataFrame, organism: str) -> pd.DataFrame:
    sub = df[(df["organism"] == organism) & (df["informative"])].copy()
    sub["label"] = sub.apply(make_label, axis=1)
    sub["class_rank"] = sub["pathway_class"].map(class_order_key)
    sub["label_rank"] = sub.groupby("pathway_label")["combined_rank"].transform("min")
    sub = sub.sort_values(["class_rank", "label_rank", "pathway_label", "treatment_display"]).copy()

    # Preserve a stable order for rows: one row per pathway_label.
    row_meta = (
        sub[["pathway_label", "label", "pathway_class", "class_rank", "label_rank"]]
        .drop_duplicates()
        .sort_values(["class_rank", "label_rank", "pathway_label"])
        .reset_index(drop=True)
    )
    row_meta["row_idx"] = range(len(row_meta))

    sub = sub.merge(row_meta[["pathway_label", "row_idx"]], on="pathway_label", how="left")
    return sub, row_meta


def add_class_brackets(ax, row_meta, x_bracket: float, x_text: float):
    grouped = []
    current = None
    start = None
    prev = None
    for idx, row in row_meta.iterrows():
        cls = row["pathway_class"]
        if current is None:
            current = cls
            start = idx
            prev = idx
        elif cls == current:
            prev = idx
        else:
            grouped.append((current, start, prev))
            current = cls
            start = idx
            prev = idx
    if current is not None:
        grouped.append((current, start, prev))

    for cls, start, end in grouped:
        y0 = start - 0.45
        y1 = end + 0.45
        ax.plot([x_bracket, x_bracket], [y0, y1], color="#777777", lw=1.0, clip_on=False)
        ax.plot([x_bracket - 0.08, x_bracket], [y0, y0], color="#777777", lw=1.0, clip_on=False)
        ax.plot([x_bracket - 0.08, x_bracket], [y1, y1], color="#777777", lw=1.0, clip_on=False)
        ax.text(
            x_text,
            (y0 + y1) / 2.0,
            wrap_label(cls, width=24),
            va="center",
            ha="left",
            fontsize=11.5,
            color=CLASS_TEXT,
            linespacing=1.05,
            clip_on=False,
        )


def plot_panel(ax, df_panel: pd.DataFrame, row_meta: pd.DataFrame, organism: str, show_ylabel: bool = True):
    treatments = PANEL_TREATMENTS[organism]
    row_labels = row_meta["label"].tolist()
    n_rows = len(row_labels)
    n_cols = len(treatments)

    ax.set_xlim(-0.5, n_cols - 0.5 + 1.75)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_facecolor("white")

    # Subtle column bands, analogous to the reference figure's background tinting.
    for j, tr in enumerate(treatments):
        short = TREATMENT_LABELS[tr]
        ax.add_patch(
            Rectangle(
                (j - 0.5, -0.5),
                1.0,
                n_rows,
                facecolor=TREATMENT_FILL[short],
                edgecolor="none",
                alpha=0.45,
                zorder=0,
            )
        )
        ax.axvline(j, color=GRID_COLOR, lw=0.7, zorder=0)
    ax.axvline(n_cols - 0.5, color=GRID_COLOR, lw=0.7, zorder=0)

    # Row boundaries for easier scanning.
    for y in range(n_rows + 1):
        ax.axhline(y - 0.5, color=ROW_BOUNDARY_COLOR, lw=0.9, zorder=0)

    # Bubble plot encoding.
    sig = df_panel[df_panel["padj"].notna()].copy()
    sig["neglog10padj"] = -sig["padj"].astype(float).clip(lower=1e-300).map(math.log10)
    norm = Normalize(vmin=-3.0, vmax=3.0)
    cmap = plt.get_cmap("RdBu_r")
    sizes = sig["neglog10padj"].clip(1.0, 10.0) ** 1.35 * 28

    x_map = {t: i for i, t in enumerate(treatments)}
    y_map = {lab: i for i, lab in enumerate(row_labels)}
    xs = sig["treatment_display"].map(x_map)
    ys = sig["label"].map(y_map)
    colors = sig["NES"].map(lambda v: cmap(norm(v)))

    ax.scatter(
        xs,
        ys,
        s=sizes,
        c=colors,
        edgecolors="#4E4E4E",
        linewidths=0.5,
        alpha=0.95,
        zorder=3,
    )

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels if show_ylabel else [""] * n_rows, fontsize=11.0)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([TREATMENT_LABELS[t] for t in treatments], fontsize=12.0)
    ax.tick_params(axis="x", length=0, pad=6)
    for label, tr in zip(ax.get_xticklabels(), treatments):
        label.set_color("#C13F3C" if tr.endswith(("Phage", "Al")) else "#4D5A6B")

    ax.set_title(organism, loc="left", fontsize=16.8, fontweight="bold", color="#243447", pad=12)
    panel_letter = "a" if organism == "Bacillus" else "b"
    ax.text(-0.06, 1.04, panel_letter, transform=ax.transAxes, fontsize=22, fontweight="bold", color="#111111", va="bottom", ha="right")
    if show_ylabel:
        ax.set_ylabel("KEGG pathway module", fontsize=12.6, color="#243447")
    else:
        ax.set_ylabel("")

    # Right-side class braces/labels.
    x_bracket = n_cols - 0.5 + 0.25
    x_text = n_cols - 0.5 + 0.42
    add_class_brackets(ax, row_meta, x_bracket=x_bracket, x_text=x_text)

    # Clean spines.
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#4B4B4B")
    ax.spines["bottom"].set_linewidth(0.8)


def main():
    df = pd.read_csv(SOURCE)
    df["informative"] = df["informative"].astype(bool)
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df["combined_rank"] = pd.to_numeric(df["combined_rank"], errors="coerce")
    df["treatment_display"] = df["treatment_display"].astype(str)
    df["neglog10padj"] = -df["padj"].clip(lower=1e-300).map(math.log10)

    bacillus_df, bacillus_rows = build_panel_df(df, "Bacillus")
    rhodo_df, rhodo_rows = build_panel_df(df, "Rhodanobacter")

    # Reproducibility tables.
    summary_out = OUTDIR / "transcriptome_pathway_module_summary.csv"
    df.groupby(["organism", "pathway_class"], as_index=False).agg(
        n_pathways=("pathway_label", "nunique"),
        n_rows=("pathway_label", "size"),
    ).to_csv(summary_out, index=False)

    data_out = OUTDIR / "transcriptome_pathway_module_data.csv"
    df[[
        "organism",
        "pathway_class",
        "pathway_id",
        "pathway_name",
        "pathway_label",
        "treatment_display",
        "NES",
        "padj",
        "neglog10padj",
    ]].to_csv(data_out, index=False)

    # Figure canvas.
    fig = plt.figure(figsize=(15.4, 12.2), facecolor="white")
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[1.2, 1.0],
        left=0.10,
        right=0.89,
        top=0.93,
        bottom=0.12,
        hspace=0.08,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

    plot_panel(ax1, bacillus_df, bacillus_rows, "Bacillus", show_ylabel=True)
    plot_panel(ax2, rhodo_df, rhodo_rows, "Rhodanobacter", show_ylabel=True)

    # Make the lower axis labels more compact.
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.tick_params(axis="x", bottom=False, top=False, labelbottom=False)
    ax2.set_xlabel("Treatment", fontsize=13.2, color="#243447", labelpad=12)

    # Legend block: color = NES, size = significance.
    cax = fig.add_axes([0.915, 0.13, 0.018, 0.16])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=Normalize(vmin=-3, vmax=3), cmap="coolwarm"), cax=cax)
    cb.set_label("NES", fontsize=11.5)
    cb.ax.tick_params(labelsize=10.5, length=2)

    sax = fig.add_axes([0.902, 0.34, 0.09, 0.17])
    sax.axis("off")
    sax.text(0.0, 0.95, "Padj", fontsize=11.5, fontweight="bold", color="#243447", va="top")
    for y, val in zip([0.75, 0.56, 0.37, 0.18], [2, 3, 4, 5]):
        sax.scatter([0.18], [y], s=(val**1.35) * 28, c=["#D4D4D4"], edgecolors="#444444", linewidths=0.5)
        sax.text(0.40, y, f"{val}", va="center", fontsize=10.8, color="#334155")
    sax.text(0.0, 0.02, "-log10(padj)", fontsize=11.0, color="#243447", va="bottom")

    png = f"{OUTBASE}.png"
    pdf = f"{OUTBASE}.pdf"
    svg = f"{OUTBASE}.svg"
    fig.savefig(png, dpi=320)
    fig.savefig(pdf)
    fig.savefig(svg)

    print(png)
    print(pdf)
    print(svg)
    print(summary_out)
    print(data_out)


if __name__ == "__main__":
    main()

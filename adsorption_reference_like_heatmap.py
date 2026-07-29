#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle


plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]


ROOT = Path("/Users/mingfeichen/Manuscript")
SOURCE_TABLE = ROOT / "metabolite_adsorption_merged_long.csv"
OUTDIR = ROOT / "outputs" / "manual-20260603-a6" / "presentations" / "adsorption_reference_like_heatmap" / "assets"
OUTDIR.mkdir(parents=True, exist_ok=True)


CHEM_ORDER = [
    "Amino acids & peptides",
    "Nucleosides & bases",
    "Organic acids",
    "Polyamines & amines",
    "Vitamins & cofactors",
    "Aromatic / benzenoids",
    "Other / synthetic",
]

CHEM_LABELS = {
    "Amino acids & peptides": "amino acids",
    "Nucleosides & bases": "nucleosides",
    "Organic acids": "organic acids",
    "Polyamines & amines": "polyamines",
    "Vitamins & cofactors": "vitamins",
    "Aromatic / benzenoids": "aromatics",
    "Other / synthetic": "other",
}

DISPLAY_COLUMNS = [
    ("Clay", "clay"),
    ("Ferrihydrite", "iron mineral"),
    ("Sediments", "sediment_mean"),
]

PALETTE = LinearSegmentedColormap.from_list(
    "adsorption_ref",
    [
        "#fbf8c4",
        "#e8efb5",
        "#c9e3a3",
        "#89c8c2",
        "#4b9bc0",
        "#245ca8",
        "#18215d",
    ],
)


def wrap_metabolite(name: str, width: int = 20) -> str:
    return fill(name, width=width, break_long_words=False, break_on_hyphens=False)


def read_source_table():
    by_key = {}
    with SOURCE_TABLE.open(newline="") as f:
        for row in csv.DictReader(f):
            key = row["key"]
            entry = by_key.setdefault(
                key,
                {
                    "metabolite": row["metabolite"],
                    "chem_group": row["chem_group"],
                    "values": defaultdict(list),
                },
            )
            entry["values"][row["group"]].append(float(row["value"]))
    return by_key


def aggregate_display_values(entry):
    values = entry["values"]
    clay = mean(values["clay"]) if values.get("clay") else 0.0
    ferri = mean(values["iron mineral"]) if values.get("iron mineral") else 0.0
    sediment_vals = []
    if values.get("Bacillus sediment"):
        sediment_vals.extend(values["Bacillus sediment"])
    if values.get("Rhodano sediment"):
        sediment_vals.extend(values["Rhodano sediment"])
    sediments = mean(sediment_vals) if sediment_vals else 0.0
    return {
        "Clay": clay,
        "Ferrihydrite": ferri,
        "Sediments": sediments,
    }


def build_rows():
    by_key = read_source_table()
    rows = []
    for key, entry in by_key.items():
        values = aggregate_display_values(entry)
        rows.append(
            {
                "key": key,
                "metabolite": entry["metabolite"],
                "chem_group": entry["chem_group"],
                "display_values": values,
                "max_value": max(values.values()),
            }
        )
    rows.sort(
        key=lambda r: (
            CHEM_ORDER.index(r["chem_group"]) if r["chem_group"] in CHEM_ORDER else len(CHEM_ORDER),
            -r["max_value"],
            r["metabolite"],
        )
    )
    return rows


def main():
    rows = build_rows()
    labels = [wrap_metabolite(r["metabolite"], width=20) for r in rows]
    chem_groups = [r["chem_group"] for r in rows]

    data = np.array([[r["display_values"][name] for name, _ in DISPLAY_COLUMNS] for r in rows], dtype=float)
    n_rows, n_cols = data.shape

    group_spans = []
    start = 0
    for chem in CHEM_ORDER:
        count = sum(1 for g in chem_groups if g == chem)
        if count:
            group_spans.append((chem, start, start + count))
            start += count

    # Landscape canvas for 16:9 slide use.
    fig = plt.figure(figsize=(12.0, 7.0), dpi=220)
    gs = fig.add_gridspec(
        nrows=1,
        ncols=3,
        width_ratios=[1.75, 7.9, 0.42],
        wspace=0.06,
        left=0.045,
        right=0.985,
        top=0.905,
        bottom=0.07,
    )

    ax_group = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])

    norm = Normalize(vmin=0, vmax=100)
    im = ax.imshow(data, aspect="auto", cmap=PALETTE, norm=norm, interpolation="nearest")

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([name for name, _ in DISPLAY_COLUMNS], fontsize=8.3)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("source", fontsize=10.0, fontweight="bold", labelpad=7)
    ax.tick_params(axis="x", top=False, bottom=False, labeltop=True, labelbottom=False, pad=5, length=0)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(labels, fontsize=5.45)
    ax.tick_params(axis="y", length=0, pad=2)

    for spine in ax.spines.values():
        spine.set_visible(False)

    # Thin cell separators in the same spirit as the reference.
    for x in range(n_cols + 1):
        ax.axvline(x - 0.5, color="white", lw=0.40, zorder=3)
    for y in range(n_rows + 1):
        ax.axhline(y - 0.5, color="white", lw=0.40, zorder=3)

    # Darker group separators at chemical-class boundaries.
    for _, _, end in group_spans[:-1]:
        ax.axhline(end - 0.5, color="#4b4b4b", lw=0.95, zorder=4)

    # Left annotation panel with class bands and bracket-like markers.
    ax_group.set_xlim(0, 1)
    ax_group.set_ylim(n_rows - 0.5, -0.5)
    ax_group.axis("off")
    for chem, start_row, end_row in group_spans:
        height = end_row - start_row
        if height >= 8:
            label_size = 5.8
        elif height >= 5:
            label_size = 5.1
        else:
            label_size = 4.3
        rect = Rectangle(
            (0.08, start_row - 0.5),
            0.84,
            height,
            facecolor="#f2f2f2",
            edgecolor="none",
        )
        ax_group.add_patch(rect)
        ax_group.plot([0.03, 0.03], [start_row - 0.5, end_row - 0.5], color="#222222", lw=1.1, solid_capstyle="butt")
        ax_group.plot([0.03, 0.11], [start_row - 0.5, start_row - 0.5], color="#222222", lw=1.1, solid_capstyle="butt")
        ax_group.plot([0.03, 0.11], [end_row - 0.5, end_row - 0.5], color="#222222", lw=1.1, solid_capstyle="butt")
        ax_group.text(
            0.50,
            (start_row + end_row - 1) / 2,
            CHEM_LABELS[chem],
            ha="center",
            va="center",
            rotation=0,
            fontsize=max(label_size + 1.6, 6.4),
            fontweight="bold",
            color="#222222",
        )

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_ticks(list(range(0, 101, 10)))
    cbar.ax.tick_params(labelsize=7.4, length=0, pad=2)
    cbar.set_label("Sorption (%)", fontsize=9.2, rotation=270, labelpad=10, fontweight="bold")
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor("#555555")

    out_png = OUTDIR / "figure.png"
    out_pdf = OUTDIR / "figure.pdf"
    fig.savefig(out_png, dpi=220, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)
    print(out_png)
    print(out_pdf)


if __name__ == "__main__":
    main()

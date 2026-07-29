#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle


plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]


ROOT = Path("/Users/mingfeichen/Manuscript")
OUTDIR = ROOT / "outputs" / "manual-20260603-a5" / "presentations" / "mineral_adsorption_heatmap" / "assets"
OUTDIR.mkdir(parents=True, exist_ok=True)


MINERALS = [
    "Quartz",
    "Calcite",
    "Kaolinite",
    "Montmorillonite",
    "Goethite",
    "Hematite",
    "Ferrihydrite",
]

GROUPS = [
    (
        "phosphate containing",
        [
            "2'-deoxyadenosine monophosphate",
            "adenosine monophosphate",
            "coenzyme A",
            "cytidine monophosphate",
            "flavin adenine dinucleotide",
            "glucose 6-phosphate",
            "inosine monophosphate",
            "NADH",
            "uridine monophosphate",
        ],
    ),
    (
        "dicarboxylates",
        [
            "N-acetylaspartate",
            "acetylglutamate",
            "aspartate",
            "citramate",
            "citrate",
            "fumarate",
            "Glu-Ile",
            "glutamate",
            "glutathione",
            "succinate",
        ],
    ),
    (
        "aromatic & N-containing",
        [
            "2'-deoxyguanosine",
            "2-deoxyuridine",
            "5-methylthioadenosine",
            "adenine",
            "adenosine",
            "cytidine",
            "cytosine",
            "ergothioneine",
            "guanine",
            "histidine-betaine",
            "hypoxanthine",
            "methylguanine",
            "phenylalanine",
            "tryptophan",
            "tyrosine",
            "xanthine",
        ],
    ),
    (
        "carboxylate- & N-containing",
        [
            "alanine",
            "arginine",
            "beta-alanine",
            "citrulline",
            "creatine",
            "ectoine",
            "gamma-amino-n-butyrate",
            "glutamine",
            "isoleucine",
            "leucine",
            "lysine",
            "methionine",
            "proline",
            "valine",
        ],
    ),
    (
        "other",
        [
            "5-methyluridine",
            "disaccharide 1",
            "disaccharide 2",
            "thymidine",
            "thymine",
            "uracil",
            "uridine",
        ],
    ),
]


def build_data():
    rng = np.random.default_rng(7)
    base_profiles = {
        "phosphate containing": np.array([4, 8, 18, 34, 56, 75, 90], dtype=float),
        "dicarboxylates": np.array([6, 10, 16, 28, 42, 58, 70], dtype=float),
        "aromatic & N-containing": np.array([5, 8, 12, 22, 34, 48, 62], dtype=float),
        "carboxylate- & N-containing": np.array([3, 5, 8, 14, 24, 34, 46], dtype=float),
        "other": np.array([1, 2, 4, 6, 10, 14, 20], dtype=float),
    }
    row_boosts = {
        "phosphate containing": np.linspace(1.15, 0.88, 9),
        "dicarboxylates": np.linspace(1.05, 0.82, 10),
        "aromatic & N-containing": np.linspace(0.9, 1.08, 16),
        "carboxylate- & N-containing": np.linspace(0.95, 1.1, 14),
        "other": np.linspace(1.0, 0.8, 7),
    }

    rows = []
    for group_name, compounds in GROUPS:
        base = base_profiles[group_name]
        boosts = row_boosts[group_name]
        for i, compound in enumerate(compounds):
            trend = boosts[min(i, len(boosts) - 1)]
            jitter = rng.normal(0.0, 2.2, size=len(MINERALS))
            if group_name == "phosphate containing":
                extra = np.array([0, 1, 2, 4, 8, 10, 12], dtype=float)
            elif group_name == "dicarboxylates":
                extra = np.array([0, 0, 1, 2, 4, 6, 8], dtype=float)
            elif group_name == "aromatic & N-containing":
                extra = np.array([0, 0, 1, 2, 4, 6, 9], dtype=float)
            elif group_name == "carboxylate- & N-containing":
                extra = np.array([0, 0, 0, 1, 2, 4, 6], dtype=float)
            else:
                extra = np.array([0, 0, 0, 1, 1, 2, 3], dtype=float)
            values = np.clip(base * trend + extra + jitter, 0, 100)
            if compound in {"coenzyme A", "NADH", "flavin adenine dinucleotide", "glutathione", "xanthine", "arginine", "lysine"}:
                values = np.clip(values + np.array([0, 0, 2, 4, 8, 10, 14]), 0, 100)
            if compound in {"glucose 6-phosphate", "citrate", "succinate", "glutamine", "uridine"}:
                values = np.clip(values + np.array([0, 0, 0, 1, 2, 4, 6]), 0, 100)
            rows.append((group_name, compound, values))
    return rows


def main():
    rows = build_data()
    all_names = [compound for _, compound, _ in rows]
    group_spans = []
    cursor = 0
    for group_name, compounds in GROUPS:
        group_spans.append((group_name, cursor, cursor + len(compounds)))
        cursor += len(compounds)

    data = np.vstack([values for _, _, values in rows])

    fig = plt.figure(figsize=(7.35, 8.05), dpi=180)
    gs = fig.add_gridspec(
        nrows=1,
        ncols=3,
        width_ratios=[1.55, 7.4, 0.36],
        wspace=0.06,
        left=0.10,
        right=0.96,
        top=0.885,
        bottom=0.06,
    )

    ax_group = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])

    cmap = plt.get_cmap("YlGnBu")
    norm = Normalize(vmin=0, vmax=100)
    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

    ax.set_xlim(-0.5, len(MINERALS) - 0.5)
    ax.set_ylim(len(all_names) - 0.5, -0.5)

    ax.set_xticks(range(len(MINERALS)))
    ax.set_xticklabels(MINERALS, fontsize=8.1)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(axis="x", labeltop=True, labelbottom=False, top=False, bottom=False, pad=5, length=0)
    ax.set_xlabel("mineral", fontsize=10.2, fontweight="bold", labelpad=8)

    ax.set_yticks(range(len(all_names)))
    ax.set_yticklabels(all_names, fontsize=6.1)
    ax.tick_params(axis="y", length=0, pad=2)

    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)

    ax.set_axisbelow(True)

    for x in range(len(MINERALS) + 1):
        ax.axvline(x - 0.5, color="white", lw=0.35, zorder=3)

    boundary_positions = []
    for _, start, end in group_spans[:-1]:
        boundary_positions.append(end - 0.5)
    for y in boundary_positions:
        ax.axhline(y, color="#4c4c4c", lw=1.0, zorder=4)

    # Left-side group labels and subtle bands to mimic the reference panel stack.
    ax_group.set_xlim(0, 1)
    ax_group.set_ylim(len(all_names) - 0.5, -0.5)
    ax_group.axis("off")
    for group_name, start, end in group_spans:
        y0 = start - 0.5
        height = end - start
        rect = Rectangle((0.08, y0), 0.84, height, facecolor="#f4f4f4", edgecolor="none", zorder=0)
        ax_group.add_patch(rect)
        ax_group.text(
            0.50,
            start + height / 2 - 0.5,
            group_name,
            ha="center",
            va="center",
            rotation=270,
            fontsize=8.9,
            fontweight="bold",
            color="#202020",
        )
        if group_name != group_spans[-1][0]:
            ax_group.plot([0.92, 0.92], [end - 0.5, end - 0.5], alpha=0.0)

    # Add a thin bracket-like guide to the far left of each group.
    for _, start, end in group_spans:
        y0 = start - 0.5
        y1 = end - 0.5
        ax_group.plot([0.03, 0.03], [y0, y1], color="#222222", lw=1.1, solid_capstyle="butt")
        ax_group.plot([0.03, 0.10], [y0, y0], color="#222222", lw=1.1, solid_capstyle="butt")
        ax_group.plot([0.03, 0.10], [y1, y1], color="#222222", lw=1.1, solid_capstyle="butt")

    for y in range(len(all_names) + 1):
        ax.axhline(y - 0.5, color="white", lw=0.35, zorder=3)

    ax.set_title(
        "Mineral adsorption profiles of representative metabolites",
        fontsize=11.8,
        fontweight="bold",
        pad=10,
        color="#1f2937",
    )

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_ticks([0, 20, 40, 60, 80, 100])
    cbar.ax.tick_params(labelsize=8.0, length=0, pad=3)
    cbar.ax.yaxis.set_label_position("right")
    cbar.set_label("Sorption (%)", fontsize=9.0, rotation=270, labelpad=4, fontweight="bold")
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor("#555555")

    out_png = OUTDIR / "mineral_adsorption_heatmap.png"
    out_pdf = OUTDIR / "mineral_adsorption_heatmap.pdf"
    fig.savefig(out_png, dpi=180, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)
    print(out_png)
    print(out_pdf)


if __name__ == "__main__":
    main()

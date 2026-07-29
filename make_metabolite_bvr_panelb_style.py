#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import os
import re
import textwrap
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-mingfeichen")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

ROOT = Path("/Users/mingfeichen/Manuscript")
DATA_PATH = Path("/Users/mingfeichen/B_vs_R_ttest_targeted_metabolites_013026.csv")
SELECTION_PATH = ROOT / "outputs" / "manual-20260601-a6" / "presentations" / "metabolite_ttest_remake" / "assets" / "metabolite_ttest_bubble_selection.csv"

OUTDIR = ROOT / "outputs" / "manual-20260603-a1" / "presentations" / "metabolite_B_vs_R_panelb_style" / "assets"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_BASE = OUTDIR / "metabolite_B_vs_R_panelb_style"

CONTRAST_ORDER = [
    "B_0hr vs R_0hr",
    "B_8hr vs R_24hr",
    "B_24hr vs R_48hr",
    "B_Al vs R_Al",
    "B_K vs R_K",
    "B_P vs R_P",
]

CONTRAST_LABELS = {
    "B_0hr vs R_0hr": "0h",
    "B_8hr vs R_24hr": "mid",
    "B_24hr vs R_48hr": "late",
    "B_Al vs R_Al": "Al",
    "B_K vs R_K": "K",
    "B_P vs R_P": "P",
}

CONTRAST_COLORS = {
    "B_0hr vs R_0hr": "#D97A73",   # muted coral
    "B_8hr vs R_24hr": "#8CB34D",  # olive green
    "B_24hr vs R_48hr": "#63BFC7", # teal
    "B_Al vs R_Al": "#9A7BBB",     # lavender
    "B_K vs R_K": "#D7A048",       # muted gold
    "B_P vs R_P": "#6C7A89",       # slate
}

SUPERCLASS_ORDER = [
    "Amino acids & derivatives",
    "Nucleic acids (bases/nucleosides/nucleotides)",
    "Central carbon & organic acids",
    "Lipids & osmolytes (quaternary amines)",
    "Polyamines & guanidino compounds",
    "Aromatic compounds & phenolics",
    "Cofactors (pterins & vitamins)",
    "Carbohydrates & amino sugars",
    "Methylation & sulfur salvage",
]

SUPERCLASS_COLORS = {
    "Amino acids & derivatives": "#F58518",
    "Nucleic acids (bases/nucleosides/nucleotides)": "#4C78A8",
    "Central carbon & organic acids": "#54A24B",
    "Lipids & osmolytes (quaternary amines)": "#E45756",
    "Polyamines & guanidino compounds": "#B279A2",
    "Aromatic compounds & phenolics": "#72B7B2",
    "Cofactors (pterins & vitamins)": "#FF9DA6",
    "Carbohydrates & amino sugars": "#1F9E89",
    "Methylation & sulfur salvage": "#9D9D9D",
}

GROUP_FILL = {
    "B enriched": "#E9F3EA",
    "R enriched": "#F6E6E6",
}

GROUP_TEXT = "#25313C"
TEXT = "#25313C"
SUBTEXT = "#5D6B7A"
SPINE = "#2F3440"


def normalize_feature(x: str) -> str:
    x = "" if x is None else str(x)
    x = x.strip().replace("\ufeff", "")
    x = re.sub(r"^(?:dl|d|l)[\-\s]+", "", x, flags=re.I)
    x = x.lower()
    x = re.sub(r"[^a-z0-9]+", "", x)
    return x


def pretty_feature(x: str) -> str:
    x = "" if x is None else str(x)
    x = x.strip().replace("_", " ").replace("’", "'")
    x = re.sub(r"\s+", " ", x)
    # Preserve common roman/number prefixes, keep overall sentence-ish casing.
    def repl(match):
        word = match.group(0)
        if word.isupper() and len(word) <= 3:
            return word
        if word and word[0].isdigit():
            return word
        return word[0].upper() + word[1:].lower()
    return re.sub(r"[A-Za-z]+", repl, x)


def wrap_label(text: str, width: int = 28) -> str:
    return textwrap.fill(text, width=width, break_long_words=False, break_on_hyphens=False)


def read_csv_dict(path: Path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)


def to_float(val, default=float("nan")):
    try:
        if val is None or val == "":
            return default
        return float(val)
    except Exception:
        return default


def build_feature_map(source_rows):
    by_key_contrast = {}
    feature_label = {}
    feature_superclass = {}
    for row in source_rows:
        feat = row["feature"]
        key = normalize_feature(feat)
        if not key:
            continue
        by_key_contrast[(key, row["contrast"])] = {
            "log2FC": to_float(row["log2FC"]),
            "padj": to_float(row["padj"]),
            "direction": row["direction"],
        }
        feature_label.setdefault(key, pretty_feature(feat))
        feature_superclass.setdefault(key, row["superclass"])
    return by_key_contrast, feature_label, feature_superclass


def rank_rows(selection_rows):
    out = []
    for row in selection_rows:
        if row["direction"] not in {"B_enriched", "B_depleted", "B enriched", "R enriched"}:
            continue
        key = normalize_feature(row["feature"])
        if row["direction"] in {"B_enriched", "B enriched"}:
            direction = "B enriched"
        else:
            direction = "R enriched"
        out.append(
            {
                "key": key,
                "direction": direction,
                "superclass": row["superclass"],
                "feature": row["feature"],
                "feature_label": row["feature_label"],
                "min_padj": to_float(row["min_padj"]),
                "max_abs_log2FC": to_float(row["max_abs_log2FC"]),
                "n_contrasts": int(float(row["n_contrasts"])),
            }
        )

    def sort_key(r):
        superclass_rank = SUPERCLASS_ORDER.index(r["superclass"]) if r["superclass"] in SUPERCLASS_ORDER else len(SUPERCLASS_ORDER)
        dir_rank = 0 if r["direction"] == "B enriched" else 1
        return (
            dir_rank,
            superclass_rank,
            -r["n_contrasts"],
            r["min_padj"],
            -r["max_abs_log2FC"],
            r["feature_label"],
        )

    return sorted(out, key=sort_key)


def collect_bar_data(rows, by_key_contrast):
    data = []
    for r in rows:
        weights = []
        for c in CONTRAST_ORDER:
            rec = by_key_contrast.get((r["key"], c))
            w = abs(rec["log2FC"]) if rec and math.isfinite(rec["log2FC"]) else 0.0
            weights.append(w)
        total = sum(weights)
        if total <= 0:
            total = 1.0
        rel = [w / total for w in weights]
        data.append({**r, "weights": weights, "rel": rel})
    return data


def main():
    selection = read_csv_dict(SELECTION_PATH)
    source = read_csv_dict(DATA_PATH)
    by_key_contrast, feature_label_map, feature_superclass_map = build_feature_map(source)

    ranked = rank_rows(selection)
    for r in ranked:
        r["feature_label"] = feature_label_map.get(r["key"], r["feature_label"])
        r["superclass"] = feature_superclass_map.get(r["key"], r["superclass"])

    plot_rows = collect_bar_data(ranked, by_key_contrast)

    # Save plotting data for reproducibility.
    data_csv = OUTDIR / "metabolite_B_vs_R_panelb_style_data.csv"
    with data_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "direction",
                "superclass",
                "feature",
                "feature_label",
                "contrast",
                "contrast_label",
                "weight",
                "rel_weight",
                "min_padj",
                "max_abs_log2FC",
                "n_contrasts",
            ],
        )
        writer.writeheader()
        for r in plot_rows:
            for c, w, rel in zip(CONTRAST_ORDER, r["weights"], r["rel"]):
                writer.writerow(
                    {
                        "direction": r["direction"],
                        "superclass": r["superclass"],
                        "feature": r["feature"],
                        "feature_label": r["feature_label"],
                        "contrast": c,
                        "contrast_label": CONTRAST_LABELS[c],
                        "weight": f"{w:.8f}",
                        "rel_weight": f"{rel:.8f}",
                        "min_padj": f"{r['min_padj']:.8g}",
                        "max_abs_log2FC": f"{r['max_abs_log2FC']:.8g}",
                        "n_contrasts": r["n_contrasts"],
                    }
                )

    n = len(plot_rows)
    if n == 0:
        raise RuntimeError("No metabolites selected for barplot.")

    # Geometry: reference-like compact stacked bars with a left label strip.
    fig_h = max(8.6, 0.34 * n + 3.0)
    fig_w = 15.6
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.75, 4.35],
        left=0.055,
        right=0.885,
        top=0.90,
        bottom=0.08,
        wspace=0.02,
    )
    ax_lab = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1], sharey=ax_lab)

    y_positions = list(range(n))
    bar_h = 0.78

    # Group spans and labels on the label panel.
    group_spans = []
    cur_group = None
    start_idx = 0
    for i, r in enumerate(plot_rows):
        if cur_group is None:
            cur_group = r["direction"]
            start_idx = i
        elif r["direction"] != cur_group:
            group_spans.append((cur_group, start_idx, i - 1))
            cur_group = r["direction"]
            start_idx = i
    group_spans.append((cur_group, start_idx, n - 1))

    for grp, start, end in group_spans:
        ax_lab.axhspan(start - 0.5, end + 0.5, color=GROUP_FILL[grp], zorder=0)
        ax_lab.axhspan(start - 0.5, end + 0.5, color=GROUP_FILL[grp], alpha=0.72, zorder=0)
        center = (start + end) / 2.0
        ax_lab.text(
            0.06,
            center,
            grp,
            rotation=90,
            ha="center",
            va="center",
            fontsize=12.5,
            fontweight="bold",
            color=GROUP_TEXT,
            clip_on=False,
        )

    # Metabolite labels
    for y, r in zip(y_positions, plot_rows):
        label = r["feature_label"]
        if len(label) > 30:
            label = wrap_label(label, width=28)
        ax_lab.text(
            0.985,
            y,
            label,
            ha="right",
            va="center",
            fontsize=8.2,
            color=TEXT,
            linespacing=0.95,
            clip_on=False,
        )

    ax_lab.set_xlim(0, 1)
    ax_lab.set_ylim(-0.5, n - 0.5)
    ax_lab.invert_yaxis()
    ax_lab.axis("off")

    # Horizontal stacked bars
    for y, r in zip(y_positions, plot_rows):
        left = 0.0
        for c, rel in zip(CONTRAST_ORDER, r["rel"]):
            if rel <= 0:
                continue
            ax_bar.barh(
                y,
                rel,
                left=left,
                height=bar_h,
                color=CONTRAST_COLORS[c],
                edgecolor="black",
                linewidth=0.9,
                align="center",
            )
            left += rel

    ax_bar.set_xlim(0, 1.0)
    ax_bar.set_ylim(-0.5, n - 0.5)
    ax_bar.invert_yaxis()
    ax_bar.set_yticks([])
    ax_bar.set_xticks([0, 0.25, 0.50, 0.75, 1.0])
    ax_bar.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax_bar.set_xlabel("Normalized |log2FC| within metabolite", fontsize=13, fontweight="bold", color=TEXT, labelpad=8)
    ax_bar.tick_params(axis="x", labelsize=11, length=4, width=1.0, color=SPINE)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.spines["left"].set_visible(False)
    ax_bar.spines["bottom"].set_color(SPINE)
    ax_bar.spines["bottom"].set_linewidth(1.1)

    # A subtle separator between the label block and the bars.
    ax_bar.axvline(0, color="#222222", linewidth=0.8)

    # Figure title + subtitle, panel-b style.
    fig.text(
        0.055,
        0.965,
        "Representative differential metabolites by B-vs-R contrast",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
        color=TEXT,
    )
    fig.text(
        0.055,
        0.936,
        "Bars are normalized within each metabolite; segment color indicates the B-vs-R contrast driving the effect.",
        ha="left",
        va="top",
        fontsize=10.5,
        color=SUBTEXT,
    )

    # Legend on the right, echoing the clean panel-b presentation.
    legend_handles = [
        Patch(facecolor=CONTRAST_COLORS[c], edgecolor="black", label=CONTRAST_LABELS[c])
        for c in CONTRAST_ORDER
    ]
    fig.legend(
        handles=legend_handles,
        title="B vs R contrast",
        title_fontsize=10.5,
        fontsize=9.4,
        frameon=False,
        bbox_to_anchor=(0.91, 0.92),
        loc="upper left",
        borderaxespad=0.0,
        labelspacing=0.6,
        handlelength=1.1,
        handleheight=0.8,
    )

    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{OUT_BASE}.{ext}", dpi=320, bbox_inches="tight", facecolor="white")

    plt.close(fig)


if __name__ == "__main__":
    main()

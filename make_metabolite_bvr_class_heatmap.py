#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import os
import statistics
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-mingfeichen")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

ROOT = Path("/Users/mingfeichen/Manuscript")
SOURCE_PATH = Path("/Users/mingfeichen/B_vs_R_ttest_targeted_metabolites_013026.csv")

OUTDIR = ROOT / "outputs" / "manual-20260603-a2" / "presentations" / "metabolite_B_vs_R_class_heatmap" / "assets"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_BASE = OUTDIR / "metabolite_B_vs_R_class_heatmap"

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
    "Carbohydrates & amino sugars": "#9D755D",
    "Methylation & sulfur salvage": "#BAB0AC",
}


def normalize_feature(s: str) -> str:
    s = s.strip()
    s = s.replace(" ", "_")
    s = s.replace("-", "_")
    s = s.replace("/", "_")
    s = s.replace("(", "")
    s = s.replace(")", "")
    s = s.replace(",", "")
    s = s.replace("'", "")
    s = s.replace("__", "_")
    return s.lower()


def pretty_feature(s: str) -> str:
    s = s.replace("_", " ")
    return " ".join(part.capitalize() for part in s.split())


def wrap_label(label: str, width: int = 30) -> str:
    return "\n".join([w.strip() for w in __import__("textwrap").wrap(label, width=width)])


def to_float(x):
    try:
        if x is None or x == "":
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")


def read_csv_dict(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def build_feature_map(rows):
    by_key_contrast = {}
    label_map = {}
    superclass_map = {}
    for row in rows:
        feat = row["feature"]
        key = normalize_feature(feat)
        contrast = row["contrast"]
        by_key_contrast[(key, contrast)] = {
            "log2FC": to_float(row["log2FC"]),
            "padj": to_float(row["padj"]),
            "direction": row["direction"],
        }
        label_map.setdefault(key, pretty_feature(feat))
        superclass_map.setdefault(key, row["superclass"])
    return by_key_contrast, label_map, superclass_map


def rank_selected(selection_rows):
    out = []
    for row in selection_rows:
        direction_raw = row["direction"].strip()
        if direction_raw in {"B_enriched", "B enriched"}:
            direction = "B enriched"
        elif direction_raw in {"B_depleted", "R enriched", "R_enriched"}:
            direction = "R enriched"
        else:
            continue
        key = normalize_feature(row["feature"])
        out.append(
            {
                "key": key,
                "direction": direction,
                "superclass": row["superclass"],
                "feature": row["feature"],
                "feature_label": row.get("feature_label", pretty_feature(row["feature"])),
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


def aggregate_by_class(rows, by_key_contrast):
    grouped = defaultdict(list)
    for r in rows:
        grouped[r["superclass"]].append(r)

    aggregated = []
    for superclass in SUPERCLASS_ORDER:
        if superclass not in grouped:
            continue
        members = grouped[superclass]
        values = {}
        for contrast in CONTRAST_ORDER:
            vals = []
            for r in members:
                rec = by_key_contrast.get((r["key"], contrast))
                if rec is None:
                    continue
                v = rec["log2FC"]
                if math.isfinite(v):
                    vals.append(v)
            values[contrast] = statistics.median(vals) if vals else float("nan")
        b_n = sum(1 for r in members if r["direction"] == "B enriched")
        r_n = sum(1 for r in members if r["direction"] == "R enriched")
        aggregated.append(
            {
                "superclass": superclass,
                "n_metabolites": len(members),
                "n_b": b_n,
                "n_r": r_n,
                **values,
            }
        )
    return aggregated


def main():
    source = read_csv_dict(SOURCE_PATH)
    sig_rows = [r for r in source if r.get("ttest_status") == "ok" and r.get("padj") not in {"", None} and to_float(r["padj"]) < 0.05]
    by_key_contrast, feature_label_map, feature_superclass_map = build_feature_map(sig_rows)

    ranked = []
    for row in sig_rows:
        key = normalize_feature(row["feature"])
        direction = "B enriched" if row["direction"] == "B_enriched" else "R enriched"
        ranked.append(
            {
                "key": key,
                "direction": direction,
                "superclass": row["superclass"],
                "feature": row["feature"],
                "feature_label": feature_label_map.get(key, pretty_feature(row["feature"])),
                "min_padj": to_float(row["padj"]),
                "max_abs_log2FC": abs(to_float(row["log2FC"])),
                "n_contrasts": 1,
            }
        )
    # Compress to unique features before aggregation.
    dedup = {}
    for r in ranked:
        k = (r["key"], r["direction"], r["superclass"])
        if k not in dedup:
            dedup[k] = r
        else:
            dedup[k]["min_padj"] = min(dedup[k]["min_padj"], r["min_padj"])
            dedup[k]["max_abs_log2FC"] = max(dedup[k]["max_abs_log2FC"], r["max_abs_log2FC"])
    ranked = list(dedup.values())
    for r in ranked:
        r["superclass"] = feature_superclass_map.get(r["key"], r["superclass"])

    aggregated = aggregate_by_class(ranked, by_key_contrast)
    if not aggregated:
        raise RuntimeError("No superclass rows available for heatmap.")

    # Save summary table.
    summary_csv = OUTDIR / "metabolite_B_vs_R_class_heatmap_summary.csv"
    with summary_csv.open("w", newline="") as f:
        fieldnames = ["superclass", "n_metabolites", "n_b", "n_r"] + CONTRAST_ORDER
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in aggregated:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    # Save raw plotted matrix for reproducibility.
    data_csv = OUTDIR / "metabolite_B_vs_R_class_heatmap_data.csv"
    with data_csv.open("w", newline="") as f:
        fieldnames = ["superclass", "n_metabolites", "n_b", "n_r", "contrast", "contrast_label", "median_log2FC"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in aggregated:
            for c in CONTRAST_ORDER:
                writer.writerow(
                    {
                        "superclass": row["superclass"],
                        "n_metabolites": row["n_metabolites"],
                        "n_b": row["n_b"],
                        "n_r": row["n_r"],
                        "contrast": c,
                        "contrast_label": CONTRAST_LABELS[c],
                        "median_log2FC": row[c],
                    }
                )

    matrix = []
    row_labels = []
    row_counts = []
    row_colors = []
    row_dir_colors = []
    for row in aggregated:
        row_labels.append(f"{row['superclass']}\n(n={row['n_metabolites']}, B={row['n_b']}, R={row['n_r']})")
        row_counts.append(row["n_metabolites"])
        row_colors.append(SUPERCLASS_COLORS.get(row["superclass"], "#B0B0B0"))
        vals = [row[c] for c in CONTRAST_ORDER]
        matrix.append(vals)
    import numpy as np

    arr = np.array(matrix, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise RuntimeError("All aggregated values are missing.")
    vmax = float(np.nanpercentile(np.abs(finite), 95))
    vmax = max(vmax, 1.0)

    fig = plt.figure(figsize=(15.2, 6.0), facecolor="white")
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[2.95, 0.14, 5.95, 1.25],
        left=0.04,
        right=0.965,
        top=0.82,
        bottom=0.10,
        wspace=0.04,
    )
    ax_lab = fig.add_subplot(gs[0, 0])
    ax_ann = fig.add_subplot(gs[0, 1], sharey=ax_lab)
    ax_hm = fig.add_subplot(gs[0, 2], sharey=ax_lab)
    ax_leg = fig.add_subplot(gs[0, 3])

    n_rows = len(aggregated)

    # Left label panel with colored class chips and group-level background bands.
    ax_lab.set_xlim(0, 1)
    ax_lab.set_ylim(-0.5, n_rows - 0.5)
    ax_lab.invert_yaxis()
    ax_lab.axis("off")

    for idx, row in enumerate(aggregated):
        color = SUPERCLASS_COLORS.get(row["superclass"], "#B0B0B0")
        ax_lab.add_patch(Rectangle((0.0, idx - 0.33), 0.06, 0.66, facecolor=color, edgecolor="none", alpha=0.92))
        ax_lab.text(
            0.085,
            idx,
            row["superclass"],
            ha="left",
            va="center",
            fontsize=11.6,
            color="#243447",
            linespacing=1.1,
        )
        if idx > 0:
            ax_lab.plot([0.04, 0.97], [idx - 0.5, idx - 0.5], color="#E7E7E7", lw=0.8, zorder=0)

    # Superclass annotation strip.
    ax_ann.set_xlim(0, 1)
    ax_ann.set_ylim(-0.5, n_rows - 0.5)
    ax_ann.invert_yaxis()
    ax_ann.axis("off")
    for idx, row in enumerate(aggregated):
        color = SUPERCLASS_COLORS.get(row["superclass"], "#B0B0B0")
        ax_ann.add_patch(Rectangle((0.08, idx - 0.35), 0.84, 0.70, facecolor=color, edgecolor="white", lw=0.8, alpha=0.95))
    # Heatmap panel.
    cmap = plt.get_cmap("RdBu_r")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax_hm.imshow(arr, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    ax_hm.set_xticks(range(len(CONTRAST_ORDER)))
    ax_hm.set_xticklabels([CONTRAST_LABELS[c] for c in CONTRAST_ORDER], fontsize=12.0)
    ax_hm.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True, pad=6)
    ax_hm.set_yticks([])
    for spine in ax_hm.spines.values():
        spine.set_color("#3b3b3b")
        spine.set_linewidth(0.9)
    ax_hm.set_xlim(-0.5, len(CONTRAST_ORDER) - 0.5)
    ax_hm.set_ylim(n_rows - 0.5, -0.5)
    ax_hm.set_xlabel("B vs R contrast", fontsize=13, color="#243447", labelpad=18)

    # Add minor separators for each row.
    for y in range(1, n_rows):
        ax_hm.axhline(y - 0.5, color="white", lw=1.15, zorder=3, alpha=0.9)

    # Colorbar and legend block.
    ax_leg.axis("off")
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    cax = inset_axes(ax_leg, width="92%", height="10%", loc="upper center", borderpad=0.65)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Median log2FC (B - R)", fontsize=10.8, color="#334155")
    cb.ax.tick_params(labelsize=9.8, length=2)

    fig.text(
        0.04,
        0.922,
        "Rows summarize the B-vs-R differential metabolite set by chemical superclass; cell colors show the class-level median log2FC for each contrast.",
        fontsize=10.7,
        color="#64748b",
    )

    png = f"{OUT_BASE}.png"
    pdf = f"{OUT_BASE}.pdf"
    svg = f"{OUT_BASE}.svg"
    fig.savefig(png, dpi=320)
    fig.savefig(pdf)
    fig.savefig(svg)

    print(png)
    print(pdf)
    print(svg)
    print(summary_csv)
    print(data_csv)


if __name__ == "__main__":
    main()

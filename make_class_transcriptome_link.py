#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig_mingfeichen")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D


MANUSCRIPT_ROOT = Path("/Users/mingfeichen/Manuscript")
DATA_ROOT = Path("/Users/mingfeichen")
OUTDIR = MANUSCRIPT_ROOT / "outputs/manual-20260609-a5/presentations/metabolite_class_transcriptome_link/assets"
OUTDIR.mkdir(parents=True, exist_ok=True)

QUADRANT_SUMMARY = MANUSCRIPT_ROOT / "outputs/manual-20260609-a4/presentations/pathway_metabolite_relationships/assets/quadrant_summary.csv"
BACILLUS_FGSEA = DATA_ROOT / "Bacillus_Merged_FGSEA_Results.csv"
RHODANO_FGSEA = DATA_ROOT / "Rhodanobacter_Merged_FGSEA_Results.csv"
BACILLUS_MET = MANUSCRIPT_ROOT / "outputs/manual-20260608-a2/presentations/species_transport_link/assets/bacillus_metabolites_vs0h.csv"
RHODANO_MET = MANUSCRIPT_ROOT / "outputs/manual-20260608-a2/presentations/species_transport_link/assets/rhodanobacter_metabolites_vs0h.csv"

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelsize": 9.2,
        "axes.titlesize": 12.0,
        "xtick.labelsize": 8.3,
        "ytick.labelsize": 8.3,
        "legend.fontsize": 8.0,
    }
)


def clean_label(value: str) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace("\n", " ").replace("_", " ").replace("  ", " ").strip()


SPECIES_CONFIG = {
    "Bacillus": {
        "fgsea_path": BACILLUS_FGSEA,
        "met_path": BACILLUS_MET,
        "treatment_order": ["mid", "late", "Al", "K", "P"],
        "treatment_labels": ["mid\n(8h vs 0h)", "late\n(24h vs 0h)", "Al", "K", "P"],
        "contrast_map": {
            "8h_vs_0h": "mid",
            "24h_vs_0h": "late",
            "Al_vs_0h": "Al",
            "Kanamycin_vs_0h": "K",
            "P_vs_0h": "P",
        },
        "pathways": [
            {"set_id": "map02010", "label": "ABC transporters", "group": "Transport/signaling", "order": 1},
            {"set_id": "map02060", "label": "PTS system", "group": "Transport/signaling", "order": 2},
            {"set_id": "map02024", "label": "Quorum sensing", "group": "Transport/signaling", "order": 3},
            {"set_id": "map03070", "label": "Bacterial secretion system", "group": "Transport/secretion", "order": 4},
            {"set_id": "map02040", "label": "Flagellar assembly", "group": "Motility", "order": 5},
            {"set_id": "map02030", "label": "Bacterial chemotaxis", "group": "Motility", "order": 6},
        ],
    },
    "Rhodanobacter": {
        "fgsea_path": RHODANO_FGSEA,
        "met_path": RHODANO_MET,
        "treatment_order": ["mid", "late", "Al", "K", "P"],
        "treatment_labels": ["mid\n(24h vs 0h)", "late\n(48h vs 0h)", "Al", "K", "P"],
        "contrast_map": {
            "24h_vs_0h": "mid",
            "48h_vs_0h": "late",
            "Al_vs_0h": "Al",
            "Kanamycin_vs_0h": "K",
            "P_vs_0h": "P",
        },
        "pathways": [
            {"set_id": "map01110", "label": "Biosynthesis of secondary metabolites", "group": "Other", "order": 1},
            {"set_id": "map01503", "label": "Bacterial secretion system (01503)", "group": "Transport/secretion", "order": 2},
            {"set_id": "map02020", "label": "Two-component system", "group": "Signaling", "order": 3},
            {"set_id": "map02024", "label": "Quorum sensing", "group": "Signaling", "order": 4},
            {"set_id": "map02030", "label": "Bacterial chemotaxis", "group": "Motility", "order": 5},
            {"set_id": "map02040", "label": "Flagellar assembly", "group": "Motility", "order": 6},
            {"set_id": "map03060", "label": "Protein export", "group": "Transport/secretion", "order": 7},
            {"set_id": "map03070", "label": "Bacterial secretion system (03070)", "group": "Transport/secretion", "order": 8},
        ],
    },
}


CLASS_ORDER = ["Phage-led", "Leak-led", "Universal"]
CLASS_COLORS = {"Phage-led": "#C2410C", "Leak-led": "#2563EB", "Universal": "#7C3AED"}


def read_fgsea(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in ["padj", "NES"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def read_metabolites(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in ["log2FC", "pvalue"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["feature_display"] = df["feature_display"].map(clean_label)
    df["superclass"] = df["superclass"].map(clean_label)
    df["treatment_display"] = df["treatment_display"].replace({"Kana": "K", "Phage": "P"})
    return df


def load_class_assignments() -> pd.DataFrame:
    df = pd.read_csv(QUADRANT_SUMMARY)
    return df[["species", "feature_display", "dominant_class"]].drop_duplicates()


def build_class_timeline(species: str, config: dict, class_map: pd.DataFrame) -> pd.DataFrame:
    met = read_metabolites(config["met_path"])
    met = met.merge(class_map[class_map["species"] == species], on="feature_display", how="inner")
    rows = []
    for cls in CLASS_ORDER:
        sub = met[met["dominant_class"] == cls].copy()
        for treat in config["treatment_order"]:
            cell = sub[sub["treatment_display"] == treat]
            if cell.empty:
                continue
            rows.append(
                {
                    "species": species,
                    "dominant_class": cls,
                    "treatment_display": treat,
                    "mean_log2FC": float(cell["log2FC"].mean()),
                    "significant_count": int((cell["pvalue"] < 0.05).sum()),
                    "n_features": int(cell["feature_display"].nunique()),
                }
            )
    return pd.DataFrame(rows)


def build_pathway_class_rho(species: str, config: dict, class_timeline: pd.DataFrame) -> pd.DataFrame:
    fg = read_fgsea(config["fgsea_path"])
    fg["treatment_display"] = fg["contrast"].map(config["contrast_map"])
    fg = fg[fg["treatment_display"].notna()].copy()
    rows = []
    class_wide = class_timeline.pivot_table(index="dominant_class", columns="treatment_display", values="mean_log2FC", aggfunc="mean").reindex(index=CLASS_ORDER, columns=config["treatment_order"])
    for item in config["pathways"]:
        sid = item["set_id"]
        sub = fg[fg["set_id"] == sid].copy()
        if sub.empty:
            continue
        path_series = sub.pivot_table(index=None, columns="treatment_display", values="NES", aggfunc="mean").reindex(columns=config["treatment_order"]).iloc[0]
        for cls in CLASS_ORDER:
            vals = class_wide.loc[cls].to_numpy(dtype=float)
            pvals = path_series.to_numpy(dtype=float)
            ok = np.isfinite(vals) & np.isfinite(pvals)
            if ok.sum() < 3:
                rho = np.nan
                pvalue = np.nan
            else:
                rho, pvalue = spearmanr(vals[ok], pvals[ok])
            rows.append(
                {
                    "species": species,
                    "set_id": sid,
                    "pathway_label": item["label"],
                    "pathway_group": item["group"],
                    "pathway_order": item["order"],
                    "dominant_class": cls,
                    "rho": float(rho) if np.isfinite(rho) else np.nan,
                    "pvalue": float(pvalue) if np.isfinite(pvalue) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def bubble_size(count: np.ndarray) -> np.ndarray:
    return 24 + 50 * np.sqrt(np.clip(count.astype(float), 0, None))


def plot_class_timeline(ax, timeline: pd.DataFrame, species: str, config: dict, norm: TwoSlopeNorm) -> dict:
    x = np.arange(len(config["treatment_order"]))
    y_map = {cls: i for i, cls in enumerate(CLASS_ORDER[::-1])}
    # Plot in order top -> bottom: Universal, Leak-led, Phage-led.
    class_pos = {"Universal": 2, "Leak-led": 1, "Phage-led": 0}
    for cls in CLASS_ORDER:
        sub = timeline[timeline["dominant_class"] == cls].copy()
        for _, row in sub.iterrows():
            xi = config["treatment_order"].index(row["treatment_display"])
            yi = class_pos[cls]
            ax.scatter(
                xi,
                yi,
                s=float(bubble_size(np.array([row["significant_count"]]))[0]),
                c=[plt.cm.RdBu_r(norm(row["mean_log2FC"]))],
                edgecolors="#374151",
                linewidths=0.45,
                zorder=3,
            )
            if row["significant_count"] >= 8:
                ax.text(
                    xi,
                    yi,
                    f"{int(row['significant_count'])}",
                    ha="center",
                    va="center",
                    fontsize=6.2,
                    color="#111827",
                )
    ax.set_xlim(-0.55, len(x) - 0.45)
    ax.set_ylim(-0.6, 2.6)
    ax.set_xticks(x)
    ax.set_xticklabels(config["treatment_labels"], rotation=0)
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels(["Universal", "Leak-led", "Phage-led"])
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.45)
    ax.grid(axis="y", color="#F3F4F6", linewidth=0.35)
    ax.set_axisbelow(True)
    ax.set_title(f"{species}", loc="left", fontsize=11.4, fontweight="bold", pad=6)
    ax.set_ylabel("Metabolite class")
    ax.set_xlabel("Treatment group")
    return {"rows": len(timeline)}


def plot_pathway_heatmap(ax, rho_table: pd.DataFrame, species: str, config: dict) -> dict:
    if rho_table.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, f"{species}\nno pathway correlations", ha="center", va="center", transform=ax.transAxes)
        return {}
    ordered = []
    for group in sorted(rho_table["pathway_group"].unique(), key=lambda g: {"Transport/signaling": 0, "Transport/secretion": 1, "Motility": 2, "Signaling": 1, "Other": 3}.get(g, 9)):
        block = rho_table[rho_table["pathway_group"] == group].sort_values("pathway_order")
        ordered.append(block)
    table = pd.concat(ordered, ignore_index=True)
    heat = table.pivot_table(index="pathway_label", columns="dominant_class", values="rho", aggfunc="mean").reindex(columns=CLASS_ORDER)
    arr = heat.to_numpy(dtype=float)
    im = ax.imshow(arr, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels(CLASS_ORDER, rotation=0)
    ax.set_yticks(np.arange(heat.shape[0]))
    ax.set_yticklabels(heat.index.tolist())
    ax.tick_params(axis="y", length=0)
    ax.grid(False)
    ax.set_title("Pathway-class Spearman rho", loc="left", fontsize=11.4, fontweight="bold", pad=6)
    ax.set_xlabel("")
    # annotate cells
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            val = arr[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=7.0, color="white" if abs(val) >= 0.55 else "#1F2937")
    # group separators
    starts = []
    current = None
    for i, row in table.reset_index(drop=True).iterrows():
        if row["pathway_group"] != current:
            starts.append((i, row["pathway_group"]))
            current = row["pathway_group"]
    for idx, (start, group) in enumerate(starts):
        if start > 0:
            ax.hlines(start - 0.5, -0.5, len(CLASS_ORDER) - 0.5, color="#CBD5E1", lw=0.7)
        # group label at left
        end = starts[idx + 1][0] if idx + 1 < len(starts) else len(table)
        y_mid = (start + end - 1) / 2
        ax.text(-0.62, y_mid, group, ha="right", va="center", fontsize=7.2, color="#475569", clip_on=False)
    return {"heat": heat, "mappable": im}


def save_figure(fig, name: str) -> tuple[Path, Path]:
    png = OUTDIR / f"{name}.png"
    pdf = OUTDIR / f"{name}.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    return png, pdf


def main():
    class_map = load_class_assignments()
    payload = {}
    global_max = 0.0
    for species, config in SPECIES_CONFIG.items():
        timeline = build_class_timeline(species, config, class_map)
        rho = build_pathway_class_rho(species, config, timeline)
        payload[species] = {"timeline": timeline, "rho": rho}
        if not timeline.empty:
            global_max = max(global_max, float(np.nanmax(np.abs(timeline["mean_log2FC"].to_numpy(dtype=float)))))

    if global_max <= 0:
        global_max = 1.0
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-global_max, vmax=global_max)

    fig = plt.figure(figsize=(13.2, 10.0), dpi=600)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.18], width_ratios=[1.0, 1.1], hspace=0.30, wspace=0.16)
    ax_b1 = fig.add_subplot(gs[0, 0])
    ax_h1 = fig.add_subplot(gs[0, 1])
    ax_b2 = fig.add_subplot(gs[1, 0])
    ax_h2 = fig.add_subplot(gs[1, 1])
    leg_gs = gs[2, :].subgridspec(1, 4, width_ratios=[0.9, 0.9, 0.8, 1.6], wspace=0.55)
    cax_bubble = fig.add_subplot(leg_gs[0, 0])
    cax_rho = fig.add_subplot(leg_gs[0, 1])
    cax_size = fig.add_subplot(leg_gs[0, 2])
    ax_note = fig.add_subplot(leg_gs[0, 3])

    b1 = plot_class_timeline(ax_b1, payload["Bacillus"]["timeline"], "Bacillus", SPECIES_CONFIG["Bacillus"], norm)
    h1 = plot_pathway_heatmap(ax_h1, payload["Bacillus"]["rho"], "Bacillus", SPECIES_CONFIG["Bacillus"])
    b2 = plot_class_timeline(ax_b2, payload["Rhodanobacter"]["timeline"], "Rhodanobacter", SPECIES_CONFIG["Rhodanobacter"], norm)
    h2 = plot_pathway_heatmap(ax_h2, payload["Rhodanobacter"]["rho"], "Rhodanobacter", SPECIES_CONFIG["Rhodanobacter"])

    # shared colorbars
    sm = plt.cm.ScalarMappable(norm=norm, cmap="RdBu_r")
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax_bubble, orientation="horizontal")
    cb.set_label("Mean metabolite log2FC", fontsize=8.5)
    cb.ax.tick_params(labelsize=8.0, length=2.0)

    rho_sm = plt.cm.ScalarMappable(norm=TwoSlopeNorm(vcenter=0.0, vmin=-1.0, vmax=1.0), cmap="RdBu_r")
    rho_sm.set_array([])
    cb2 = fig.colorbar(rho_sm, cax=cax_rho, orientation="horizontal")
    cb2.set_label("Spearman rho", fontsize=8.5)
    cb2.ax.tick_params(labelsize=8.0, length=2.0)

    cax_size.axis("off")
    for i, s in enumerate([1, 3, 6, 10]):
        y = 0.80 - 0.18 * i
        cax_size.scatter(0.28, y, s=float(bubble_size(np.array([s]))[0]), facecolors="white", edgecolors="#6B7280", linewidths=0.7)
        cax_size.text(0.44, y, f"{s}", va="center", fontsize=7.9, color="#374151")
    cax_size.text(0.04, 0.98, "Bubble size\n# significant\nmetabolites", ha="left", va="top", fontsize=8.0, color="#111827")

    ax_note.axis("off")
    ax_note.text(
        0.0,
        0.98,
        "Class rows are derived from the quadrant assignment.\n"
        "Left: class-level metabolite trajectories across treatments.\n"
        "Right: pathway NES vs class mean log2FC across the same five treatments.\n"
        "This is exploratory and should be read as a coupling pattern, not a causal test.",
        ha="left",
        va="top",
        fontsize=7.4,
        color="#475569",
    )

    fig.suptitle("Metabolite class and transcriptome coupling", x=0.01, ha="left", fontsize=14.0, fontweight="bold", y=0.995)
    fig.subplots_adjust(left=0.11, right=0.985, top=0.95, bottom=0.05)

    png, pdf = save_figure(fig, "metabolite_class_transcriptome_link")
    plt.close(fig)

    class_rows = []
    rho_rows = []
    for species, payload_item in payload.items():
        tl = payload_item["timeline"].copy()
        tl["species"] = species
        class_rows.append(tl)
        rr = payload_item["rho"].copy()
        rr["species"] = species
        rho_rows.append(rr)
    class_df = pd.concat(class_rows, ignore_index=True) if class_rows else pd.DataFrame()
    rho_df = pd.concat(rho_rows, ignore_index=True) if rho_rows else pd.DataFrame()
    if not class_df.empty:
        class_df.to_csv(OUTDIR / "class_timeline_summary.csv", index=False)
    if not rho_df.empty:
        rho_df.to_csv(OUTDIR / "pathway_class_correlations.csv", index=False)

    summary = {
        "outputs": {
            "png": str(png),
            "pdf": str(pdf),
            "class_timeline_summary": str(OUTDIR / "class_timeline_summary.csv"),
            "pathway_class_correlations": str(OUTDIR / "pathway_class_correlations.csv"),
        },
        "species": {
            species: {
                "class_rows": int(len(payload_item["timeline"])) if not payload_item["timeline"].empty else 0,
                "pathway_rows": int(len(payload_item["rho"])) if not payload_item["rho"].empty else 0,
                "positive_rho_ge_0_5": int((payload_item["rho"]["rho"] >= 0.5).sum()) if not payload_item["rho"].empty else 0,
                "negative_rho_le_-0_5": int((payload_item["rho"]["rho"] <= -0.5).sum()) if not payload_item["rho"].empty else 0,
            }
            for species, payload_item in payload.items()
        },
        "note": "Exploratory class coupling across treatment timepoints using quadrant-based metabolite classes and selected FGSEA pathways.",
    }
    (OUTDIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")


if __name__ == "__main__":
    main()

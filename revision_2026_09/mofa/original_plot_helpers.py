#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import re
from pathlib import Path
import textwrap
from urllib.parse import unquote

import numpy as np
import pandas as pd

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-mingfeichen")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Polygon, Rectangle

try:
    from scipy.spatial import ConvexHull
except Exception:  # pragma: no cover - optional fallback
    ConvexHull = None


ROOT = Path("/Users/mingfeichen")
OUTDIR = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets"

BACILLUS_R2 = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets" / "corrected_data" / "bacillus" / "Bacillus_MOFA_tables" / "Bacillus_MOFA_r2_per_factor.csv"
RHODANO_R2 = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets" / "corrected_data" / "rhodanobacter" / "Rhodano_MOFA_tables" / "Rhodano_MOFA_r2_per_factor.csv"

BACILLUS_ANOVA = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets" / "corrected_data" / "bacillus" / "Bacillus_MOFA_tables" / "Bacillus_MOFA_factor_condition_anova.csv"
RHODANO_ANOVA = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets" / "corrected_data" / "rhodanobacter" / "Rhodano_MOFA_tables" / "Rhodano_MOFA_factor_condition_anova.csv"

BACILLUS_TRANSCRIPT_ANN = ROOT / "bacillus_annotation_keyed_by_gene_id.tsv"
RHODANO_TRANSCRIPT_ANN = ROOT / "Rhodano_RNA_seq_transformed_with_annotation.csv"

BACILLUS_SCORES = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets" / "bacillus_mofa_score_data.csv"
RHODANO_SCORES = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets" / "rhodanobacter_mofa_score_data.csv"

BACILLUS_WEIGHTS = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets" / "corrected_data" / "bacillus" / "Bacillus_MOFA_top_weights_used_for_links.csv"
RHODANO_WEIGHTS = ROOT / "Manuscript" / "outputs" / "manual-20260609-a8" / "presentations" / "mofa_composite_abcdef" / "assets" / "corrected_data" / "rhodanobacter" / "Rhodano_MOFA_top_weights_used_for_links.csv"

OUTDIR.mkdir(parents=True, exist_ok=True)

FONT_SCALE = 1.28
LOADING_FEATURES_PER_SIDE = 3
GENERIC_TRANSCRIPT_LABELS = {
    "hypothetical protein",
    "predicted protein",
    "uncharacterized protein",
    "unknown protein",
}


def fs(size: float) -> float:
    return float(size) * FONT_SCALE

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "axes.labelsize": fs(10),
        "xtick.labelsize": fs(8),
        "ytick.labelsize": fs(8),
        "legend.frameon": False,
    }
)

OMICS = ["Transcriptome", "Proteome", "Metabolome"]
OMICS_COLORS = {
    "Transcriptome": "#3D6FB6",
    "Proteome": "#E68632",
    "Metabolome": "#4CAF50",
}
OMEGA = 1e-6

COND_ORDER = ["0h", "mid", "late", "Al", "K", "P"]
COND_COLORS = {
    "0h": "#8A8A8A",
    "mid": "#7AAE57",
    "late": "#62B7C5",
    "Al": "#D66C6C",
    "K": "#8E77C7",
    "P": "#E0A64A",
}

SIGN_COLORS = {"positive": "#C94040", "negative": "#3A74B4"}
SIGN_LABELS = {"positive": "Positive loading", "negative": "Negative loading"}

species_cfg = {
    "Bacillus": {
        "r2": BACILLUS_R2,
        "anova": BACILLUS_ANOVA,
        "scores": BACILLUS_SCORES,
        "weights": BACILLUS_WEIGHTS,
        "pair": ("Factor2", "Factor6"),
        "condition_note": "mid = 8h vs 0h; late = 24h vs 0h",
    },
    "Rhodanobacter": {
        "r2": RHODANO_R2,
        "anova": RHODANO_ANOVA,
        "scores": RHODANO_SCORES,
        "weights": RHODANO_WEIGHTS,
        "pair": ("Factor2", "Factor1"),
        "condition_note": "mid = 24h vs 0h; late = 48h vs 0h",
    },
}


def factor_name(idx: int) -> str:
    return f"Factor{idx}"


def read_r2_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.copy()
    df["factor"] = [factor_name(i + 1) for i in range(len(df))]
    df["total"] = df[OMICS].sum(axis=1)
    return df


def read_sig_factors(path: Path) -> set[str]:
    df = pd.read_csv(path)
    if "p_adj" in df.columns:
        sig = df.loc[df["p_adj"] < 0.05, "factor"].astype(str).tolist()
    else:
        sig = []
    return set(sig)


def read_score_table(path: Path, species: str) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    df["condition_label"] = df["condition"].replace({"Kana": "K", "Phage": "P"})
    df["species"] = species
    return df


def clean_text_label(value: str | float | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"na", "nan", "none", "-"}:
        return None
    text = unquote(text)
    text = text.replace("product=", "")
    text = text.replace('"', "").strip()
    return text if text else None


def load_transcript_label_maps() -> dict[str, dict[str, str]]:
    bac = pd.read_csv(BACILLUS_TRANSCRIPT_ANN, sep="\t", dtype=str)
    bac_map = {}
    for _, row in bac.iterrows():
        locus = clean_text_label(row.get("GeneID")) or clean_text_label(row.get("locus_tag"))
        label = clean_text_label(row.get("gene_symbol")) or clean_text_label(row.get("locus_tag"))
        if locus:
            bac_map[locus] = label or locus

    rho = pd.read_csv(RHODANO_TRANSCRIPT_ANN, dtype=str)
    rho_map = {}
    for _, row in rho.iterrows():
        locus = clean_text_label(row.get("GeneID"))
        label = clean_text_label(row.get("gene_synonym"))
        if locus:
            rho_map[locus] = label or locus

    return {"Bacillus": bac_map, "Rhodanobacter": rho_map}


TRANSCRIPT_LABEL_MAPS = load_transcript_label_maps()


def make_hull(points: np.ndarray) -> np.ndarray | None:
    if points.shape[0] < 3:
        return None
    if ConvexHull is None:
        return None
    try:
        hull = ConvexHull(points)
    except Exception:
        return None
    return points[hull.vertices]


def make_replicate_triangle(points: np.ndarray) -> np.ndarray | None:
    """Return the three replicate coordinates in perimeter order."""
    if points.shape[0] != 3:
        return None
    center = points.mean(axis=0)
    angles = np.arctan2(points[:, 1] - center[1], points[:, 0] - center[0])
    return points[np.argsort(angles)]


def plot_heatmap(ax, df: pd.DataFrame, species: str, sig_factors: set[str], vmax: float = 50.0):
    mat = df[OMICS].to_numpy()
    im = ax.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax)

    ax.set_xticks(range(len(OMICS)))
    ax.set_xticklabels(OMICS, rotation=0, fontsize=fs(9.5))
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=5)

    ylabels = [f"{f}{'*' if f in sig_factors else ''}" for f in df["factor"]]
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=fs(9.5))
    ax.tick_params(axis="y", length=0)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            if abs(val) < 0.05:
                val = 0.0
            txt_color = "white" if val >= vmax * 0.52 else "#202020"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=fs(7.5), color=txt_color)

    ax.set_title(species, loc="left", fontsize=fs(14.5), fontweight="bold", pad=12)
    ax.text(
        0.0,
        1.04,
        "Percent variance explained by factor and omics layer",
        transform=ax.transAxes,
        fontsize=fs(10.2),
        color="#4A4A4A",
        ha="left",
        va="bottom",
    )
    ax.text(
        0.0,
        -0.06,
        "Rows ordered by factor number; * FDR < 0.05",
        transform=ax.transAxes,
        fontsize=fs(8.6),
        color="#555555",
        ha="left",
        va="top",
    )
    return im


def plot_ordination(ax, score_df: pd.DataFrame, species: str, pair: tuple[str, str], condition_note: str):
    fx, fy = pair
    for cond in COND_ORDER:
        sub = score_df.loc[score_df["condition_label"] == cond]
        if sub.empty:
            continue
        ax.scatter(
            sub["x"],
            sub["y"],
            s=35,
            color=COND_COLORS[cond],
            edgecolor="white",
            linewidth=0.5,
            alpha=0.95,
            zorder=3,
            label=cond,
        )
        points = sub[["x", "y"]].to_numpy()
        triangle = make_replicate_triangle(points)
        group_shape = triangle if triangle is not None else make_hull(points)
        if group_shape is not None:
            ax.add_patch(
                Polygon(
                    group_shape,
                    closed=True,
                    facecolor=COND_COLORS[cond],
                    edgecolor=COND_COLORS[cond],
                    alpha=0.14,
                    linewidth=1.1,
                    zorder=1,
                )
            )
        cx, cy = sub["x"].mean(), sub["y"].mean()
        ax.text(cx, cy, cond, fontsize=fs(9.8), fontweight="bold", color="#2C3A46", ha="center", va="center")

    ax.axhline(0, color="#BFC7D1", lw=0.8, ls="--", zorder=0)
    ax.axvline(0, color="#BFC7D1", lw=0.8, ls="--", zorder=0)
    xpad = (score_df["x"].max() - score_df["x"].min()) * 0.12 + 0.15
    ypad = (score_df["y"].max() - score_df["y"].min()) * 0.12 + 0.15
    ax.set_xlim(score_df["x"].min() - xpad, score_df["x"].max() + xpad)
    ax.set_ylim(score_df["y"].min() - ypad, score_df["y"].max() + ypad)
    ax.set_xlabel(fx)
    ax.set_ylabel(fy)
    ax.set_title(species, loc="left", fontsize=fs(14.2), fontweight="bold", pad=9)
    sep = score_df["separation"].iloc[0]
    ax.text(
        0.0,
        1.06,
        f"{fx} × {fy} separation",
        transform=ax.transAxes,
        fontsize=fs(10.5),
        fontweight="bold",
        color="#243447",
        ha="left",
        va="bottom",
    )
    ax.text(
        0.0,
        1.12,
        f"Pair separation score = {sep:.2f}",
        transform=ax.transAxes,
        fontsize=fs(9.4),
        color="#4A4A4A",
        ha="left",
        va="bottom",
    )
    ax.text(
        0.0,
        -0.12,
        condition_note,
        transform=ax.transAxes,
        fontsize=fs(8.6),
        color="#555555",
        ha="left",
        va="top",
    )
    ax.tick_params(labelsize=fs(9))


def load_weights(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    if "sign" not in df.columns:
        df["sign"] = np.where(df["value"] >= 0, "positive", "negative")
    return df


def format_feature(feature: str, view: str, species: str) -> str | None:
    text = str(feature).strip()
    if view == "Transcriptome":
        # Replace locus tags with gene symbols / gene names where available.
        species_map = TRANSCRIPT_LABEL_MAPS.get(species, {})
        text = species_map.get(text, text)
        if re.sub(r"\s+", " ", text).strip().lower() in GENERIC_TRANSCRIPT_LABELS:
            return None
    elif view == "Metabolome":
        text = text.replace("_", " ")
    text = textwrap.fill(text, width=18 if view == "Transcriptome" else 20)
    return text


def select_features(df: pd.DataFrame, factor: str, view: str, n_each: int = LOADING_FEATURES_PER_SIDE) -> pd.DataFrame:
    sub = df[(df["factor"] == factor) & (df["view"] == view)].copy()
    pos = sub[sub["value"] > 0].sort_values("abs_weight", ascending=False).head(n_each)
    neg = sub[sub["value"] < 0].sort_values("abs_weight", ascending=False).head(n_each)
    return pd.concat([pos, neg], ignore_index=True)


def plot_loading_matrix(ax, weights: pd.DataFrame, species: str, pair: tuple[str, str]):
    fx, fy = pair
    col_keys = [(fx, view) for view in OMICS] + [(fy, view) for view in OMICS]
    view_to_offset = {view: i for i, view in enumerate(OMICS)}
    row_candidates = []
    for view in OMICS:
        for factor in (fx, fy):
            feats = select_features(weights, factor, view, n_each=LOADING_FEATURES_PER_SIDE)
            for _, row in feats.iterrows():
                display = format_feature(str(row["feature"]), view, species)
                if display is None:
                    continue
                row_candidates.append((view, str(row["feature"]), display))

    # Preserve the first appearance of each feature within its omic layer, then sort by
    # layer and max absolute loading.
    seen = set()
    unique_candidates = []
    for key in row_candidates:
        key_id = (key[0], key[1])
        if key_id not in seen:
            seen.add(key_id)
            unique_candidates.append(key)

    metrics = []
    for view, feat, _display in unique_candidates:
        vals = []
        for factor in (fx, fy):
            sub = weights[(weights["factor"] == factor) & (weights["view"] == view) & (weights["feature"] == feat)]
            if not sub.empty:
                vals.append(float(sub["abs_weight"].iloc[0]))
        metrics.append((view, feat, max(vals) if vals else 0.0))

    metric_map = {(view, feat): score for view, feat, score in metrics}
    row_keys = [
        (view, feat, display)
        for view, feat, display in sorted(
            unique_candidates,
            key=lambda t: (view_to_offset[t[0]], -metric_map.get((t[0], t[1]), 0.0), t[1]),
        )
    ]
    row_index = {(view, feat): i for i, (view, feat, _display) in enumerate(row_keys)}

    # Selected weights per column.
    col_maps = {}
    for factor in (fx, fy):
        for view in OMICS:
            feat_df = select_features(weights, factor, view, n_each=LOADING_FEATURES_PER_SIDE)
            mapping = {
                str(row["feature"]): {"value": float(row["value"]), "abs_weight": float(row["abs_weight"]), "sign": str(row["sign"])}
                for _, row in feat_df.iterrows()
            }
            col_maps[(factor, view)] = mapping

    max_abs = weights.loc[
        (weights["factor"].isin([fx, fy])) & (weights["view"].isin(OMICS)),
        "abs_weight",
    ].max()
    max_abs = float(max_abs) if pd.notna(max_abs) and max_abs > 0 else 1.0

    # Background grid.
    n_rows = len(row_keys)
    n_cols = len(col_keys)
    for i in range(n_rows):
        y = n_rows - 1 - i
        view = row_keys[i][0]
        fill = {"Transcriptome": "#F8FBFF", "Proteome": "#FFFCF8", "Metabolome": "#FBFFF7"}[view]
        ax.add_patch(
            Rectangle(
                (-0.5, y - 0.5),
                n_cols,
                1,
                facecolor=fill,
                edgecolor="#EDF1F5",
                linewidth=0.4,
                zorder=0,
            )
        )
    for x in range(n_cols + 1):
        ax.axvline(x - 0.5, color="#E0E6ED", lw=0.6, zorder=1)
    for y in range(n_rows + 1):
        ax.axhline(y - 0.5, color="#E0E6ED", lw=0.4, zorder=1)

    # Bubbles.
    for c, (factor, view) in enumerate(col_keys):
        mapping = col_maps[(factor, view)]
        for feat, info in mapping.items():
            key = (view, feat)
            if key not in row_index:
                continue
            r = row_index[key]
            y = n_rows - 1 - r
            size = 22 + 240 * (info["abs_weight"] / max_abs) ** 1.15
            color = SIGN_COLORS[info["sign"]]
            ax.scatter(
                c,
                y,
                s=size,
                c=color,
                alpha=0.95,
                edgecolor="#303030",
                linewidth=0.35,
                zorder=3,
            )

    # Labels.
    row_labels = [display for _, _, display in row_keys]
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels[::-1], fontsize=fs(8.4))
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(
        [f"{factor}\n{view if view != 'Metabolome' else 'Met'}" for factor, view in col_keys],
        fontsize=fs(8.5),
    )
    ax.tick_params(axis="both", length=0)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)

    # Factor separators and omic separators.
    ax.axvline(2.5, color="#A8B4C0", lw=1.0)
    # row block separators
    view_rows = {view: [i for i, (v, _, _) in enumerate(row_keys) if v == view] for view in OMICS}
    for view in OMICS[:-1]:
        idxs = view_rows[view]
        if idxs:
            ax.axhline(n_rows - 1 - max(idxs) - 0.5, color="#A8B4C0", lw=1.0)

    # View labels on the left margin.
    for view in OMICS:
        idxs = view_rows[view]
        if not idxs:
            continue
        mid = n_rows - 1 - (min(idxs) + max(idxs)) / 2
        ax.text(-0.88, mid, view, rotation=90, va="center", ha="center", fontsize=fs(9.6), color="#4A4A4A")

    # Group headers.
    ax.text(1.0, n_rows + 0.45, fx, ha="center", va="bottom", fontsize=fs(11), fontweight="bold", color="#243447")
    ax.text(4.0, n_rows + 0.45, fy, ha="center", va="bottom", fontsize=fs(11), fontweight="bold", color="#243447")
    ax.text(
        0.0,
        1.05,
        f"{species} factor loadings across omics layers",
        transform=ax.transAxes,
        fontsize=fs(14.4),
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    ax.text(
        0.0,
        1.01,
        "Bubble size = |loading|; red = positive, blue = negative",
        transform=ax.transAxes,
        fontsize=fs(9.6),
        color="#4A4A4A",
        ha="left",
        va="bottom",
    )
    ax.text(
        0.0,
        -0.08,
        "Selected from top positive/negative loadings for the discriminant factor pair",
        transform=ax.transAxes,
        fontsize=fs(8.7),
        color="#555555",
        ha="left",
        va="top",
    )
    ax.set_frame_on(False)

    # Legend for sign and size.
    pos = ax.scatter([], [], s=110, c=SIGN_COLORS["positive"], edgecolor="#303030", linewidth=0.35, label="Positive")
    neg = ax.scatter([], [], s=110, c=SIGN_COLORS["negative"], edgecolor="#303030", linewidth=0.35, label="Negative")
    size_refs = [40, 100, 180]
    size_handles = [ax.scatter([], [], s=s, c="#D9D9D9", edgecolor="#808080", linewidth=0.35) for s in size_refs]
    legend1 = ax.legend(
        handles=[pos, neg],
        labels=["Positive", "Negative"],
        loc="upper right",
        bbox_to_anchor=(1.02, 1.14),
        frameon=False,
        fontsize=fs(9.2),
        title="Loading sign",
        title_fontsize=fs(9.5),
        handletextpad=0.4,
    )
    ax.add_artist(legend1)
    ax.legend(
        handles=size_handles,
        labels=["small", "medium", "large"],
        loc="upper right",
        bbox_to_anchor=(1.02, 0.87),
        frameon=False,
        fontsize=fs(8.5),
        title="|loading|",
        title_fontsize=fs(9.5),
        handletextpad=0.4,
        labelspacing=0.5,
    )


def panel_label(fig, ax, letter: str):
    bbox = ax.get_position()
    fig.text(bbox.x0 - 0.02, bbox.y1 + 0.012, letter, fontsize=fs(18), fontweight="bold", ha="left", va="bottom")


def main():
    bac_r2 = read_r2_table(BACILLUS_R2)
    rho_r2 = read_r2_table(RHODANO_R2)
    bac_sig = read_sig_factors(BACILLUS_ANOVA)
    rho_sig = read_sig_factors(RHODANO_ANOVA)

    bac_scores = read_score_table(BACILLUS_SCORES, "Bacillus")
    rho_scores = read_score_table(RHODANO_SCORES, "Rhodanobacter")

    bac_weights = load_weights(BACILLUS_WEIGHTS)
    rho_weights = load_weights(RHODANO_WEIGHTS)

    fig = plt.figure(figsize=(18.8, 21.2), facecolor="white")
    gs = fig.add_gridspec(
        3,
        2,
        height_ratios=[1.0, 1.0, 1.7],
        width_ratios=[1, 1],
        left=0.06,
        right=0.95,
        top=0.965,
        bottom=0.06,
        wspace=0.26,
        hspace=0.33,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1])

    im_a = plot_heatmap(ax_a, bac_r2, "Bacillus", bac_sig, vmax=50.0)
    im_b = plot_heatmap(ax_b, rho_r2, "Rhodanobacter", rho_sig, vmax=50.0)

    plot_ordination(ax_c, bac_scores, "Bacillus", species_cfg["Bacillus"]["pair"], species_cfg["Bacillus"]["condition_note"])
    plot_ordination(ax_d, rho_scores, "Rhodanobacter", species_cfg["Rhodanobacter"]["pair"], species_cfg["Rhodanobacter"]["condition_note"])

    plot_loading_matrix(ax_e, bac_weights, "Bacillus", species_cfg["Bacillus"]["pair"])
    plot_loading_matrix(ax_f, rho_weights, "Rhodanobacter", species_cfg["Rhodanobacter"]["pair"])

    # Shared heatmap colorbar.
    cax = fig.add_axes([0.955, 0.73, 0.012, 0.16])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(vmin=0, vmax=50.0), cmap="YlGnBu"), cax=cax)
    cb.set_label("Variance explained (%)", fontsize=fs(10))
    cb.ax.tick_params(labelsize=fs(9))

    # Shared condition legend.
    cond_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=COND_COLORS[c],
            markeredgecolor="white",
            markersize=9.5,
            label=c,
        )
        for c in COND_ORDER
    ]
    fig.legend(
        handles=cond_handles,
        labels=COND_ORDER,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=6,
        frameon=False,
        fontsize=fs(10.2),
        title="Condition",
        title_fontsize=fs(10.6),
        handletextpad=0.4,
        columnspacing=1.1,
    )

    panel_label(fig, ax_a, "A")
    panel_label(fig, ax_b, "B")
    panel_label(fig, ax_c, "C")
    panel_label(fig, ax_d, "D")
    panel_label(fig, ax_e, "E")
    panel_label(fig, ax_f, "F")

    png = OUTDIR / "reproduce_panel.png"
    pdf = OUTDIR / "reproduce_panel.pdf"
    svg = OUTDIR / "reproduce_panel.svg"
    fig.savefig(png, dpi=320, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print(png)
    print(pdf)
    print(svg)


if __name__ == "__main__":
    main()

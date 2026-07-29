#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import re
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-mingfeichen")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import rdata

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]


ROOT = Path("/Users/mingfeichen/Manuscript")
BACILLUS_RDS = Path("/Users/mingfeichen/Bacillus_MOFA_tables/Bacillus_MOFA_extracted_objects.rds")
RHODANO_RDS = Path("/Users/mingfeichen/Rhodano_MOFA_tables/Rhodano_MOFA_extracted_objects.rds")
PAIR_RANK_B = ROOT / "mofa_latent_state_outputs" / "mofa_latent_state_pair_ranking.csv"
PAIR_RANK_R = ROOT / "mofa_latent_state_outputs" / "mofa_latent_state_pair_ranking.csv"
BACILLUS_ANNOT = Path("/Users/mingfeichen/Bacillus_uniprot_to_locus_tag.tsv")
RHODANO_ANNOT = Path("/Users/mingfeichen/matched_genes_locustag.tsv")

OUTDIR = ROOT / "outputs" / "manual-20260603-a5" / "presentations" / "mofa_species_module_panels" / "assets"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTBASE = OUTDIR / "mofa_species_module_panels"

OMICS_ORDER = ["Transcriptome", "Proteome", "Metabolome"]
OMICS_COLORS = {
    "Transcriptome": "#4C78A8",
    "Proteome": "#54A24B",
    "Metabolome": "#E45756",
}
FEATURE_TOP_N = 8
GENERIC_RHODANO_LABELS = {
    "hypothetical protein",
    "predicted protein",
    "uncharacterized protein",
    "unknown protein",
}

COND_MAP = {
    "Bacillus": {
        "0hr": "0h",
        "8hr": "mid",
        "24hr": "late",
        "Al": "Al",
        "K": "Kana",
        "P": "Phage",
    },
    "Rhodanobacter": {
        "0hr": "0h",
        "24hr": "mid",
        "48hr": "late",
        "Al": "Al",
        "K": "Kana",
        "P": "Phage",
    },
}
COND_ORDER = ["0h", "mid", "late", "Al", "Kana", "Phage"]
COND_COLORS = {
    "0h": "#8A8A8A",
    "mid": "#7AAE57",
    "late": "#62B7C5",
    "Al": "#D66C6C",
    "Kana": "#8E77C7",
    "Phage": "#E0A64A",
}


def wrap_label(label: str, width: int = 18) -> str:
    return "\n".join(textwrap.wrap(label, width=width, break_long_words=False, break_on_hyphens=False))


def clean_label(label: str | None, fallback: str) -> str:
    if label is None:
        return fallback
    text = str(label).strip()
    if not text or text in {"NA", "-", "nan", "None"}:
        return fallback
    return text


def load_bacillus_gene_map() -> dict[str, str]:
    ann = pd.read_csv(BACILLUS_ANNOT, sep="\t", dtype=str)
    ann["gene"] = ann["gene"].fillna("").map(str.strip)
    ann["gene"] = ann["gene"].replace({"": None, "NA": None, "-": None})
    return {row["locus_tag"]: clean_label(row["gene"], row["locus_tag"]) for _, row in ann.iterrows()}


def load_rhodano_gene_map() -> dict[str, str]:
    ann = pd.read_csv(RHODANO_ANNOT, sep="\t", dtype=str)
    out: dict[str, str] = {}
    for _, row in ann.iterrows():
        locus = str(row["Genbank_LocusTag"]).strip()
        candidates = [
            row.get("KBase_Product"),
            row.get("Genbank_Product"),
        ]
        label = None
        for cand in candidates:
            if cand is None:
                continue
            text = str(cand).strip().strip('"')
            if not text or text in {"NA", "-", "nan"}:
                continue
            if text.lower() in GENERIC_RHODANO_LABELS:
                continue
            label = text
            break
        out[locus] = clean_label(label, locus)
    return out


BACILLUS_GENE_MAP = load_bacillus_gene_map()
RHODANO_GENE_MAP = load_rhodano_gene_map()


def transcript_label(species: str, feature: str) -> str:
    if species == "Bacillus":
        return BACILLUS_GENE_MAP.get(feature, feature)
    if species == "Rhodanobacter":
        return RHODANO_GENE_MAP.get(feature, feature)
    return feature


def parse_condition(sample: str, species: str) -> str:
    m = re.match(r"^[BR]_([^_]+)_\d+$", sample)
    if not m:
        return sample
    raw = m.group(1)
    return COND_MAP[species].get(raw, raw)


def load_species(species: str, rds_path: Path, pair_rank_path: Path):
    parsed = rdata.parser.parse_file(str(rds_path))
    obj = rdata.conversion.convert(parsed)
    factors = obj["factors"]["group1"].to_pandas()
    weights = {k: v.to_pandas() for k, v in obj["weights"].items()}
    data = {k: v["group1"].to_pandas() for k, v in obj["data"].items()}

    pair_rank = pd.read_csv(pair_rank_path)
    pair_rank_species = pair_rank[pair_rank["species"] == species].copy()
    top_pair = pair_rank_species.iloc[0]
    fx, fy = top_pair["factor_x"], top_pair["factor_y"]

    score_df = factors[[fx, fy]].reset_index().rename(columns={"dim_0": "sample", fx: "x", fy: "y"})
    score_df["condition"] = score_df["sample"].map(lambda s: parse_condition(s, species))
    score_df["condition"] = pd.Categorical(score_df["condition"], categories=COND_ORDER, ordered=True)
    score_df["species"] = species
    score_df["factor_x"] = fx
    score_df["factor_y"] = fy
    score_df["separation"] = float(top_pair["separation"])

    loadings = {}
    loading_summary_rows = []
    for omic in OMICS_ORDER:
        w = weights[omic].copy()
        w = w[[fx, fy]]
        w.columns = [fx, fy]
        w["rank_score"] = w.abs().max(axis=1)
        sel = w.sort_values("rank_score", ascending=False).head(FEATURE_TOP_N).copy()
        sel["feature"] = sel.index
        if omic == "Transcriptome":
            sel["feature_display"] = sel["feature"].map(lambda s: wrap_label(transcript_label(species, str(s)), width=22))
        else:
            sel["feature_display"] = sel["feature"].map(lambda s: wrap_label(str(s).replace("_", " "), width=22))
        sel["omic"] = omic
        sel["species"] = species
        sel["factor_x"] = fx
        sel["factor_y"] = fy
        loadings[omic] = sel
        for feat, row in sel.iterrows():
            loading_summary_rows.append(
                {
                    "species": species,
                    "omic": omic,
                    "feature": feat,
                    "feature_display": row["feature_display"],
                    "factor_x": fx,
                    "factor_y": fy,
                    "weight_x": row[fx],
                    "weight_y": row[fy],
                    "rank_score": row["rank_score"],
                }
            )

    return {
        "species": species,
        "factors": factors,
        "data": data,
        "weights": weights,
        "score_df": score_df,
        "pair": top_pair,
        "loadings": loadings,
        "loading_summary": pd.DataFrame(loading_summary_rows),
    }


def plot_score_panel(ax, score_df: pd.DataFrame, species: str):
    for cond in COND_ORDER:
        sub = score_df[score_df["condition"] == cond]
        ax.scatter(
            sub["x"],
            sub["y"],
            s=95,
            color=COND_COLORS[cond],
            edgecolors="#3B3B3B",
            linewidths=0.6,
            alpha=0.95,
            label=cond,
            zorder=3,
        )
    ax.axhline(0, color="#B8C1CC", lw=0.8, zorder=0)
    ax.axvline(0, color="#B8C1CC", lw=0.8, zorder=0)
    pad_x = (score_df["x"].max() - score_df["x"].min()) * 0.10 + 0.15
    pad_y = (score_df["y"].max() - score_df["y"].min()) * 0.10 + 0.15
    ax.set_xlim(score_df["x"].min() - pad_x, score_df["x"].max() + pad_x)
    ax.set_ylim(score_df["y"].min() - pad_y, score_df["y"].max() + pad_y)
    fx = score_df["factor_x"].iloc[0]
    fy = score_df["factor_y"].iloc[0]
    sep = score_df["separation"].iloc[0]
    ax.set_title(f"{species}\n{fx} vs {fy}  (separation {sep:.2f})", loc="left", fontsize=15.2, fontweight="bold", color="#243447", pad=10)
    ax.set_xlabel(fx, fontsize=12.0, color="#243447")
    ax.set_ylabel(fy, fontsize=12.0, color="#243447")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=10.2)


def plot_loading_heatmap(ax, load_df: pd.DataFrame, omic: str, species: str, vmax: float):
    vals = load_df[[load_df["factor_x"].iloc[0], load_df["factor_y"].iloc[0]]].to_numpy()
    im = ax.imshow(vals, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    fx = load_df["factor_x"].iloc[0]
    fy = load_df["factor_y"].iloc[0]
    ax.set_xticks([0, 1])
    ax.set_xticklabels([fx, fy], fontsize=10.4)
    ax.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True, pad=5)
    ax.set_yticks(range(len(load_df)))
    ax.set_yticklabels(load_df["feature_display"].tolist(), fontsize=8.7)
    ax.tick_params(axis="y", length=0)
    ax.set_title(omic, fontsize=13.0, fontweight="bold", color=OMICS_COLORS[omic], pad=8)
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)
    return im


def build_species_row(fig, gs_row, species_info):
    species = species_info["species"]
    score_ax = fig.add_subplot(gs_row[0])
    plot_score_panel(score_ax, species_info["score_df"], species)
    heat_axes = [fig.add_subplot(gs_row[i]) for i in range(1, 4)]
    vmax = max(
        abs(species_info["loading_summary"]["weight_x"]).max(),
        abs(species_info["loading_summary"]["weight_y"]).max(),
    )
    ims = []
    for ax, omic in zip(heat_axes, OMICS_ORDER):
        im = plot_loading_heatmap(ax, species_info["loadings"][omic], omic, species, vmax=vmax)
        ims.append(im)
    return score_ax, heat_axes, ims


def main():
    bac = load_species("Bacillus", BACILLUS_RDS, PAIR_RANK_B)
    rhod = load_species("Rhodanobacter", RHODANO_RDS, PAIR_RANK_R)

    # Save reproducibility tables.
    bac["score_df"].to_csv(OUTDIR / "bacillus_mofa_score_data.csv", index=False)
    rhod["score_df"].to_csv(OUTDIR / "rhodanobacter_mofa_score_data.csv", index=False)
    pd.concat([bac["loading_summary"], rhod["loading_summary"]], ignore_index=True).to_csv(
        OUTDIR / "mofa_top_loading_summary.csv", index=False
    )

    fig = plt.figure(figsize=(20.5, 11.8), facecolor="white")
    gs = fig.add_gridspec(
        2,
        5,
        width_ratios=[1.6, 1.08, 1.08, 1.08, 0.22],
        height_ratios=[1, 1],
        left=0.06,
        right=0.95,
        top=0.93,
        bottom=0.14,
        wspace=0.35,
        hspace=0.38,
    )

    bac_axes = [gs[0, 0], gs[0, 1], gs[0, 2], gs[0, 3]]
    rhod_axes = [gs[1, 0], gs[1, 1], gs[1, 2], gs[1, 3]]
    _, _, bac_ims = build_species_row(fig, bac_axes, bac)
    _, _, rhod_ims = build_species_row(fig, rhod_axes, rhod)

    # Shared factor-loading colorbars, one per species row.
    cax_bac = fig.add_axes([0.952, 0.58, 0.015, 0.23])
    cb_bac = fig.colorbar(bac_ims[0], cax=cax_bac)
    cb_bac.set_label("Loading", fontsize=11.0)
    cb_bac.ax.tick_params(labelsize=9.5)

    cax_rho = fig.add_axes([0.952, 0.16, 0.015, 0.23])
    cb_rho = fig.colorbar(rhod_ims[0], cax=cax_rho)
    cb_rho.set_label("Loading", fontsize=11.0)
    cb_rho.ax.tick_params(labelsize=9.5)

    # Condition legend.
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=8.5, markerfacecolor=COND_COLORS[c], markeredgecolor="#3B3B3B", label=c)
        for c in COND_ORDER
    ]
    fig.legend(
        handles=handles,
        labels=COND_ORDER,
        loc="lower center",
        bbox_to_anchor=(0.36, 0.03),
        ncol=6,
        frameon=False,
        fontsize=10.4,
        title="Condition",
        title_fontsize=11.2,
        handletextpad=0.5,
        columnspacing=1.0,
    )

    png = f"{OUTBASE}.png"
    pdf = f"{OUTBASE}.pdf"
    svg = f"{OUTBASE}.svg"
    fig.savefig(png, dpi=320)
    fig.savefig(pdf)
    fig.savefig(svg)

    print(png)
    print(pdf)
    print(svg)
    print(OUTDIR / "bacillus_mofa_score_data.csv")
    print(OUTDIR / "rhodanobacter_mofa_score_data.csv")
    print(OUTDIR / "mofa_top_loading_summary.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import os
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig_mingfeichen")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


MANUSCRIPT_ROOT = Path("/Users/mingfeichen/Manuscript")
DATA_ROOT = Path("/Users/mingfeichen")
OUTDIR = MANUSCRIPT_ROOT / "outputs/manual-20260609-a2/presentations/transporter_metabolite_link/assets"
OUTDIR.mkdir(parents=True, exist_ok=True)

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
        "axes.titlesize": 11.0,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
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
        "transport_candidates": [
            "map02010",
            "map02060",
            "map02024",
            "map03070",
            "map02040",
            "map02030",
            "map02020",
            "map03060",
            "map01503",
        ],
        "transport_labels": {
            "map02010": ("ABC transporters", "Transport / signaling", 1),
            "map02060": ("PTS system", "Transport / signaling", 2),
            "map02024": ("Quorum sensing", "Transport / signaling", 3),
            "map03070": ("Bacterial secretion system", "Transport / secretion", 4),
            "map02040": ("Flagellar assembly", "Signal transduction / motility", 5),
            "map02030": ("Bacterial chemotaxis", "Signal transduction / motility", 6),
            "map02020": ("Two-component system", "Transport / signaling", 7),
            "map03060": ("Protein export", "Transport / secretion", 8),
            "map01503": ("Bacterial secretion system", "Transport / secretion", 9),
        },
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
        "transport_candidates": [
            "map02020",
            "map03070",
            "map03060",
            "map02040",
            "map02030",
            "map02024",
            "map02010",
            "map01503",
            "map01110",
        ],
        "transport_labels": {
            "map02020": ("Two-component system", "Transport / signaling", 1),
            "map03070": ("Bacterial secretion system", "Transport / secretion", 2),
            "map03060": ("Protein export", "Transport / secretion", 3),
            "map02040": ("Flagellar assembly", "Signal transduction / motility", 4),
            "map02030": ("Bacterial chemotaxis", "Signal transduction / motility", 5),
            "map02024": ("Quorum sensing", "Transport / signaling", 6),
            "map02010": ("ABC transporters", "Transport / signaling", 7),
            "map01503": ("Bacterial secretion system", "Transport / secretion", 8),
            "map01110": ("Biosynthesis of secondary metabolites", "Secondary metabolism / transport", 9),
        },
    },
}


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


def scale_sizes(padj: pd.Series) -> np.ndarray:
    neglog = -np.log10(np.clip(padj.to_numpy(dtype=float), 1e-300, None))
    neglog = np.nan_to_num(neglog, nan=0.0, posinf=10.0, neginf=0.0)
    return 28.0 + np.clip(neglog, 0.0, 6.5) * 28.0


def section_bounds(items: list[dict]) -> tuple[float, float]:
    ys = [item["y"] for item in items]
    return min(ys), max(ys)


def build_transport_table(species: str, config: dict, df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for set_id in config["transport_candidates"]:
        if set_id not in config["transport_labels"]:
            continue
        label, pathway_class, order = config["transport_labels"][set_id]
        sub = df[df["set_id"] == set_id].copy()
        if sub.empty:
            continue
        sub["treatment_display"] = sub["contrast"].map(config["contrast_map"])
        sub = sub[sub["treatment_display"].notna()]
        if sub.empty:
            continue
        any_sig = bool((sub["padj"] < 0.05).any())
        if not any_sig:
            continue
        for treatment in config["treatment_order"]:
            cell = sub[sub["treatment_display"] == treatment]
            if cell.empty:
                rows.append(
                    {
                        "species": species,
                        "section": "Transporter pathways",
                        "subsection": pathway_class,
                        "row_type": "transport",
                        "row_label": label,
                        "row_key": set_id,
                        "row_order": order,
                        "treatment_display": treatment,
                        "effect_value": np.nan,
                        "padj": np.nan,
                        "significant": False,
                    }
                )
                continue
            cell = cell.sort_values(["padj", "NES"], ascending=[True, False]).iloc[0]
            rows.append(
                {
                    "species": species,
                    "section": "Transporter pathways",
                    "subsection": pathway_class,
                    "row_type": "transport",
                    "row_label": label,
                    "row_key": set_id,
                    "row_order": order,
                    "treatment_display": treatment,
                    "effect_value": float(cell["NES"]),
                    "padj": float(cell["padj"]),
                    "significant": bool(cell["padj"] < 0.05),
                }
            )
    return pd.DataFrame(rows)


def select_metabolites(
    df: pd.DataFrame,
    n_universal: int = 4,
    n_phage: int = 5,
    n_leak: int = 5,
) -> pd.DataFrame:
    sig = df[df["pvalue"] < 0.05].copy()
    if sig.empty:
        return pd.DataFrame()
    leak_treatments = {"late", "Al", "K"}
    phage_p = (
        sig[sig["treatment_display"] == "Phage"]
        .groupby("feature_display")["pvalue"]
        .min()
        .rename("phage_p")
    )
    leak_p = (
        sig[sig["treatment_display"].isin(leak_treatments)]
        .groupby("feature_display")["pvalue"]
        .min()
        .rename("leak_p")
    )
    summary = (
        sig.groupby(["feature_display", "superclass"], as_index=False)
        .agg(
            n_treatments=("treatment_display", "nunique"),
            min_p=("pvalue", "min"),
            max_abs=("log2FC", lambda s: float(np.abs(s).max())),
            order=("order", "min"),
        )
        .merge(phage_p.reset_index(), on="feature_display", how="left")
        .merge(leak_p.reset_index(), on="feature_display", how="left")
        .sort_values(["min_p", "feature_display"], ascending=[True, True])
    )
    max_n = int(summary["n_treatments"].max())
    universal = summary[summary["n_treatments"] == max_n].sort_values(["min_p", "feature_display"]).head(n_universal)
    phage = summary[
        summary["phage_p"].notna() & ~summary["feature_display"].isin(universal["feature_display"])
    ]
    phage = phage.sort_values(["phage_p", "min_p", "feature_display"], ascending=[True, True, True]).head(n_phage)
    leak = summary[
        summary["leak_p"].notna()
        & ~summary["feature_display"].isin(universal["feature_display"])
        & ~summary["feature_display"].isin(phage["feature_display"])
    ]
    leak = leak.sort_values(["leak_p", "min_p", "feature_display"], ascending=[True, True, True]).head(n_leak)
    selected = pd.concat([universal, phage, leak], ignore_index=True)
    selected["recurrence_class"] = np.where(
        selected["feature_display"].isin(universal["feature_display"]),
        "Universal core",
        np.where(
            selected["feature_display"].isin(phage["feature_display"]),
            "Phage-enriched",
            "Leak-associated",
        ),
    )
    return selected


def build_metabolite_table(species: str, config: dict, df: pd.DataFrame) -> pd.DataFrame:
    selected = select_metabolites(df)
    if selected.empty:
        return pd.DataFrame()
    selected_labels = selected["feature_display"].tolist()
    selected = selected.copy()
    selected["species"] = species
    selected["section"] = "Metabolites"
    selected["row_type"] = "metabolite"
    selected["row_key"] = selected["feature_display"]
    selected["row_order"] = np.arange(1, len(selected) + 1)

    sig = df[df["pvalue"] < 0.05].copy()
    sig["feature_display"] = sig["feature_display"].map(clean_label)
    sig = sig[sig["feature_display"].isin(selected_labels)]
    sig["recurrence_class"] = sig["feature_display"].map(
        selected.set_index("feature_display")["recurrence_class"].to_dict()
    )
    sig["species"] = species
    sig["section"] = "Metabolites"
    sig["row_type"] = "metabolite"
    sig["row_key"] = sig["feature_display"]
    sig["row_order"] = sig["feature_display"].map(
        selected.set_index("feature_display")["row_order"].to_dict()
    )
    sig["significant"] = True
    sig["treatment_display"] = sig["treatment_display"].map(lambda x: clean_label(x))
    return sig


def layout_sections(rows: list[dict], section_label: str, gap_after: float = 1.2) -> list[dict]:
    y = 0.0
    sectioned: list[dict] = []
    for idx, row in enumerate(rows):
        if idx > 0:
            y -= gap_after
        row_copy = dict(row)
        row_copy["y"] = y
        sectioned.append(row_copy)
        y -= 1.0
    return sectioned


def prepare_species(species: str) -> dict:
    config = SPECIES_CONFIG[species]
    fgsea = read_fgsea(config["fgsea_path"])
    metabolite = read_metabolites(config["met_path"])

    transport_table = build_transport_table(species, config, fgsea)
    met_sig = build_metabolite_table(species, config, metabolite)

    # Convert the selected metabolite summary back to a row-level table with each treatment column.
    met_selected = select_metabolites(metabolite)
    met_rows = []
    if not met_selected.empty:
        for _, item in met_selected.sort_values(["recurrence_class", "n_treatments", "min_p", "feature_display"], ascending=[True, True, True, True]).iterrows():
            feat = item["feature_display"]
            superc = item["superclass"]
            rec_class = item["recurrence_class"]
            order = int(item["order"])
            feat_df = metabolite[metabolite["feature_display"] == feat].copy()
            for treatment in config["treatment_order"]:
                cell = feat_df[feat_df["treatment_display"] == treatment]
                if cell.empty:
                    met_rows.append(
                        {
                            "species": species,
                            "section": "Metabolites",
                            "subsection": rec_class,
                            "superclass": superc,
                            "row_type": "metabolite",
                            "row_label": feat,
                            "row_key": feat,
                            "row_order": order,
                            "treatment_display": treatment,
                            "effect_value": np.nan,
                            "padj": np.nan,
                            "significant": False,
                            "recurrence_class": rec_class,
                        }
                    )
                    continue
                cell = cell.sort_values(["pvalue", "log2FC"], ascending=[True, False]).iloc[0]
                met_rows.append(
                    {
                        "species": species,
                        "section": "Metabolites",
                        "subsection": rec_class,
                        "superclass": superc,
                        "row_type": "metabolite",
                        "row_label": feat,
                        "row_key": feat,
                        "row_order": order,
                        "treatment_display": treatment,
                        "effect_value": float(cell["log2FC"]),
                        "padj": float(cell["pvalue"]),
                        "significant": bool(cell["pvalue"] < 0.05),
                        "recurrence_class": rec_class,
                    }
                )

    meta_row_table = pd.DataFrame(met_rows)
    if not meta_row_table.empty:
        meta_row_table["row_rank"] = meta_row_table["row_key"].map(
            {k: i for i, k in enumerate(met_selected.sort_values(["recurrence_class", "n_treatments", "min_p", "feature_display"], ascending=[True, True, True, True])["feature_display"].tolist(), start=1)}
        )
    else:
        meta_row_table["row_rank"] = []

    transport_rows = []
    if not transport_table.empty:
        transport_rows = transport_table.copy()
        transport_rows["row_rank"] = transport_rows["row_key"].map(
            {k: i for i, k in enumerate(transport_table.drop_duplicates("row_key").sort_values(["row_order", "row_label"])["row_key"].tolist(), start=1)}
        )
    return {
        "transport": transport_table,
        "metabolites": meta_row_table,
        "metabolite_summary": met_selected,
        "transport_selection": transport_rows,
    }


def rows_with_positions(df: pd.DataFrame, section_name: str, gap_after: float = 1.2) -> list[dict]:
    base_rows: list[dict] = []
    seen = []
    for row_key, sub in df.drop_duplicates("row_key").sort_values(["row_order", "row_label"]).groupby("row_key", sort=False):
        first = sub.iloc[0]
        base_rows.append(
            {
                "row_key": row_key,
                "row_label": first["row_label"],
                "subsection": first.get("subsection", ""),
                "row_type": first["row_type"],
                "superclass": first.get("superclass", ""),
                "recurrence_class": first.get("recurrence_class", ""),
            }
        )
        seen.append(row_key)
    positioned = layout_sections(base_rows, section_name, gap_after=gap_after)
    return positioned


def plot_species(ax, species: str, spec_data: dict) -> dict:
    config = SPECIES_CONFIG[species]
    transport = spec_data["transport"].copy()
    metabolites = spec_data["metabolites"].copy()

    transport_base = (
        transport.drop_duplicates("row_key")
        .sort_values(["row_order", "row_label"])
        .loc[:, ["row_key", "row_label", "subsection", "row_type"]]
        .to_dict("records")
    )
    metabolite_base = (
        metabolites.drop_duplicates("row_key")
        .sort_values(["recurrence_class", "row_order", "row_label"])
        .loc[:, ["row_key", "row_label", "subsection", "row_type", "superclass", "recurrence_class"]]
        .to_dict("records")
    )

    transport_rows = layout_sections(transport_base, "Transporter pathways", gap_after=0.35)
    transport_label_map = {row["row_key"]: row["y"] for row in transport_rows}
    transport["y"] = transport["row_key"].map(transport_label_map)

    metabolite_groups = []
    for rec_class in ["Universal core", "Phage-enriched", "Leak-associated"]:
        subset = [row for row in metabolite_base if row.get("recurrence_class", "") == rec_class]
        if subset:
            positioned = layout_sections(subset, rec_class, gap_after=0.35)
            metabolite_groups.extend(positioned)
    if not metabolite_groups:
        metabolite_groups = layout_sections(metabolite_base, "Metabolites", gap_after=0.35)

    if transport_rows and metabolite_groups:
        transport_min = min(row["y"] for row in transport_rows)
        metabolite_offset = transport_min - 1.25
        adjusted_groups = []
        for rec_class in ["Universal core", "Phage-enriched", "Leak-associated"]:
            class_rows = [row for row in metabolite_groups if row.get("recurrence_class") == rec_class]
            if not class_rows:
                continue
            current = []
            for row in class_rows:
                row_copy = dict(row)
                row_copy["y"] = row_copy["y"] + metabolite_offset
                current.append(row_copy)
            adjusted_groups.extend(current)
            metabolite_offset = min(row["y"] for row in current) - 0.95
        if adjusted_groups:
            metabolite_groups = adjusted_groups

    metabolite_label_map = {row["row_key"]: row["y"] for row in metabolite_groups}
    metabolites["y"] = metabolites["row_key"].map(metabolite_label_map)

    combined = pd.concat([transport, metabolites], ignore_index=True, sort=False)
    combined["neglog10padj"] = -np.log10(np.clip(combined["padj"].astype(float), 1e-300, None))
    combined["row_label"] = combined["row_label"].map(clean_label)
    combined["section"] = combined["section"].fillna("")

    transport_vals = combined[(combined["row_type"] == "transport") & combined["significant"]]["effect_value"].dropna()
    met_vals = combined[(combined["row_type"] == "metabolite") & combined["significant"]]["effect_value"].dropna()
    if len(transport_vals):
        tmax = float(np.nanmax(np.abs(transport_vals)))
    else:
        tmax = 1.0
    if len(met_vals):
        mmax = float(np.nanmax(np.abs(met_vals)))
    else:
        mmax = 1.0
    transport_norm = TwoSlopeNorm(vcenter=0.0, vmin=-tmax, vmax=tmax)
    metabolite_norm = TwoSlopeNorm(vcenter=0.0, vmin=-mmax, vmax=mmax)

    x = np.arange(len(config["treatment_order"]))

    # Transport section.
    transport_plot = combined[combined["row_type"] == "transport"].copy()
    if not transport_plot.empty:
        for _, row in transport_plot.iterrows():
            tx = config["treatment_order"].index(row["treatment_display"])
            if bool(row["significant"]):
                ax.scatter(
                    tx,
                    row["y"],
                    s=float(scale_sizes(pd.Series([row["padj"]]))[0]),
                    c=[plt.cm.RdBu_r(transport_norm(row["effect_value"]))],
                    edgecolors="#3F3F46",
                    linewidths=0.4,
                    zorder=3,
                )
            else:
                ax.scatter(
                    tx,
                    row["y"],
                    s=18,
                    facecolors="none",
                    edgecolors="#CBD5E1",
                    linewidths=0.45,
                    zorder=2,
                )

    # Metabolite section.
    metabolite_plot = combined[combined["row_type"] == "metabolite"].copy()
    if not metabolite_plot.empty:
        for _, row in metabolite_plot.iterrows():
            tx = config["treatment_order"].index(row["treatment_display"])
            if bool(row["significant"]):
                ax.scatter(
                    tx,
                    row["y"],
                    s=float(scale_sizes(pd.Series([row["padj"]]))[0]),
                    c=[plt.cm.RdBu_r(metabolite_norm(row["effect_value"]))],
                    edgecolors="#3F3F46",
                    linewidths=0.4,
                    zorder=3,
                )
            else:
                ax.scatter(
                    tx,
                    row["y"],
                    s=18,
                    facecolors="none",
                    edgecolors="#CBD5E1",
                    linewidths=0.45,
                    zorder=2,
                )

    # Background bands and section separators.
    all_transport = transport_rows
    all_metabolites = metabolite_groups
    if all_transport:
        ymin_t, ymax_t = section_bounds(all_transport)
        ax.axhspan(ymin_t - 0.45, ymax_t + 0.45, color="#F8FAFC", zorder=0)
        ax.hlines(ymin_t - 0.95, -0.45, len(x) - 0.55, color="#D1D5DB", linewidth=0.6, zorder=1)
    if all_metabolites:
        ymin_m, ymax_m = section_bounds(all_metabolites)
        ax.axhspan(ymin_m - 0.45, ymax_m + 0.45, color="#FCFCFD", zorder=0)
        ax.hlines(ymin_m - 0.95, -0.45, len(x) - 0.55, color="#D1D5DB", linewidth=0.6, zorder=1)

    # Row labels.
    yticks = [row["y"] for row in all_transport] + [row["y"] for row in all_metabolites]
    ylabels = [clean_label(row["row_label"]) for row in all_transport] + [clean_label(row["row_label"]) for row in all_metabolites]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)

    # Group annotations on the right.
    def annotate_group(rows: list[dict], text_map: dict[str, str], x_text: float):
        grouped = {}
        for row in rows:
            grouped.setdefault(row["subsection"], []).append(row["y"])
        for key, ys in grouped.items():
            if not ys:
                continue
            label = text_map.get(key, key)
            y_mid = (min(ys) + max(ys)) / 2
            ax.text(
                x_text,
                y_mid,
                label,
                ha="left",
                va="center",
                fontsize=6.9,
                color="#374151",
                clip_on=False,
            )
            ax.plot([x_text - 0.10, x_text - 0.10], [min(ys) - 0.35, max(ys) + 0.35], color="#CBD5E1", lw=0.7, clip_on=False)

    annotate_group(
        all_transport,
        {
            "Transport / signaling": "Transport/signaling",
            "Transport / secretion": "Transport/secretion",
            "Signal transduction / motility": "Motility",
            "Secondary metabolism / transport": "Secondary metab./transport",
        },
        4.65,
    )
    if all_metabolites:
        universal_rows = [row for row in all_metabolites if row.get("recurrence_class") == "Universal core"]
        phage_rows = [row for row in all_metabolites if row.get("recurrence_class") == "Phage-enriched"]
        leak_rows = [row for row in all_metabolites if row.get("recurrence_class") == "Leak-associated"]
        if universal_rows:
            ymin_u, ymax_u = section_bounds(universal_rows)
            ax.axhspan(ymin_u - 0.40, ymax_u + 0.40, color="#F7FAFF", zorder=0)
        if phage_rows:
            ymin_r, ymax_r = section_bounds(phage_rows)
            ax.axhspan(ymin_r - 0.40, ymax_r + 0.40, color="#FFF8F2", zorder=0)
        if leak_rows:
            ymin_l, ymax_l = section_bounds(leak_rows)
            ax.axhspan(ymin_l - 0.40, ymax_l + 0.40, color="#F4FAF5", zorder=0)

    if any(row.get("recurrence_class") == "Universal core" for row in all_metabolites):
        ys = [row["y"] for row in all_metabolites if row.get("recurrence_class") == "Universal core"]
        ax.text(
            -0.82,
            (min(ys) + max(ys)) / 2,
            "Universal core",
            ha="left",
            va="center",
            fontsize=7.7,
            color="#475569",
            fontstyle="italic",
            clip_on=False,
        )
    if any(row.get("recurrence_class") == "Phage-enriched" for row in all_metabolites):
        ys = [row["y"] for row in all_metabolites if row.get("recurrence_class") == "Phage-enriched"]
        ax.text(
            -0.82,
            (min(ys) + max(ys)) / 2,
            "Phage-enriched",
            ha="left",
            va="center",
            fontsize=7.7,
            color="#475569",
            fontstyle="italic",
            clip_on=False,
        )
    if any(row.get("recurrence_class") == "Leak-associated" for row in all_metabolites):
        ys = [row["y"] for row in all_metabolites if row.get("recurrence_class") == "Leak-associated"]
        ax.text(
            -0.82,
            (min(ys) + max(ys)) / 2,
            "Leak-associated",
            ha="left",
            va="center",
            fontsize=7.7,
            color="#475569",
            fontstyle="italic",
            clip_on=False,
        )

    ax.set_xlim(-0.6, len(x) - 0.05)
    ax.set_ylim((min(yticks) - 1.0) if yticks else -1, (max(yticks) + 1.0) if yticks else 1)
    ax.set_xticks(x)
    ax.set_xticklabels(config["treatment_labels"], rotation=0)
    ax.tick_params(axis="x", length=0)
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.45)
    ax.grid(axis="y", color="#F3F4F6", linewidth=0.35)
    ax.set_axisbelow(True)
    ax.set_title(species, loc="left", fontsize=11.0, fontweight="bold", pad=7)
    ax.set_xlabel("Treatment group")
    ax.set_ylabel("")
    return {
        "transport_norm": transport_norm,
        "metabolite_norm": metabolite_norm,
        "transport_max": tmax,
        "metabolite_max": mmax,
        "transport_rows": transport_plot,
        "metabolite_rows": metabolite_plot,
        "all_rows": combined,
    }


def make_colorbar(fig, cax, norm, title, cmap=plt.cm.RdBu_r):
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.55)
    cb.ax.tick_params(labelsize=8.0, length=2.0, colors="#374151")
    cb.set_label(title, fontsize=8.6, color="#111827")
    return cb


def main():
    species_data = {species: prepare_species(species) for species in SPECIES_CONFIG}

    fig = plt.figure(figsize=(12.2, 9.0), dpi=600)
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.0, 1.0, 0.28], hspace=0.24)
    ax_b = fig.add_subplot(gs[0, 0])
    ax_r = fig.add_subplot(gs[1, 0])
    leg_gs = gs[2, 0].subgridspec(1, 3, width_ratios=[1.05, 1.05, 0.92], wspace=0.55)
    cax_transport = fig.add_subplot(leg_gs[0, 0])
    cax_metabolite = fig.add_subplot(leg_gs[0, 1])
    ax_legend = fig.add_subplot(leg_gs[0, 2])

    b_plot = plot_species(ax_b, "Bacillus", species_data["Bacillus"])
    r_plot = plot_species(ax_r, "Rhodanobacter", species_data["Rhodanobacter"])

    make_colorbar(
        fig,
        cax_transport,
        b_plot["transport_norm"],
        "Transporter pathways: NES",
    )
    make_colorbar(
        fig,
        cax_metabolite,
        r_plot["metabolite_norm"],
        "Metabolites: log2FC",
    )

    ax_legend.axis("off")
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#CBD5E1", markersize=7.5, label="padj >= 0.05"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#111827", markeredgecolor="#111827", markersize=4.2, label="padj < 0.05"),
    ]
    ax_legend.legend(
        handles=legend_handles,
        loc="upper left",
        frameon=False,
        handletextpad=0.6,
        borderaxespad=0.0,
        title="-log10(padj)",
        title_fontsize=8.4,
    )
    ax_legend.text(
        0.0,
        0.05,
        "Shared treatment axis links transporter FGSEA and metabolite t-tests.\nK = Kanamycin, P = Phage.\nHollow circles indicate non-significant cells.",
        ha="left",
        va="bottom",
        fontsize=7.4,
        color="#475569",
        transform=ax_legend.transAxes,
    )

    fig.subplots_adjust(left=0.18, right=0.965, top=0.98, bottom=0.06)

    png_path = OUTDIR / "reproduce_panel.png"
    pdf_path = OUTDIR / "reproduce_panel.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    # Persist a transparent audit table for the selected rows.
    audit_rows = []
    for species, bundle in species_data.items():
        for _, row in bundle["transport"].iterrows():
            audit_rows.append(
                {
                    "species": species,
                    "section": "Transporter pathways",
                    "row_label": row["row_label"],
                    "subsection": row["subsection"],
                    "treatment_display": row["treatment_display"],
                    "effect_value": row["effect_value"],
                    "padj": row["padj"],
                    "significant": row["significant"],
                    "kind": "transport",
                }
            )
        for _, row in bundle["metabolites"].iterrows():
            audit_rows.append(
                {
                    "species": species,
                    "section": "Metabolites",
                    "row_label": row["row_label"],
                    "subsection": row["subsection"],
                    "treatment_display": row["treatment_display"],
                    "effect_value": row["effect_value"],
                    "padj": row["padj"],
                    "significant": row["significant"],
                    "kind": "metabolite",
                }
            )
    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(OUTDIR / "reproduce_panel_selected_data.csv", index=False)

    summary = {
        "species": {
            species: {
                "transport_rows": int(bundle["transport"]["row_key"].nunique()) if not bundle["transport"].empty else 0,
                "metabolite_rows": int(bundle["metabolites"]["row_key"].nunique()) if not bundle["metabolites"].empty else 0,
                "transport_significant_cells": int(bundle["transport"]["significant"].sum()) if not bundle["transport"].empty else 0,
                "metabolite_significant_cells": int(bundle["metabolites"]["significant"].sum()) if not bundle["metabolites"].empty else 0,
            }
            for species, bundle in species_data.items()
        },
        "outputs": {
            "png": str(png_path),
            "pdf": str(pdf_path),
            "selected_data": str(OUTDIR / "reproduce_panel_selected_data.csv"),
        },
    }
    (OUTDIR / "reproduce_panel_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import sqlite3
import textwrap
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-mingfeichen")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec


ROOT = Path("/Users/mingfeichen")
OUTDIR = Path(
    "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a1/presentations/fgsea_panel_a_reproduction/assets"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

BACILLUS_PATH = ROOT / "Bacillus_Merged_FGSEA_Results.csv"
RHODANO_PATH = ROOT / "Rhodanobacter_Merged_FGSEA_Results.csv"
ANNOTATION_PATH = Path("/Users/mingfeichen/Downloads/kegg_pathways_rhodanobacter_bacillus.csv")

KEGG_MODULE_DB = Path(
    "/Users/mingfeichen/mambaforge/envs/anvio-8/lib/python3.10/site-packages/anvio/data/misc/KEGG/MODULES.db"
)

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

SPECIES_SPECS = {
    "Bacillus": {
        "path": BACILLUS_PATH,
        "row_order": ["8h_vs_0h", "24h_vs_0h", "Al_vs_0h", "Kanamycin_vs_0h", "P_vs_0h"],
        "xlabels": ["mid", "late", "Al", "K", "P"],
        "xcolors": ["#DCCF7A", "#A9D9C6", "#B9C5E5", "#E7B9B9", "#E6C9A6"],
    },
    "Rhodanobacter": {
        "path": RHODANO_PATH,
        "row_order": ["24h_vs_0h", "48h_vs_0h", "Al_vs_0h", "Kanamycin_vs_0h", "P_vs_0h"],
        "xlabels": ["mid", "late", "Al", "K", "P"],
        "xcolors": ["#A9D9C6", "#D8C7E8", "#E7B9B9", "#E6C9A6", "#E3D48B"],
    },
}

PATHWAY_LOOKUP = {
    "map00020": "Citrate cycle",
    "map00051": "Fructose and mannose metabolism",
    "map00190": "Oxidative phosphorylation",
    "map00230": "Purine metabolism",
    "map00240": "Pyrimidine metabolism",
    "map00250": "Alanine, aspartate and glutamate metabolism",
    "map00270": "Cysteine and methionine metabolism",
    "map00290": "Valine, leucine and isoleucine biosynthesis",
    "map00340": "Histidine metabolism",
    "map00500": "Starch and sucrose metabolism",
    "map00520": "Amino sugar and nucleotide sugar metabolism",
    "map00550": "Peptidoglycan biosynthesis",
    "map00620": "Pyruvate metabolism",
    "map00630": "Glyoxylate and dicarboxylate metabolism",
    "map00650": "Butanoate metabolism",
    "map00660": "C5-branched dibasic acid metabolism",
    "map00720": "Carbon fixation pathways",
    "map00740": "Riboflavin metabolism",
    "map00780": "Biotin metabolism",
    "map00790": "Folate biosynthesis",
    "map00920": "Sulfur metabolism",
    "map00970": "Aminoacyl-tRNA biosynthesis",
    "map01053": "Biosynthesis of siderophore group nonribosomal peptides",
    "map01100": "Metabolic pathways",
    "map01110": "Biosynthesis of secondary metabolites",
    "map01120": "Microbial metabolism in diverse environments",
    "map01130": "Biosynthesis of secondary metabolites",
    "map01200": "Carbon metabolism",
    "map01210": "2-Oxocarboxylic acid metabolism",
    "map01212": "Fatty acid metabolism",
    "map01230": "Biosynthesis of amino acids",
    "map01240": "Biosynthesis of cofactors",
    "map01501": "Beta-lactam resistance",
    "map01502": "Vancomycin resistance",
    "map01503": "Bacterial secretion system",
    "map02010": "ABC transporters",
    "map02020": "Two-component system",
    "map02024": "Quorum sensing",
    "map02025": "Biofilm formation",
    "map02030": "Bacterial chemotaxis",
    "map02040": "Flagellar assembly",
    "map02060": "Phosphotransferase system",
    "map03010": "Ribosome",
    "map03060": "Protein export",
    "map03070": "Bacterial secretion system",
    "map03420": "Nucleotide excision repair",
    "map04112": "Cell cycle",
}

GROUP_ORDER = [
    "Signal transduction / motility",
    "Transport / secretion",
    "Translation / ribosome",
    "Nucleotide metabolism",
    "Amino acid / carbon metabolism",
    "Cofactor / vitamin metabolism",
    "Cell envelope / resistance",
    "Energy / respiration",
    "Other",
]

GROUP_COLORS = {
    "Signal transduction / motility": "#D94E4E",
    "Transport / secretion": "#E39F3A",
    "Translation / ribosome": "#6C8EC7",
    "Nucleotide metabolism": "#8E6BBE",
    "Amino acid / carbon metabolism": "#48A36C",
    "Cofactor / vitamin metabolism": "#D07BB2",
    "Cell envelope / resistance": "#C77D5A",
    "Energy / respiration": "#5E9CBE",
    "Other": "#8D8D8D",
}

TOP_LEVEL_ORDER = [
    "Metabolism",
    "Genetic Information Processing",
    "Environmental Information Processing",
    "Cellular Processes",
    "Organismal Systems",
    "Human Diseases",
]

COMBINED_CONDITIONS = [
    {
        "label": "Bacillus\n8h vs 0h",
        "species": ["Bacillus"],
        "contrast": "8h_vs_0h",
        "color": "#DCCF7A",
        "offsets": {"Bacillus": 0.0},
        "markers": {"Bacillus": "o"},
    },
    {
        "label": "Bacillus\n24h vs 0h",
        "species": ["Bacillus"],
        "contrast": "24h_vs_0h",
        "color": "#A9D9C6",
        "offsets": {"Bacillus": 0.0},
        "markers": {"Bacillus": "o"},
    },
    {
        "label": "Bacillus\nAl vs 0h",
        "species": ["Bacillus"],
        "contrast": "Al_vs_0h",
        "color": "#E7B9B9",
        "offsets": {"Bacillus": 0.0},
        "markers": {"Bacillus": "o"},
    },
    {
        "label": "Bacillus\nKana vs 0h",
        "species": ["Bacillus"],
        "contrast": "Kanamycin_vs_0h",
        "color": "#E6C9A6",
        "offsets": {"Bacillus": 0.0},
        "markers": {"Bacillus": "o"},
    },
    {
        "label": "Bacillus\nPhage vs 0h",
        "species": ["Bacillus"],
        "contrast": "P_vs_0h",
        "color": "#E3D48B",
        "offsets": {"Bacillus": 0.0},
        "markers": {"Bacillus": "o"},
    },
    {
        "label": "Rhodanobacter\n24h vs 0h",
        "species": ["Rhodanobacter"],
        "contrast": "24h_vs_0h",
        "color": "#CDBFE9",
        "offsets": {"Rhodanobacter": 0.0},
        "markers": {"Rhodanobacter": "D"},
    },
    {
        "label": "Rhodanobacter\n48h vs 0h",
        "species": ["Rhodanobacter"],
        "contrast": "48h_vs_0h",
        "color": "#D8C7E8",
        "offsets": {"Rhodanobacter": 0.0},
        "markers": {"Rhodanobacter": "D"},
    },
    {
        "label": "Rhodanobacter\nAl vs 0h",
        "species": ["Rhodanobacter"],
        "contrast": "Al_vs_0h",
        "color": "#E7B9B9",
        "offsets": {"Rhodanobacter": 0.0},
        "markers": {"Rhodanobacter": "D"},
    },
    {
        "label": "Rhodanobacter\nKana vs 0h",
        "species": ["Rhodanobacter"],
        "contrast": "Kanamycin_vs_0h",
        "color": "#E6C9A6",
        "offsets": {"Rhodanobacter": 0.0},
        "markers": {"Rhodanobacter": "D"},
    },
    {
        "label": "Rhodanobacter\nPhage vs 0h",
        "species": ["Rhodanobacter"],
        "contrast": "P_vs_0h",
        "color": "#E3D48B",
        "offsets": {"Rhodanobacter": 0.0},
        "markers": {"Rhodanobacter": "D"},
    },
]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def load_annotations(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    annotations: dict[str, dict[str, object]] = {}
    for row in read_csv_rows(path):
        set_id = f"map{row['pathway_id']}"
        key = f"{row['genus']}|{set_id}"
        annotations[key] = {
            "set_id": set_id,
            "label": row["pathway_name"],
            "pathways": row["pathway_name"],
            "top_level_category": row["top_level_category"],
            "subcategory": row["subcategory"],
            "comparison_status": row["comparison_status"],
            "genus": row["genus"],
            "organism": row["organism"],
            "organism_pathway_id": row["organism_pathway_id"],
            "annotation_order": len(annotations),
        }
    return annotations


def load_module_lookup(db_path: Path) -> dict[str, dict[str, str]]:
    if not db_path.exists():
        return {}
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute("SELECT DISTINCT module FROM modules")
    modules = [row[0] for row in cur.fetchall()]
    out: dict[str, dict[str, str]] = {}
    for mod in modules:
        cur.execute("SELECT data_name, data_value FROM modules WHERE module=?", (mod,))
        name = None
        pathways: list[str] = []
        for data_name, data_value in cur.fetchall():
            if data_name == "NAME" and name is None:
                name = data_value
            elif data_name == "PATHWAY" and data_value:
                pathways.append(data_value)
        if name or pathways:
            out[mod] = {"name": name or mod, "pathways": ",".join(pathways)}
    con.close()
    return out


def pretty_label(text: str) -> str:
    text = (text or "").strip()
    text = text.replace("_", " ")
    text = " ".join(text.split())
    return text


def classify_group(label: str, pathways: str = "") -> str:
    txt = f"{label} {pathways}".lower()
    if any(k in txt for k in ["flagellar", "chemotaxis", "two-component", "biofilm", "quorum"]):
        return "Signal transduction / motility"
    if any(k in txt for k in ["transport", "efflux", "secretion", "abc", "export", "pts", "permease"]):
        return "Transport / secretion"
    if any(k in txt for k in ["ribosome", "translation", "ribosomal"]):
        return "Translation / ribosome"
    if any(k in txt for k in ["purine", "pyrimidine", "nucleotide", "dna repair"]):
        return "Nucleotide metabolism"
    if any(
        k in txt
        for k in [
            "amino acid",
            "histidine",
            "cysteine",
            "methionine",
            "valine",
            "leucine",
            "isoleucine",
            "alanine",
            "aspartate",
            "glutamate",
            "pyruvate",
            "carbon metabolism",
            "citrate cycle",
            "glyoxylate",
            "starch",
            "sucrose",
            "butanoate",
            "fatty acid",
        ]
    ):
        return "Amino acid / carbon metabolism"
    if any(k in txt for k in ["biotin", "riboflavin", "folate", "cofactor", "vitamin", "sulfur"]):
        return "Cofactor / vitamin metabolism"
    if any(k in txt for k in ["peptidoglycan", "vancomycin", "beta-lactam", "bacitracin", "cell cycle"]):
        return "Cell envelope / resistance"
    if any(k in txt for k in ["respiration", "oxidative phosphorylation", "nadh", "tca", "acetyl-coa"]):
        return "Energy / respiration"
    return "Other"


def label_for_set_id(set_id: str, module_lookup: dict[str, dict[str, str]]) -> tuple[str, str]:
    if set_id.startswith("map"):
        return PATHWAY_LOOKUP.get(set_id, set_id), PATHWAY_LOOKUP.get(set_id, set_id)
    info = module_lookup.get(set_id)
    if info:
        return info["name"], info["pathways"]
    return set_id, ""


def summarize_species(
    rows: list[dict[str, str]],
    annotations: dict[str, dict[str, object]],
    species: str,
) -> list[dict[str, object]]:
    by_id: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        set_id = row["set_id"]
        if not set_id.startswith("map"):
            continue
        padj = float(row["padj"])
        nes = float(row["NES"])
        if not math.isfinite(padj) or not math.isfinite(nes):
            continue
        by_id[set_id].append({"contrast": row["contrast"], "padj": padj, "nes": nes})

    ordered_ann = [
        ann
        for set_id, ann in sorted(
            ((sid, ann) for sid, ann in annotations.items() if ann["genus"] == species),
            key=lambda item: (
                TOP_LEVEL_ORDER.index(item[1]["top_level_category"]) if item[1]["top_level_category"] in TOP_LEVEL_ORDER else len(TOP_LEVEL_ORDER),
                item[1]["subcategory"],
                item[1]["annotation_order"],
                item[0],
            ),
        )
    ]

    summaries: list[dict[str, object]] = []
    seen: set[str] = set()
    for ann in ordered_ann:
        set_id = ann["set_id"]
        recs = by_id.get(set_id, [])
        if not recs:
            continue
        sig = [r for r in recs if r["padj"] < 0.05]
        best = min(recs, key=lambda r: r["padj"])
        group = f"{ann['top_level_category']} / {ann['subcategory']}"
        summaries.append(
            {
                "species": species,
                "set_id": set_id,
                "label": pretty_label(ann["label"]),
                "pathways": pretty_label(ann["pathways"]),
                "group": pretty_label(group),
                "top_level_category": ann["top_level_category"],
                "subcategory": ann["subcategory"],
                "comparison_status": ann["comparison_status"],
                "n_sig": len(sig),
                "best_padj": best["padj"],
                "best_nes": best["nes"],
                "max_abs_nes": max(abs(r["nes"]) for r in recs),
                "signed_sum": sum(r["nes"] for r in recs),
                "annotation_order": ann["annotation_order"],
            }
        )
        seen.add(set_id)

    # Keep any significant pathway hits that are not present in the curated annotation CSV.
    for set_id in sorted(by_id):
        if set_id in seen:
            continue
        recs = by_id[set_id]
        if not recs:
            continue
        sig = [r for r in recs if r["padj"] < 0.05]
        best = min(recs, key=lambda r: r["padj"])
        label, pathways = label_for_set_id(set_id, {})
        summaries.append(
            {
                "species": species,
                "set_id": set_id,
                "label": pretty_label(label),
                "pathways": pretty_label(pathways),
                "group": "Unannotated / other",
                "top_level_category": "Unannotated",
                "subcategory": "other",
                "comparison_status": "unannotated",
                "n_sig": len(sig),
                "best_padj": best["padj"],
                "best_nes": best["nes"],
                "max_abs_nes": max(abs(r["nes"]) for r in recs),
                "signed_sum": sum(r["nes"] for r in recs),
                "annotation_order": 10**9,
            }
        )

    for item in summaries:
        if item["signed_sum"] > 0:
            item["direction"] = "up"
        elif item["signed_sum"] < 0:
            item["direction"] = "down"
        else:
            item["direction"] = "up" if item["best_nes"] >= 0 else "down"

    def rank_key(item: dict[str, object]):
        return (
            TOP_LEVEL_ORDER.index(item["top_level_category"]) if item["top_level_category"] in TOP_LEVEL_ORDER else len(TOP_LEVEL_ORDER),
            item["subcategory"],
            item["annotation_order"],
            0 if item["direction"] == "up" else 1,
            -int(item["n_sig"]),
            -float(item["max_abs_nes"]),
            float(item["best_padj"]),
            item["set_id"],
        )

    ordered = sorted(summaries, key=rank_key)
    return ordered


def build_matrix(selected_summaries: list[dict[str, object]], row_order: list[str]):
    out = []
    for item in selected_summaries:
        row = {
            "set_id": item["set_id"],
            "label": item["label"],
            "group": item["group"],
            "top_level_category": item["top_level_category"],
            "subcategory": item["subcategory"],
            "comparison_status": item["comparison_status"],
            "pathways": item["pathways"],
            "values": {contrast: None for contrast in row_order},
        }
        out.append(row)
    return out


def estimate_bubble_size(padj: float, max_size: float = 11.0) -> float:
    value = max(0.0, min(-math.log10(padj), max_size))
    return 16.0 + 12.0 * value ** 1.18


def build_combined_rows(
    annotations: dict[str, dict[str, object]],
    raw_rows_by_species: dict[str, list[dict[str, str]]],
):
    fgsea_lookup: dict[tuple[str, str, str], dict[str, float]] = {}
    for species, rows in raw_rows_by_species.items():
        for row in rows:
            set_id = row["set_id"]
            if not set_id.startswith("map"):
                continue
            try:
                padj = float(row["padj"])
                nes = float(row["NES"])
            except ValueError:
                continue
            if not math.isfinite(padj) or not math.isfinite(nes):
                continue
            fgsea_lookup[(species, set_id, row["contrast"])] = {"padj": padj, "nes": nes}

    row_meta: dict[tuple[str, str], dict[str, object]] = {}
    for ann in annotations.values():
        set_id = ann["set_id"]
        species = ann["genus"]
        item = row_meta.setdefault(
            (species, set_id),
            {
                "set_id": set_id,
                "species": species,
                "labels": set(),
                "top_level_categories": set(),
                "subcategories": set(),
                "comparison_statuses": set(),
                "annotation_order": ann["annotation_order"],
            },
        )
        item["labels"].add(ann["label"])
        item["top_level_categories"].add(ann["top_level_category"])
        item["subcategories"].add(ann["subcategory"])
        item["comparison_statuses"].add(ann["comparison_status"])
        item["annotation_order"] = min(int(item["annotation_order"]), int(ann["annotation_order"]))

    for species, rows in raw_rows_by_species.items():
        for row in rows:
            set_id = row["set_id"]
            if not set_id.startswith("map"):
                continue
            if (species, set_id) in row_meta:
                continue
            try:
                padj = float(row["padj"])
                nes = float(row["NES"])
            except ValueError:
                continue
            if not math.isfinite(padj) or not math.isfinite(nes):
                continue
            row_meta[(species, set_id)] = {
                "set_id": set_id,
                "species": species,
                "labels": {PATHWAY_LOOKUP.get(set_id, set_id)},
                "top_level_categories": {"Unannotated"},
                "subcategories": {"other"},
                "comparison_statuses": {"unannotated"},
                "annotation_order": 10**9,
            }

    rows: list[dict[str, object]] = []
    for (_, set_id), meta in row_meta.items():
        labels = sorted(meta["labels"])
        top_levels = sorted(meta["top_level_categories"], key=lambda x: TOP_LEVEL_ORDER.index(x) if x in TOP_LEVEL_ORDER else len(TOP_LEVEL_ORDER))
        subcategories = sorted(meta["subcategories"])
        comparison_statuses = sorted(meta["comparison_statuses"])
        if comparison_statuses:
            comparison_status = comparison_statuses[0]
        else:
            comparison_status = "unannotated"

        rows.append(
            {
                "set_id": set_id,
                "species": meta["species"],
                "label": labels[0],
                "top_level_category": top_levels[0] if top_levels else "Unannotated",
                "subcategory": subcategories[0] if subcategories else "other",
                "comparison_status": comparison_status,
                "annotation_order": int(meta["annotation_order"]),
                "cells": {idx: [] for idx in range(len(COMBINED_CONDITIONS))},
            }
        )

    rows.sort(
        key=lambda item: (
            0 if item["species"] == "Bacillus" else 1,
            TOP_LEVEL_ORDER.index(item["top_level_category"]) if item["top_level_category"] in TOP_LEVEL_ORDER else len(TOP_LEVEL_ORDER),
            item["subcategory"],
            item["annotation_order"],
            item["label"],
            item["set_id"],
        )
    )

    for row in rows:
        for idx, cond in enumerate(COMBINED_CONDITIONS):
            if row["species"] not in cond["species"]:
                continue
            for species in cond["species"]:
                rec = fgsea_lookup.get((species, row["set_id"], cond["contrast"]))
                if rec is None:
                    continue
                row["cells"][idx].append(
                    {
                        "species": species,
                        "contrast": cond["contrast"],
                        "padj": rec["padj"],
                        "nes": rec["nes"],
                        "offset": cond["offsets"][species],
                        "marker": cond["markers"][species],
                    }
                )

    return rows, None


def draw_combined_panel(ax, ann_ax, rows: list[dict[str, object]]):
    n_rows = len(rows)
    n_cols = len(COMBINED_CONDITIONS)
    cmap = plt.get_cmap("RdBu_r")
    norm = TwoSlopeNorm(vmin=-3.6, vcenter=0.0, vmax=3.6)

    for idx, cond in enumerate(COMBINED_CONDITIONS):
        ax.axvspan(idx - 0.5, idx + 0.5, color=cond["color"], alpha=0.18, zorder=0)
        ax.axvline(idx - 0.5, color="#E5E7EB", lw=0.8, zorder=1)
    ax.axvline(n_cols - 0.5, color="#E5E7EB", lw=0.8, zorder=1)

    for y, row in enumerate(rows):
        for x, cell_items in row["cells"].items():
            for item in cell_items:
                ax.scatter(
                    x + item["offset"],
                    y,
                    s=estimate_bubble_size(item["padj"]),
                    facecolor=cmap(norm(item["nes"])),
                    edgecolor="#4B5563" if item["species"] == "Bacillus" else "#1F2937",
                    linewidth=0.55,
                    marker=item["marker"],
                    alpha=0.96,
                    zorder=3,
                )
        if y > 0 and rows[y - 1]["species"] != row["species"]:
            ax.axhline(y - 0.5, color="#9CA3AF", lw=1.2, zorder=2)

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([cond["label"] for cond in COMBINED_CONDITIONS], rotation=45, ha="right", fontsize=9)
    for tick, cond in zip(ax.get_xticklabels(), COMBINED_CONDITIONS):
        tick.set_color(cond["color"])
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(
        [f"{'Bac' if row['species'] == 'Bacillus' else 'Rho'} {row['set_id']}" for row in rows],
        fontsize=7.2,
    )
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", color="#F3F4F6", linewidth=0.8, zorder=0)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#D1D5DB")
    ax.spines["bottom"].set_color("#D1D5DB")
    ax.set_ylabel("KEGG pathway / module", fontsize=11)

    ann_ax.set_xlim(0, 1)
    ann_ax.set_ylim(n_rows - 0.5, -0.5)
    ann_ax.axis("off")
    groups = []
    current = None
    for idx, row in enumerate(rows):
        group = row["top_level_category"]
        if current is None or group != current["group"]:
            if current is not None:
                current["end"] = idx - 1
                groups.append(current)
            current = {"group": group, "start": idx}
    if current is not None:
        current["end"] = n_rows - 1
        groups.append(current)

    for g in groups:
        start = g["start"]
        end = g["end"]
        mid = (start + end) / 2
        ann_ax.plot([0.03, 0.03], [start, end], color="#7C7C7C", lw=1.0, clip_on=False)
        ann_ax.plot([0.03, 0.085], [start, start], color="#7C7C7C", lw=1.0, clip_on=False)
        ann_ax.plot([0.03, 0.085], [end, end], color="#7C7C7C", lw=1.0, clip_on=False)
        ann_ax.text(
            0.12,
            mid,
            textwrap.fill(g["group"], width=28),
            ha="left",
            va="center",
            fontsize=8.4,
            color="#404040",
        )
        if g["start"] > 0 and rows[g["start"] - 1]["species"] != rows[g["start"]]["species"]:
            ann_ax.plot([0.03, 0.98], [g["start"] - 0.5, g["start"] - 0.5], color="#9CA3AF", lw=1.1, clip_on=False)

    return rows, None


def draw_panel(
    ax,
    ann_ax,
    species: str,
    selected: list[dict[str, object]],
    raw_rows: list[dict[str, str]],
    module_lookup: dict[str, dict[str, str]],
):
    spec = SPECIES_SPECS[species]
    row_order = spec["row_order"]
    xlabels = spec["xlabels"]
    xcolors = spec["xcolors"]

    matrix = build_matrix(selected, row_order)
    selected_ids = {item["set_id"] for item in selected}
    raw_by_set: dict[str, list[dict[str, str]]] = defaultdict(list)
    for rec in raw_rows:
        if rec["set_id"] in selected_ids:
            raw_by_set[rec["set_id"]].append(rec)

    n_rows = len(matrix)
    if n_rows == 0:
        ax.text(0.5, 0.5, f"No significant features for {species}", ha="center", va="center")
        return []

    # Recompute row positions.
    y_positions = list(range(n_rows))
    row_lookup = {r["set_id"]: y for r, y in zip(matrix, y_positions)}

    # Background treatment bands.
    for idx, (label, color) in enumerate(zip(xlabels, xcolors)):
        ax.axvspan(idx - 0.5, idx + 0.5, color=color, alpha=0.18, zorder=0)
        ax.axvline(idx - 0.5, color="#E5E7EB", lw=0.8, zorder=1)
    ax.axvline(len(xlabels) - 0.5, color="#E5E7EB", lw=0.8, zorder=1)

    cmap = plt.get_cmap("RdBu_r")
    norm = TwoSlopeNorm(vmin=-3.6, vcenter=0.0, vmax=3.6)

    for y, row in enumerate(matrix):
        recs = raw_by_set.get(row["set_id"], [])
        for x, contrast in enumerate(row_order):
            rec = next((rr for rr in recs if rr["contrast"] == contrast), None)
            if not rec:
                continue
            nes = float(rec["NES"])
            padj = float(rec["padj"])
            ax.scatter(
                x,
                y,
                s=estimate_bubble_size(padj),
                facecolor=cmap(norm(nes)),
                edgecolor="#5A4B4B",
                linewidth=0.45,
                alpha=0.96,
                zorder=3,
            )

    ax.set_xlim(-0.5, len(row_order) - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks(range(len(row_order)))
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=9)
    for tick, color in zip(ax.get_xticklabels(), xcolors):
        tick.set_color(color)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([row["set_id"] for row in matrix], fontsize=8)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", color="#F3F4F6", linewidth=0.8, zorder=0)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#D1D5DB")
    ax.spines["bottom"].set_color("#D1D5DB")
    ax.set_title(species, loc="left", fontsize=13, fontweight="bold", pad=5)
    ax.set_ylabel("KEGG module / pathway", fontsize=11)

    # Right-side pathway annotations.
    ann_ax.set_xlim(0, 1)
    ann_ax.set_ylim(n_rows - 0.5, -0.5)
    ann_ax.axis("off")
    groups = []
    current = None
    for idx, row in enumerate(matrix):
        group = row["subcategory"]
        if current is None or group != current["group"]:
            if current is not None:
                current["end"] = idx - 1
                groups.append(current)
            current = {"group": group, "start": idx}
    if current is not None:
        current["end"] = n_rows - 1
        groups.append(current)

    for g in groups:
        start = g["start"]
        end = g["end"]
        mid = (start + end) / 2
        label = g["group"]
        ann_ax.plot([0.03, 0.03], [start, end], color="#7C7C7C", lw=1.0, clip_on=False)
        ann_ax.plot([0.03, 0.085], [start, start], color="#7C7C7C", lw=1.0, clip_on=False)
        ann_ax.plot([0.03, 0.085], [end, end], color="#7C7C7C", lw=1.0, clip_on=False)
        ann_ax.text(
            0.12,
            mid,
            textwrap.fill(label, width=30),
            ha="left",
            va="center",
            fontsize=8.2,
            color="#404040",
        )

    return matrix


def save_outputs(fig):
    out_png = OUTDIR / "reproduce_panel.png"
    out_pdf = OUTDIR / "reproduce_panel.pdf"
    fig.savefig(out_png, dpi=320, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    return out_png, out_pdf


def main():
    annotations = load_annotations(ANNOTATION_PATH)
    raw_rows: dict[str, list[dict[str, str]]] = {}
    selected_all: dict[str, list[dict[str, object]]] = {}

    for species, spec in SPECIES_SPECS.items():
        raw = read_csv_rows(spec["path"])
        raw_rows[species] = raw
        selected_all[species] = summarize_species(raw, annotations, species)

    bac_rows = selected_all["Bacillus"]
    rho_rows = selected_all["Rhodanobacter"]
    max_rows = max(len(bac_rows), len(rho_rows))
    main_h = max(14.5, 0.095 * max_rows + 4.2)
    legend_h = 1.7
    fig = plt.figure(figsize=(18.0, main_h + legend_h))
    gs = GridSpec(
        2,
        5,
        figure=fig,
        height_ratios=[main_h, legend_h],
        width_ratios=[1.0, 0.30, 0.10, 1.0, 0.30],
        hspace=0.10,
        wspace=0.14,
    )

    ax_bac = fig.add_subplot(gs[0, 0])
    ann_bac = fig.add_subplot(gs[0, 1])
    spacer = fig.add_subplot(gs[0, 2])
    spacer.axis("off")
    ax_rho = fig.add_subplot(gs[0, 3])
    ann_rho = fig.add_subplot(gs[0, 4])
    legend_ax = fig.add_subplot(gs[1, :])
    legend_ax.axis("off")

    bac_matrix = draw_panel(ax_bac, ann_bac, "Bacillus", bac_rows, raw_rows["Bacillus"], None)
    rho_matrix = draw_panel(ax_rho, ann_rho, "Rhodanobacter", rho_rows, raw_rows["Rhodanobacter"], None)

    cmap = plt.get_cmap("RdBu_r")
    sm = plt.cm.ScalarMappable(norm=TwoSlopeNorm(vmin=-3.6, vcenter=0.0, vmax=3.6), cmap=cmap)
    cax = fig.add_axes([0.27, 0.05, 0.22, 0.02])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("Enrichment (NES)", fontsize=10)
    cb.ax.tick_params(labelsize=9, length=2)

    size_ax = fig.add_axes([0.53, 0.03, 0.15, 0.08])
    size_ax.axis("off")
    size_ax.text(0.0, 0.97, "-log10(FDR)", fontsize=10, fontweight="bold", va="top")
    for idx, val in enumerate([2, 3, 4, 5]):
        size = estimate_bubble_size(10 ** (-val))
        size_ax.scatter(0.15 + idx * 0.18, 0.42, s=size, facecolor="white", edgecolor="#808080", linewidth=0.9)
        size_ax.text(0.15 + idx * 0.18, 0.08, str(val), ha="center", va="center", fontsize=9, color="#404040")
    size_ax.set_xlim(0, 1)
    size_ax.set_ylim(0, 1)

    note_ax = fig.add_axes([0.70, 0.02, 0.28, 0.09])
    note_ax.axis("off")
    note_ax.text(
        0.0,
        0.92,
        "Abbreviations:\nP, phage; K, kanamycin; Al, aluminum.\nAll conditions are shown within each species panel.\nAnnotation: curated KEGG pathway CSV.",
        fontsize=9,
        color="#4B5563",
        va="top",
    )

    out_png, out_pdf = save_outputs(fig)

    selected_rows = []
    for species, items in selected_all.items():
        for item in items:
            selected_rows.append(
                {
                    "species": species,
                    "set_id": item["set_id"],
                    "label": item["label"],
                    "group": item["group"],
                    "top_level_category": item["top_level_category"],
                    "subcategory": item["subcategory"],
                    "comparison_status": item["comparison_status"],
                    "pathways": item["pathways"],
                    "direction": item["direction"],
                    "n_sig": item["n_sig"],
                    "best_padj": item["best_padj"],
                    "best_nes": item["best_nes"],
                    "max_abs_nes": item["max_abs_nes"],
                }
            )

    with (OUTDIR / "reproduce_panel_selected_data.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "species",
                "set_id",
                "label",
                "group",
                "top_level_category",
                "subcategory",
                "comparison_status",
                "pathways",
                "direction",
                "n_sig",
                "best_padj",
                "best_nes",
                "max_abs_nes",
            ],
        )
        writer.writeheader()
        for item in selected_rows:
            writer.writerow(item)

    meta = {
        "output_png": str(out_png),
        "output_pdf": str(out_pdf),
        "bacillus_rows": len(bac_matrix),
        "rhodanobacter_rows": len(rho_matrix),
        "selected_counts": {species: len(items) for species, items in selected_all.items()},
        "annotation_source": str(ANNOTATION_PATH),
        "plot_scope": "all annotated pathways within each species panel",
    }
    (OUTDIR / "reproduce_panel_metadata.json").write_text(json.dumps(meta, indent=2))

    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")
    print(f"Selected features written to {OUTDIR / 'reproduce_panel_selected_data.csv'}")


if __name__ == "__main__":
    main()

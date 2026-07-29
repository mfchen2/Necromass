#!/usr/bin/env python3
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


ROOT = Path("/Users/mingfeichen")
OUTDIR = ROOT / "Manuscript"
FILES = [
    {
        "source": "Necromass heatmap pool",
        "path": OUTDIR / "Targeted_adsorption_necromass_heatmap_041326_grouped.csv",
        "value_scale": 100.0,
        "groups": [
            ("Bacillus sediment", "Bacillus ratio"),
            ("Rhodano sediment", "Rhodano ratio"),
        ],
        "name_field": "lib_Compound_Name",
    },
    {
        "source": "Clay pH 7",
        "path": OUTDIR / "clay_adsorption_pH7_experiment.csv",
        "value_scale": 100.0,
        "groups": [
            ("clay", "Absorbance (clay)"),
        ],
        "name_field": "Metabolites",
    },
    {
        "source": "Ferrihydrite pH 7",
        "path": OUTDIR / "ferrihydrite_adsorption_pH7.csv",
        "groups": [
            ("iron mineral", "adsorption_percent"),
        ],
        "name_field": "metabolite",
    },
]

COLORS = {
    "clay": "#f4a261",
    "iron mineral": "#2a9d8f",
    "Bacillus sediment": "#577590",
    "Rhodano sediment": "#e76f51",
}

CHEM_GROUP_ORDER = [
    "Amino acids & peptides",
    "Nucleosides & bases",
    "Organic acids",
    "Polyamines & amines",
    "Vitamins & cofactors",
    "Polyols & glycans",
    "Aromatic / benzenoids",
    "Carbohydrates & sugars",
    "Other / synthetic",
]

CHEM_GROUP_COLORS = {
    "Amino acids & peptides": "#f58518",
    "Nucleosides & bases": "#4c78a8",
    "Organic acids": "#e45756",
    "Polyamines & amines": "#54a24b",
    "Vitamins & cofactors": "#b279a2",
    "Polyols & glycans": "#72b7b2",
    "Aromatic / benzenoids": "#e377c2",
    "Carbohydrates & sugars": "#9c755f",
    "Other / synthetic": "#8c8c8c",
}


def parse_percent(text):
    if text is None:
        return None
    text = str(text).strip().replace("%", "")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def normalize_name(name):
    """
    Normalize metabolite names for matching.
    Handles capitalization differences, punctuation, and common leading
    stereochemical prefixes like L-, D-, and DL-.
    """
    if name is None:
        return ""
    name = str(name).strip().replace("’", "'")
    name = re.sub(r"^(?:dl|d|l)[\-\s]+", "", name, flags=re.IGNORECASE)
    name = re.sub(r"^\ufeff", "", name)
    name = re.sub(r"\s+", " ", name)
    name = name.lower()
    name = re.sub(r"[^a-z0-9]+", "", name)
    return name


def prettify_name(name):
    """
    Convert names like '5-METHYLTHIOADENOSINE' to '5-Methylthioadenosine'
    while keeping punctuation and numbers intact.
    """
    name = str(name).strip().replace("’", "'")
    name = re.sub(r"^(?:dl|d|l)[\-\s]+", "", name, flags=re.IGNORECASE)
    name = name.replace("_", " ")
    name = re.sub(r"\s+", " ", name)

    def repl(match):
        word = match.group(0)
        if len(word) <= 2 and word.isupper():
            return word
        return word[0].upper() + word[1:].lower()

    return re.sub(r"[A-Za-z]+", repl, name)


def classify_chemical_group(name):
    name = normalize_name(name)
    if any(token in name for token in ["glucose", "fructose", "galactose", "mannose", "ribose", "maltose", "sucrose"]):
        return "Carbohydrates & sugars"
    if any(token in name for token in ["xylitol", "sorbitol", "mannitol", "ribitol", "glycerol", "glycol"]):
        return "Polyols & glycans"
    if any(token in name for token in ["phenyl", "benzene", "benzoic", "pyridin", "cinnam", "phenethyl", "benzoyl"]):
        return "Aromatic / benzenoids"
    if any(
        token in name
        for token in [
            "adenosine",
            "adenine",
            "hypoxanthine",
            "xanthine",
            "cytosine",
            "cytidine",
            "guanosine",
            "guanine",
            "uridine",
            "uracil",
            "thymine",
            "deoxyadenosine",
            "deoxyguanosine",
            "methyluridine",
        ]
    ):
        return "Nucleosides & bases"
    if any(
        token in name
        for token in [
            "glutamicacid",
            "arginine",
            "lysine",
            "phenylalanine",
            "isoleucine",
            "alanine",
            "leucine",
            "histidine",
            "proline",
            "ornithine",
            "valine",
        ]
    ):
        return "Amino acids & peptides"
    if any(
        token in name
        for token in [
            "choline",
            "betaine",
            "spermidine",
            "tyramine",
            "diethanolamine",
            "trometamol",
        ]
    ):
        return "Polyamines & amines"
    if any(
        token in name
        for token in [
            "pyruvicacid",
            "malicacid",
            "methylsuccinicacid",
            "citramalate",
            "hydroxybutyricacid",
            "hydroxymethylpentanoicacid",
            "urocanicacid",
            "pyroglutamicacid",
            "glutamicacid",
            "adipicacid",
        ]
    ):
        return "Organic acids"
    if any(token in name for token in ["pyridoxine", "pyridoxate", "nicotinamide", "pterin"]):
        return "Vitamins & cofactors"
    return "Other / synthetic"


def read_source(cfg):
    by_key = {}
    display_name = {}
    order = []
    seen_order = set()
    value_scale = float(cfg.get("value_scale", 1.0))

    with cfg["path"].open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        header = [h.lstrip("\ufeff") for h in header]
        col_idx = {name: idx for idx, name in enumerate(header)}

        for row in reader:
            raw_name = row[col_idx[cfg["name_field"]]]
            key = normalize_name(raw_name)
            if not key:
                continue

            if key not in display_name:
                display_name[key] = prettify_name(raw_name)
            if key not in seen_order:
                seen_order.add(key)
                order.append(key)

            for group_name, value_col in cfg["groups"]:
                value = parse_percent(row[col_idx[value_col]])
                if value is None or value <= 0:
                    continue
                value *= value_scale
                by_key[(key, group_name)] = max(by_key.get((key, group_name), -math.inf), value)

    return by_key, display_name, order


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Aggregate each source separately, then merge into one long table.
    merged = []
    display_lookup = {}
    source_orders = {}
    for cfg in FILES:
        agg, names, order = read_source(cfg)
        for key, value in names.items():
            display_lookup.setdefault(key, value)
        source_orders[cfg["source"]] = order
        for (key, group_name), value in agg.items():
            merged.append(
                {
                    "source": cfg["source"],
                    "group": group_name,
                    "key": key,
                    "metabolite": names[key],
                    "chem_group": classify_chemical_group(names[key]),
                    "value": value,
                }
            )

    # Keep only metabolites observed in at least two distinct groups.
    group_presence = defaultdict(set)
    for row in merged:
        if row["value"] > 0:
            group_presence[row["key"]].add(row["group"])
    keep_keys = {k for k, groups in group_presence.items() if len(groups) >= 2}
    merged = [row for row in merged if row["key"] in keep_keys]

    # Order metabolites by the updated necromass pool first, then append any
    # additional metabolites from the mineral datasets.
    max_by_key = defaultdict(float)
    chem_group_by_key = {}
    for row in merged:
        max_by_key[row["key"]] = max(max_by_key[row["key"]], row["value"])
        chem_group_by_key[row["key"]] = row["chem_group"]

    order = []
    seen = set()
    for source_name in [FILES[0]["source"], FILES[1]["source"], FILES[2]["source"]]:
        for key in source_orders.get(source_name, []):
            if key in max_by_key and key not in seen:
                seen.add(key)
                order.append(key)

    remaining = [k for k in max_by_key.keys() if k not in seen]
    remaining.sort(
        key=lambda k: (
            CHEM_GROUP_ORDER.index(chem_group_by_key.get(k, "Other / synthetic"))
            if chem_group_by_key.get(k, "Other / synthetic") in CHEM_GROUP_ORDER
            else len(CHEM_GROUP_ORDER),
            -max_by_key[k],
            display_lookup.get(k, k),
        ),
    )
    order.extend(remaining)
    order_map = {k: i for i, k in enumerate(order)}

    # Save merged table for reproducibility.
    csv_out = OUTDIR / "metabolite_adsorption_merged_long.csv"
    with csv_out.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["source", "group", "chem_group", "metabolite", "value", "key"]
        )
        writer.writeheader()
        for row in sorted(merged, key=lambda r: (order_map[r["key"]], r["source"], r["group"])):
            writer.writerow(row)

    # Build a single-panel lollipop-style plot spanning the updated metabolite pool.
    fig, ax = plt.subplots(figsize=(15.5, max(7, 0.30 * len(order) + 1.8)))
    y_positions = list(range(len(order)))
    bar_h = 0.15

    offsets = {
        "clay": -1.6 * bar_h,
        "iron mineral": -0.55 * bar_h,
        "Bacillus sediment": 0.55 * bar_h,
        "Rhodano sediment": 1.6 * bar_h,
    }

    for group_name in ["clay", "iron mineral", "Bacillus sediment", "Rhodano sediment"]:
        xs = [None] * len(order)
        for row in merged:
            if row["group"] == group_name:
                xs[order_map[row["key"]]] = row["value"]
        y = [yy + offsets[group_name] for yy in y_positions]
        vals = [0 if v is None else v for v in xs]
        ax.hlines(
            y,
            0,
            vals,
            linewidth=4.5,
            color=COLORS[group_name],
            alpha=0.95,
        )
        ax.scatter(
            vals,
            y,
            s=16,
            color=COLORS[group_name],
            zorder=3,
            label=group_name,
        )

    ax.set_title(
        "Merged metabolite adsorption across mineral and necromass particle types",
        fontsize=18,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlabel("Adsorption (%)", fontsize=14, fontweight="bold")
    ax.set_yticks(y_positions)
    ax.set_yticklabels([display_lookup[k] for k in order], fontsize=9.5)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    ax.tick_params(axis="x", labelsize=11)
    ax.set_xlim(0, 105)

    # Add a colored strip to indicate chemical group.
    for y, key in zip(y_positions, order):
        chem_group = chem_group_by_key.get(key, "Other / synthetic")
        ax.scatter(
            [-3.0],
            [y],
            s=90,
            marker="s",
            color=CHEM_GROUP_COLORS[chem_group],
            clip_on=False,
            zorder=4,
        )
    ax.text(
        -0.11,
        1.015,
        "Chemical group",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    legend_handles = [Patch(color=COLORS[g], label=g) for g in [
        "clay",
        "iron mineral",
        "Bacillus sediment",
        "Rhodano sediment",
    ]]
    chem_handles = [Patch(color=CHEM_GROUP_COLORS[g], label=g) for g in CHEM_GROUP_ORDER]
    leg1 = ax.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=False,
        fontsize=10,
        ncol=2,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=chem_handles,
        loc="lower right",
        frameon=False,
        fontsize=9,
        ncol=2,
        title="Chemical group",
        title_fontsize=10,
    )

    fig.tight_layout()

    png_out = OUTDIR / "metabolite_adsorption_merged_barplot.png"
    fig.savefig(png_out, dpi=220, bbox_inches="tight")

    print(f"Wrote {png_out}")
    print(f"Wrote {csv_out}")


if __name__ == "__main__":
    main()

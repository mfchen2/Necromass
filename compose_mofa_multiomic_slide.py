#!/usr/bin/env python3

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np


ROOT = Path("/Users/mingfeichen/Manuscript")
COMPONENT_DIR = ROOT / "outputs" / "manual-20260601-a8" / "presentations" / "mofa_multiomic_slide" / "assets"
OUT_BASE = COMPONENT_DIR / "mofa_multiomic_slide_ready"

IMG_RHO = COMPONENT_DIR / "rhodanobacter_latent_state_slide.png"
IMG_BAC = COMPONENT_DIR / "bacillus_latent_state_slide.png"
IMG_BRIDGE = ROOT / "combined_bacillus_rhodano_bridge_heatmaps.png"


def crop_white(arr, thresh=245, pad=12):
    gray = arr[..., :3].mean(axis=2)
    mask = gray < thresh
    if not mask.any():
        return arr
    ys, xs = np.where(mask)
    y0 = max(int(ys.min()) - pad, 0)
    y1 = min(int(ys.max()) + pad + 1, arr.shape[0])
    x0 = max(int(xs.min()) - pad, 0)
    x1 = min(int(xs.max()) + pad + 1, arr.shape[1])
    return arr[y0:y1, x0:x1]


def load_trim(path, top_crop=0):
    arr = mpimg.imread(path)
    if arr.dtype.kind == "f":
        # Keep consistent alpha handling for PNGs loaded as floats.
        arr = (np.clip(arr, 0, 1) * 255).astype(np.uint8)
    else:
        arr = arr.copy()
    arr = crop_white(arr)
    if top_crop > 0 and arr.shape[0] > top_crop:
        arr = arr[top_crop:, :, :]
    return arr


def add_panel_label(ax, label):
    ax.text(
        0.01,
        0.985,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        color="#111827",
        bbox=dict(facecolor="white", edgecolor="none", pad=0.2, alpha=0.65),
    )


def main():
    img_rho = load_trim(IMG_RHO)
    img_bac = load_trim(IMG_BAC)
    img_bridge = load_trim(IMG_BRIDGE, top_crop=95)

    fig = plt.figure(figsize=(16, 9), facecolor="white")
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.95, 1.35],
        hspace=0.10,
        wspace=0.06,
        left=0.03,
        right=0.985,
        bottom=0.035,
        top=0.90,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    for ax, img in [(ax1, img_rho), (ax2, img_bac), (ax3, img_bridge)]:
        ax.imshow(img, aspect="auto")
        ax.axis("off")

    add_panel_label(ax1, "A")
    add_panel_label(ax2, "B")
    add_panel_label(ax3, "C")

    fig.text(
        0.03,
        0.975,
        "MOFA links shared latent states to sparse transcript-protein-metabolite bridge modules",
        ha="left",
        va="top",
        fontsize=19,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.03,
        0.945,
        "Top row: best factor pairs by species. Bottom panel: cross-omic bridge modules linking transcripts, proteins, and metabolites.",
        ha="left",
        va="top",
        fontsize=10.5,
        color="#4B5563",
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{OUT_BASE}.{ext}", dpi=300, facecolor="white")


if __name__ == "__main__":
    main()

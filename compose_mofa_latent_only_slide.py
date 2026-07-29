#!/usr/bin/env python3

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np


ROOT = Path("/Users/mingfeichen/Manuscript")
COMPONENT_DIR = ROOT / "outputs" / "manual-20260601-a8" / "presentations" / "mofa_multiomic_slide" / "assets"
OUT_BASE = ROOT / "outputs" / "manual-20260601-a9" / "presentations" / "mofa_latent_only_slide" / "assets" / "mofa_latent_only_slide"

IMG_RHO = COMPONENT_DIR / "rhodanobacter_latent_state_slide.png"
IMG_BAC = COMPONENT_DIR / "bacillus_latent_state_slide.png"


def crop_white(arr, thresh=245, pad=10):
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
        arr = (np.clip(arr, 0, 1) * 255).astype(np.uint8)
    else:
        arr = arr.copy()
    arr = crop_white(arr, pad=12)
    if top_crop > 0 and arr.shape[0] > top_crop:
        arr = arr[top_crop:, :, :]
    return arr


def add_panel_label(ax, label):
    ax.text(
        0.012,
        0.985,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=20,
        fontweight="bold",
        color="#111827",
        bbox=dict(facecolor="white", edgecolor="none", pad=0.15, alpha=0.60),
    )


def main():
    out_parent = OUT_BASE.parent
    out_parent.mkdir(parents=True, exist_ok=True)

    img_rho = load_trim(IMG_RHO)
    img_bac = load_trim(IMG_BAC)

    fig = plt.figure(figsize=(16, 9), facecolor="white")
    gs = fig.add_gridspec(
        1,
        2,
        left=0.035,
        right=0.985,
        top=0.90,
        bottom=0.04,
        wspace=0.07,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    ax1.imshow(img_rho, aspect="auto")
    ax1.axis("off")
    add_panel_label(ax1, "A")

    ax2.imshow(img_bac, aspect="auto")
    ax2.axis("off")
    add_panel_label(ax2, "B")

    fig.text(
        0.035,
        0.975,
        "MOFA latent ordinations are enough to summarize the multi-omics link",
        ha="left",
        va="top",
        fontsize=19,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.035,
        0.945,
        "The same treatment structure separates the omics in shared latent space; bridge modules can stay in the supplement.",
        ha="left",
        va="top",
        fontsize=10.5,
        color="#4B5563",
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{OUT_BASE}.{ext}", dpi=300, facecolor="white")


if __name__ == "__main__":
    main()

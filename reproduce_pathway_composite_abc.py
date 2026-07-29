#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig_mingfeichen_composite_abc")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


MANUSCRIPT_ROOT = Path("/Users/mingfeichen/Manuscript")
OUTDIR = MANUSCRIPT_ROOT / "outputs/manual-20260609-a8/presentations/pathway_composite_abc/assets"
OUTDIR.mkdir(parents=True, exist_ok=True)

SOURCE_A = MANUSCRIPT_ROOT / "outputs/manual-20260609-a7/presentations/pathway_metabolite_relationships_ab/assets/figure_a_pathway_concordance_scatter.png"
SOURCE_B = MANUSCRIPT_ROOT / "outputs/manual-20260609-a7/presentations/pathway_metabolite_relationships_ab/assets/figure_b_pathway_examples.png"
SOURCE_C = MANUSCRIPT_ROOT / "outputs/manual-20260609-a1/presentations/fgsea_panel_a_reproduction/assets/reproduce_panel.png"


def load_image(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    return mpimg.imread(path)


def trim_whitespace(img, pad: int = 18, threshold: int = 250):
    if img.ndim == 3 and img.shape[2] >= 3:
        rgb = img[:, :, :3]
    else:
        rgb = img
    mask = (rgb < (threshold / 255.0)).any(axis=2)
    if not mask.any():
        return img
    ys, xs = mask.nonzero()
    y0 = max(int(ys.min()) - pad, 0)
    y1 = min(int(ys.max()) + pad + 1, img.shape[0])
    x0 = max(int(xs.min()) - pad, 0)
    x1 = min(int(xs.max()) + pad + 1, img.shape[1])
    return img[y0:y1, x0:x1]


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.set_title(label, loc="left", fontsize=18, fontweight="bold", pad=8)


def main() -> None:
    img_a = load_image(SOURCE_A)
    img_b = load_image(SOURCE_B)
    img_c = load_image(SOURCE_C)

    img_a = trim_whitespace(img_a)
    img_b = trim_whitespace(img_b)
    img_c = trim_whitespace(img_c)

    fig = plt.figure(figsize=(18.5, 34.0), facecolor="white")
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[0.8, 1.0, 1.9],
        left=0.015,
        right=0.985,
        top=0.985,
        bottom=0.02,
        wspace=0.0,
        hspace=0.11,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[2, 0])

    for ax, img, label in [
        (ax_a, img_a, "A"),
        (ax_b, img_b, "B"),
        (ax_c, img_c, "C"),
    ]:
        ax.imshow(img)
        ax.set_axis_off()
        add_panel_label(ax, label)

    fig.savefig(OUTDIR / "reproduce_panel.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(OUTDIR / "reproduce_panel.pdf", bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)

    print(f"Wrote {OUTDIR / 'reproduce_panel.png'}")
    print(f"Wrote {OUTDIR / 'reproduce_panel.pdf'}")


if __name__ == "__main__":
    main()

# Necromass Figure Code

Current figure and supporting-analysis scripts for the September 2026 manuscript revision are organized in [`revision_2026_09/`](revision_2026_09/). The detailed README there lists entry points, required inputs, methods and reproducibility limits.

## Figure index

| Manuscript item | Code |
|---|---|
| Figure 1: untargeted metabolite overlap, volcano plots and PCoA | [`figure1/`](revision_2026_09/figure1/) |
| Figure 2: corrected targeted-metabolite profiles | [`analysis/make_corrected_necromass_fig2.R`](revision_2026_09/analysis/make_corrected_necromass_fig2.R) |
| Figure 3: endpoint pathway and proteome–metabolome analyses | [`analysis/`](revision_2026_09/analysis/) |
| Figure 4: corrected MOFA2 analysis and factor-loading panels | [`mofa/`](revision_2026_09/mofa/) |
| Figure 5: sediment retention | [`analysis/make_corrected_mineral_fig.R`](revision_2026_09/analysis/make_corrected_mineral_fig.R) |
| Supplementary figures and retention predictions | [`figures/`](revision_2026_09/figures/), [`retention/`](revision_2026_09/retention/) |

This is a code release, not a data release or turnkey pipeline. Input workbooks, model objects and other study data are not included; configure paths and inspect the detailed workflow notes before running scripts. Figure 1 plots can be regenerated from the processed feature matrix, but raw-to-matrix LC–MS/MS processing is incomplete because the raw files and authoritative final filtering/ion-mode merge workflow were not available. See [`RAW_PROCESSING_PROVENANCE.md`](revision_2026_09/figure1/RAW_PROCESSING_PROVENANCE.md).

The repository contains only the current figure and supporting-analysis bundle; superseded exploratory scripts have been removed. No license has been added.

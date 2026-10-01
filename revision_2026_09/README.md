# September 2026 source snapshot

## Provenance and execution

`SOURCE_MANIFEST.csv` records source locations and SHA-256 hashes for bundled files. Most analysis scripts are retained verbatim; Figure 1 utilities were parameterized for portable input/output paths, and the Gower PCoA script was reconstructed from the processed matrix. The original workspace root was `/Users/mingfeichen/Manuscript`; many inputs reside under `/Users/mingfeichen`.

Before running, inspect each script's input/output constants, replace them with local paths, and use a new output directory. Do not execute the entire folder indiscriminately: some scripts are sequential revisions or produce overlapping outputs. Some MOFA scripts automatically install packages and contain exploratory postprocessing; review them before execution. Scripts called through `source()` or imports must also be available at their configured paths. This is not a turnkey workflow or a locked historical environment.

## Analysis map

| Analysis | Entry points | Required inputs |
|---|---|---|
| Figure 1 and untargeted PCoA | `figure1/` | Processed 2,787-feature, 36-sample Supplementary File S1 workbook; see `figure1/RAW_PROCESSING_PROVENANCE.md` |
| Endpoint metabolite contrasts | `analysis/recompute_endpoint_reference_comparison.R` or `.py` | Corrected targeted and processed untargeted matrices; inspect each implementation before choosing |
| Figure 2 targeted profiles | `analysis/make_corrected_necromass_fig2.R` | Corrected necromass workbook and contrast outputs; verify reference settings against endpoint script |
| Targeted PERMANOVA | `analysis/recompute_corrected_targeted_permanova.R` | Corrected targeted workbook and sample grouping |
| Figure 3 pathway contrasts | `analysis/recompute_endpoint_transcriptome_fgsea.R`, then `analysis/rebuild_endpoint_metabolite_pathway_summary.R` | RNA counts, sample metadata, KEGG gene sets, metabolite mapping and endpoint contrasts |
| Proteome integration | `analysis/make_proteome_metabolome_concordance_endpoint.R` | Proteome enrichment and endpoint metabolite pathway summaries |
| MOFA training | `mofa/MOFA2_Bacillus_corrected.R`, `mofa/MOFA2_corrected.R` (Rhodanobacter) | Matched RNA/protein/corrected metabolite matrices |
| Figure 4 verification | `mofa/audit_mofa_figure4.R`, then `mofa/render_verified.py` | Saved corrected HDF5/RDS models, audit CSVs, transcript annotations and plotting helpers |
| Figure 5 sediment retention | `analysis/make_corrected_mineral_fig.R` | Corrected workbook's separate mineral sheet |
| Sediment model | `retention/rebuild_corrected_mineral_prediction.R` | `New_metabolites_corrected.xlsx`, sheet `Supplememntary Table Y Mineral` |
| Cross-sorbent S10 | `retention/build_all_sorbent_predictions.R` | Sediment output CSV, `metabolite_adsorption_merged_long.csv`, descriptors in `Necromass Supplementary Table.xlsx` |
| Supplementary assembly | `figures/` | Original figure PDFs/images and caption source script |

The Figure 1 raw-processing audit distinguishes recovered MZmine batch XMLs from the unavailable raw `.mzML` files and final filtering/merge workflow. Downstream plots are reproducible from processed Supplementary File S1, but raw-to-matrix regeneration is incomplete. The available matrix also yields PCoA variance labels different from the prior artwork; see the audit before reusing its labels.

The verified MOFA renderer uses the bundled `mofa/original_plot_helpers.py`; configure its helper and input/output paths for the local workspace. The renderer output directory must contain the audit CSV exports. The helper reads annotations at import time.

Untreated mid/late contrasts use 0 h; stress endpoints use the corresponding late untreated reference in the endpoint scripts. Direct organism contrasts match conditions. MOFA pairs are selected for treatment separation, not total variance; factor numbers are organism-specific.

## Retention models

Natural sediment: `(mean NA - mean sediment)/mean NA` per metabolite/organism, clipped to [0,1], then [0.0001,0.9999] for logit transformation. `lmer(logit_response ~ charge_class + source + (1|metabolite))` uses REML and Wald intervals. Both organism observations are held out together for each metabolite. Predictions for the unseen metabolite use zero random intercepts.

Pooled ridge: `model.matrix(~ sorbent + chemical_group + charge_class)`; unstandardized treatment-coded predictors, no interactions/source/random effects. Lambda is fixed at 1, not tuned; the intercept is unpenalized and the objective uses unscaled residual sums of squares. All observations with the same normalized metabolite name are held out across organisms and sorbents. Normalization lowercases and removes non-alphanumeric characters; it is not chemical synonym resolution. Predictions are inverse-logit transformed and performance is reported by sorbent. These models test unseen metabolites, not unseen sorbents. Operational annotations and heterogeneous assay sources limit mechanistic interpretation. Retention measures aqueous-phase removal, not uniquely surface adsorption.

## Dependencies

R packages include readxl, readr, dplyr, tidyr, stringr, tibble, forcats, ggplot2, patchwork, scales, vegan, lme4, openxlsx, ragg, svglite, ComplexHeatmap, circlize, DESeq2, fgsea, MOFA2, matrixStats, data.table, reticulate and GGally. MOFA2 also needs its supported Python backend, mofapy2. Inspect each script for its complete package list and installation behavior.

Python dependencies are listed in `requirements.txt` without version pinning. This is not a reconstructed training environment.

## Verification and limitations

Packaging checks cover Python/R syntax, copied-source checksums and candidate credential patterns. The scientific analyses were not rerun end-to-end during packaging. Absolute paths, external data and intermediate-file dependencies remain. The release excludes workbooks, models, private reviewer credentials, cover letters and manuscript files. The included scripts cover the current manuscript figures and selected supporting analyses; they do not replace raw LC–MS/MS, RNA-seq or proteomics processing pipelines.

References: Bates et al. (2015), doi:10.18637/jss.v067.i01; Hoerl & Kennard (1970), doi:10.1080/00401706.1970.10488634; Roberts et al. (2017), doi:10.1111/ecog.02881.

# Manuscript Analysis Scripts

This folder collects the manuscript-specific analysis and figure-generation scripts located in `/Users/mingfeichen/Manuscript`.

## Organization

- `01_targeted_metabolomics/`: targeted metabolite summaries, tests, ordination, and composite plots.
- `02_pathway_integration/`: transcriptome pathway and exometabolome integration figures.
- `03_mofa_multiomics/`: MOFA2 latent-state panels, cross-omic loading plots, and composite figures.
- `04_adsorption/`: sediment or mineral adsorption summaries, heatmaps, bar plots, and predicted-versus-observed plots.
- `05_supporting_figures/`: additional figure-reproduction scripts found under the workspace `outputs/` directory.

## Manuscript mapping

- Figure 2, targeted metabolomics: scripts in `01_targeted_metabolomics/`.
- Figure 3, pathway-level transcriptome and exometabolome comparison: scripts in `02_pathway_integration/`.
- Figure 4, MOFA2 and cross-omic modules: scripts in `03_mofa_multiomics/`.
- Figure 5, targeted metabolite retention or adsorption: scripts in `04_adsorption/`.
- Figure 1, untargeted LC-MS/MS feature processing and ordination: the final raw-processing script was not found in the local manuscript workspace.

## Transcriptome-specific scripts

The transcriptome-related downstream scripts are included in `02_pathway_integration/`, `03_mofa_multiomics/`, and `05_supporting_figures/`. These cover pathway-level integration, transcriptome–metabolite specificity, transporter links, cross-omic modules, and figure reproduction. The upstream RNA-seq workflow described in the manuscript, including trimming, Bowtie2 alignment, featureCounts, DESeq2 contrasts, and the eggNOG/KEGG crosswalk, was not found as a complete manuscript-specific script in the local workspace.

## Reproducibility notes

These scripts were copied without changing their contents. Several scripts use absolute paths such as `/Users/mingfeichen/Manuscript` or `/Users/mingfeichen`, and some expect Excel workbooks, CSV/TSV files, RDS objects, or intermediate outputs that are not included here. They may therefore require path edits and the corresponding input files before execution.

The local workspace did not contain the complete upstream scripts for raw LC-MS/MS processing, RNA-seq preprocessing, proteomics processing, or the full MOFA2 model-fitting step. The scripts here primarily reproduce downstream analyses and figure generation from already processed inputs.

The MOFA2 model-fitting script and saved model object were also not found. The available MOFA2 scripts are downstream visualization and loading-analysis scripts. See `03_mofa_multiomics/MOFA2_FITTING_STATUS.md` for the missing reproducibility components.

## Suggested execution order

1. Prepare the required input tables and intermediate analysis objects.
2. Run targeted metabolomics scripts.
3. Run pathway integration scripts.
4. Run MOFA2 and cross-omic figure scripts.
5. Run adsorption scripts.
6. Compare generated outputs with the manuscript figures and supplementary tables.

## Provenance

Source directory: `/Users/mingfeichen/Manuscript`

Assembly date: 2026-07-29

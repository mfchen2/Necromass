# Necromass Manuscript Analysis Code

## September 2026 revision

Start with [revision_2026_09/README.md](revision_2026_09/README.md) for corrected-metabolite, endpoint-reference, MOFA and retention analyses. Existing root-level scripts are preserved as historical workflows. The notes below describe the earlier July collection and are superseded where the revision documentation differs.

This is a source-code release, not a self-contained data release or a validated one-command pipeline. Inputs and absolute paths must be configured before execution. No experimental data, manuscript documents or private reviewer credentials are included in the revision bundle. See the manuscript data-availability statement for deposited data.

The corrected MOFA fitting scripts have now been located and are included; the historical statement below that they were unavailable is no longer current. The repository also contains previously deposited RNA-seq preprocessing scripts at its root. No license has been added; licensing remains an author decision.

## Historical July collection

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
- Figure 1, untargeted downstream figures and ordination: see `revision_2026_09/figure1/`; raw-to-matrix processing remains incomplete.

## Transcriptome-specific scripts

The transcriptome-related downstream scripts are included in `02_pathway_integration/`, `03_mofa_multiomics/`, and `05_supporting_figures/`. These cover pathway-level integration, transcriptome–metabolite specificity, transporter links, cross-omic modules, and figure reproduction. The upstream RNA-seq workflow described in the manuscript, including trimming, Bowtie2 alignment, featureCounts, DESeq2 contrasts, and the eggNOG/KEGG crosswalk, was not found as a complete manuscript-specific script in the local workspace.

## Reproducibility notes

Most historical scripts were copied without changing their contents. Figure 1 utilities were parameterized for portable input/output paths, and the Gower PCoA script was reconstructed from the processed matrix. Several older scripts use absolute paths such as `/Users/mingfeichen/Manuscript` or `/Users/mingfeichen`, and some expect workbooks, CSV/TSV files, RDS objects, or intermediate outputs not included here.

The Figure 1 audit recovered MZmine 3.7.2 positive/negative batch XMLs from the matching NECPATH3 acquisition and added a Gower PCoA generator for the processed matrix. Raw `.mzML` data and the authoritative final filtering/ion-mode merge script remain unavailable, so raw-to-matrix regeneration is incomplete. See `revision_2026_09/figure1/RAW_PROCESSING_PROVENANCE.md`. Other upstream LC-MS/MS, RNA-seq, and proteomics components may also remain incomplete.

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

# Figure 1 Untargeted LC-MS/MS provenance

## What is reproducible here

The repository contains downstream Figure 1 scripts for the species-specific feature-overlap panels, endpoint-reference volcano plots, and the untargeted Gower PCoA. They run from the 2,787-row, 36-sample processed matrix distributed as Supplementary File S1; the workbook is not duplicated in this code repository. Set `UNTARGETED_INPUT_PATH` to that workbook and `UNTARGETED_OUTPUT_DIR` to a writable output directory before running each script.

```sh
export UNTARGETED_INPUT_PATH="/path/to/Supplementary files S1.xlsx"
export UNTARGETED_OUTPUT_DIR="figures/figure1"
export UNTARGETED_REFERENCE_MODE=endpoint
cd /path/to/Necromass/revision_2026_09/figure1
Rscript make_species_specific_untargeted_overlap.R
Rscript make_untargeted_metabolites_volcano.R
Rscript make_untargeted_pcoa.R
```

The overlap script writes separate *Bacillus* and *Rhodanobacter* panels and overlap-count tables. The volcano script performs replicate-level Welch tests on log2-transformed intensities with a pseudocount of 1, BH correction within contrast, and the current endpoint-reference scheme. The PCoA script uses the sample-by-feature intensity matrix, Gower dissimilarities (`vegan::vegdist(method = "gower")`) and classical multidimensional scaling (`stats::cmdscale`). It writes coordinates and axis variance alongside the figure.

## Upstream feature processing: located files and remaining gap

Two MZmine batch parameter files are included under `mzmine/`. They are the positive- and negative-mode MZmine 3.7.2 XMLs from the 22 October 2025 NECPATH3 acquisition. The raw-file naming in that acquisition contains the 36 `supern` biological samples represented in the processed matrix, plus treatment and extraction controls. The XMLs preserve the saved MZmine processing parameters and export settings, but retain facility-specific `/global/cfs/...` file paths; those paths must be remapped to the raw-file and output locations before execution.

The raw `.mzML` files are not included. I did not find an authoritative final script that applies the manuscript's post-export retention-time/peak-height/extraction-control filters, reconciles positive and negative ion tables, and produces this exact 2,787-row workbook. The MZmine XMLs alone therefore do **not** regenerate Supplementary File S1. The older manuscript text says MZmine 2.0, whereas the recovered batch files identify MZmine 3.7.2; reconcile that version with the historical instrument-processing record before changing the manuscript methods. Do not describe raw-to-matrix processing as fully reproducible from this repository until the raw files and final filter/merge script are deposited.

## Ordination audit

Running the saved Gower/PCoA calculation on the available Supplementary File S1 workbook produces PCoA1 = 60.4% and PCoA2 = 16.7% of positive eigenvalue variance. The earlier PCoA artwork is labelled 60.9% and 16.9%. The script reports the calculated values and warns on this discrepancy; the original labels should not be carried forward without verifying the exact historical matrix or preprocessing transformation used for that artwork.

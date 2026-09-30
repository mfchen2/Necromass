# MOFA2 Fitting Status

## Update: 30 September 2026

The historical status below is superseded. Corrected fitting scripts for both organisms are included in `revision_2026_09/mofa/`. Corrected HDF5 models were located and audited locally but are not distributed in this code-only release. The verified pairs are Bacillus Factor2/Factor6 and Rhodanobacter Factor2/Factor1. The revised renderer replaces the genuine Factor3 loadings previously shown in Rhodanobacter panel F, as well as correcting the label on panel D. Input matrices, annotations and saved model objects are required to rerun the audit; see the revision README.

## Historical status (superseded)

No manuscript-specific MOFA2 model-fitting script or saved MOFA2 model object was found in the local workspace.

The scripts in this directory perform downstream analyses and figure generation from precomputed MOFA2 outputs, including latent-factor rankings, factor scores, feature loadings, cross-omic modules, and composite panels.

To make the analysis fully reproducible, the following missing items should be added when available:

- The script that imports and aligns transcriptome, proteome, and metabolome matrices.
- The normalization, filtering, and feature-selection code.
- The MOFA2 model and training options, including likelihoods, view scaling, group centering, convergence settings, and iteration limits.
- The saved MOFA2 model object or a reproducible model-fitting command.
- The exact input files used to generate the factor scores and loadings.

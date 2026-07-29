# MOFA2 Fitting Status

No manuscript-specific MOFA2 model-fitting script or saved MOFA2 model object was found in the local workspace.

The scripts in this directory perform downstream analyses and figure generation from precomputed MOFA2 outputs, including latent-factor rankings, factor scores, feature loadings, cross-omic modules, and composite panels.

To make the analysis fully reproducible, the following missing items should be added when available:

- The script that imports and aligns transcriptome, proteome, and metabolome matrices.
- The normalization, filtering, and feature-selection code.
- The MOFA2 model and training options, including likelihoods, view scaling, group centering, convergence settings, and iteration limits.
- The saved MOFA2 model object or a reproducible model-fitting command.
- The exact input files used to generate the factor scores and loadings.

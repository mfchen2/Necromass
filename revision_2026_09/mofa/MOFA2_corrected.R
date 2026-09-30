############################################################
# Full MOFA2 pipeline for integrating transcriptomics,
# proteomics, and metabolomics from CSV files
#
# 1) Install & load packages
# 2) Read CSVs
# 3) Align samples across omics
# 4) Normalize / transform:
#       - RNA: DESeq2 + VST
#       - Proteomics / metabolomics: log2 + row scaling
# 5) Feature selection (most variable)
# 6) Build and run MOFA2 model
# 7) Basic plots
############################################################

##############################
# 1. Install and load packages
##############################

required_pkgs <- c(
  "MOFA2",      # Multi-omics factor analysis
  "DESeq2",     # RNA normalization / VST
  "matrixStats",# rowVars, rowSds
  "data.table"  # fast reading / writing (optional)
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " ...")
    if (pkg == "MOFA2" || pkg == "DESeq2") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("GGally", quietly = TRUE)) {
  BiocManager::install("GGally", ask = FALSE, update = FALSE)
}
library(reticulate)
library(GGally)

# Force reticulate to use the conda env you created
use_condaenv("mofa-env", required = TRUE)

# Verify this points to .../envs/mofa-env/bin/python
py_config()

# Optional sanity check
py_module_available("mofapy2")

library(MOFA2)
library(DESeq2)
library(matrixStats)
library(data.table)

########################
# 2. User file locations
########################

# >>>>>>> EDIT THESE THREE LINES <<<<<<<
rna_file   <- "~/Rhodano_mRNAseq/counts/gene_counts_prokka.csv"        # transcriptomics (counts)
prot_file  <- "~/Desktop/Jupyter/Full_list_proteins_matrix_counts_sum_Rhodano_MOFA.csv"        # proteomics intensities
metab_file <- "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a8/presentations/mofa_composite_abcdef/assets/corrected_data/rhodanobacter/Rhodano_metabolites_corrected_for_mofa.csv"      # ion-suppression-corrected metabolomics intensities

########################################
# 3. Helper to read CSV with rownames
########################################

read_omics_csv <- function(file) {
  # Assumes first column = feature ID, remaining columns = samples
  df <- read.csv(file, header = TRUE, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  # Drop any completely NA rows
  keep <- rowSums(is.na(mat)) < ncol(mat)
  mat[keep, , drop = FALSE]
}

message("Reading input CSV files...")

rna_counts <- read_omics_csv(rna_file)
prot_mat   <- read_omics_csv(prot_file)
metab_mat  <- read_omics_csv(metab_file)

# Match metabolite timepoint notation to the RNA/protein sample IDs.
colnames(metab_mat) <- sub("_(0|8|24|48)h_", "_\\1hr_", colnames(metab_mat))

cat("Dimensions:\n")
cat("  RNA counts:        ", dim(rna_counts)[1], "features x", dim(rna_counts)[2], "samples\n")
cat("  Proteomics:        ", dim(prot_mat)[1],   "features x", dim(prot_mat)[2],   "samples\n")
cat("  Metabolomics:      ", dim(metab_mat)[1],  "features x", dim(metab_mat)[2],  "samples\n\n")

######################################################
# 4. Align samples across the three omics matrices
######################################################

common_samples <- Reduce(intersect, list(
  colnames(rna_counts),
  colnames(prot_mat),
  colnames(metab_mat)
))

if (length(common_samples) == 0) {
  stop("No overlapping sample names across the three omics. ",
       "Please check your CSV column names.")
}

cat("Number of common samples across all 3 omics:", length(common_samples), "\n")
cat("Common sample IDs:\n")
print(common_samples)

# Subset and reorder columns to the same sample order
rna_counts <- rna_counts[, common_samples, drop = FALSE]
prot_mat   <- prot_mat[, common_samples, drop = FALSE]
metab_mat  <- metab_mat[, common_samples, drop = FALSE]

# metab_mat <- read.csv("~/All_targeted_metabolites_011426.csv",row.names=1)
# metab_mat <- as.matrix(metab_mat[,-c(1,2)])

########################################
# 5. Normalization / transformations
########################################

###########
# 5.1 RNA-seq
###########

# Extract sample names from the RNA count matrix
samples <- colnames(rna_counts)

# Derive condition by removing trailing "_1", "_2", "_3" from names
# e.g. "R_0hr_1" -> "R_0hr", "R_Al_3" -> "R_Al"
condition <- sub("_[0-9]+$", "", samples)

# Build colData for DESeq2
rna_coldata <- data.frame(
  row.names = samples,
  condition = factor(condition)
)

# Check it looks right
table(rna_coldata$condition)
# Should show:
# R_0hr  R_24hr  R_48hr  R_P  R_K  R_Al 
#    3      3      3     3    3    3

# Use this colData in DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = round(rna_counts),
  colData   = rna_coldata,
  design    = ~ 1   # still ~1 since we’re only doing VST
)

vsd <- vst(dds, blind = TRUE)
rna_norm <- assay(vsd)


#############################
# 5.2 Proteomics & metabolomics
#############################

# For continuous intensities: log2(x + 1) and row-wise scaling
log2_plus1 <- function(mat) {
  # Replace negative values with NA before log
  mat[mat < 0] <- NA
  log2(mat + 1)
}

scale_rows <- function(mat) {
  # Row-wise z-score: (x - mean) / sd
  m <- rowMeans(mat, na.rm = TRUE)
  s <- matrixStats::rowSds(mat, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- 1
  mat_scaled <- sweep(mat, 1, m, "-")
  mat_scaled <- sweep(mat_scaled, 1, s, "/")
  rownames(mat_scaled) <- rownames(mat)
  colnames(mat_scaled) <- colnames(mat)
  mat_scaled
}

# rna_B <- read.csv("~/exo_pellet_harmonized_gene_hits_Bacillus_only_010826.csv",row.names=1)


prot_log  <- log2_plus1(prot_mat)
metab_log <- log2_plus1(metab_mat)

prot_norm  <- scale_rows(prot_log)
metab_norm <- scale_rows(metab_log)

cat("Proteomics normalization done. Resulting matrix:",
    dim(prot_norm)[1], "features x", dim(prot_norm)[2], "samples\n")
cat("Metabolomics normalization done. Resulting matrix:",
    dim(metab_norm)[1], "features x", dim(metab_norm)[2], "samples\n")

#######################################
# 6. Feature selection (most variable)
#######################################

select_most_variable <- function(mat, n_keep) {
  n_keep <- min(n_keep, nrow(mat))
  vars <- matrixStats::rowVars(mat, na.rm = TRUE)
  idx <- order(vars, decreasing = TRUE)[seq_len(n_keep)]
  mat[idx, , drop = FALSE]
}

# You can tweak these numbers depending on how many features you have
rna_norm_filt   <- select_most_variable(rna_norm,   n_keep = 3000)
prot_norm_filt  <- select_most_variable(prot_norm,  n_keep = 2000)
metab_norm_filt <- select_most_variable(metab_norm, n_keep = 1000)

cat("\nAfter feature selection:\n")
cat("  RNA:        ", dim(rna_norm_filt)[1],   "features x", dim(rna_norm_filt)[2],   "samples\n")
cat("  Proteomics: ", dim(prot_norm_filt)[1],  "features x", dim(prot_norm_filt)[2],  "samples\n")
cat("  Metabolome: ", dim(metab_norm_filt)[1], "features x", dim(metab_norm_filt)[2], "samples\n")

############################################
# 7. Build MOFA object from three omics views
############################################

data_list <- list(
  Transcriptome = rna_norm_filt,
  Proteome      = prot_norm_filt,
  Metabolome    = metab_norm_filt
)

mofa_obj <- create_mofa(data_list)

# Quick data overview plot (saved to file)
pdf("Rhodano_MOFA_data_overview.pdf", width = 7, height = 5)
plot_data_overview(mofa_obj)
dev.off()

######################################################
# 8. Set data, model, and training options for MOFA2
######################################################

# Data options
data_opts <- get_default_data_options(mofa_obj)
data_opts$scale_views   <- TRUE   # scale each view internally
data_opts$scale_groups  <- FALSE  # single group
data_opts$center_groups <- TRUE

# Model options
model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- 10      # you can tune this
model_opts$likelihoods <- c(
  Transcriptome = "gaussian",
  Proteome      = "gaussian",
  Metabolome    = "gaussian"
)

# Training options
train_opts <- get_default_training_options(mofa_obj)
train_opts$convergence_mode <- "medium"
train_opts$maxiter          <- 1000
train_opts$verbose          <- TRUE

########################################
# 9. Prepare and run the MOFA2 model
########################################

mofa_prep <- prepare_mofa(
  object           = mofa_obj,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

outfile <- "Rhodano_MOFA_multiomics_model.hdf5"

mofa_trained <- run_mofa(
  object       = mofa_prep,
  outfile      = outfile,
  use_basilisk = FALSE   # let MOFA2 manage Python env
)

cat("\nMOFA training finished. Model saved to:", outfile, "\n")

########################
# 10. Downstream analyses
########################

# 10.1 Variance explained
pdf("Rhodano_MOFA_variance_explained.pdf", width = 7, height = 5)
plot_variance_explained(mofa_trained)
dev.off()

# 10.2 If you have sample metadata (conditions/time/etc.), add here:
# >>>>>>> EDIT THIS TO YOUR REAL METADATA <<<<<<<
# Make sure we assign the same condition factor to the MOFA object
# colnames(rna_norm) should be in the same order as the samples
samples_metadata(mofa_trained)$condition <- rna_coldata[colnames(rna_norm), "condition"]

# Example: plot MOFA factors colored by your conditions (R_0hr, R_24hr, etc.)
pdf("Rhodano_MOFA_factors_by_condition.pdf", width = 7, height = 5)
plot_factors(
  mofa_trained,
  factors  = c(1:5),
  color_by = "condition"
)
dev.off()

# 10.3 Top features per factor for each view
dir.create("Rhodano_MOFA_top_features", showWarnings = FALSE)

views <- c("Transcriptome", "Proteome", "Metabolome")
for (v in views) {
  for (f in 1:3) {
    pdf(sprintf("Rhodano_MOFA_top_features/%s_Factor%d.pdf", v, f),
        width = 7, height = 5)
    
    p <- plot_weights(
      mofa_trained,
      view      = v,
      factor    = f,
      nfeatures = 20
    )
    print(p)  # <- this actually draws the plot into the PDF
    
    dev.off()
  }
}


cat("\nBasic plots written to: \n",
    "  MOFA_data_overview.pdf\n",
    "  MOFA_variance_explained.pdf\n",
    "  MOFA_factors_by_condition.pdf\n",
    "  MOFA_top_features/ (per view x factor)\n")

cat("\nDone.\n")

####analyzing the results got from MOFA####
ve <- get_variance_explained(mofa_trained)
# Total R2 explained per view:
ve$r2_total[[1]]

# R2 per factor per view:
head(ve$r2_per_factor[[1]])

# One-factor view is often clearer than scatterplots
pdf("Rhodano_MOFA_factor_annotation.pdf", 9, 10)
for (k in 1:5) {
  print(plot_factor(mofa_trained, factors = k, color_by = "condition"))
}
dev.off()

library(MOFA2)

# 1) Factors (long format)
Z  <- get_factors(mofa_trained, as.data.frame = TRUE)  # sample, factor, value :contentReference[oaicite:1]{index=1}
md <- samples_metadata(mofa_trained)[, c("sample", "condition")]

df <- merge(Z, md, by = "sample")

# 2) Make sure condition is a factor and has >1 level overall
df$condition <- factor(df$condition)
stopifnot(nlevels(df$condition) > 1)

# 3) Compute ANOVA p-values per factor safely
factor_ids <- sort(unique(as.character(df$factor)))

p_anova <- setNames(rep(NA_real_, length(factor_ids)), factor_ids)

for (fk in factor_ids) {
  d <- df[df$factor == fk, , drop = FALSE]
  
  # If this factor's samples have only one condition level, ANOVA is undefined
  if (nlevels(droplevels(d$condition)) < 2) next
  
  p_anova[fk] <- tryCatch(
    summary(aov(value ~ condition, data = d))[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
}

p_adj <- p.adjust(p_anova, method = "BH")

res <- data.frame(
  factor = names(p_anova),
  p      = as.numeric(p_anova),
  p_adj  = as.numeric(p_adj),
  row.names = NULL
)

# Drop NAs (factors that couldn't be tested)
res <- res[is.finite(res$p) & !is.na(res$p), ]
res <- res[order(res$p_adj), ]

res

library(dplyr)

Z  <- get_factors(mofa_trained, as.data.frame = TRUE)
md <- samples_metadata(mofa_trained)[, c("sample","condition")]
df <- left_join(Z, md, by = "sample")

# For the key factors
key_factors <- c("Factor1","Factor3","Factor5","Factor6","Factor4")

df %>%
  filter(factor %in% key_factors) %>%
  group_by(factor, condition) %>%
  summarise(mean = mean(value), sd = sd(value), .groups="drop") %>%
  arrange(factor, desc(mean))

pdf("MOFA_key_factors_by_condition.pdf", 8, 10)
for (k in c(1,3,5,6,4)) {
  print(plot_factor(mofa_trained, factors = k, color_by = "condition",
                    dodge = TRUE, add_violin = TRUE))
}
dev.off()

# Positive end
plot_top_weights(mofa_trained, view="Proteome", factor=5, sign="positive", nfeatures=25)
# Negative end
plot_top_weights(mofa_trained, view="Transcriptome", factor=1, sign="negative", nfeatures=25)


  

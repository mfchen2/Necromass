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
  "data.table",  # fast reading / writing (optional)
  "dplyr"
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
rna_file   <- "~/Bacillus_RNAseq/counts/gene_counts.csv"        # transcriptomics (counts)
prot_file  <- "~/Desktop/Jupyter/Full_list_proteins_matrix_counts_sum_Bacillus_MOFA.csv"        # proteomics intensities
metab_file <- "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a8/presentations/mofa_composite_abcdef/assets/corrected_data/bacillus/Bacillus_metabolites_corrected_for_mofa.csv"      # ion-suppression-corrected metabolomics intensities

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
pdf("Bacillus_MOFA_data_overview.pdf", width = 7, height = 5)
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

outfile <- "Bacillus_MOFA_multiomics_model.hdf5"

mofa_trained <- run_mofa(
  object       = mofa_prep,
  outfile      = outfile,
  use_basilisk = FALSE   # let MOFA2 manage Python env
)

cat("\nMOFA training finished. Model saved to:", outfile, "\n")

########################
# 10. Downstream analyses
########################

library(ggplot2)
library(dplyr)

# Optional: give factors human-readable names
factors_names(mofa_trained) <- paste0("Factor", seq_len(get_dimensions(mofa_trained)$K))

# 10.0 Add sample metadata to the model following the MOFA2 downstream vignette
sample_metadata <- data.frame(
  sample    = common_samples,
  condition = as.character(rna_coldata[common_samples, "condition"]),
  stringsAsFactors = FALSE
)
samples_metadata(mofa_trained) <- sample_metadata

# Helper: save ggplot objects or lists of ggplots
save_plot_object <- function(plot_obj, file, width = 8, height = 6) {
  pdf(file, width = width, height = height)
  if (inherits(plot_obj, "ggplot") || inherits(plot_obj, "gg")) {
    print(plot_obj)
  } else if (is.list(plot_obj)) {
    for (p in plot_obj) {
      if (inherits(p, "ggplot") || inherits(p, "gg")) print(p)
    }
  } else {
    print(plot_obj)
  }
  dev.off()
}

# Helper: extract top weights as a table
get_top_weights_table <- function(model, view, factor, n = 50) {
  w <- get_weights(model, views = view, factors = factor)[[view]]
  if (is.matrix(w)) {
    w <- w[, 1]
  }
  out <- data.frame(
    feature = names(w),
    weight  = as.numeric(w),
    stringsAsFactors = FALSE
  )
  out[order(abs(out$weight), decreasing = TRUE), , drop = FALSE][seq_len(min(n, nrow(out))), ]
}

dir.create("Bacillus_MOFA_top_features", showWarnings = FALSE)
dir.create("Bacillus_MOFA_data_patterns", showWarnings = FALSE)
dir.create("Bacillus_MOFA_tables", showWarnings = FALSE)

# 10.1 Variance decomposition
ve <- get_variance_explained(mofa_trained)

write.csv(ve$r2_total[[1]],
          file = "Bacillus_MOFA_tables/Bacillus_MOFA_r2_total.csv",
          row.names = FALSE)
write.csv(ve$r2_per_factor[[1]],
          file = "Bacillus_MOFA_tables/Bacillus_MOFA_r2_per_factor.csv",
          row.names = FALSE)

save_plot_object(
  plot_variance_explained(mofa_trained, x = "view", y = "factor"),
  "Bacillus_MOFA_variance_explained.pdf",
  width = 8, height = 6
)

save_plot_object(
  plot_variance_explained(mofa_trained, x = "view", y = "factor", plot_total = TRUE),
  "Bacillus_MOFA_variance_explained_with_total.pdf",
  width = 8, height = 6
)

# 10.2 Visualisation of factors
n_factors_to_plot <- min(5, get_dimensions(mofa_trained)$K)

save_plot_object(
  plot_factor(
    mofa_trained,
    factors      = seq_len(n_factors_to_plot),
    color_by     = "condition",
    dodge        = TRUE,
    add_violin   = TRUE,
    violin_alpha = 0.25
  ),
  "Bacillus_MOFA_factor_annotation.pdf",
  width = 9, height = 10
)

save_plot_object(
  plot_factors(
    mofa_trained,
    factors  = seq_len(min(3, n_factors_to_plot)),
    color_by = "condition"
  ),
  "Bacillus_MOFA_factors_by_condition.pdf",
  width = 8, height = 6
)

# 10.3 Visualisation of feature weights
views <- c("Transcriptome", "Proteome", "Metabolome")

for (v in views) {
  for (f in seq_len(min(3, get_dimensions(mofa_trained)$K))) {
    
    save_plot_object(
      plot_weights(
        mofa_trained,
        view      = v,
        factor    = f,
        nfeatures = 20,
        scale     = TRUE,
        abs       = FALSE
      ),
      sprintf("Bacillus_MOFA_top_features/%s_Factor%d_weights.pdf", v, f),
      width = 8, height = 6
    )
    
    save_plot_object(
      plot_top_weights(
        mofa_trained,
        view      = v,
        factor    = f,
        nfeatures = 20
      ),
      sprintf("Bacillus_MOFA_top_features/%s_Factor%d_top_weights.pdf", v, f),
      width = 8, height = 6
    )
    
    write.csv(
      get_top_weights_table(mofa_trained, view = v, factor = f, n = 50),
      file = sprintf("Bacillus_MOFA_tables/%s_Factor%d_top_weights.csv", v, f),
      row.names = FALSE
    )
  }
}

# 10.4 Visualisation of covariation patterns in the input data
for (v in views) {
  for (f in seq_len(min(3, get_dimensions(mofa_trained)$K))) {
    
    pdf(sprintf("Bacillus_MOFA_data_patterns/%s_Factor%d_heatmap.pdf", v, f),
        width = 8, height = 6)
    plot_data_heatmap(
      mofa_trained,
      view = v,
      factor = f,
      features = 20,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE
    )
    dev.off()
    
    save_plot_object(
      plot_data_scatter(
        mofa_trained,
        view = v,
        factor = f,
        features = 5,
        add_lm = TRUE,
        color_by = "condition"
      ),
      sprintf("Bacillus_MOFA_data_patterns/%s_Factor%d_scatter.pdf", v, f),
      width = 8, height = 6
    )
  }
}

# 10.5 Extract factors, weights and data for custom downstream analysis
factors_list <- get_factors(mofa_trained, factors = "all")
weights_list <- get_weights(mofa_trained, views = "all", factors = "all")
data_list_mofa <- get_data(mofa_trained)

saveRDS(
  list(
    factors = factors_list,
    weights = weights_list,
    data    = data_list_mofa
  ),
  file = "Bacillus_MOFA_tables/Bacillus_MOFA_extracted_objects.rds"
)

# Long-format factors for simple statistics / modeling
Z_long <- get_factors(mofa_trained, as.data.frame = TRUE)
df_factors <- left_join(Z_long, sample_metadata, by = "sample")
df_factors$condition <- factor(df_factors$condition)

# 10.6 ANOVA: which factors are associated with condition?
factor_ids <- unique(as.character(df_factors$factor))

anova_results <- lapply(factor_ids, function(fk) {
  d <- df_factors[df_factors$factor == fk, , drop = FALSE]
  
  if (nlevels(droplevels(d$condition)) < 2) {
    return(data.frame(factor = fk, p = NA_real_))
  }
  
  pval <- tryCatch(
    summary(aov(value ~ condition, data = d))[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
  
  data.frame(factor = fk, p = pval)
}) %>%
  bind_rows() %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

write.csv(anova_results,
          file = "Bacillus_MOFA_tables/Bacillus_MOFA_factor_condition_anova.csv",
          row.names = FALSE)

print(anova_results)

# 10.7 Plot the most condition-associated factors
key_factors <- anova_results %>%
  filter(!is.na(p_adj)) %>%
  slice_head(n = 5) %>%
  pull(factor)

key_factor_nums <- unique(as.integer(gsub("[^0-9]", "", key_factors)))
key_factor_nums <- key_factor_nums[!is.na(key_factor_nums)]

if (length(key_factor_nums) > 0) {
  save_plot_object(
    plot_factor(
      mofa_trained,
      factors = key_factor_nums,
      color_by = "condition",
      dodge = TRUE,
      add_violin = TRUE,
      violin_alpha = 0.25
    ),
    "MOFA_key_factors_by_condition.pdf",
    width = 8, height = 10
  )
  
  factor_summary <- df_factors %>%
    filter(factor %in% key_factors) %>%
    group_by(factor, condition) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor, desc(mean))
  
  write.csv(
    factor_summary,
    file = "Bacillus_MOFA_tables/Bacillus_MOFA_key_factor_condition_summary.csv",
    row.names = FALSE
  )
}

cat("
Updated downstream outputs written to:
",
    "  Bacillus_MOFA_variance_explained.pdf
",
    "  Bacillus_MOFA_variance_explained_with_total.pdf
",
    "  Bacillus_MOFA_factor_annotation.pdf
",
    "  Bacillus_MOFA_factors_by_condition.pdf
",
    "  Bacillus_MOFA_top_features/
",
    "  Bacillus_MOFA_data_patterns/
",
    "  Bacillus_MOFA_tables/
")

## =========================================================
## MOFA2 cross-omics link table
## - Builds pairwise and triplet links across views
## - Optional raw-data validation by Spearman correlation
## - No dplyr needed
## =========================================================

library(MOFA2)

make_mofa_link_table <- function(
    model,
    factors,
    views = c("Transcriptome", "Proteome", "Metabolome"),
    top_n = 15,
    min_abs_weight = NULL,
    same_sign_only = TRUE,
    add_raw_cor = TRUE,
    cor_method = "spearman",
    denoise = FALSE,
    max_pairs_per_block = 200,
    max_triplets_per_block = 200
) {
  
  ## ----------------------------
  ## 0. Normalize factor names
  ## ----------------------------
  if (is.numeric(factors)) {
    factors <- paste0("Factor", factors)
  }
  factors <- as.character(factors)
  
  ## ----------------------------
  ## 1. Pull MOFA weights
  ## ----------------------------
  W <- get_weights(model, as.data.frame = TRUE)
  W$factor <- as.character(W$factor)
  W$view   <- as.character(W$view)
  W$feature <- as.character(W$feature)
  W$abs_weight <- abs(W$value)
  W$sign <- ifelse(W$value >= 0, "positive", "negative")
  
  W <- W[W$factor %in% factors & W$view %in% views, , drop = FALSE]
  if (!is.null(min_abs_weight)) {
    W <- W[W$abs_weight >= min_abs_weight, , drop = FALSE]
  }
  
  if (nrow(W) == 0) {
    stop("No weights left after filtering. Check factor names / view names / thresholds.")
  }
  
  ## ----------------------------
  ## 2. Keep top features per factor x view x sign
  ## ----------------------------
  top_blocks <- list()
  idx <- 1
  
  for (f in unique(W$factor)) {
    for (v in views) {
      for (s in c("positive", "negative")) {
        tmp <- W[W$factor == f & W$view == v & W$sign == s, 
                 c("factor", "view", "feature", "value", "abs_weight", "sign"),
                 drop = FALSE]
        if (nrow(tmp) == 0) next
        tmp <- tmp[order(-tmp$abs_weight), , drop = FALSE]
        tmp <- head(tmp, top_n)
        top_blocks[[idx]] <- tmp
        idx <- idx + 1
      }
    }
  }
  
  topW <- do.call(rbind, top_blocks)
  rownames(topW) <- NULL
  
  ## ----------------------------
  ## 3. Build pairwise cross-view links
  ## ----------------------------
  pair_list <- list()
  idx <- 1
  
  for (f in unique(topW$factor)) {
    signs_to_use <- if (same_sign_only) c("positive", "negative") else unique(topW$sign)
    
    for (s in signs_to_use) {
      block <- topW[topW$factor == f & topW$sign == s, , drop = FALSE]
      present_views <- intersect(views, unique(block$view))
      if (length(present_views) < 2) next
      
      view_pairs <- combn(present_views, 2, simplify = FALSE)
      
      for (vp in view_pairs) {
        a <- block[block$view == vp[1], , drop = FALSE]
        b <- block[block$view == vp[2], , drop = FALSE]
        if (nrow(a) == 0 || nrow(b) == 0) next
        
        ab <- merge(a, b, by = NULL, suffixes = c("_1", "_2"))
        ab$pair_score <- ab$abs_weight_1 * ab$abs_weight_2
        ab <- ab[order(-ab$pair_score), , drop = FALSE]
        ab <- head(ab, max_pairs_per_block)
        
        out <- data.frame(
          factor     = ab$factor_1,
          sign       = ab$sign_1,
          view1      = ab$view_1,
          feature1   = ab$feature_1,
          weight1    = ab$value_1,
          abs_w1     = ab$abs_weight_1,
          view2      = ab$view_2,
          feature2   = ab$feature_2,
          weight2    = ab$value_2,
          abs_w2     = ab$abs_weight_2,
          pair_score = ab$pair_score,
          stringsAsFactors = FALSE
        )
        
        pair_list[[idx]] <- out
        idx <- idx + 1
      }
    }
  }
  
  pair_links <- if (length(pair_list) > 0) do.call(rbind, pair_list) else data.frame()
  if (nrow(pair_links) > 0) {
    pair_links <- pair_links[order(pair_links$factor, -pair_links$pair_score), , drop = FALSE]
    rownames(pair_links) <- NULL
  }
  
  ## ----------------------------
  ## 4. Build 3-view triplets
  ## ----------------------------
  triplet_links <- data.frame()
  
  if (length(views) >= 3) {
    triplet_list <- list()
    idx <- 1
    
    for (f in unique(topW$factor)) {
      signs_to_use <- if (same_sign_only) c("positive", "negative") else unique(topW$sign)
      
      for (s in signs_to_use) {
        block <- topW[topW$factor == f & topW$sign == s, , drop = FALSE]
        present_views <- intersect(views, unique(block$view))
        if (length(present_views) < 3) next
        
        a <- block[block$view == present_views[1], , drop = FALSE]
        b <- block[block$view == present_views[2], , drop = FALSE]
        c <- block[block$view == present_views[3], , drop = FALSE]
        if (nrow(a) == 0 || nrow(b) == 0 || nrow(c) == 0) next
        
        ab  <- merge(a, b, by = NULL, suffixes = c("_1", "_2"))
        abc <- merge(ab, c, by = NULL)
        
        abc$triplet_score <- abc$abs_weight_1 * abc$abs_weight_2 * abc$abs_weight
        abc <- abc[order(-abc$triplet_score), , drop = FALSE]
        abc <- head(abc, max_triplets_per_block)
        
        out <- data.frame(
          factor        = abc$factor_1,
          sign          = abc$sign_1,
          view1         = abc$view_1,
          feature1      = abc$feature_1,
          weight1       = abc$value_1,
          abs_w1        = abc$abs_weight_1,
          view2         = abc$view_2,
          feature2      = abc$feature_2,
          weight2       = abc$value_2,
          abs_w2        = abc$abs_weight_2,
          view3         = abc$view,
          feature3      = abc$feature,
          weight3       = abc$value,
          abs_w3        = abc$abs_weight,
          triplet_score = abc$triplet_score,
          stringsAsFactors = FALSE
        )
        
        triplet_list[[idx]] <- out
        idx <- idx + 1
      }
    }
    
    if (length(triplet_list) > 0) {
      triplet_links <- do.call(rbind, triplet_list)
      triplet_links <- triplet_links[order(triplet_links$factor, -triplet_links$triplet_score), , drop = FALSE]
      rownames(triplet_links) <- NULL
    }
  }
  
  ## ----------------------------
  ## 5. Optional: validate pair links in the original data
  ## get_data(..., as.data.frame=TRUE) returns:
  ## (view, feature, sample, value)
  ## ----------------------------
  if (add_raw_cor && nrow(pair_links) > 0) {
    
    feat_list <- list()
    all_views <- unique(c(pair_links$view1, pair_links$view2))
    
    for (v in all_views) {
      f1 <- pair_links$feature1[pair_links$view1 == v]
      f2 <- pair_links$feature2[pair_links$view2 == v]
      feat_list[[v]] <- unique(c(f1, f2))
    }
    
    ## fetch each view separately to avoid MOFA2 multi-view get_data() bug
    feat_list <- feat_list[vapply(feat_list, length, integer(1)) > 0]
    
    D_list <- lapply(names(feat_list), function(v) {
      get_data(
        model,
        views = v,
        features = setNames(list(feat_list[[v]]), v),
        as.data.frame = TRUE,
        add_intercept = FALSE,
        denoise = denoise,
        na.rm = TRUE
      )
    })
    
    D <- do.call(rbind, D_list)
    rownames(D) <- NULL
    
    D$key <- paste(D$view, D$feature, sep = "||")
    splitD <- split(D[, c("sample", "value"), drop = FALSE], D$key)
    
    raw_cor   <- rep(NA_real_, nrow(pair_links))
    raw_p     <- rep(NA_real_, nrow(pair_links))
    n_shared  <- rep(0L, nrow(pair_links))
    
    for (i in seq_len(nrow(pair_links))) {
      k1 <- paste(pair_links$view1[i], pair_links$feature1[i], sep = "||")
      k2 <- paste(pair_links$view2[i], pair_links$feature2[i], sep = "||")
      
      d1 <- splitD[[k1]]
      d2 <- splitD[[k2]]
      if (is.null(d1) || is.null(d2)) next
      
      m <- merge(d1, d2, by = "sample", suffixes = c("_1", "_2"))
      n_shared[i] <- nrow(m)
      
      if (nrow(m) >= 3 &&
          stats::sd(m$value_1, na.rm = TRUE) > 0 &&
          stats::sd(m$value_2, na.rm = TRUE) > 0) {
        ct <- suppressWarnings(stats::cor.test(m$value_1, m$value_2, method = cor_method))
        raw_cor[i] <- unname(ct$estimate)
        raw_p[i]   <- ct$p.value
      }
    }
    
    pair_links$n_shared_samples <- n_shared
    pair_links$raw_cor <- raw_cor
    pair_links$raw_cor_p <- raw_p
    
    ## combined score: MOFA linkage + raw correlation support
    pair_links$combined_score <- with(
      pair_links,
      pair_score * ifelse(is.na(raw_cor), 0, abs(raw_cor))
    )
    
    pair_links <- pair_links[order(pair_links$factor, -pair_links$combined_score, -pair_links$pair_score), , drop = FALSE]
    rownames(pair_links) <- NULL
  }
  
  ## ----------------------------
  ## 6. Return results
  ## ----------------------------
  list(
    top_weights   = topW,
    pair_links    = pair_links,
    triplet_links = triplet_links
  )
}

## choose factors you care about
## examples:
# factors_to_use <- c("Factor1", "Factor3", "Factor5")
## or
factors_to_use <- key_factors

res_links <- make_mofa_link_table(
  model = mofa_trained,
  factors = key_factors,   # or c(1,3,5)
  views = c("Transcriptome", "Proteome", "Metabolome"),
  top_n = 15,
  same_sign_only = TRUE,
  add_raw_cor = TRUE,
  cor_method = "spearman",
  denoise = FALSE
)

## top weighted features used to build the links
head(res_links$top_weights, 20)

## cross-omics pairs
head(res_links$pair_links, 30)

## 3-way links
head(res_links$triplet_links, 30)

## save
write.csv(res_links$top_weights,   "Bacillus_MOFA_top_weights_used_for_links.csv", row.names = FALSE)
write.csv(res_links$pair_links,    "Bacillus_MOFA_pair_links.csv",                 row.names = FALSE)
write.csv(res_links$triplet_links, "Bacillus_MOFA_triplet_links.csv",              row.names = FALSE)

library(ggplot2)

tw <- read.csv("Bacillus_MOFA_top_weights_used_for_links.csv", stringsAsFactors = FALSE)
tw$signed_weight <- ifelse(tw$sign == "negative", -tw$abs_weight, tw$abs_weight)

plot_one_factor <- function(df, fac = "Factor5") {
  sub <- df[df$factor == fac, , drop = FALSE]
  
  pieces <- lapply(split(sub, sub$view), function(x) {
    x <- x[order(x$signed_weight), , drop = FALSE]
    x$feature_plot <- factor(x$feature, levels = x$feature)
    x
  })
  sub2 <- do.call(rbind, pieces)
  
  ggplot(sub2, aes(x = signed_weight, y = feature_plot)) +
    geom_segment(aes(x = 0, xend = signed_weight, yend = feature_plot)) +
    geom_point(aes(color = view), size = 2.5) +
    geom_vline(xintercept = 0, linetype = 2) +
    facet_wrap(~ view, scales = "free_y", ncol = 1) +
    labs(
      x = "Signed MOFA weight",
      y = NULL,
      title = paste("Top weights for", fac),
      subtitle = "Positive side = features higher in samples with positive factor values"
    ) +
    theme_bw()
}

p <- plot_one_factor(tw, "Factor3")
print(p)
ggsave("Factor3_top_weights_lollipop.pdf", p, width = 8, height = 10)

#!/usr/bin/env Rscript

# ============================================================
# Full end-to-end WGCNA (RNA-seq transcriptome) + metabolite linking
# - Inputs:
#   1) Expression matrix: genes x samples (normalized log-like values: VST/rlog/logCPM)
#   2) Metabolite matrix: metabolites x samples (abundance)
#   3) Metadata table: samples x covariates (must include a sample ID column)
#
# - Outputs:
#   * WGCNA modules, eigengenes
#   * Module↔metabolite correlation tables (raw + optional residualized)
#   * Heatmaps
#   * Hub genes for top module-metabolite links
#   * Optional Cytoscape exports for selected modules
#
# Notes:
# - WGCNA should NOT be run on raw counts. Use DESeq2::vst / rlog or limma-voom logCPM.
# - If you only have raw counts, see the OPTIONAL section near the bottom.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(matrixStats)
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ---------------------------
# User-configurable parameters
# ---------------------------
CONFIG <- list(
  expr_file      = "Rhodano_mRNA_seq_vst_011826.csv",       # genes x samples (first column = gene IDs)
  met_file       = "Rhodano_metabolites_011826.csv",     # metabolites x samples (first column = metabolite IDs)
  meta_file      = "~/Rhodano_mRNAseq/data/meta/sample_info.csv",        # samples x covariates (must include sample ID column)
  meta_sample_id = "sample",            # name of sample ID column in metadata
  
  out_dir        = "Rhodano_wgcna_out",
  
  # Transcriptome preprocessing
  gene_var_quantile = 0.75,   # keep top (1 - quantile) variable genes. 0.75 => top 25%
  max_missing_frac  = 0.2,    # remove genes/samples with too many NAs (WGCNA will also check)
  sample_outlier_cut_height = NA,  # set numeric to drop outlier samples; NA disables auto-drop
  
  # Metabolomics preprocessing
  met_log_transform = TRUE,  # log10(x + 1)
  met_zscore        = TRUE,  # z-score each metabolite across samples
  
  # WGCNA settings
  networkType     = "signed", # "signed" is often preferred for biology
  corType         = "bicor",  # robust correlation
  minModuleSize   = 30,
  mergeCutHeight  = 0.25,
  power_override  = NA,       # if NA, will pick based on scale-free fit; else use integer
  
  # Module↔metabolite association thresholds (for reporting)
  assoc_fdr_thresh = 0.05,
  assoc_cor_thresh = 0.5,     # for raw correlations
  assoc_cor_thresh_resid = 0.4, # for residualized correlations
  
  # Residualization (covariate adjustment) using limma
  do_residualize = TRUE,
  # Put a model formula here using metadata column names.
  # Example: "~ 0 + Condition + Timepoint + Batch"
  # If you want an intercept model: "~ Condition + Timepoint + Batch"
  resid_formula = "~ 0 + condition + replicate",
  
  # Hub gene extraction
  # For each significant module-met pair, take top N hub genes (high |kME| and high |GS|)
  top_hubs_per_pair = 30,
  
  # Cytoscape export (optional)
  do_cytoscape = TRUE,
  cytoscape_modules = NULL,          # NULL = export all modules
  cytoscape_min_module_size = 20,    # skip tiny modules
  cytoscape_tom_threshold = 0.10
)

# ---------------------------
# Helper functions
# ---------------------------

dir.create(CONFIG$out_dir, recursive = TRUE, showWarnings = FALSE)

logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf(...), "\n")
}

stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing file(s): ", paste(missing, collapse=", "))
}

read_matrix_tsv <- function(path) {
  df <- fread(path) |> as.data.frame()
  if (ncol(df) < 2) stop("File has <2 columns: ", path)
  rownames(df) <- df[[1]]
  df <- df[,-1, drop=FALSE]
  return(df)
}

write_tsv <- function(df, path) {
  fwrite(as.data.table(df), file = path, sep = "\t")
}

save_plot_pdf <- function(path, width=8, height=6, expr) {
  pdf(path, width=width, height=height)
  on.exit(dev.off(), add=TRUE)
  force(expr)
}

scale_free_pick_power <- function(datExpr, networkType="signed", corType="bicor") {
  powers <- c(1:20, seq(22, 30, 2))
  sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = networkType,
    corFnc = if (corType == "bicor") "bicor" else "cor",
    verbose = 5
  )
  
  # Choose the smallest power with signed R^2 >= 0.85 if possible,
  # else choose power at the max R^2 (or elbow-ish).
  fit <- sft$fitIndices
  # fit columns: Power, SFT.R.sq, slope, truncated R^2, mean.k, median.k, max.k
  # WGCNA uses: -sign(slope)*SFT.R.sq (commonly plotted). We'll follow that convention.
  signedR2 <- -sign(fit[,3]) * fit[,2]
  ok <- which(signedR2 >= 0.85)
  if (length(ok)) {
    chosen <- fit[min(ok), 1]
  } else {
    chosen <- fit[which.max(signedR2), 1]
  }
  list(sft=sft, chosen_power=chosen, signedR2=signedR2, powers=powers)
}

residualize_matrix <- function(mat_samples_x_features, meta, formula_str) {
  # Residualize each feature (column) using limma::lmFit for speed.
  # mat: samples x features
  mm <- model.matrix(as.formula(formula_str), data = meta)
  if (nrow(mm) != nrow(mat_samples_x_features)) {
    stop("Model matrix rows do not match sample count.")
  }
  fit <- lmFit(t(mat_samples_x_features), mm)  # limma expects features x samples, so transpose
  coef <- fit$coefficients                     # features x covariates
  fitted <- coef %*% t(mm)                     # features x samples
  resid <- t(mat_samples_x_features) - fitted  # features x samples
  t(resid)                                     # back to samples x features
}

safe_bicor <- function(x, y) {
  # x: samples x p, y: samples x q
  bicor(x, y, use="pairwise.complete.obs")
}

cor_and_fdr <- function(corMat, nSamples) {
  pMat <- corPvalueStudent(corMat, nSamples = nSamples)
  fdr <- matrix(p.adjust(as.vector(pMat), method="fdr"), nrow=nrow(pMat),
                dimnames=dimnames(pMat))
  list(pMat=pMat, fdr=fdr)
}

# ---------------------------
# Main
# ---------------------------

stop_if_missing(c(CONFIG$expr_file, CONFIG$met_file, CONFIG$meta_file))

logmsg("Reading expression: %s", CONFIG$expr_file)
expr <- read_matrix_tsv(CONFIG$expr_file)   # genes x samples

logmsg("Reading metabolites: %s", CONFIG$met_file)
met  <- read_matrix_tsv(CONFIG$met_file)    # metabolites x samples

logmsg("Reading metadata: %s", CONFIG$meta_file)
meta <- fread(CONFIG$meta_file) |> as.data.frame()

if (!CONFIG$meta_sample_id %in% colnames(meta)) {
  stop("Metadata is missing sample ID column: ", CONFIG$meta_sample_id)
}
rownames(meta) <- meta[[CONFIG$meta_sample_id]]

# Align samples across expr/met/meta
common_samples <- Reduce(intersect, list(colnames(expr), colnames(met), rownames(meta)))
if (length(common_samples) < 6) {
  stop("Too few common samples after alignment (need >= 6 ideally). Found: ", length(common_samples))
}
expr <- expr[, common_samples, drop=FALSE]
met  <- met[,  common_samples, drop=FALSE]
meta <- meta[ common_samples, , drop=FALSE]

stopifnot(all(colnames(expr) == colnames(met)))
stopifnot(all(colnames(expr) == rownames(meta)))

logmsg("Aligned samples: %d", length(common_samples))
logmsg("Expr dims (genes x samples): %d x %d", nrow(expr), ncol(expr))
logmsg("Met dims  (mets  x samples): %d x %d", nrow(met),  ncol(met))

# ---------------------------
# Preprocess expression
# ---------------------------
logmsg("Preprocessing expression for WGCNA (samples x genes)")
datExpr0 <- t(as.matrix(expr))  # samples x genes

# Missingness filter (simple)
gene_na_frac <- colMeans(is.na(datExpr0))
samp_na_frac <- rowMeans(is.na(datExpr0))
keep_genes <- gene_na_frac <= CONFIG$max_missing_frac
keep_samps <- samp_na_frac <= CONFIG$max_missing_frac
datExpr0 <- datExpr0[keep_samps, keep_genes, drop=FALSE]

# WGCNA QC
gsg <- goodSamplesGenes(datExpr0, verbose=3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop=FALSE]
}

# Variance filter
geneVar <- colVars(datExpr0)
cutoff <- quantile(geneVar, CONFIG$gene_var_quantile, na.rm=TRUE)
datExpr <- datExpr0[, geneVar >= cutoff, drop=FALSE]

logmsg("After QC + variance filter: %d samples x %d genes", nrow(datExpr), ncol(datExpr))

# Update met/meta to match any sample removal
keep_samples <- rownames(datExpr)
met  <- met[, keep_samples, drop=FALSE]
meta <- meta[keep_samples, , drop=FALSE]

# Sample clustering for outliers
sampleTree <- hclust(dist(datExpr), method="average")
save_plot_pdf(file.path(CONFIG$out_dir, "sample_clustering.pdf"), 10, 6, {
  plot(sampleTree, main="Sample clustering (expression)", xlab="", sub="", cex=0.8)
  if (!is.na(CONFIG$sample_outlier_cut_height)) {
    abline(h=CONFIG$sample_outlier_cut_height, col="red", lty=2)
  }
})

if (!is.na(CONFIG$sample_outlier_cut_height)) {
  clust <- cutreeStatic(sampleTree, cutHeight=CONFIG$sample_outlier_cut_height, minSize=10)
  keep <- (clust == 1)
  datExpr <- datExpr[keep, , drop=FALSE]
  keep_samples <- rownames(datExpr)
  met  <- met[, keep_samples, drop=FALSE]
  meta <- meta[keep_samples, , drop=FALSE]
  logmsg("After outlier removal: %d samples x %d genes", nrow(datExpr), ncol(datExpr))
}

# ---------------------------
# Preprocess metabolomics
# ---------------------------
logmsg("Preprocessing metabolites (samples x metabolites)")
datMet0 <- t(as.matrix(met))  # samples x metabolites

if (CONFIG$met_log_transform) {
  datMet0 <- log10(datMet0 + 1)
}

datMet <- datMet0
if (CONFIG$met_zscore) {
  datMet <- scale(datMet)
}

# ---------------------------
# Pick soft threshold power (or use override)
# ---------------------------
softPower <- CONFIG$power_override
if (is.na(softPower)) {
  logmsg("Picking soft-threshold power via scale-free topology fit...")
  pick <- scale_free_pick_power(datExpr, networkType=CONFIG$networkType, corType=CONFIG$corType)
  softPower <- pick$chosen_power
  
  # Save diagnostic plots
  sft <- pick$sft
  fit <- sft$fitIndices
  signedR2 <- pick$signedR2
  powers <- pick$powers
  
  save_plot_pdf(file.path(CONFIG$out_dir, "soft_threshold_diagnostics.pdf"), 11, 5, {
    par(mfrow=c(1,2))
    plot(fit[,1], signedR2, xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit, signed R^2",
         type="n", main="Scale independence")
    text(fit[,1], signedR2, labels=powers, cex=0.7, col="red")
    abline(h=0.85, col="blue", lty=2)
    abline(v=softPower, col="darkgreen", lty=2)
    
    plot(fit[,1], fit[,5], xlab="Soft Threshold (power)",
         ylab="Mean Connectivity", type="n", main="Mean connectivity")
    text(fit[,1], fit[,5], labels=powers, cex=0.7, col="red")
    abline(v=softPower, col="darkgreen", lty=2)
  })
  
} else {
  logmsg("Using user-provided softPower = %d", softPower)
}
logmsg("Selected softPower = %d", softPower)

# ---------------------------
# Build WGCNA modules
# ---------------------------
logmsg("Running blockwiseModules (this can take time)...")
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = CONFIG$networkType,
  corType = if (CONFIG$corType == "bicor") "bicor" else "pearson",
  minModuleSize = CONFIG$minModuleSize,
  mergeCutHeight = CONFIG$mergeCutHeight,
  reassignThreshold = 0,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(CONFIG$out_dir, "TOM"),
  numericLabels = TRUE,
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
MEs0 <- net$MEs
MEs <- orderMEs(MEs0)

logmsg("Modules found (incl. grey):")
logmsg("%s", paste(names(sort(table(moduleColors), decreasing=TRUE)), collapse=", "))

# Save module assignment table
gene_ids <- colnames(datExpr)
gene2module <- data.frame(Gene=gene_ids, Module=moduleColors, stringsAsFactors=FALSE)
write_tsv(gene2module, file.path(CONFIG$out_dir, "gene_module_assignments.tsv"))

# Plot dendrogram + modules for first block
save_plot_pdf(file.path(CONFIG$out_dir, "gene_dendrogram_modules.pdf"), 12, 6, {
  plotDendroAndColors(
    net$dendrograms[[1]],
    moduleColors[net$blockGenes[[1]]],
    groupLabels = "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05
  )
})

# Save eigengenes
write_tsv(data.frame(SampleID=rownames(MEs), MEs, check.names=FALSE),
          file.path(CONFIG$out_dir, "module_eigengenes.tsv"))

# ---------------------------
# Module ↔ Metabolite correlations (raw)
# ---------------------------
logmsg("Computing module eigengene ↔ metabolite correlations (raw)")
stopifnot(rownames(MEs) == rownames(datExpr))
stopifnot(rownames(datMet) == rownames(datExpr))

corMat <- safe_bicor(MEs, datMet)
stats <- cor_and_fdr(corMat, nSamples=nrow(datExpr))
pMat <- stats$pMat
fdrMat <- stats$fdr

# Save matrices
write_tsv(as.data.frame(corMat, check.names=FALSE), file.path(CONFIG$out_dir, "ME_met_cor_raw.tsv"))
write_tsv(as.data.frame(pMat,  check.names=FALSE), file.path(CONFIG$out_dir, "ME_met_p_raw.tsv"))
write_tsv(as.data.frame(fdrMat,check.names=FALSE), file.path(CONFIG$out_dir, "ME_met_fdr_raw.tsv"))

# ------------------------------------------------------------
# Heatmap (ROBUST v2): drop bad metabolite columns first, then rows
# ------------------------------------------------------------

# 1) Remove metabolite columns that contain any non-finite correlations
bad_cols <- colSums(!is.finite(corMat)) > 0

if (any(bad_cols)) {
  message(sprintf(
    "Heatmap: removing %d metabolites (columns) with NA/Inf correlations.",
    sum(bad_cols)
  ))
}

corMat_hm <- corMat[, !bad_cols, drop = FALSE]

# 2) Now re-check rows AFTER removing bad columns
bad_rows <- rowSums(!is.finite(corMat_hm)) > 0
if (any(bad_rows)) {
  message(sprintf(
    "Heatmap: removing %d modules (rows) still containing NA/Inf after column filtering.",
    sum(bad_rows)
  ))
  corMat_hm <- corMat_hm[!bad_rows, , drop = FALSE]
}

# 3) Final safeguard: replace any remaining non-finite values with 0
# (should be rare after filtering)
corMat_hm[!is.finite(corMat_hm)] <- 0

# 4) Only plot if something remains
if (nrow(corMat_hm) == 0 || ncol(corMat_hm) == 0) {
  warning(sprintf("Heatmap skipped: corMat_hm is %d x %d after filtering.",
                  nrow(corMat_hm), ncol(corMat_hm)))
} else {
  save_plot_pdf(file.path(CONFIG$out_dir, "ME_met_heatmap_raw.pdf"), 12, 10, {
    pheatmap(corMat_hm,
             main = "Module eigengenes vs metabolites (bicor, raw)",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             fontsize_row = 8,
             fontsize_col = 6)
  })
}



# Significant pairs table (raw)
sig_idx <- which(fdrMat < CONFIG$assoc_fdr_thresh & abs(corMat) >= CONFIG$assoc_cor_thresh, arr.ind=TRUE)
sig_raw <- data.frame(
  Module = rownames(corMat)[sig_idx[,1]],
  Metabolite = colnames(corMat)[sig_idx[,2]],
  bicor = corMat[sig_idx],
  p = pMat[sig_idx],
  fdr = fdrMat[sig_idx],
  stringsAsFactors = FALSE
)
sig_raw <- sig_raw[order(sig_raw$fdr, -abs(sig_raw$bicor)), ]
write_tsv(sig_raw, file.path(CONFIG$out_dir, "ME_met_significant_raw.tsv"))
logmsg("Significant raw module-met pairs: %d", nrow(sig_raw))

# ---------------------------
# Optional residualization and correlations (recommended for complex designs)
# ---------------------------
# --- Residualize on condition + replicate (metadata columns) ---
# This version assumes you have BOTH:
#   meta$condition  (or whatever you name it)
#   meta$replicate  (replicate ID / batch-like factor)
#
# It will residualize MEs and metabolites using: ~ 0 + condition + replicate

sig_resid <- data.frame()

if (CONFIG$do_residualize) {
  
  meta2 <- meta[rownames(datExpr), , drop = FALSE]
  
  # ---- set your metadata column names here ----
  cond_col <- "condition"
  rep_col  <- "replicate"
  
  # ---- validate presence ----
  missing_cols <- setdiff(c(cond_col, rep_col), colnames(meta2))
  if (length(missing_cols) > 0) {
    stop("Residualization failed: metadata column(s) not found: ",
         paste(missing_cols, collapse=", "),
         "\nAvailable columns: ", paste(colnames(meta2), collapse=", "))
  }
  
  # ---- coerce to factors (important for model.matrix) ----
  meta2[[cond_col]] <- as.factor(meta2[[cond_col]])
  meta2[[rep_col]]  <- as.factor(meta2[[rep_col]])
  
  # ---- build formula ----
  CONFIG$resid_formula <- paste0("~ 0 + ", cond_col, " + ", rep_col)
  logmsg("Residualizing MEs and metabolites using formula: %s", CONFIG$resid_formula)
  
  # ---- residualize ----
  MEs_res <- residualize_matrix(MEs, meta2, CONFIG$resid_formula)
  colnames(MEs_res) <- colnames(MEs); rownames(MEs_res) <- rownames(MEs)
  
  Met_res <- residualize_matrix(datMet, meta2, CONFIG$resid_formula)
  colnames(Met_res) <- colnames(datMet); rownames(Met_res) <- rownames(datMet)
  
  # ---- correlate residuals ----
  corMat_r <- safe_bicor(MEs_res, Met_res)
  stats_r  <- cor_and_fdr(corMat_r, nSamples = nrow(datExpr))
  pMat_r   <- stats_r$pMat
  fdrMat_r <- stats_r$fdr
  
  write_tsv(as.data.frame(corMat_r, check.names=FALSE),
            file.path(CONFIG$out_dir, "ME_met_cor_resid.tsv"))
  write_tsv(as.data.frame(pMat_r, check.names=FALSE),
            file.path(CONFIG$out_dir, "ME_met_p_resid.tsv"))
  write_tsv(as.data.frame(fdrMat_r, check.names=FALSE),
            file.path(CONFIG$out_dir, "ME_met_fdr_resid.tsv"))
  
  # ---- robust heatmap (drop bad metabolite columns first) ----
  bad_cols <- colSums(!is.finite(corMat_r)) > 0
  corMat_r_hm <- corMat_r[, !bad_cols, drop=FALSE]
  bad_rows <- rowSums(!is.finite(corMat_r_hm)) > 0
  corMat_r_hm <- corMat_r_hm[!bad_rows, , drop=FALSE]
  corMat_r_hm[!is.finite(corMat_r_hm)] <- 0
  
  if (nrow(corMat_r_hm) > 0 && ncol(corMat_r_hm) > 0) {
    save_plot_pdf(file.path(CONFIG$out_dir, "ME_met_heatmap_resid.pdf"), 12, 10, {
      pheatmap(corMat_r_hm,
               main="Module eigengenes vs metabolites (bicor, residualized on condition+replicate)",
               clustering_distance_rows="correlation",
               clustering_distance_cols="correlation",
               fontsize_row=8, fontsize_col=6)
    })
  } else {
    warning(sprintf("Residualized heatmap skipped: matrix is %d x %d after filtering.",
                    nrow(corMat_r_hm), ncol(corMat_r_hm)))
  }
  
  # ---- significant pairs ----
  sig_idx_r <- which(fdrMat_r < CONFIG$assoc_fdr_thresh &
                       abs(corMat_r) >= CONFIG$assoc_cor_thresh_resid, arr.ind=TRUE)
  
  sig_resid <- data.frame(
    Module = rownames(corMat_r)[sig_idx_r[,1]],
    Metabolite = colnames(corMat_r)[sig_idx_r[,2]],
    bicor = corMat_r[sig_idx_r],
    p = pMat_r[sig_idx_r],
    fdr = fdrMat_r[sig_idx_r],
    stringsAsFactors = FALSE
  )
  sig_resid <- sig_resid[order(sig_resid$fdr, -abs(sig_resid$bicor)), ]
  write_tsv(sig_resid, file.path(CONFIG$out_dir, "ME_met_significant_resid.tsv"))
  logmsg("Significant residualized module-met pairs: %d", nrow(sig_resid))
}


# ---------------------------
# Hub gene extraction for each significant pair
# ---------------------------
# ---------------------------
# Hub gene extraction (FIXED for numeric MEs like ME0/ME1)
# ---------------------------

logmsg("Computing module membership (kME) for all genes")
kME  <- as.data.frame(bicor(datExpr, MEs, use="pairwise.complete.obs"))
kMEp <- as.data.frame(corPvalueStudent(as.matrix(kME), nSamples=nrow(datExpr)))

colnames(kME)  <- paste0("kME_",  colnames(MEs))   # e.g., kME_ME0, kME_MEturquoise
colnames(kMEp) <- paste0("kMEp_", colnames(MEs))

geneInfo_base <- data.frame(
  Gene = colnames(datExpr),
  ModuleColor = moduleColors,   # color assignment per gene
  stringsAsFactors = FALSE
)
geneInfo_base <- cbind(geneInfo_base, kME, kMEp)

# Choose which significant list to use
pairs_for_hubs <- sig_resid
pairs_label <- "resid"
if (nrow(pairs_for_hubs) == 0) {
  pairs_for_hubs <- sig_raw
  pairs_label <- "raw"
}

if (nrow(pairs_for_hubs) == 0) {
  logmsg("No significant module-metabolite pairs (raw or resid); skipping hub gene tables.")
} else {
  
  hub_dir <- file.path(CONFIG$out_dir, paste0("hub_genes_", pairs_label))
  dir.create(hub_dir, showWarnings=FALSE, recursive=TRUE)
  
  # Helper: convert ME name -> module color string used in moduleColors
  ME_to_color <- function(me_name) {
    # me_name can be "ME0" or "MEturquoise"
    if (grepl("^ME\\d+$", me_name)) {
      lab <- as.numeric(sub("^ME", "", me_name))
      return(labels2colors(lab))   # 0->grey, 1->turquoise, etc.
    }
    if (grepl("^ME", me_name)) {
      return(sub("^ME", "", me_name))  # "MEturquoise" -> "turquoise"
    }
    # fallback: already a color?
    return(me_name)
  }
  
  # If you used residualized pairs, use residualized metabolite values for GS
  Met_for_GS <- datMet
  if (pairs_label == "resid" && CONFIG$do_residualize) {
    meta2 <- meta[rownames(datExpr), , drop=FALSE]
    # Use your existing CONFIG$resid_formula (condition+replicate, etc.)
    Met_for_GS <- residualize_matrix(datMet, meta2, CONFIG$resid_formula)
    colnames(Met_for_GS) <- colnames(datMet)
    rownames(Met_for_GS) <- rownames(datMet)
  }
  
  logmsg("Extracting hub genes for %d %s module-met pairs", nrow(pairs_for_hubs), pairs_label)
  
  for (i in seq_len(nrow(pairs_for_hubs))) {
    
    me_name <- pairs_for_hubs$Module[i]        # e.g., "ME0" or "MEturquoise"
    metName <- pairs_for_hubs$Metabolite[i]
    if (!metName %in% colnames(Met_for_GS)) next
    
    mod_color <- ME_to_color(me_name)          # convert ME0 -> turquoise, etc.
    if (!mod_color %in% unique(moduleColors)) next
    
    # Genes in this module (by color)
    inMod <- (geneInfo_base$ModuleColor == mod_color)
    if (sum(inMod) == 0) next
    
    # Gene significance for this metabolite
    y <- as.numeric(Met_for_GS[, metName])
    GS  <- as.numeric(bicor(datExpr, y, use="pairwise.complete.obs"))
    GSp <- as.numeric(corPvalueStudent(GS, nSamples=nrow(datExpr)))
    
    geneInfo <- geneInfo_base
    geneInfo$GS <- GS
    geneInfo$GSp <- GSp
    
    # Correct kME column name (based on ME column naming)
    kme_col <- paste0("kME_", me_name)         # e.g., kME_ME0
    if (!kme_col %in% colnames(geneInfo)) next
    
    hubs <- geneInfo[inMod, , drop=FALSE]
    hubs <- hubs[order(-abs(hubs[[kme_col]]), -abs(hubs$GS)), , drop=FALSE]
    hubs_top <- head(hubs, CONFIG$top_hubs_per_pair)
    
    safe_met <- gsub("[^A-Za-z0-9_.-]+", "_", metName)
    out_path <- file.path(hub_dir, paste0("hubs_", mod_color, "__", safe_met, ".tsv"))
    write_tsv(hubs_top, out_path)
    
    # GS vs kME plot
    plot_path <- file.path(hub_dir, paste0("GS_vs_kME_", mod_color, "__", safe_met, ".pdf"))
    save_plot_pdf(plot_path, 7, 6, {
      plot(abs(hubs[[kme_col]]), abs(hubs$GS),
           xlab=paste0("Module membership |", kme_col, "|"),
           ylab=paste0("Gene significance |bicor(gene, ", metName, ")|"),
           main=paste0("Hub structure: ", mod_color, " vs ", metName),
           pch=16)
    })
  }
  
  logmsg("Hub gene extraction finished. Output: %s", hub_dir)
}


# ---------------------------
# Optional Cytoscape export
# ---------------------------
if (CONFIG$do_cytoscape) {
  
  out_dir <- file.path(CONFIG$out_dir, "cytoscape")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Decide which modules to export
  if (is.null(CONFIG$cytoscape_modules) || length(CONFIG$cytoscape_modules) == 0) {
    mods <- sort(unique(moduleColors))
    mods <- setdiff(mods, c("grey", "gold"))  # grey = unassigned; 'gold' sometimes appears, optional
  } else if (length(CONFIG$cytoscape_modules) == 1 && toupper(CONFIG$cytoscape_modules[1]) == "ALL") {
    mods <- sort(unique(moduleColors))
    mods <- setdiff(mods, "grey")
  } else {
    mods <- CONFIG$cytoscape_modules
    mods <- mods[mods != ""]
  }
  
  # Module sizes
  mod_sizes <- table(moduleColors)
  mods <- mods[mods %in% names(mod_sizes) & mod_sizes[mods] >= CONFIG$cytoscape_min_module_size]
  
  logmsg("Cytoscape export: exporting %d modules (min size %d, TOM threshold %.2f)",
         length(mods), CONFIG$cytoscape_min_module_size, CONFIG$cytoscape_tom_threshold)
  
  # Iterate modules and export TOM edges/nodes
  for (mod_col in mods) {
    
    inModule <- (moduleColors == mod_col)
    modGenes <- colnames(datExpr)[inModule]
    nGenes <- length(modGenes)
    
    if (nGenes < CONFIG$cytoscape_min_module_size) next
    
    logmsg("Cytoscape export: module %s (n=%d)", mod_col, nGenes)
    
    # Compute TOM for this module only (keeps memory manageable)
    datExpr_mod <- datExpr[, inModule, drop = FALSE]
    
    TOM_mod <- TOMsimilarityFromExpr(datExpr_mod,
                                     power = softPower,
                                     networkType = "signed",
                                     corType = "bicor")
    dimnames(TOM_mod) <- list(modGenes, modGenes)
    
    edgeFile <- file.path(out_dir, paste0("CytoscapeInput-edges-", mod_col, ".txt"))
    nodeFile <- file.path(out_dir, paste0("CytoscapeInput-nodes-", mod_col, ".txt"))
    
    exportNetworkToCytoscape(
      TOM_mod,
      edgeFile = edgeFile,
      nodeFile = nodeFile,
      weighted = TRUE,
      threshold = CONFIG$cytoscape_tom_threshold,
      nodeNames = modGenes,
      nodeAttr = rep(mod_col, length(modGenes))
    )
  }
  
  logmsg("Cytoscape export finished: %s", out_dir)
}


# ---------------------------
# Save workspace objects for reproducibility
# ---------------------------
saveRDS(list(
  config=CONFIG,
  datExpr=datExpr,
  datMet=datMet,
  meta=meta,
  net=net,
  moduleColors=moduleColors,
  MEs=MEs,
  sig_raw=sig_raw,
  sig_resid=sig_resid
), file = file.path(CONFIG$out_dir, "wgcna_results.rds"))

logmsg("Done. Outputs written to: %s", CONFIG$out_dir)

# ============================================================
# OPTIONAL: If you only have raw counts (genes x samples), normalize first:
#
# 1) DESeq2 VST (recommended)
#   library(DESeq2)
#   counts <- read_matrix_tsv("counts.tsv")  # genes x samples
#   coldata <- fread("metadata.tsv") |> as.data.frame()
#   rownames(coldata) <- coldata$SampleID
#   coldata <- coldata[colnames(counts), , drop=FALSE]
#   dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(counts)),
#                                colData=coldata,
#                                design=~ 1)
#   dds <- estimateSizeFactors(dds)
#   vsd <- vst(dds, blind=TRUE)
#   expr_vst <- assay(vsd)  # genes x samples
#   # write out as TSV with gene IDs as first column if needed
# ============================================================

# ============================================================
# FULL WORKING CODE: Build Cytoscape network (modules + genes + metabolites)
# Fixes rbindlist errors caused by duplicate kME columns / weird kME types
# Outputs:
#   <out_prefix>_nodes.tsv
#   <out_prefix>_edges.tsv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -------------------------
# 1) Helper: make a clean kME table (Gene + one kME_<ME> col per module)
# -------------------------
make_kME_clean <- function(kME_df, modules) {
  stopifnot("Gene" %in% colnames(kME_df))
  modules <- unique(as.character(modules))
  kcols <- paste0("kME_", modules)
  
  # For each kME column name, take the FIRST matching column index
  idx <- sapply(kcols, function(nm) {
    w <- which(colnames(kME_df) == nm)
    if (length(w) == 0) NA_integer_ else w[1]
  })
  
  if (any(is.na(idx))) {
    missing <- kcols[is.na(idx)]
    stop("Missing kME columns: ", paste(missing, collapse = ", "),
         "\nAvailable kME columns (first 50): ",
         paste(head(colnames(kME_df)[grepl("^kME_", colnames(kME_df))], 50), collapse = ", "))
  }
  
  out <- kME_df[, c("Gene", colnames(kME_df)[idx]), drop = FALSE]
  colnames(out) <- c("Gene", kcols)   # enforce expected names
  out$Gene <- as.character(out$Gene)
  
  # Force numeric kME columns (protect against accidental data.frame columns)
  for (kc in kcols) out[[kc]] <- as.numeric(out[[kc]])
  
  out
}

# -------------------------
# 2) Function: build module-gene-metabolite network for Cytoscape
# -------------------------
make_module_gene_met_network <- function(datExpr, MEs, moduleColors, kME,
                                         sig_pairs,
                                         out_prefix = "module_gene_met",
                                         top_genes_per_module = 30,
                                         fdr_max = 0.05,
                                         cor_min = 0.6,
                                         modules_keep = NULL,
                                         gene_anno = NULL) {
  
  # Filter module-met pairs
  sig_pairs <- sig_pairs[sig_pairs$fdr <= fdr_max & abs(sig_pairs$bicor) >= cor_min, , drop = FALSE]
  if (nrow(sig_pairs) == 0) stop("No module-metabolite pairs after filtering (try relaxing thresholds).")
  
  mods <- unique(as.character(sig_pairs$Module))
  if (!is.null(modules_keep)) mods <- intersect(mods, modules_keep)
  
  # Ensure kME has required columns
  need <- c("Gene", paste0("kME_", mods))
  miss <- setdiff(need, colnames(kME))
  if (length(miss) > 0) {
    stop("kME missing required columns: ", paste(miss, collapse = ", "),
         "\nTip: run kME_clean <- make_kME_clean(geneInfo_base, unique(sig_pairs$Module)) and pass kME_clean.")
  }
  
  # ---- Build module -> gene edges and gene nodes ----
  gene_edges_list <- list()
  gene_nodes_list <- list()
  
  for (m in mods) {
    kcol <- paste0("kME_", m)
    
    kval_all <- as.numeric(kME[[kcol]])
    ord <- order(-abs(kval_all))
    top_genes <- kME$Gene[ord][1:min(top_genes_per_module, length(ord))]
    top_genes <- as.character(top_genes)
    
    idx <- match(top_genes, kME$Gene)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) next
    
    tg   <- as.character(kME$Gene[idx])
    kval <- as.numeric(kME[[kcol]][idx])
    
    gene_edges_list[[m]] <- data.frame(
      source = rep(m, length(tg)),
      target = tg,
      edge_type = rep("module_gene", length(tg)),
      weight = abs(kval),
      sign = ifelse(kval >= 0, "+", "-"),
      stringsAsFactors = FALSE
    )
    
    gene_nodes_list[[m]] <- data.frame(
      id = tg,
      node_type = rep("gene", length(tg)),
      module = rep(m, length(tg)),
      stringsAsFactors = FALSE
    )
  }
  
  # Convert lists to tables
  gene_edges <- rbindlist(gene_edges_list, fill = TRUE)
  gene_nodes <- unique(rbindlist(gene_nodes_list, fill = TRUE))
  
  # ---- Module -> metabolite edges and metabolite nodes ----
  met_edges <- data.frame(
    source = as.character(sig_pairs$Module),
    target = as.character(sig_pairs$Metabolite),
    edge_type = "module_metabolite",
    weight = abs(as.numeric(sig_pairs$bicor)),
    sign = ifelse(as.numeric(sig_pairs$bicor) >= 0, "+", "-"),
    fdr = as.numeric(sig_pairs$fdr),  # <--- This caused the mismatch (Column 6)
    stringsAsFactors = FALSE
  )
  
  met_nodes <- unique(data.frame(
    id = as.character(sig_pairs$Metabolite),
    node_type = "metabolite",
    module = NA_character_,
    stringsAsFactors = FALSE
  ))
  
  # ---- Module nodes ----
  mod_nodes <- data.frame(
    id = mods,
    node_type = "module",
    module = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # ---- Combine nodes (consistent schema) ----
  # Combine nodes using rbindlist with fill=TRUE for safety
  nodes <- unique(rbindlist(list(mod_nodes, gene_nodes, met_nodes), fill = TRUE))
  
  # Optional annotation join (expects gene_anno with Gene column)
  if (!is.null(gene_anno)) {
    nodes <- merge(nodes, gene_anno, by.x = "id", by.y = "Gene", all.x = TRUE)
  }
  
  # ---- Combine edges (THE FIX) ----
  # Use rbindlist(list(...), fill = TRUE) to handle the missing 'fdr' column in gene_edges
  edges <- rbindlist(list(gene_edges, met_edges), fill = TRUE)
  
  # ---- Write ----
  fwrite(nodes, paste0(out_prefix, "_nodes.tsv"), sep = "\t")
  fwrite(edges, paste0(out_prefix, "_edges.tsv"), sep = "\t")
  
  invisible(list(nodes = nodes, edges = edges))
}

# ============================================================
# 3) HOW TO CALL IT (EXACTLY)
# ============================================================
# Assumes you already have these objects in your environment:
#   datExpr, MEs, moduleColors, geneInfo_base, sig_raw, CONFIG$out_dir

# (A) Build a clean kME table (prevents duplicate-name issues)
mods_needed <- unique(as.character(sig_raw$Module))
kME_clean <- make_kME_clean(geneInfo_base, mods_needed)

# (B) Make Cytoscape tables
res <- make_module_gene_met_network(
  datExpr = datExpr,
  MEs = MEs,
  moduleColors = moduleColors,
  kME = kME_clean,  # <-- use the clean kME table
  sig_pairs = sig_raw,
  out_prefix = file.path(CONFIG$out_dir, "Cytoscape_module_gene_met"),
  top_genes_per_module = 30,
  fdr_max = 0.05,
  cor_min = 0.6
)

# (C) Confirm output
list.files(CONFIG$out_dir, pattern = "Cytoscape_module_gene_met_.*\\.tsv$", full.names = TRUE)

# Optional: peek
# head(res$nodes); head(res$edges)

# ============================================================
# FULL UPDATED CODE (R): ggraph network with:
# - nodes/edges from Cytoscape TSV exports
# - igraph build (deduplicated vertices)
# - highlight selected modules + neighborhood/edges
# - node colors by type (module/gene/metabolite)
# - edge colors by sign (+/-)
# - labels: highlighted modules + metabolites in highlight neighborhood + top hub genes
# - TRUE repel labels using ggrepel + ggraph::get_nodes() (fixes x/y error)
# - export high-res PDF + 600 dpi PNG
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(ggrepel)
})

# -------------------------
# USER SETTINGS
# -------------------------
OUT_PREFIX <- file.path(CONFIG$out_dir, "Cytoscape_module_gene_met")  # change if needed

modules_to_highlight <- c("ME0","ME1","ME2")   # <-- set to c("ME12","ME4") etc.
N_GENE_LABELS <- 15                # top genes to label (by degree, within highlight neighborhood)

LABEL_METABOLITES_IN_HIGHLIGHT <- TRUE
LABEL_ALL_MODULES <- FALSE         # TRUE labels all module nodes; FALSE labels only highlighted module nodes

# Node colors
node_cols <- c(
  module = "#1f77b4",
  gene = "#2ca02c",
  metabolite = "#ff7f0e"
)

# Edge sign colors
edge_cols <- c(
  `+` = "#d62728",
  `-` = "#1f77b4"
)

# -------------------------
# 1) READ TSVs
# -------------------------
nodes <- fread(paste0(OUT_PREFIX, "_nodes.tsv"))
edges <- fread(paste0(OUT_PREFIX, "_edges.tsv"))

# types
nodes[, id := as.character(id)]
nodes[, node_type := as.character(node_type)]
if (!("module" %in% colnames(nodes))) nodes[, module := NA_character_]
nodes[, module := as.character(module)]

edges[, source := as.character(source)]
edges[, target := as.character(target)]
edges[, sign := as.character(sign)]
edges[, edge_type := as.character(edge_type)]
edges[, weight := as.numeric(weight)]
if (!("touches_highlight_module" %in% colnames(edges))) {
  # will be computed on the graph later; placeholder not needed here
}

# -------------------------
# 2) DEDUPLICATE NODES (igraph requires unique vertex names)
# -------------------------
nodes_u <- nodes[, .(
  node_type = node_type[1],
  module = {
    m <- module[which(!is.na(module) & module != "")]
    if (length(m) > 0) m[1] else NA_character_
  }
), by = id]

keep_ids <- unique(c(edges$source, edges$target))
nodes_u <- nodes_u[id %in% keep_ids]

setnames(nodes_u, "id", "name")  # igraph vertex id column

# -------------------------
# 3) BUILD GRAPH (igraph)
# -------------------------
g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes_u)

# -------------------------
# 4) HIGHLIGHT TAGS ON GRAPH
# -------------------------
V(g)$is_highlight_module <- (V(g)$node_type == "module" & V(g)$name %in% modules_to_highlight)

highlight_neighbors <- function(g, module_names) {
  mods <- V(g)[V(g)$node_type == "module" & V(g)$name %in% module_names]
  if (length(mods) == 0) return(rep(FALSE, vcount(g)))
  nb <- unique(unlist(neighborhood(g, order = 1, nodes = mods)))
  seq_len(vcount(g)) %in% nb
}
V(g)$is_in_highlight_neighborhood <- highlight_neighbors(g, modules_to_highlight)

edge_ends <- ends(g, E(g))
E(g)$touches_highlight_module <- (edge_ends[,1] %in% modules_to_highlight) |
  (edge_ends[,2] %in% modules_to_highlight)

# -------------------------
# 5) LABEL SELECTION
# -------------------------
V(g)$deg <- degree(g)

# Modules
if (LABEL_ALL_MODULES) {
  V(g)$label_module <- (V(g)$node_type == "module")
} else {
  V(g)$label_module <- (V(g)$node_type == "module" & V(g)$name %in% modules_to_highlight)
}

# Metabolites (in highlight neighborhood)
V(g)$label_metabolite <- FALSE
if (LABEL_METABOLITES_IN_HIGHLIGHT) {
  V(g)$label_metabolite <- (V(g)$node_type == "metabolite" & V(g)$is_in_highlight_neighborhood)
}

# Genes: top by degree, within highlight neighborhood (fallback to all genes if empty)
gene_candidates <- which(V(g)$node_type == "gene" & V(g)$is_in_highlight_neighborhood)
if (length(gene_candidates) == 0) gene_candidates <- which(V(g)$node_type == "gene")

top_gene_vids <- gene_candidates[order(-V(g)$deg[gene_candidates])][1:min(N_GENE_LABELS, length(gene_candidates))]
V(g)$label_gene <- FALSE
V(g)$label_gene[top_gene_vids] <- TRUE

# Final label string
V(g)$label <- ""
V(g)$label[V(g)$label_module] <- V(g)$name[V(g)$label_module]
V(g)$label[V(g)$label_metabolite] <- V(g)$name[V(g)$label_metabolite]
V(g)$label[V(g)$label_gene] <- V(g)$name[V(g)$label_gene]

# -------------------------
# 6) PLOT WITH ggraph + ggrepel (FIXED x/y)
# -------------------------
tg <- as_tbl_graph(g)

p <- ggraph(tg, layout = "fr") +
  
  # edges
  geom_edge_link(aes(edge_colour = sign,
                     edge_width  = weight,
                     edge_alpha  = ifelse(touches_highlight_module, 0.9, 0.12))) +
  scale_edge_colour_manual(values = edge_cols, name = "Interaction sign") +
  scale_edge_width_continuous(range = c(0.2, 2.8), name = "Weight") +
  scale_edge_alpha(range = c(0.05, 1), guide = "none") +
  
  # nodes
  geom_node_point(aes(color = node_type,
                      alpha = ifelse(is_in_highlight_neighborhood, 0.95, 0.18)),
                  size = 2.4) +
  scale_color_manual(values = node_cols, name = "Node type") +
  scale_alpha(range = c(0.1, 1), guide = "none") +
  
  # repel labels: MUST use get_nodes() so x/y are available
  ggrepel::geom_text_repel(
    data = get_nodes(),
    aes(x = x, y = y, label = label),
    size = 3,
    min.segment.length = 0,
    max.overlaps = 200,
    box.padding = 0.35,
    point.padding = 0.2
  ) +
  
  theme_void()

print(p)

# -------------------------
# 7) EXPORT HIGH-RES FIGURES
# -------------------------
out_pdf <- file.path(CONFIG$out_dir, "ggraph_module_gene_met_network.pdf")
out_png <- file.path(CONFIG$out_dir, "ggraph_module_gene_met_network_600dpi.png")

ggsave(out_pdf, p, width = 10, height = 8, device = cairo_pdf)
ggsave(out_png, p, width = 10, height = 8, dpi = 600)

# -------------------------
# 8) OPTIONAL: Export just the highlighted subnetwork (clean manuscript panel)
# -------------------------
keep_v <- which(V(g)$is_in_highlight_neighborhood)
g_sub <- induced_subgraph(g, vids = keep_v)
tg_sub <- as_tbl_graph(g_sub)

# Label again for subgraph: label all modules + metabolites, plus top genes
V(g_sub)$deg <- degree(g_sub)
V(g_sub)$label <- ""
V(g_sub)$label[V(g_sub)$node_type == "module"] <- V(g_sub)$name[V(g_sub)$node_type == "module"]
V(g_sub)$label[V(g_sub)$node_type == "metabolite"] <- V(g_sub)$name[V(g_sub)$node_type == "metabolite"]

gene_vids <- which(V(g_sub)$node_type == "gene")
top_gene_vids_sub <- gene_vids[order(-V(g_sub)$deg[gene_vids])][1:min(N_GENE_LABELS, length(gene_vids))]
V(g_sub)$label[top_gene_vids_sub] <- V(g_sub)$name[top_gene_vids_sub]

p_sub <- ggraph(tg_sub, layout = "fr") +
  geom_edge_link(aes(edge_colour = sign,
                     edge_width  = weight,
                     edge_alpha  = weight)) +
  scale_edge_colour_manual(values = edge_cols, name = "Interaction sign") +
  scale_edge_width_continuous(range = c(0.2, 2.8), name = "Weight") +
  scale_edge_alpha(range = c(0.1, 1), guide = "none") +
  geom_node_point(aes(color = node_type), size = 2.6) +
  scale_color_manual(values = node_cols, name = "Node type") +
  ggrepel::geom_text_repel(
    data = get_nodes(),
    aes(x = x, y = y, label = label),
    size = 3,
    max.overlaps = 200
  ) +
  theme_void()

print(p_sub)

ggsave(file.path(CONFIG$out_dir, "ggraph_subnetwork_highlight.pdf"),
       p_sub, width = 9, height = 7, device = cairo_pdf)
ggsave(file.path(CONFIG$out_dir, "ggraph_subnetwork_highlight_600dpi.png"),
       p_sub, width = 9, height = 7, dpi = 600)


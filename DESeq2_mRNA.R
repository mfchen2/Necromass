#!/usr/bin/env Rscript
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("apeglm")
BiocManager::install("DESeq2")
BiocManager::install("tidyverse")
BiocManager::install("pheatmap")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
})

####DRAM Annotation####
## -------------------- PATHS & BASIC SETTINGS --------------------

count_file       <- "~/Rhodano_mRNAseq/counts/gene_counts.txt"                    # featureCounts output (SAF-based)
sample_info_file <- "~/Rhodano_mRNAseq/data/meta/sample_info.csv"                 # metadata: sample, condition, ...
# annot_file       <- "~/Rhodano_mRNAseq/data/ref/FW104_10B01_gene_annotation.tsv"  # GeneID / product / ko mapping

out_dir   <- "~/Rhodano_mRNAseq/results_deseq2"
plots_dir <- file.path(out_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

## EDIT THESE: reference and comparison level (must match sample_info$condition)
ref_level   <- "0h"      # control / baseline
contrast_to <- "48h"     # e.g. "24h", "48h", "Al", "K", "P"
contrast    <- c("condition", contrast_to, ref_level)

## -------------------- READ COUNTS --------------------

raw_counts <- read.delim(count_file,
                         comment.char = "#",
                         check.names = FALSE)

message("Columns in count file:")
print(colnames(raw_counts))

# featureCounts: cols 1–6 annotation, 7:n samples
count_mat <- raw_counts[, 7:ncol(raw_counts)]
rownames(count_mat) <- raw_counts$Geneid

# Clean column names: "bam/R_0h_1.sorted.bam" -> "R_0h_1"
colnames(count_mat) <- basename(colnames(count_mat))
colnames(count_mat) <- sub("\\.sorted\\.bam$", "", colnames(count_mat))
colnames(count_mat) <- sub("\\.bam$", "", colnames(count_mat))

message("Cleaned count matrix column names:")
print(colnames(count_mat))

## -------------------- READ SAMPLE METADATA --------------------

sample_info <- read.csv(sample_info_file, stringsAsFactors = TRUE)
message("Samples in sample_info:")
print(sample_info$sample)

# Check & align order
if (!all(colnames(count_mat) %in% sample_info$sample)) {
  missing <- setdiff(colnames(count_mat), sample_info$sample)
  stop("Some columns in count matrix are missing from sample_info$sample:\n",
       paste(missing, collapse = ", "))
}

sample_info <- sample_info[match(colnames(count_mat), sample_info$sample), ]
stopifnot(all(colnames(count_mat) == sample_info$sample))

## -------------------- DESeq2 OBJECT --------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = ~ condition
)

# Set reference condition
dds$condition <- relevel(dds$condition, ref = ref_level)

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

## -------------------- RUN DESeq2 --------------------

dds <- DESeq(dds)

message("Available coefficients:")
print(resultsNames(dds))

res <- results(dds, contrast = contrast)

message("DESeq2 results summary:")
print(summary(res))

# Save raw (unannotated) DESeq2 results
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)

raw_out <- file.path(out_dir,
                     paste0("DE_", contrast_to, "_vs_", ref_level, "_raw.csv"))
write.csv(res_df, file = raw_out, row.names = FALSE)
message("Raw DE results saved to: ", raw_out)

## -------------------- READ GENE ANNOTATION & JOIN --------------------

# Expect columns: GeneID, product, ko
annot <- read.delim(annot_file,
                    header = TRUE,
                    stringsAsFactors = FALSE)

if (!"GeneID" %in% colnames(annot)) {
  stop("Annotation file must have a 'GeneID' column.")
}

res_annot <- res_df %>%
  left_join(annot, by = "GeneID")

message("Annotated results (head):")
print(head(res_annot[, c("GeneID", "product", "ko", "log2FoldChange", "padj")]))

annot_out <- file.path(
  out_dir,
  paste0("DE_", contrast_to, "_vs_", ref_level, "_annotated.csv")
)
write.csv(res_annot, file = annot_out, row.names = FALSE)
message("Annotated DE results saved to: ", annot_out)

## -------------------- VST TRANSFORMATION (FOR PCA & HEATMAPS) --------------------

vsd <- vst(dds, blind = FALSE)   # variance-stabilizing transform, standard for RNA-seq QC/plots :contentReference[oaicite:1]{index=1}

## -------------------- VOLCANO PLOT (ONLY ANNOTATED GENES LABELED) --------------------

volc_df <- res_annot %>%
  mutate(
    neg_log10_padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "NS"
    ),
    # label ONLY uses product; unannotated genes get NA, so they won't be labeled
    label = ifelse(!is.na(product) & product != "", product, NA_character_)
  )

# all dots will still be plotted from volc_df
# but top_labs for labels is restricted to annotated, significant genes
top_labs <- volc_df %>%
  filter(sig != "NS",
         !is.na(label),
         !is.na(padj)) %>%
  arrange(padj) %>%
  head(20)

p_volcano <- ggplot(volc_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(
    values = c("Down" = "#3b4cc0", "NS" = "grey80", "Up" = "#b40426")
  ) +
  ggrepel::geom_text_repel(
    data = top_labs,
    aes(label = label),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.4,
    segment.size = 0.2
  ) +
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    color = "Status",
    title = paste("Volcano:", contrast_to, "vs", ref_level),
    subtitle = "All genes shown; only annotated Up/Down labeled"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

pdf(file.path(plots_dir, paste0("volcano_", contrast_to, "_vs_", ref_level, ".pdf")),
    width = 4, height = 4)
print(p_volcano)
dev.off()

## -------------------- PCA PLOT (SAMPLES) --------------------

pca_dat <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_dat, "percentVar"))

p_pca <- ggplot(pca_dat,
                aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    x = paste0("PC1 (", percentVar[1], "% variance)"),
    y = paste0("PC2 (", percentVar[2], "% variance)"),
    color = "Condition",
    title = "PCA of samples (vst)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

pdf(file.path(plots_dir, "PCA_samples.pdf"), width = 4.5, height = 4.5)
print(p_pca)
dev.off()

## -------------------- SAMPLE–SAMPLE CORRELATION HEATMAP --------------------

vsd_mat <- assay(vsd)
cor_mat <- cor(vsd_mat)  # Pearson is standard

ann <- as.data.frame(colData(vsd)[, c("condition"), drop = FALSE])

pdf(file.path(plots_dir, "heatmap_sample_correlation.pdf"),
    width = 4.5, height = 4)
pheatmap(
  cor_mat,
  annotation_col = ann,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  clustering_method = "ward.D2",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 9,
  main = "Sample–sample correlation (vst)"
)
dev.off()

## -------------------- TOP DE GENES HEATMAP --------------------

# top 50 by smallest padj
top_genes <- res_annot %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  pull(GeneID) %>%
  head(50)

mat_top <- assay(vsd)[top_genes, ]

# Optional: shorter labels from product
short_labels <- res_annot %>%
  filter(GeneID %in% top_genes) %>%
  select(GeneID, product) %>%
  mutate(
    short_product = ifelse(
      !is.na(product) & product != "",
      sub("\\s*\\[.*$", "", product),  # strip bracketed EC/KO
      GeneID
    )
  )

# match row order
short_labels <- short_labels[match(rownames(mat_top), short_labels$GeneID), ]
rownames(mat_top) <- short_labels$short_product

# row-wise z-score
mat_scaled <- t(scale(t(mat_top)))   # common practice for RNA-seq heatmaps :contentReference[oaicite:3]{index=3}

pdf(file.path(plots_dir,
              paste0("heatmap_top50_", contrast_to, "_vs_", ref_level, ".pdf")),
    width = 5, height = 6)
pheatmap(
  mat_scaled,
  annotation_col = ann,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  clustering_method = "ward.D2",
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6.5,
  main = paste("Top 50 DE genes:", contrast_to, "vs", ref_level)
)
dev.off()

message("All DESeq2 + annotated plots complete.")

####prokka####
# -------------------- Paths --------------------
count_file <- "~/Rhodano_mRNAseq/counts/gene_counts_prokka.txt"
sample_info_file <- "~/Rhodano_mRNAseq/data/meta/sample_info.csv"
out_dir   <- "~/Rhodano_mRNAseq/results_deseq2_prokka"
plots_dir <- file.path(out_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

## EDIT THESE: reference and comparison level (must match sample_info$condition)
ref_level   <- "48h"      # control / baseline
contrast_to <- "Kanamycin"     # e.g. "24h", "48h", "Al", "K", "P"
contrast    <- c("condition", contrast_to, ref_level)

# -------------------- Read counts --------------------
raw_counts <- read.delim(count_file,
                         comment.char = "#",
                         check.names = FALSE)

message("Columns in count file:")
print(colnames(raw_counts))

count_mat <- raw_counts[, 2:ncol(raw_counts)]
rownames(count_mat) <- raw_counts$Geneid

# Clean column names: bam/R_0h_1.sorted.bam -> R_0h_1
# colnames(count_mat) <- basename(colnames(count_mat))
# colnames(count_mat) <- sub("\\.sorted\\.bam$", "", colnames(count_mat))
# colnames(count_mat) <- sub("\\.bam$", "", colnames(count_mat))

message("Cleaned count matrix column names:")
print(colnames(count_mat))

# -------------------- Read sample metadata --------------------
sample_info <- read.csv(sample_info_file, stringsAsFactors = TRUE)
message("Samples in sample_info:")
print(sample_info$sample)

# Check name consistency
if (!all(colnames(count_mat) %in% sample_info$sample)) {
  missing <- setdiff(colnames(count_mat), sample_info$sample)
  stop("Some columns in count matrix are missing from sample_info$sample:\n",
       paste(missing, collapse = ", "))
}

sample_info <- sample_info[match(colnames(count_mat), sample_info$sample), ]
stopifnot(all(colnames(count_mat) == sample_info$sample))

# -------------------- DESeq2 object --------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = ~ condition
)

## Set reference level for condition if desired
## dds$condition <- relevel(dds$condition, ref = "0h")

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# -------------------- Run DESeq2 --------------------
dds <- DESeq(dds)

message("Available coefficients:")
print(resultsNames(dds))

res <- results(dds, contrast = contrast)

message("DESeq2 results summary:")
print(summary(res))

# Save raw (unannotated) DESeq2 results
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)

raw_out <- file.path(out_dir,
                     paste0("DE_", contrast_to, "_vs_", ref_level, "_raw.csv"))
write.csv(res_df, file = raw_out, row.names = FALSE)
message("Raw DE results saved to: ", raw_out)

## -------------------- VST TRANSFORMATION (FOR PCA & HEATMAPS) --------------------

vsd <- vst(dds, blind = FALSE)   # variance-stabilizing transform, standard for RNA-seq QC/plots :contentReference[oaicite:1]{index=1}

## -------------------- VOLCANO PLOT (ONLY ANNOTATED GENES LABELED) --------------------

# assume `res` is a DESeq2 results object converted to data.frame
res_df <- as.data.frame(res) %>%
  mutate(
    # strip the FHOMJJCM_ prefix from the rownames
    gene = sub("^FHOMJJCM_", "", rownames(res)),
    neg_log10_padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "NS"
    )
  )

# choose top genes to label
top_labs <- res_df %>%
  filter(sig != "NS") %>%
  arrange(padj) %>%
  head(20)

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(
    values = c("Down" = "#3b4cc0", "NS" = "grey80", "Up" = "#b40426")
  ) +
  geom_text_repel(
    data = top_labs,
    aes(label = gene),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.4,
    segment.size = 0.2
  ) +
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    color = "Status",
    title = "Volcano plot",
    subtitle = "|log2FC| > 1, FDR < 0.05"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

pdf(file.path(plots_dir, paste0("volcano_", contrast_to, "_vs_", ref_level, ".pdf")),
    width = 4, height = 4)
print(p_volcano)
dev.off()

## -------------------- PCA PLOT (SAMPLES) --------------------

pca_dat <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_dat, "percentVar"))

p_pca <- ggplot(pca_dat,
                aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    x = paste0("PC1 (", percentVar[1], "% variance)"),
    y = paste0("PC2 (", percentVar[2], "% variance)"),
    color = "Condition",
    title = "PCA of samples (vst)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

pdf(file.path(plots_dir, "PCA_samples_ref.pdf"), width = 4.5, height = 4.5)
print(p_pca)
dev.off()


####prokka####
# -------------------- Paths --------------------
count_file <- "~/Rhodano_mRNAseq/counts/gene_counts_refseq.txt"
sample_info_file <- "~/Rhodano_mRNAseq/data/meta/sample_info.csv"
out_dir   <- "~/Rhodano_mRNAseq/results_deseq2_refseq"
plots_dir <- file.path(out_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

## EDIT THESE: reference and comparison level (must match sample_info$condition)
#ref_level   <- "48h"      # control / baseline
# contrast_to <- "Kanamycin"     # e.g. "24h", "48h", "Al", "K", "P"
# contrast    <- c("condition", contrast_to, ref_level)

# -------------------- Read counts --------------------
raw_counts <- read.delim(count_file,
                         comment.char = "#",
                         check.names = FALSE)

message("Columns in count file:")
print(colnames(raw_counts))

count_mat <- raw_counts[, 2:ncol(raw_counts)]
rownames(count_mat) <- raw_counts$Geneid

# Clean column names: bam/R_0h_1.sorted.bam -> R_0h_1
# colnames(count_mat) <- basename(colnames(count_mat))
# colnames(count_mat) <- sub("\\.sorted\\.bam$", "", colnames(count_mat))
# colnames(count_mat) <- sub("\\.bam$", "", colnames(count_mat))

message("Cleaned count matrix column names:")
print(colnames(count_mat))

# -------------------- Read sample metadata --------------------
sample_info <- read.csv(sample_info_file, stringsAsFactors = TRUE)
message("Samples in sample_info:")
print(sample_info$sample)

# Check name consistency
if (!all(colnames(count_mat) %in% sample_info$sample)) {
  missing <- setdiff(colnames(count_mat), sample_info$sample)
  stop("Some columns in count matrix are missing from sample_info$sample:\n",
       paste(missing, collapse = ", "))
}

sample_info <- sample_info[match(colnames(count_mat), sample_info$sample), ]
stopifnot(all(colnames(count_mat) == sample_info$sample))

# -------------------- DESeq2 object --------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = ~ condition
)

## Set reference level for condition if desired
## dds$condition <- relevel(dds$condition, ref = "0h")

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# -------------------- Run DESeq2 --------------------
dds <- DESeq(dds)
# ---- after you create dds and run dds <- DESeq(dds) ----

alpha <- 0.05        # significance threshold
min_lfc <- 0         # set e.g. 1 for |log2FC| >= 1 requirement (optional)

rp_level <- "P"
others <- c("0h", "24h", "48h", "Al", "Kanamycin")  # edit to match your levels exactly

# sanity checks
message("Condition levels:")
print(levels(dds$condition))
stopifnot(rp_level %in% levels(dds$condition))
stopifnot(all(others %in% levels(dds$condition)))

# store results per contrast
res_list <- list()

for (other in others) {
  res <- results(dds, contrast = c("condition", rp_level, other))
  df <- as.data.frame(res)
  df$GeneID <- rownames(df)
  df$contrast <- paste0(rp_level, "_vs_", other)
  
  # keep just what we need
  df <- df[, c("GeneID", "log2FoldChange", "padj", "contrast")]
  res_list[[other]] <- df
}

# build a wide table so you can export the shared genes with per-contrast stats
wide <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE),
               lapply(names(res_list), function(other) {
                 df <- res_list[[other]]
                 colnames(df)[colnames(df) == "log2FoldChange"] <- paste0("log2FC_", other)
                 colnames(df)[colnames(df) == "padj"] <- paste0("padj_", other)
                 df[, c("GeneID", paste0("log2FC_", other), paste0("padj_", other))]
               }))

# helper: logical vector for "significant in all" + direction
padj_cols <- paste0("padj_", others)
lfc_cols  <- paste0("log2FC_", others)

sig_all <- apply(wide[, padj_cols, drop = FALSE], 1, function(p) all(!is.na(p) & p < alpha))

enriched_all <- sig_all &
  apply(wide[, lfc_cols, drop = FALSE], 1, function(lfc) all(!is.na(lfc) & lfc > 0 & abs(lfc) >= min_lfc))

depleted_all <- sig_all &
  apply(wide[, lfc_cols, drop = FALSE], 1, function(lfc) all(!is.na(lfc) & lfc < 0 & abs(lfc) >= min_lfc))

shared_enriched <- wide[enriched_all, ]
shared_depleted <- wide[depleted_all, ]

# sort by "worst" padj across contrasts (smallest max padj = strongest consistently)
shared_enriched$max_padj <- apply(shared_enriched[, padj_cols, drop = FALSE], 1, max, na.rm = TRUE)
shared_depleted$max_padj <- apply(shared_depleted[, padj_cols, drop = FALSE], 1, max, na.rm = TRUE)

shared_enriched <- shared_enriched[order(shared_enriched$max_padj), ]
shared_depleted <- shared_depleted[order(shared_depleted$max_padj), ]

# export ONLY the shared sets
out_enriched <- file.path(out_dir, paste0("SHARED_enriched_in_", rp_level, "_across_", length(others), "_contrasts.csv"))
out_depleted <- file.path(out_dir, paste0("SHARED_depleted_in_", rp_level, "_across_", length(others), "_contrasts.csv"))

write.csv(shared_enriched, out_enriched, row.names = FALSE)
write.csv(shared_depleted, out_depleted, row.names = FALSE)

message("Exported shared enriched: ", out_enriched, " (n=", nrow(shared_enriched), ")")
message("Exported shared depleted: ", out_depleted, " (n=", nrow(shared_depleted), ")")

# -------------------- DESeq2 object --------------------
# If you have a separate time column in sample_info (recommended), use: design = ~ time + condition
# Otherwise stick to ~ condition
design_formula <- ~ condition

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = design_formula
)

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

message("Condition levels:")
print(levels(dds$condition))

# ---- set this to your phage group label exactly ----
P_level <- "P"

# ---- set others to your other 5 labels exactly ----
others <- c("0h", "24h", "48h", "Al", "Kanamycin")

stopifnot(P_level %in% levels(dds$condition))
stopifnot(all(others %in% levels(dds$condition)))

# -------------------- get results per contrast --------------------
res_list <- lapply(others, function(other) {
  res <- results(dds, contrast = c("condition", P_level, other))
  df <- as.data.frame(res)
  df$GeneID <- rownames(df)
  df$other  <- other
  df
})
names(res_list) <- others

# -------------------- build wide table of log2FC + padj for each comparator ----
wide <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE),
               lapply(others, function(other) {
                 df <- res_list[[other]] %>%
                   dplyr::select(GeneID, log2FoldChange, padj)
                 colnames(df)[2:3] <- c(paste0("log2FC_", other), paste0("padj_", other))
                 df
               }))

padj_cols <- paste0("padj_", others)
lfc_cols  <- paste0("log2FC_", others)

sig_all <- apply(wide[, padj_cols, drop = FALSE], 1, function(p)
  all(!is.na(p) & p < alpha)
)

up_in_P_all <- sig_all &
  apply(wide[, lfc_cols, drop = FALSE], 1, function(lfc)
    all(!is.na(lfc) & lfc >= min_lfc)
  )

down_in_P_all <- sig_all &
  apply(wide[, lfc_cols, drop = FALSE], 1, function(lfc)
    all(!is.na(lfc) & lfc <= -min_lfc)
  )

P_up   <- wide[up_in_P_all, ]
P_down <- wide[down_in_P_all, ]

# rank by "worst" (largest) padj among comparisons, then by mean LFC
P_up$max_padj   <- apply(P_up[, padj_cols, drop = FALSE], 1, max, na.rm = TRUE)
P_up$mean_lfc   <- rowMeans(P_up[, lfc_cols, drop = FALSE], na.rm = TRUE)

P_down$max_padj <- apply(P_down[, padj_cols, drop = FALSE], 1, max, na.rm = TRUE)
P_down$mean_lfc <- rowMeans(P_down[, lfc_cols, drop = FALSE], na.rm = TRUE)

P_up   <- P_up[order(P_up$max_padj, -P_up$mean_lfc), ]
P_down <- P_down[order(P_down$max_padj,  P_down$mean_lfc), ]

# -------------------- write outputs --------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(wide,
          file.path(out_dir, paste0("P_vs_all_wide_alpha", alpha, "_lfc", min_lfc, ".csv")),
          row.names = FALSE)

write.csv(P_up,
          file.path(out_dir, paste0("P_UP_consistent_alpha", alpha, "_lfc", min_lfc, ".csv")),
          row.names = FALSE)

write.csv(P_down,
          file.path(out_dir, paste0("P_DOWN_consistent_alpha", alpha, "_lfc", min_lfc, ".csv")),
          row.names = FALSE)

message("P-up (consistent across all comparisons): n=", nrow(P_up))
message("P-down (consistent across all comparisons): n=", nrow(P_down))



#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
})

## ------------------------------------------------------------
## 1. Define your DE files and stress names
## ------------------------------------------------------------

# EDIT THESE to match your files / stresses
de_files <- c(
  "DE_8h_vs_0h_merged_with_annotation.csv",
  "DE_24h_vs_0h_merged_with_annotation.csv"
)

stress_names <- c("8h", "24h")  # must match order of files above

stopifnot(length(de_files) == length(stress_names))

## ------------------------------------------------------------
## 2. Read and harmonize DE tables
## ------------------------------------------------------------

read_de <- function(file, stress_label) {
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  # GeneID
  if ("GeneID" %in% names(df)) {
    df <- df %>% mutate(GeneID = sub("^ID=", "", as.character(GeneID)))
  } else {
    stop(paste("No 'GeneID' column found in", file))
  }
  
  # log2FC
  if ("log2FoldChange" %in% names(df)) {
    df <- df %>% rename(log2FC = log2FoldChange)
  } else if (!"log2FC" %in% names(df)) {
    stop(paste("No 'log2FoldChange' or 'log2FC' column found in", file))
  }
  
  # padj
  if ("padj" %in% names(df)) {
    # ok
  } else if ("pvalue" %in% names(df)) {
    df <- df %>% rename(padj = pvalue)
  } else {
    stop(paste("No 'padj' or 'pvalue' column found in", file))
  }
  
  # require gene_synonym (from your merged annotation)
  if (!"gene_synonym" %in% names(df)) {
    stop(paste("No 'gene_synonym' column found in", file,
               "- make sure you merged annotation first."))
  }
  
  df %>%
    mutate(
      stress = stress_label
    )
}

de_list <- purrr::map2(de_files, stress_names, read_de)

de_all <- bind_rows(de_list)

## ------------------------------------------------------------
## 3. Restrict to annotated genes & flag significance
## ------------------------------------------------------------
# de_all <- read.csv("~/proteome_annotated_aging_112625.csv")

de_all <- de_all %>%
  # keep only genes with annotation
  mutate(
    gene_synonym = as.character(gene_synonym),
    gene_synonym = ifelse(gene_synonym == "", NA_character_, gene_synonym)
  ) %>%
  filter(!is.na(gene_synonym)) %>%
  # significance flags
  mutate(
    neg_log10_padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FC >  1 ~ "Up",
      padj < 0.05 & log2FC < -1 ~ "Down",
      TRUE                      ~ "NS"
    ),
    sig = factor(sig, levels = c("Down", "NS", "Up"))
  )

de_sig <- de_all %>%
  filter(sig != "NS")

## ------------------------------------------------------------
## 4. Volcano plots (facet by stress, free y-scale)
##     Only annotated genes included; label top 20 per stress
## ------------------------------------------------------------

top_labs <- de_all %>%
  filter(sig != "NS", !is.na(padj)) %>%
  group_by(stress) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = 20) %>%   # top 20 annotated per stress
  ungroup()

p_volcano <- ggplot(de_all, aes(x = log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(
    values = c("Down" = "#3b4cc0", "NS" = "grey80", "Up" = "#b40426"),
    name   = "Status"
  ) +
  ggrepel::geom_text_repel(
    data = top_labs,
    aes(label = gene_synonym),
    size          = 3,
    max.overlaps  = 40,
    box.padding   = 0.4,
    segment.size  = 0.2,
    show.legend   = FALSE
  ) +
  facet_wrap(~ stress, nrow = 1, scales = "free_y") +  # <--- free y per stress
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    title = "Annotated DE genes across stresses",
    subtitle = "padj < 0.05, |log2FC| > 1 highlighted; labels = gene_synonym"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold")
  )

print(p_volcano)

ggsave(
  "Bacillustranscriptome_volcano_annotated_faceted_freeY_aging.pdf",
  p_volcano,
  width = 10,
  height = 5
)

## ------------------------------------------------------------
## 5. (Optional) Heatmap of annotated top genes across stresses
##     Uses only annotated genes; union of top 20 per stress
## ------------------------------------------------------------

TOP_PER_STRESS <- 20

top_genes <- de_sig %>%
  group_by(stress) %>%
  arrange(padj) %>%
  slice_head(n = TOP_PER_STRESS) %>%
  ungroup() %>%
  pull(GeneID) %>%
  unique()

length(top_genes)

mat_long <- de_all %>%
  filter(GeneID %in% top_genes) %>%
  select(GeneID, gene_synonym, stress, log2FC)

mat_wide <- mat_long %>%
  distinct() %>%
  tidyr::pivot_wider(
    names_from  = stress,
    values_from = log2FC
  ) %>%
  as.data.frame()

rownames(mat_wide) <- mat_wide$gene_synonym  # rownames = annotated name
mat_wide$GeneID <- NULL
mat_wide$gene_synonym <- NULL

mat <- as.matrix(mat_wide)

# scale per gene
mat_scaled <- t(scale(t(mat)))  # center & scale rows

pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  main = "Annotated top DE genes across stresses (scaled log2FC)",
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(50)
)

pdf("transcriptome_annotated_top_DE_heatmap_across_stresses.pdf",
    width = 6, height = 8)
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  main = "Annotated top DE genes across stresses (scaled log2FC)",
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(50)
)
dev.off()

## ------------------------------------------------------------
## 3. Restrict to annotated genes & flag significance
## ------------------------------------------------------------

de_all <- read.csv("~/Downloads/proteome_annotated_sup_aging_121525.csv")
de_all <- de_all %>%
  # keep only genes with annotation
  mutate(
    gene_synonym = as.character(gene_synonym),
    gene_synonym = ifelse(gene_synonym == "", NA_character_, gene_synonym)
  ) %>%
  filter(!is.na(gene_synonym)) %>%
  # significance flags
  mutate(
    neg_log10_padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FC >  1 ~ "Up",
      padj < 0.05 & log2FC < -1 ~ "Down",
      TRUE                      ~ "NS"
    ),
    sig = factor(sig, levels = c("Down", "NS", "Up"))
  )

de_sig <- de_all %>%
  filter(sig != "NS")

## ------------------------------------------------------------
## 4. Volcano plots (facet by stress, free y-scale)
##     Only annotated genes included; label top 20 per stress
## ------------------------------------------------------------

top_labs <- de_all %>%
  filter(sig != "NS", !is.na(padj)) %>%
  group_by(stress) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = 20) %>%   # top 20 annotated per stress
  ungroup()

p_volcano <- ggplot(de_all, aes(x = log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(
    values = c("Down" = "#3b4cc0", "NS" = "grey80", "Up" = "#b40426"),
    name   = "Status"
  ) +
  ggrepel::geom_text_repel(
    data = top_labs,
    aes(label = gene_synonym),
    size          = 3,
    max.overlaps  = 40,
    box.padding   = 0.4,
    segment.size  = 0.2,
    show.legend   = FALSE
  ) +
  facet_wrap(~ stress, nrow = 1, scales = "fixed") +  # <--- free y per stress
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    title = "Annotated DE genes across stresses",
    subtitle = "padj < 0.05, |log2FC| > 1 highlighted; labels = gene_synonym"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold")
  )

print(p_volcano)

ggsave(
  "proteomics_volcano_sup_annotated_aging_121525.pdf",
  p_volcano,
  width = 10,
  height = 4
)

## ------------------------------------------------------------
## 5. (Optional) Heatmap of annotated top genes across stresses
##     Uses only annotated genes; union of top 20 per stress
## ------------------------------------------------------------

TOP_PER_STRESS <- 20

top_genes <- de_sig %>%
  group_by(stress) %>%
  arrange(padj) %>%
  slice_head(n = TOP_PER_STRESS) %>%
  ungroup() %>%
  pull(GeneID) %>%
  unique()

length(top_genes)

mat_long <- de_all %>%
  filter(GeneID %in% top_genes) %>%
  select(GeneID, gene_synonym, stress, log2FC)

mat_wide <- mat_long %>%
  distinct() %>%
  tidyr::pivot_wider(
    names_from  = stress,
    values_from = log2FC
  ) %>%
  as.data.frame()

rownames(mat_wide) <- mat_wide$gene_synonym  # rownames = annotated name
mat_wide$GeneID <- NULL
mat_wide$gene_synonym <- NULL

mat <- as.matrix(mat_wide)

# scale per gene
mat_scaled <- t(scale(t(mat)))  # center & scale rows

pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  main = "Annotated top DE genes across stresses (scaled log2FC)",
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(50)
)

pdf("transcriptome_annotated_top_DE_heatmap_across_stresses.pdf",
    width = 6, height = 8)
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  main = "Annotated top DE genes across stresses (scaled log2FC)",
  color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(50)
)
dev.off()

####compare prot and trans####
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

## ------------------------------------------------------------
## 0. Starting point: two data frames in memory
## de_trans     : transcriptomics DE
## de_all_prot  : proteomics DE
## ------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

## Starting point: two data frames already in memory:
##   de_trans      = transcriptomics DE
##   de_all_prot   = proteomics DE

de_rna  <- de_trans
de_prot <- de_all

## 1. Helper to standardize DE tables ------------------------------------

standardize_de <- function(df, source_label) {
  df <- df %>%
    mutate(across(everything(), ~ .x))  # no-op, just explicit
  
  # log2FC
  if ("log2FoldChange" %in% names(df)) {
    df <- df %>% rename(log2FC = log2FoldChange)
  }
  if (!"log2FC" %in% names(df)) {
    stop(paste("No 'log2FC' or 'log2FoldChange' in", source_label))
  }
  
  # padj
  if (!"padj" %in% names(df)) {
    if ("pvalue" %in% names(df)) {
      df <- df %>% rename(padj = pvalue)
    } else {
      stop(paste("No 'padj' or 'pvalue' in", source_label))
    }
  }
  
  # gene_synonym
  if (!"gene_synonym" %in% names(df)) {
    stop(paste("No 'gene_synonym' in", source_label,
               "- make sure annotation was merged."))
  }
  
  df %>%
    mutate(
      gene_synonym = as.character(gene_synonym),
      gene_synonym = ifelse(gene_synonym == "", NA_character_, gene_synonym)
    )
}

de_rna  <- standardize_de(de_rna,  "RNA")
de_prot <- standardize_de(de_prot, "Protein")

## 2. Build a case-insensitive join key ----------------------------------

# function to normalize synonyms: lower-case, trim, optionally strip punctuation
normalize_syn <- function(x) {
  x %>%
    tolower() %>%
    stringr::str_trim() %>%
    # if you want to be even more aggressive:
    # stringr::str_replace_all("[^a-z0-9]", "")
    identity()
}

de_rna <- de_rna %>%
  mutate(
    gene_key = normalize_syn(gene_synonym)
  )

de_prot <- de_prot %>%
  mutate(
    gene_key = normalize_syn(gene_synonym)
  )

## 3. (Optional) restrict to one stress if you have 'stress' ----------------
## If both tables have 'stress' and you want to compare only one condition:
focus_stress <- "Al"
de_rna  <- de_rna  %>% filter(stress == focus_stress)
de_prot <- de_prot %>% filter(stress == focus_stress)

## 4. Significance flags -------------------------------------------------

de_rna <- de_rna %>%
  mutate(
    sig_rna = case_when(
      padj < 0.05 & log2FC >  1 ~ "Up",
      padj < 0.05 & log2FC < -1 ~ "Down",
      TRUE                      ~ "NS"
    )
  )

de_prot <- de_prot %>%
  mutate(
    sig_prot = case_when(
      padj < 0.05 & log2FC >  0.58 ~ "Up",   # ~1.5-fold if you like
      padj < 0.05 & log2FC < -0.58 ~ "Down",
      TRUE                         ~ "NS"
    )
  )

## 5. Join by normalized gene_key (and stress if present) ----------------

join_cols <- c("gene_key")
if ("stress" %in% names(de_rna) && "stress" %in% names(de_prot)) {
  join_cols <- c("gene_key", "stress")
}

overlap <- de_rna %>%
  select(
    all_of(join_cols),
    gene_synonym_rna  = gene_synonym,
    log2FC_rna        = log2FC,
    padj_rna          = padj,
    sig_rna
  ) %>%
  inner_join(
    de_prot %>%
      select(
        all_of(join_cols),
        gene_synonym_prot = gene_synonym,
        log2FC_prot       = log2FC,
        padj_prot         = padj,
        sig_prot
      ),
    by = join_cols
  )

nrow(overlap)
head(overlap)

## 6. Concordance classification ----------------------------------------

overlap <- overlap %>%
  mutate(
    direction_rna  = case_when(
      sig_rna == "Up"   ~ "Up",
      sig_rna == "Down" ~ "Down",
      TRUE              ~ "NS"
    ),
    direction_prot = case_when(
      sig_prot == "Up"   ~ "Up",
      sig_prot == "Down" ~ "Down",
      TRUE               ~ "NS"
    ),
    concordance = case_when(
      direction_rna %in% c("Up", "Down") &
        direction_rna == direction_prot                ~ "Concordant",
      direction_rna %in% c("Up", "Down") &
        direction_prot %in% c("Up", "Down") &
        direction_rna != direction_prot                ~ "Discordant",
      TRUE                                              ~ "RNA-only or Prot-only"
    ),
    concordance = factor(
      concordance,
      levels = c("Concordant", "Discordant", "RNA-only or Prot-only")
    )
  )

overlap %>%
  count(concordance) %>%
  print()

## 7. Scatter: log2FC RNA vs Protein ------------------------------------

# choose a single label to show on the plot; prefer RNA name if different
overlap <- overlap %>%
  mutate(
    label = dplyr::coalesce(gene_synonym_rna, gene_synonym_prot)
  )

overlap_labeled <- overlap %>%
  mutate(
    combined_score = -log10(padj_rna + 1e-16) + -log10(padj_prot + 1e-16)
  ) %>%
  arrange(desc(combined_score), desc(abs(log2FC_rna) + abs(log2FC_prot))) %>%
  slice_head(n = 30)   # label top 30

p_scatter <- ggplot(overlap,
                    aes(x = log2FC_rna,
                        y = log2FC_prot,
                        color = concordance)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(alpha = 0.7, size = 2) +
  ggrepel::geom_text_repel(
    data         = overlap_labeled,
    aes(label    = label),
    size         = 3,
    max.overlaps = 40,
    box.padding  = 0.4,
    segment.size = 0.2,
    show.legend  = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Concordant"            = "#1b9e77",
      "Discordant"            = "#d95f02",
      "RNA-only or Prot-only" = "grey60"
    ),
    name = "RNA vs Protein"
  ) +
  labs(
    x = "Transcript log2FC",
    y = "Protein log2FC",
    title = "Concordance between transcriptome and proteome",
    subtitle = "Joined case-insensitively by gene_synonym"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_scatter)

ggsave(
  "RNA_vs_protein_log2FC_scatter_concordance_case_insensitive.pdf",
  p_scatter,
  width = 5,
  height = 5
)

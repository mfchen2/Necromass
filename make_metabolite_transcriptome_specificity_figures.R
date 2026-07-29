#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

base_dir <- "/Users/mingfeichen"
out_dir <- file.path(
  getwd(),
  "outputs",
  "manual-20260601-a5",
  "presentations",
  "metabolite_transcriptome_specificity",
  "assets"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.02), margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(0.92), margin = margin(b = 6)),
      axis.title = element_text(face = "plain"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.border = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", size = rel(0.85))
    )
}

save_all <- function(plot, basename, width, height) {
  ggsave(file.path(out_dir, paste0(basename, ".png")), plot, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(basename, ".pdf")), plot, width = width, height = height)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(out_dir, paste0(basename, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  }
}

species_from_treatment <- function(tr) {
  if_else(str_starts(tr, "B_"), "Bacillus", "Rhodanobacter")
}

treatment_display_met <- function(tr) {
  recode(
    tr,
    B_0hr = "B 0h",
    B_8hr = "B mid",
    B_24hr = "B late",
    B_Al = "B Al",
    B_K = "B Kana",
    B_P = "B Phage",
    R_0hr = "R 0h",
    R_24hr = "R mid",
    R_48hr = "R late",
    R_Al = "R Al",
    R_K = "R Kana",
    R_P = "R Phage"
  )
}

treatment_display_tx <- function(organism, tr) {
  case_when(
    organism == "Bacillus" & tr == "8hr" ~ "B mid",
    organism == "Bacillus" & tr == "24hr" ~ "B late",
    organism == "Bacillus" & tr == "Al" ~ "B Al",
    organism == "Bacillus" & tr == "K" ~ "B Kana",
    organism == "Bacillus" & tr == "P" ~ "B Phage",
    organism == "Rhodanobacter" & tr == "24hr" ~ "R mid",
    organism == "Rhodanobacter" & tr == "48hr" ~ "R late",
    organism == "Rhodanobacter" & tr == "Al" ~ "R Al",
    organism == "Rhodanobacter" & tr == "K" ~ "R Kana",
    organism == "Rhodanobacter" & tr == "P" ~ "R Phage",
    TRUE ~ NA_character_
  )
}

met_superclass_order <- c(
  "Amino acids & derivatives",
  "Nucleic acids (bases/nucleosides/nucleotides)",
  "Central carbon & organic acids",
  "Lipids & osmolytes (quaternary amines)",
  "Polyamines & guanidino compounds",
  "Aromatic compounds & phenolics",
  "Methylation & sulfur salvage",
  "Cofactors (pterins)",
  "Sugars & phosphorylated sugars"
)

met_palette <- c(
  "Amino acids & derivatives" = "#F58518",
  "Nucleic acids (bases/nucleosides/nucleotides)" = "#4C78A8",
  "Central carbon & organic acids" = "#54A24B",
  "Lipids & osmolytes (quaternary amines)" = "#E45756",
  "Polyamines & guanidino compounds" = "#B279A2",
  "Aromatic compounds & phenolics" = "#72B7B2",
  "Methylation & sulfur salvage" = "#9D9D9D",
  "Cofactors (pterins)" = "#FF9DA6",
  "Sugars & phosphorylated sugars" = "#1F9E89"
)

pathway_lookup <- tribble(
  ~pathway_id, ~pathway_name, ~pathway_class,
  "map00230", "Purine metabolism", "Nucleotide metabolism",
  "map00240", "Pyrimidine metabolism", "Nucleotide metabolism",
  "map00220", "Arginine biosynthesis", "Amino acid / nitrogen metabolism",
  "map00250", "Alanine, aspartate and glutamate metabolism", "Amino acid / nitrogen metabolism",
  "map00270", "Cysteine and methionine metabolism", "Amino acid / nitrogen metabolism",
  "map00280", "Valine, leucine and isoleucine degradation", "Amino acid / nitrogen metabolism",
  "map00290", "Valine, leucine and isoleucine biosynthesis", "Amino acid / nitrogen metabolism",
  "map00650", "Butanoate metabolism", "Central carbon / energy metabolism",
  "map00660", "C5-branched dibasic acid metabolism", "Central carbon / energy metabolism",
  "map00680", "Methane metabolism", "Central carbon / energy metabolism",
  "map00620", "Pyruvate metabolism", "Central carbon / energy metabolism",
  "map00020", "Citrate cycle (TCA cycle)", "Central carbon / energy metabolism",
  "map00190", "Oxidative phosphorylation", "Central carbon / energy metabolism",
  "map00790", "Folate biosynthesis", "Cofactor / one-carbon metabolism",
  "map00970", "Aminoacyl-tRNA biosynthesis", "Translation / genetic processing",
  "map01110", "Biosynthesis of secondary metabolites", "Secondary metabolism / transport",
  "map01210", "2-Oxocarboxylic acid metabolism", "Amino acid / nitrogen metabolism",
  "map01230", "Biosynthesis of amino acids", "Amino acid / nitrogen metabolism",
  "map02010", "ABC transporters", "Transport / signaling",
  "map02020", "Two-component system", "Transport / signaling",
  "map00340", "Histidine metabolism", "Amino acid / nitrogen metabolism"
)

pathway_class_order <- c(
  "Amino acid / nitrogen metabolism",
  "Nucleotide metabolism",
  "Central carbon / energy metabolism",
  "Cofactor / one-carbon metabolism",
  "Transport / signaling",
  "Translation / genetic processing",
  "Secondary metabolism / transport",
  "Other"
)

pathway_palette <- c(
  "Amino acid / nitrogen metabolism" = "#F58518",
  "Nucleotide metabolism" = "#4C78A8",
  "Central carbon / energy metabolism" = "#54A24B",
  "Cofactor / one-carbon metabolism" = "#B279A2",
  "Transport / signaling" = "#E45756",
  "Translation / genetic processing" = "#72B7B2",
  "Secondary metabolism / transport" = "#9D9D9D",
  "Other" = "#C7C7C7"
)

pathway_class_order_main <- pathway_class_order[pathway_class_order != "Other"]

make_metabolite_plot <- function() {
  met <- read_csv(file.path(base_dir, "All_Metabolite_one_vs_rest_significant_markers.csv"), show_col_types = FALSE) %>%
    mutate(
      species = species_from_treatment(treatment),
      treatment_display = treatment_display_met(treatment),
      direction = as.character(direction),
      superclass = factor(superclass, levels = met_superclass_order)
    ) %>%
    filter(direction == "enriched", !is.na(padj), padj < 0.05) %>%
    mutate(
      treatment_display = factor(
        treatment_display,
        levels = c("B 0h", "B mid", "B late", "B Phage", "B Kana", "B Al", "R 0h", "R mid", "R late", "R Phage", "R Kana", "R Al")
      )
    )

  summary_counts <- met %>%
    count(species, treatment_display, superclass, name = "n") %>%
    mutate(
      treatment_display = fct_relevel(
        treatment_display,
        "B 0h", "B mid", "B late", "B Phage", "B Kana", "B Al",
        "R 0h", "R mid", "R late", "R Phage", "R Kana", "R Al"
      ),
      superclass = fct_drop(superclass)
    )

  exemplar_rows <- met %>%
    group_by(species, superclass) %>%
    slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(species, superclass, Gene)

  exemplar_points <- met %>%
    semi_join(exemplar_rows, by = c("species", "superclass", "Gene")) %>%
    mutate(point_size = pmin(8, 1.3 + -log10(padj)))

  exemplar_b <- exemplar_points %>%
    filter(species == "Bacillus") %>%
    mutate(
      treatment_display = factor(
        treatment_display,
        levels = c("B 0h", "B mid", "B late", "B Phage", "B Kana", "B Al")
      ),
      Gene = fct_reorder2(Gene, as.numeric(treatment_display), -log10(padj))
    )

  exemplar_r <- exemplar_points %>%
    filter(species == "Rhodanobacter") %>%
    mutate(
      treatment_display = factor(
        treatment_display,
        levels = c("R 0h", "R mid", "R late", "R Phage", "R Kana", "R Al")
      ),
      Gene = fct_reorder2(Gene, as.numeric(treatment_display), -log10(padj))
    )

  row_order_b <- exemplar_b %>%
    distinct(superclass, Gene, padj) %>%
    group_by(superclass, Gene) %>%
    summarise(score = min(padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(superclass, score) %>%
    pull(Gene)

  row_order_r <- exemplar_r %>%
    distinct(superclass, Gene, padj) %>%
    group_by(superclass, Gene) %>%
    summarise(score = min(padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(superclass, score) %>%
    pull(Gene)

  exemplar_b <- exemplar_b %>%
    mutate(Gene = factor(Gene, levels = rev(unique(row_order_b))))
  exemplar_r <- exemplar_r %>%
    mutate(Gene = factor(Gene, levels = rev(unique(row_order_r))))

  p_summary <- ggplot(summary_counts, aes(x = treatment_display, y = n, fill = superclass)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2) +
    facet_wrap(~species, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = met_palette, drop = FALSE) +
    labs(
      title = "Enriched metabolites are concentrated in distinct chemical superclasses",
      x = NULL,
      y = "Number of significant enriched metabolites",
      fill = "Superclass"
    ) +
    theme_pub(base_size = 10) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold", size = 10)
    )

  p_exemplar_b <- ggplot(exemplar_b, aes(x = treatment_display, y = Gene)) +
    geom_point(aes(size = point_size, fill = superclass), shape = 21, color = "black", stroke = 0.25, alpha = 0.95) +
    scale_fill_manual(values = met_palette, drop = FALSE) +
    scale_size_identity() +
    labs(
      title = "Bacillus",
      x = NULL,
      y = NULL,
      size = "-log10(padj)",
      fill = "Superclass"
    ) +
    theme_pub(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "right"
    )

  p_exemplar_r <- ggplot(exemplar_r, aes(x = treatment_display, y = Gene)) +
    geom_point(aes(size = point_size, fill = superclass), shape = 21, color = "black", stroke = 0.25, alpha = 0.95) +
    scale_fill_manual(values = met_palette, drop = FALSE) +
    scale_size_identity() +
    labs(
      title = "Rhodanobacter",
      x = NULL,
      y = NULL,
      size = "-log10(padj)",
      fill = "Superclass"
    ) +
    theme_pub(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "right"
    )

  p_exemplar <- (p_exemplar_b | p_exemplar_r) + plot_layout(guides = "collect")

  combined <- p_summary / p_exemplar + plot_annotation(tag_levels = "A")
  list(
    plot = combined,
    summary = summary_counts,
    exemplar = bind_rows(exemplar_b, exemplar_r)
  )
}

make_transcriptome_plot <- function() {
  joined <- read_csv(file.path(getwd(), "pathway_prediction_vs_exometabolome_outputs", "joined_gene_pathway_vs_exomet_pathway_all.csv"), show_col_types = FALSE) %>%
    left_join(pathway_lookup, by = "pathway_id") %>%
    mutate(
      pathway_class = if_else(is.na(pathway_class), "Other", pathway_class),
      pathway_name = if_else(is.na(pathway_name), pathway_id, pathway_name),
      treatment_display = treatment_display_tx(organism, treatment)
    )

  joined_specific <- joined %>%
    filter(informative, pathway_class != "Other") %>%
    mutate(
      pathway_class = factor(pathway_class, levels = pathway_class_order_main),
      treatment_display = factor(
        treatment_display,
        levels = c("B mid", "B late", "B Phage", "B Kana", "B Al", "R mid", "R late", "R Phage", "R Kana", "R Al")
      )
    )

  exemplar_paths <- joined_specific %>%
    group_by(organism, pathway_class) %>%
    slice_min(order_by = combined_rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(pathway_class, pathway_id, pathway_name) %>%
    distinct()

  heatmap_df <- joined_specific %>%
    semi_join(exemplar_paths, by = c("pathway_class", "pathway_id", "pathway_name")) %>%
    mutate(
      treatment_group = as.character(treatment_display),
      pathway_label = paste0(pathway_class, " | ", pathway_name, " (", pathway_id, ")")
    ) %>%
    group_by(pathway_class, pathway_label) %>%
    mutate(row_rank = min(combined_rank, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(pathway_class, row_rank, pathway_label) %>%
    mutate(
      treatment_group = factor(
        treatment_group,
        levels = c("B mid", "B late", "B Phage", "B Kana", "B Al", "R mid", "R late", "R Phage", "R Kana", "R Al")
      ),
      pathway_label = factor(pathway_label, levels = rev(unique(pathway_label)))
    )

  p_heatmap <- ggplot(heatmap_df, aes(x = treatment_group, y = pathway_label, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      name = "NES"
    ) +
    labs(
      title = "Representative informative pathways show treatment and phylogeny specificity",
      x = NULL,
      y = NULL
    ) +
    theme_pub(base_size = 12) +
    theme(
      plot.title = element_text(size = rel(1.05), face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.position = "right"
    )
    
  list(
    plot = p_heatmap,
    summary = joined_specific %>%
      distinct(organism, treatment, pathway_id, pathway_class) %>%
      count(organism, treatment, pathway_class, name = "n") %>%
      mutate(
        treatment_display = treatment_display_tx(organism, treatment),
        treatment_display = factor(
          treatment_display,
          levels = c("B mid", "B late", "B Phage", "B Kana", "B Al", "R mid", "R late", "R Phage", "R Kana", "R Al")
        ),
        pathway_class = factor(pathway_class, levels = pathway_class_order_main)
      ),
    exemplar = heatmap_df
  )
}

met_res <- make_metabolite_plot()
tx_res <- make_transcriptome_plot()

save_all(met_res$plot, "metabolite_specificity_panels", width = 14.5, height = 8.75)
save_all(tx_res$plot, "transcriptome_specificity_panels", width = 14.5, height = 8.75)

write_csv(met_res$summary, file.path(out_dir, "metabolite_superclass_summary.csv"))
write_csv(met_res$exemplar, file.path(out_dir, "metabolite_exemplar_points.csv"))
write_csv(tx_res$summary, file.path(out_dir, "transcriptome_pathway_class_summary.csv"))
write_csv(tx_res$exemplar, file.path(out_dir, "transcriptome_exemplar_heatmap.csv"))

message("Wrote figure outputs to: ", out_dir)

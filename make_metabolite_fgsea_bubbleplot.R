#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(face = "bold", colour = "#1F2937"),
      axis.text = element_text(colour = "#1F2937"),
      plot.title = element_text(face = "bold", colour = "#111827"),
      plot.subtitle = element_text(colour = "#4B5563"),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text = element_text(face = "bold", colour = "#111827"),
      legend.title = element_text(face = "bold", colour = "#111827"),
      legend.text = element_text(colour = "#1F2937"),
      plot.margin = margin(5, 6, 5, 6)
    )
}

save_pub <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white"
  )
}

input_files <- c(
  Rhodanobacter = "/Users/mingfeichen/Rhodanobacter_Merged_FGSEA_Results.csv",
  Bacillus = "/Users/mingfeichen/Bacillus_Merged_FGSEA_Results.csv"
)

out_dir <- "/Users/mingfeichen/Manuscript/outputs/manual-20260606-a3/presentations/metabolite_fgsea_bubbleplot/assets"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

treatment_order <- c("mid", "late", "Al", "Kana", "Phage")
treatment_map <- c(
  "Rhodanobacter" = "R",
  "Bacillus" = "B"
)

pathway_lookup <- tribble(
  ~pathway_id, ~pathway_name, ~pathway_class,
  "map00230", "Purine metabolism", "Nucleotide metabolism",
  "map00240", "Pyrimidine metabolism", "Nucleotide metabolism",
  "map00250", "Alanine, aspartate and glutamate metabolism", "Amino acid / nitrogen metabolism",
  "map00290", "Valine, leucine and isoleucine biosynthesis", "Amino acid / nitrogen metabolism",
  "map00051", "Fructose and mannose metabolism", "Central carbon / energy metabolism",
  "map00500", "Starch and sucrose metabolism", "Central carbon / energy metabolism",
  "map00550", "Peptidoglycan biosynthesis", "Secondary metabolism / transport",
  "map00620", "Pyruvate metabolism", "Central carbon / energy metabolism",
  "map00630", "Glyoxylate and dicarboxylate metabolism", "Central carbon / energy metabolism",
  "map00650", "Butanoate metabolism", "Central carbon / energy metabolism",
  "map00680", "Methane metabolism", "Central carbon / energy metabolism",
  "map00790", "Folate biosynthesis", "Cofactor / one-carbon metabolism",
  "map01130", "Biosynthesis of secondary metabolites", "Secondary metabolism / transport",
  "map01200", "Carbon metabolism", "Central carbon / energy metabolism",
  "map01210", "2-Oxocarboxylic acid metabolism", "Amino acid / nitrogen metabolism",
  "map01230", "Biosynthesis of amino acids", "Amino acid / nitrogen metabolism",
  "map02010", "ABC transporters", "Transport / signaling",
  "map02020", "Two-component system", "Transport / signaling",
  "map02025", "Biofilm formation", "Transport / signaling",
  "map02030", "Bacterial chemotaxis", "Transport / signaling",
  "map02040", "Flagellar assembly", "Transport / signaling",
  "map02060", "Phosphotransferase system", "Transport / signaling",
  "map03010", "Ribosome", "Translation / genetic processing",
  "map03060", "Protein export", "Translation / genetic processing",
  "map03070", "Bacterial secretion system", "Translation / genetic processing",
  "map04112", "Cell cycle", "Transport / signaling"
)

pathway_class_order <- c(
  "Transport / signaling",
  "Amino acid / nitrogen metabolism",
  "Nucleotide metabolism",
  "Central carbon / energy metabolism",
  "Cofactor / one-carbon metabolism",
  "Translation / genetic processing",
  "Secondary metabolism / transport",
  "Other"
)

pretty_name <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim() %>%
    str_replace_all("([a-z])([A-Z])", "\\1 \\2") %>%
    str_to_sentence()
}

load_fgsea <- function(path, organism) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      organism = organism,
      treatment = case_when(
        contrast == "24h_vs_0h" & organism == "Bacillus" ~ "late",
        contrast == "24h_vs_0h" & organism == "Rhodanobacter" ~ "mid",
        contrast == "48h_vs_0h" ~ "late",
        contrast == "8h_vs_0h" ~ "mid",
        contrast == "P_vs_0h" ~ "Phage",
        contrast == "Kanamycin_vs_0h" ~ "Kana",
        contrast == "Al_vs_0h" ~ "Al",
        TRUE ~ NA_character_
      ),
      pathway_id = set_id,
      padj = as.numeric(padj),
      NES = as.numeric(NES)
    ) %>%
    filter(!is.na(treatment), treatment %in% treatment_order)
}

select_top_pathways <- function(df, n_keep = 8) {
  df %>%
    filter(str_starts(pathway_id, "map"), is.finite(padj), is.finite(NES), padj < 0.05) %>%
    left_join(pathway_lookup, by = "pathway_id") %>%
    mutate(
      pathway_class = if_else(is.na(pathway_class), "Other", pathway_class),
      pathway_name = if_else(is.na(pathway_name), pathway_id, pathway_name),
      pathway_label = if_else(pathway_name == pathway_id, pathway_id, paste0(pathway_name, " (", pathway_id, ")")),
      pathway_class = factor(pathway_class, levels = pathway_class_order)
    ) %>%
    group_by(organism, pathway_id, pathway_label, pathway_class) %>%
    summarise(
      n_contrasts = n_distinct(treatment),
      min_padj = min(padj, na.rm = TRUE),
      max_abs_nes = max(abs(NES), na.rm = TRUE),
      best_nes = NES[which.min(padj)][1],
      .groups = "drop"
    ) %>%
    group_by(organism) %>%
    arrange(desc(n_contrasts), desc(max_abs_nes), min_padj, pathway_label, .by_group = TRUE) %>%
    slice_head(n = n_keep) %>%
    ungroup()
}

make_panel <- function(df, organism_name, show_x = TRUE) {
  selected_paths <- df %>%
    filter(organism == organism_name) %>%
    arrange(pathway_class, desc(max_abs_nes), min_padj, pathway_label)

  panel_df <- fgsea_all %>%
    filter(organism == organism_name, pathway_id %in% selected_paths$pathway_id, padj < 0.05, str_starts(pathway_id, "map")) %>%
    left_join(
      selected_paths %>% select(pathway_id, pathway_label, pathway_class, min_padj, max_abs_nes),
      by = "pathway_id"
    ) %>%
    mutate(
      treatment = factor(treatment, levels = treatment_order),
      pathway_class = factor(pathway_class, levels = pathway_class_order),
      pathway_label = factor(pathway_label, levels = rev(selected_paths$pathway_label))
    ) %>%
    arrange(pathway_class, desc(max_abs_nes), min_padj)

  ggplot(panel_df, aes(x = treatment, y = pathway_label)) +
    geom_point(
      aes(fill = NES, size = pmin(-log10(padj), 10)),
      shape = 21,
      colour = "#2E2E2E",
      stroke = 0.18,
      alpha = 0.92
    ) +
    facet_grid(pathway_class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_discrete(drop = FALSE) +
    scale_fill_gradient2(
      low = "#3B73B9",
      mid = "white",
      high = "#D1495B",
      midpoint = 0,
      name = "NES"
    ) +
    scale_size_continuous(
      range = c(1.6, 5.4),
      name = "-log10(padj)"
    ) +
    labs(
      title = organism_name,
      x = if (show_x) "Treatment" else NULL,
      y = "KEGG pathway (selected)"
    ) +
    theme_pub(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7.7),
      legend.position = "bottom",
      legend.box = "vertical",
      strip.placement = "outside",
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5),
      panel.grid.major = element_line(color = "#EEF2F7", linewidth = 0.25),
      panel.spacing.y = unit(0.12, "in"),
      plot.margin = margin(5, 8, 4, 8)
    )
}

fgsea_all <- bind_rows(
  load_fgsea(input_files[["Rhodanobacter"]], "Rhodanobacter"),
  load_fgsea(input_files[["Bacillus"]], "Bacillus")
)

selected <- select_top_pathways(fgsea_all, n_keep = 8)

if (nrow(selected) == 0) {
  stop("No significant map pathways passed the selection filters.")
}

plot_bacillus <- make_panel(selected, "Bacillus", show_x = FALSE)
plot_rhodano <- make_panel(selected, "Rhodanobacter", show_x = TRUE) +
  theme(legend.position = "bottom")
plot_bacillus <- plot_bacillus + theme(legend.position = "none")

figure <- plot_bacillus / plot_rhodano +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Metabolite FGSEA bubble plot highlights stress-dependent pathway programs",
    subtitle = "Only the strongest KEGG map pathways are shown. Bubble color encodes NES and bubble size encodes significance (-log10 adjusted p-value).",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, colour = "#111827", hjust = 0.5),
      plot.subtitle = element_text(size = 9.3, colour = "#4B5563", hjust = 0.5)
    )
  )

write_csv(selected, file.path(out_dir, "metabolite_fgsea_bubbleplot_selected.csv"))
save_pub(figure, "metabolite_fgsea_bubbleplot.png", width = 13.8, height = 8.8)
save_pub(figure, "metabolite_fgsea_bubbleplot.pdf", width = 13.8, height = 8.8)
save_pub(figure, "metabolite_fgsea_bubbleplot.svg", width = 13.8, height = 8.8)

message("Wrote: ", file.path(out_dir, "metabolite_fgsea_bubbleplot.png"))
message("Wrote: ", file.path(out_dir, "metabolite_fgsea_bubbleplot.pdf"))
message("Wrote: ", file.path(out_dir, "metabolite_fgsea_bubbleplot.svg"))

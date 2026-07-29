#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ---------------------------
# User-adjustable parameters
# ---------------------------
home_dir <- path.expand("~")

rhodano_fgsea_file <- file.path(home_dir, "Rhodanobacter_Merged_FGSEA_Results.csv")
bacillus_fgsea_file <- file.path(home_dir, "Bacillus_Merged_FGSEA_Results.csv")
metabolite_summary_file <- file.path(home_dir, "metabolite_pathway_summary_vs0h.tsv")

out_dir <- file.path(getwd(), "pathway_prediction_vs_exometabolome_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

padj_thresh <- 0.05
abs_nes_thresh <- 1.0
min_mets <- 2
min_sig_mets <- 1
n_show_heatmap <- 10
label_top_n <- 2

allowed_pathway_ids_file <- NA_character_

treatment_map <- c(
  "P_vs_0h" = "P",
  "Kanamycin_vs_0h" = "K",
  "Al_vs_0h" = "Al",
  "24h_vs_0h" = "24hr",
  "48h_vs_0h" = "48hr",
  "8h_vs_0h" = "8hr"
)

treatment_order <- c("P", "K", "Al", "8hr", "24hr", "48hr")
treatment_display <- c(
  "P" = "P",
  "K" = "K",
  "Al" = "Al",
  "8hr" = "8 h",
  "24hr" = "24 h",
  "48hr" = "48 h"
)

organism_order <- c("Rhodanobacter", "Bacillus")
panel_levels <- c(
  paste("Rhodanobacter", treatment_order, sep = " | "),
  paste("Bacillus", treatment_order, sep = " | ")
)

resolve_existing_path <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input file: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

safe_concordance <- function(x, y) {
  ok <- is.finite(x) & is.finite(y) & (x != 0) & (y != 0)
  if (sum(ok) < 3) return(NA_real_)
  mean(sign(x[ok]) == sign(y[ok]))
}

panel_title_from_label <- function(x) {
  parts <- str_split_fixed(as.character(x), " \\| ", 2)
  treatment_display[parts[, 2]]
}

strip_pathway_prefix <- function(x) {
  str_replace(as.character(x), ".*\\|", "")
}

base_plot_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_blank(),
    axis.title = element_text(size = 9.5),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7.8)
  )

load_fgsea <- function(path, organism, phylogeny) {
  read_csv(path, show_col_types = FALSE) %>%
    transmute(
      organism = organism,
      phylogeny = phylogeny,
      contrast = contrast,
      treatment = recode(contrast, !!!treatment_map),
      pathway_id = set_id,
      padj = as.numeric(padj),
      NES = as.numeric(NES)
    )
}

load_metabolite_summary <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      phylogeny = phylogeny,
      treatment = treatment,
      pathway_id = pathway_id,
      n_mets = as.numeric(n_mets),
      n_sig = as.numeric(n_sig),
      mean_log2FC = as.numeric(mean_log2FC),
      frac_up = as.numeric(frac_up),
      top_up = top_up,
      top_down = top_down
    )
}

make_scatter_plot <- function(joined, label_top_n) {
  scatter_df <- joined %>%
    filter(gene_sig | met_sig) %>%
    mutate(point_size = pmax(1, n_sig))

  label_df <- joined %>%
    filter(informative) %>%
    group_by(organism, treatment, panel_label) %>%
    arrange(desc(combined_rank), .by_group = TRUE) %>%
    slice_head(n = label_top_n) %>%
    ungroup()

  ggplot(scatter_df, aes(x = NES, y = met_score)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_point(aes(color = relationship, size = point_size), alpha = 0.82) +
    geom_text_repel(
      data = label_df,
      aes(label = pathway_id),
      size = 2.4,
      box.padding = 0.18,
      point.padding = 0.12,
      segment.alpha = 0.45,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    facet_wrap(~ panel_label, ncol = 5, labeller = as_labeller(panel_title_from_label)) +
    scale_color_manual(
      values = c(
        "same direction" = "#D1495B",
        "opposite direction" = "#2F6DB5",
        "weak/zero" = "grey75",
        "missing" = "grey88"
      ),
      breaks = c("same direction", "opposite direction", "weak/zero"),
      name = NULL
    ) +
    scale_size_continuous(name = "n significant metabolites", range = c(1.3, 3.8)) +
    labs(
      title = "Pathway-level concordance between transcriptome (FGSEA) and exometabolome",
      x = "Gene-pathway enrichment (NES)",
      y = "Observed exometabolite pathway score"
    ) +
    base_plot_theme +
    theme(
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      axis.text.x = element_text(size = 7.5, margin = margin(t = 2)),
      axis.text.y = element_text(size = 7.5, margin = margin(r = 3)),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(8, 12, 14, 12)
    )
}

make_heatmap_data <- function(joined, n_show_heatmap) {
  heatmap_top <- joined %>%
    filter(informative) %>%
    group_by(organism, treatment, panel_label) %>%
    arrange(desc(combined_rank), .by_group = TRUE) %>%
    slice_head(n = n_show_heatmap) %>%
    ungroup() %>%
    mutate(panel_pathway = paste(panel_label, pathway_id, sep = "|"))

  panel_levels_present <- heatmap_top %>%
    distinct(panel_label, panel_pathway, combined_rank) %>%
    arrange(panel_label, combined_rank) %>%
    pull(panel_pathway)

  heatmap_long <- heatmap_top %>%
    select(organism, treatment, panel_label, panel_pathway, pathway_id, NES, met_score) %>%
    pivot_longer(cols = c(NES, met_score), names_to = "metric", values_to = "value") %>%
    group_by(panel_label, metric) %>%
    mutate(
      value_scaled = {
        denom <- max(abs(value), na.rm = TRUE)
        if (is.finite(denom) && denom > 0) value / denom else rep(0, dplyr::n())
      }
    ) %>%
    ungroup() %>%
    mutate(
      panel_pathway = factor(panel_pathway, levels = rev(unique(panel_levels_present))),
      metric = recode(metric, NES = "NES", met_score = "Exomet\nscore")
    )

  list(
    heatmap_top = heatmap_top,
    heatmap_long = heatmap_long
  )
}

make_heatmap_plot <- function(heatmap_long) {
  ggplot(heatmap_long, aes(x = metric, y = panel_pathway, fill = value_scaled)) +
    geom_tile(color = "white", linewidth = 0.35) +
    facet_wrap(~ panel_label, ncol = 5, scales = "free_y", labeller = as_labeller(panel_title_from_label)) +
    scale_y_discrete(labels = strip_pathway_prefix) +
    scale_fill_gradient2(
      low = "#3B73B9",
      mid = "white",
      high = "#D1495B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Within-metric\nscaled score"
    ) +
    labs(
      title = "Top pathways by concordance within each condition",
      x = NULL,
      y = "KEGG pathway (map)"
    ) +
    base_plot_theme +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 7.2, margin = margin(t = 2)),
      axis.text.y = element_text(size = 6.8, margin = margin(r = 3)),
      legend.position = "right",
      plot.margin = margin(8, 20, 14, 12)
    )
}

make_bubble_plot <- function(fgsea_all) {
  plot_data <- fgsea_all %>%
    filter(str_starts(pathway_id, "map")) %>%
    mutate(
      neg_log_padj = pmin(-log10(padj), 10),
      is_significant = padj <= padj_thresh,
      treatment = factor(treatment, levels = treatment_order),
      treatment_display = factor(treatment_display[as.character(treatment)], levels = treatment_display[treatment_order]),
      organism = factor(organism, levels = organism_order)
    )

  sig_pathways <- plot_data %>%
    filter(is_significant) %>%
    pull(pathway_id) %>%
    unique()

  plot_data_filtered <- plot_data %>%
    filter(pathway_id %in% sig_pathways) %>%
    filter(is.finite(NES), is.finite(neg_log_padj), !is.na(treatment_display)) %>%
    group_by(pathway_id) %>%
    mutate(max_signal = max(neg_log_padj, na.rm = TRUE), max_abs_nes = max(abs(NES), na.rm = TRUE)) %>%
    ungroup()

  pathway_levels <- plot_data_filtered %>%
    distinct(pathway_id, max_signal, max_abs_nes) %>%
    arrange(desc(max_signal), desc(max_abs_nes), pathway_id) %>%
    pull(pathway_id)

  plot_data_filtered <- plot_data_filtered %>%
    mutate(pathway_id = factor(pathway_id, levels = rev(pathway_levels)))

  ggplot(plot_data_filtered, aes(x = treatment_display, y = pathway_id)) +
    geom_point(
      aes(fill = NES, size = neg_log_padj),
      shape = 21,
      color = "#2E2E2E",
      stroke = 0.15,
      alpha = 0.9
    ) +
    facet_wrap(~ organism, nrow = 1, scales = "free_x") +
    scale_fill_gradient2(
      low = "#3B73B9",
      mid = "white",
      high = "#D1495B",
      midpoint = 0,
      name = "NES\n(direction)"
    ) +
    scale_size_continuous(range = c(1.2, 5.2), name = "-log10(padj)") +
    labs(
      title = "Global pathway enrichment landscape across conditions",
      x = NULL,
      y = "KEGG pathway (map)"
    ) +
    base_plot_theme +
    theme(
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5, margin = margin(t = 4)),
      axis.text.y = element_text(size = 6.7, margin = margin(r = 3)),
      legend.position = "right",
      plot.margin = margin(8, 20, 18, 12)
    )
}

save_plot_set <- function(scatter_plot, heatmap_plot, bubble_plot, composite_plot, out_dir) {
  ggsave(
    filename = file.path(out_dir, "joined_gene_pathway_vs_exomet_scatter.png"),
    plot = scatter_plot,
    width = 16,
    height = 10,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "joined_gene_pathway_vs_exomet_heatmap.png"),
    plot = heatmap_plot,
    width = 15,
    height = 10.5,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "Combined_FGSEA_BubblePlot.png"),
    plot = bubble_plot,
    width = 14,
    height = 8.5,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "Combined_FGSEA_BubblePlot.pdf"),
    plot = bubble_plot,
    width = 14,
    height = 8.5,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "integrated_pathway_vs_exometabolome_figure.png"),
    plot = composite_plot,
    width = 18,
    height = 13.5,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "integrated_pathway_vs_exometabolome_figure.pdf"),
    plot = composite_plot,
    width = 18,
    height = 13.5,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "integrated_pathway_vs_exometabolome_figure_no_panel_subtitles.png"),
    plot = composite_plot,
    width = 20,
    height = 15.5,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "integrated_pathway_vs_exometabolome_figure_no_panel_subtitles.pdf"),
    plot = composite_plot,
    width = 20,
    height = 15.5,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "integrated_pathway_vs_exometabolome_figure_shortlabels.png"),
    plot = composite_plot,
    width = 20,
    height = 15.5,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, "integrated_pathway_vs_exometabolome_figure_shortlabels.pdf"),
    plot = composite_plot,
    width = 20,
    height = 15.5,
    bg = "white"
  )
}

write_figure_legend <- function(out_dir) {
  legend_lines <- c(
    "Figure. Integration of transcriptome pathway activity with exometabolome changes across stresses in Bacillus and Rhodanobacter.",
    "",
    "(a) Pathway-level concordance between transcriptome and exometabolome responses across conditions. Each point represents a KEGG pathway shared between FGSEA and the exometabolome pathway summary. The x-axis shows transcriptome pathway enrichment (FGSEA NES), and the y-axis shows the observed exometabolite pathway score. Red points indicate the same direction of change between transcriptome and exometabolome, blue points indicate opposite direction, and point size reflects the number of significant metabolites contributing to the pathway score.",
    "",
    "(b) Top informative pathways within each condition. Heatmaps show the leading pathways ranked by combined transcriptome and exometabolome signal. Values are scaled within each metric for visualization, with red indicating higher positive values and blue indicating lower negative values.",
    "",
    "(c) Global pathway enrichment landscape across conditions. Bubble plots summarize FGSEA results across all conditions in Rhodanobacter and Bacillus. Bubble color indicates NES direction and magnitude, and bubble size indicates significance as -log10(adjusted P value).",
    "",
    "Abbreviations: P, phage; K, kanamycin; Al, aluminum."
  )

  writeLines(legend_lines, file.path(out_dir, "integrated_pathway_vs_exometabolome_figure_legend.txt"))
}

# ---------------------------
# Load data
# ---------------------------
rhodano_fgsea <- load_fgsea(resolve_existing_path(rhodano_fgsea_file), "Rhodanobacter", "R")
bacillus_fgsea <- load_fgsea(resolve_existing_path(bacillus_fgsea_file), "Bacillus", "B")
metabolite_summary <- load_metabolite_summary(resolve_existing_path(metabolite_summary_file))

fgsea_all <- bind_rows(rhodano_fgsea, bacillus_fgsea)

if (!is.na(allowed_pathway_ids_file) && file.exists(allowed_pathway_ids_file)) {
  keep_ids <- read_lines(allowed_pathway_ids_file) %>%
    str_trim() %>%
    { .[. != ""] } %>%
    unique()
  fgsea_all <- fgsea_all %>% filter(pathway_id %in% keep_ids)
  metabolite_summary <- metabolite_summary %>% filter(pathway_id %in% keep_ids)
}

joined <- fgsea_all %>%
  inner_join(metabolite_summary, by = c("phylogeny", "treatment", "pathway_id")) %>%
  mutate(
    organism = factor(organism, levels = organism_order),
    treatment = factor(treatment, levels = treatment_order),
    panel_label = factor(paste(organism, treatment, sep = " | "), levels = panel_levels),
    frac_sig = if_else(n_mets > 0, n_sig / n_mets, 0),
    met_score = mean_log2FC * frac_sig,
    gene_sig = !is.na(padj) & padj <= padj_thresh & abs(NES) >= abs_nes_thresh,
    met_sig = !is.na(n_mets) & !is.na(n_sig) & n_mets >= min_mets & n_sig >= min_sig_mets,
    informative = gene_sig & met_sig,
    relationship = case_when(
      !is.finite(NES) | !is.finite(met_score) ~ "missing",
      abs(NES) < 1e-12 | abs(met_score) < 1e-12 ~ "weak/zero",
      sign(NES) == sign(met_score) ~ "same direction",
      TRUE ~ "opposite direction"
    ),
    combined_rank = abs(NES) + abs(met_score)
  ) %>%
  filter(!is.na(panel_label))

summary_tbl <- joined %>%
  group_by(organism, phylogeny, treatment) %>%
  summarise(
    n_pathways_overlap = n(),
    n_gene_sig = sum(gene_sig, na.rm = TRUE),
    n_met_sig = sum(met_sig, na.rm = TRUE),
    n_informative = sum(informative, na.rm = TRUE),
    spearman_all = safe_cor(NES, met_score, method = "spearman"),
    spearman_informative = safe_cor(NES[informative], met_score[informative], method = "spearman"),
    concordance_all = safe_concordance(NES, met_score),
    concordance_informative = safe_concordance(NES[informative], met_score[informative]),
    .groups = "drop"
  ) %>%
  arrange(organism, treatment)

relationship_summary <- joined %>%
  filter(informative) %>%
  count(organism, treatment, panel_label, relationship, name = "n_pathways") %>%
  group_by(organism, treatment, panel_label) %>%
  mutate(frac_pathways = n_pathways / sum(n_pathways)) %>%
  ungroup() %>%
  arrange(organism, treatment, desc(n_pathways))

write_csv(joined, file.path(out_dir, "joined_gene_pathway_vs_exomet_pathway_all.csv"))
write_csv(summary_tbl, file.path(out_dir, "joined_gene_pathway_vs_exomet_pathway_summary.csv"))
write_csv(relationship_summary, file.path(out_dir, "joined_gene_pathway_vs_exomet_relationship_summary.csv"))

heatmap_data <- make_heatmap_data(joined, n_show_heatmap)
write_csv(
  heatmap_data$heatmap_top %>%
    select(organism, phylogeny, treatment, pathway_id, padj, NES, n_mets, n_sig, frac_sig, mean_log2FC, met_score, relationship, top_up, top_down),
  file.path(out_dir, "top_informative_joined_pathways_for_heatmap.csv")
)

scatter_plot <- make_scatter_plot(joined, label_top_n)
heatmap_plot <- make_heatmap_plot(heatmap_data$heatmap_long)
bubble_plot <- make_bubble_plot(fgsea_all)

top_row <- scatter_plot + heatmap_plot + plot_layout(widths = c(2.0, 1.15))
composite_plot <- top_row / bubble_plot +
  plot_layout(heights = c(1.0, 1.0)) +
  plot_annotation(
    title = "Integration of transcriptome pathway activity with exometabolome changes across stresses in Bacillus and Rhodanobacter",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.tag = element_text(face = "bold", size = 15),
      plot.tag.position = c(0, 1)
    )
  )

save_plot_set(scatter_plot, heatmap_plot, bubble_plot, composite_plot, out_dir)
write_figure_legend(out_dir)

message("Wrote outputs to: ", out_dir)
message("- joined_gene_pathway_vs_exomet_pathway_all.csv")
message("- joined_gene_pathway_vs_exomet_pathway_summary.csv")
message("- joined_gene_pathway_vs_exomet_relationship_summary.csv")
message("- top_informative_joined_pathways_for_heatmap.csv")
message("- joined_gene_pathway_vs_exomet_scatter.png")
message("- joined_gene_pathway_vs_exomet_heatmap.png")
message("- Combined_FGSEA_BubblePlot.png")
message("- Combined_FGSEA_BubblePlot.pdf")
message("- integrated_pathway_vs_exometabolome_figure.png")
message("- integrated_pathway_vs_exometabolome_figure.pdf")
message("- integrated_pathway_vs_exometabolome_figure_shortlabels.png")
message("- integrated_pathway_vs_exometabolome_figure_shortlabels.pdf")
message("- integrated_pathway_vs_exometabolome_figure_legend.txt")

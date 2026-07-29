#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

theme_nature <- function(base_size = 7.7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(face = "bold", size = base_size, colour = "#111827"),
      axis.text = element_text(size = base_size - 0.35, colour = "#1F2937"),
      plot.title = element_text(face = "bold", size = base_size + 1.1, colour = "#111827"),
      plot.subtitle = element_text(size = base_size - 0.15, colour = "#4B5563"),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text = element_text(face = "bold", size = base_size - 0.05, colour = "#111827"),
      legend.title = element_text(face = "bold", size = base_size - 0.1, colour = "#111827"),
      legend.text = element_text(size = base_size - 0.2, colour = "#1F2937"),
      legend.key.height = unit(3.3, "mm"),
      legend.key.width = unit(8, "mm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_pub <- function(plot, filename, width_mm = 203, height_mm = 152, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot)
  dev.off()
}

clean_met_label <- function(x) {
  x %>%
    str_remove("_(positive|negative)$") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("^2 4-dihydroxypteridine$", "2,4-dihydroxypteridine") %>%
    str_replace_all("^3 4-dihydroxybenzoic acid$", "3,4-dihydroxybenzoic acid") %>%
    str_squish()
}

theme_set(theme_nature())

out_dir <- "/Users/mingfeichen/Manuscript/outputs/manual-20260608-a2/presentations/species_transport_link/assets"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

met_file <- "/Users/mingfeichen/metabolite_log2FC_ttests_vs0hr.csv"
tx_file <- "/Users/mingfeichen/Manuscript/outputs/manual-20260603-a4/presentations/transcriptome_pathway_modules/assets/transcriptome_pathway_module_data.csv"

met_superset <- tibble::tribble(
  ~feature, ~display, ~superclass, ~order,
  "Betaine_positive", "Betaine", "Lipids & osmolytes", 1,
  "Carnitine_positive", "Carnitine", "Lipids & osmolytes", 2,
  "sn-Glycero-3-Phosphocholine_positive", "sn-Glycero-3-Phosphocholine", "Lipids & osmolytes", 3,
  "N-trimethyllysine_positive", "N-trimethyllysine", "Amino acids & derivatives", 4,
  "4-guanidinobutanoic_acid_positive", "4-guanidinobutanoic acid", "Polyamines & guanidino compounds", 5,
  "Adenosine_positive", "Adenosine", "Nucleic acids", 6,
  "Guanosine_negative", "Guanosine", "Nucleic acids", 7,
  "Cytosine_positive", "Cytosine", "Nucleic acids", 8,
  "Xanthine_positive", "Xanthine", "Nucleic acids", 9,
  "Uracil_negative", "Uracil", "Nucleic acids", 10,
  "3-methyladenine_positive", "3-methyladenine", "Nucleic acids", 11,
  "5-methylcytosine_positive", "5-methylcytosine", "Nucleic acids", 12,
  "Hypoxanthine_positive", "Hypoxanthine", "Nucleic acids", 13,
  "pterin_positive", "pterin", "Cofactors & vitamins", 14,
  "2_4-dihydroxypteridine_negative", "2,4-dihydroxypteridine", "Cofactors & vitamins", 15,
  "glycerol_2-phosphoric_acid_negative", "glycerol 2-phosphoric acid", "Central carbon & organic acids", 16,
  "2-hydroxybutyric_acid_negative", "2-hydroxybutyric acid", "Central carbon & organic acids", 17,
  "4-hydroxybenzoic_acid_negative", "4-hydroxybenzoic acid", "Aromatic compounds & phenolics", 18,
  "orotic_acid_negative", "orotic acid", "Nucleic acids", 19,
  "N-Acetyl-Glutamine_negative", "N-Acetyl-Glutamine", "Amino acids & derivatives", 20,
  "trans-4-hydroxyproline_positive", "trans-4-hydroxyproline", "Amino acids & derivatives", 21,
  "5-oxo-proline_negative", "5-oxo-proline", "Amino acids & derivatives", 22,
  "Guanine_positive", "Guanine", "Nucleic acids", 23,
  "4-aminobutanoic_acid_positive", "4-aminobutanoic acid", "Amino acids & derivatives", 24,
  "agmatine_sulfuric_acid_positive", "agmatine sulfuric acid", "Polyamines & guanidino compounds", 25,
  "creatine_positive", "creatine", "Amino acids & derivatives", 26,
  "methionine_positive", "methionine", "Amino acids & derivatives", 27
)

species_configs <- list(
  Bacillus = list(
    phylogeny = "B",
    title = "Bacillus",
    subtitle = "log2FC vs 0h; 8h and 24h recoded to mid/late.",
    met_treatment_map = c("8hr" = "mid", "24hr" = "late", "Al" = "Al", "K" = "Kana", "P" = "Phage"),
    met_features = c(
      "Betaine_positive",
      "Carnitine_positive",
      "sn-Glycero-3-Phosphocholine_positive",
      "N-trimethyllysine_positive",
      "4-guanidinobutanoic_acid_positive",
      "Adenosine_positive",
      "Guanosine_negative",
      "Cytosine_positive",
      "Xanthine_positive",
      "Uracil_negative",
      "3-methyladenine_positive",
      "pterin_positive",
      "2_4-dihydroxypteridine_negative",
      "glycerol_2-phosphoric_acid_negative",
      "2-hydroxybutyric_acid_negative",
      "4-hydroxybenzoic_acid_negative"
    ),
    tx_levels = c(
      "ABC transporters (map02010)",
      "Biosynthesis of secondary metabolites (map01110)"
    ),
    tx_focus = c("map02010", "map01110"),
    tx_title = "Transporters"
  ),
  Rhodanobacter = list(
    phylogeny = "R",
    title = "Rhodanobacter",
    subtitle = "log2FC vs 0h; 24h and 48h recoded to mid/late.",
    met_treatment_map = c("24hr" = "mid", "48hr" = "late", "Al" = "Al", "K" = "Kana", "P" = "Phage"),
    met_features = c(
      "5-methylcytosine_positive",
      "Xanthine_positive",
      "Hypoxanthine_positive",
      "2_4-dihydroxypteridine_negative",
      "trans-4-hydroxyproline_positive",
      "5-oxo-proline_negative",
      "Cytosine_positive",
      "glycerol_2-phosphoric_acid_negative",
      "pterin_positive",
      "N-Acetyl-Glutamic_acid_negative",
      "orotic_acid_negative",
      "3-methyladenine_positive",
      "Guanine_positive",
      "Uracil_negative",
      "4-guanidinobutanoic_acid_positive",
      "N-Acetyl-Glutamine_negative"
    ),
    tx_levels = c(
      "Two-component system (map02020)",
      "Biosynthesis of secondary metabolites (map01110)"
    ),
    tx_focus = c("map02020", "map01110"),
    tx_title = "Sensing / transport"
  )
)

build_species_figure <- function(species_name, config, met_raw, tx_raw) {
  met_species <- met_raw %>%
    filter(phylogeny == config$phylogeny) %>%
    mutate(
      treatment_display = recode(treatment, !!!config$met_treatment_map),
      treatment_display = factor(treatment_display, levels = c("mid", "late", "Al", "Kana", "Phage")),
      feature_feature = feature,
      feature_display = clean_met_label(feature),
      feature_display = factor(
        feature_display,
        levels = rev(clean_met_label(met_superset$display[match(config$met_features, met_superset$feature)]))
      )
    ) %>%
    filter(feature %in% config$met_features)

  met_catalog <- met_superset %>%
    filter(feature %in% config$met_features) %>%
    mutate(
      feature_display = clean_met_label(display),
      feature_display = factor(feature_display, levels = rev(clean_met_label(display[order(order)]))),
      superclass = factor(
        superclass,
        levels = c(
          "Lipids & osmolytes",
          "Amino acids & derivatives",
          "Polyamines & guanidino compounds",
          "Nucleic acids",
          "Cofactors & vitamins",
          "Central carbon & organic acids",
          "Aromatic compounds & phenolics"
        )
      )
    )

  met_plot <- met_catalog %>%
    select(feature, feature_display, superclass, order) %>%
    tidyr::expand_grid(treatment_display = factor(c("mid", "late", "Al", "Kana", "Phage"), levels = c("mid", "late", "Al", "Kana", "Phage"))) %>%
    left_join(
      met_species %>%
        select(feature, treatment_display, log2FC, pvalue),
      by = c("feature", "treatment_display")
    ) %>%
    mutate(
      feature_display = factor(feature_display, levels = rev(unique(met_catalog$feature_display[order(met_catalog$order)])))
    )

  panel_met <- ggplot(met_plot, aes(x = treatment_display, y = feature_display)) +
    geom_tile(aes(fill = log2FC), colour = "white", linewidth = 0.22) +
    facet_grid(superclass ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_fill_gradient2(
      low = "#3B73B9",
      mid = "white",
      high = "#D1495B",
      midpoint = 0,
      name = "log2FC\n(vs 0h)"
    ) +
    labs(
      title = paste0(species_name, " exometabolites"),
      subtitle = config$subtitle,
      x = NULL,
      y = "Representative metabolites"
    ) +
    theme_nature(base_size = 7.2) +
    theme(
      plot.title = element_text(size = 8.9, face = "bold", hjust = 0),
      axis.text.x = element_text(size = 7.0),
      axis.text.y = element_text(size = 6.8),
      strip.text.y.left = element_text(angle = 0, size = 7.1),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major = element_line(colour = "#EEF2F7", linewidth = 0.22),
      panel.spacing.y = unit(0.10, "in"),
      plot.margin = margin(3, 7, 3, 7)
    )

  tx_species <- tx_raw %>%
    mutate(
      treatment_display = str_remove(treatment_display, "^[BR] "),
      treatment_display = factor(treatment_display, levels = c("mid", "late", "Al", "Kana", "Phage"))
    ) %>%
    filter(str_detect(pathway_id, "^map"), pathway_id %in% config$tx_focus, organism == species_name) %>%
    group_by(pathway_id, treatment_display) %>%
    arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      pathway_label = paste0(pathway_name, " (", pathway_id, ")"),
      pathway_label = factor(
        pathway_label,
        levels = rev(config$tx_levels)
      ),
      pathway_class = factor(
        pathway_class,
        levels = c("Transport / signaling", "Secondary metabolism / transport")
      )
    )

  tx_plot <- tidyr::expand_grid(
    tx_species %>% distinct(pathway_label, pathway_class),
    treatment_display = factor(c("mid", "late", "Al", "Kana", "Phage"), levels = c("mid", "late", "Al", "Kana", "Phage"))
  ) %>%
    left_join(
      tx_species %>%
        select(pathway_label, pathway_class, treatment_display, NES, padj),
      by = c("pathway_label", "pathway_class", "treatment_display")
    )

  panel_tx <- ggplot(tx_plot, aes(x = treatment_display, y = pathway_label)) +
    geom_tile(fill = "#FAFBFD", colour = "white", linewidth = 0.22) +
    geom_point(
      aes(fill = NES, size = -log10(padj)),
      shape = 21, colour = "#1F2937", stroke = 0.22, alpha = 0.96
    ) +
    facet_grid(pathway_class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_fill_gradient2(
      low = "#3B73B9",
      mid = "white",
      high = "#D1495B",
      midpoint = 0,
      name = "NES"
    ) +
    scale_size_continuous(range = c(1.8, 5.0), name = "-log10(padj)") +
    labs(
      title = paste0(species_name, " ", config$tx_title),
      x = NULL,
      y = "Selected pathways"
    ) +
    theme_nature(base_size = 7.2) +
    theme(
      plot.title = element_text(size = 8.9, face = "bold", hjust = 0),
      axis.text.x = element_text(size = 7.0),
      axis.text.y = element_text(size = 6.8),
      strip.text.y.left = element_text(angle = 0, size = 7.1),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major = element_line(colour = "#EEF2F7", linewidth = 0.22),
      panel.spacing.y = unit(0.10, "in"),
      plot.margin = margin(3, 7, 3, 7)
    )

  figure <- panel_met | panel_tx +
    plot_layout(widths = c(1.12, 0.88), guides = "collect") +
    plot_annotation(
      title = paste0(species_name, ": exometabolite release linked to transport remodeling"),
      subtitle = config$subtitle,
      tag_levels = "a",
      theme = theme(
        plot.title = element_text(face = "bold", size = 11.0, colour = "#111827", hjust = 0.5),
        plot.subtitle = element_text(size = 8.0, colour = "#4B5563", hjust = 0.5)
      )
    ) &
    theme(plot.tag = element_text(face = "bold", size = 9))

  list(
    figure = figure,
    panel_met = panel_met,
    panel_tx = panel_tx,
    met_plot = met_plot,
    tx_plot = tx_plot
  )
}

met_raw <- read_csv(met_file, show_col_types = FALSE)
tx_raw <- read_csv(tx_file, show_col_types = FALSE)

for (species_name in names(species_configs)) {
  config <- species_configs[[species_name]]
  res <- build_species_figure(species_name, config, met_raw, tx_raw)
  met_prefix <- file.path(out_dir, paste0(tolower(species_name), "_metabolites_vs0h"))
  tx_prefix <- file.path(out_dir, paste0(tolower(species_name), "_transporters"))
  write_csv(res$met_plot, paste0(met_prefix, ".csv"))
  write_csv(res$tx_plot, paste0(tx_prefix, ".csv"))
  save_pub(res$panel_met, met_prefix, width_mm = 203, height_mm = 152)
  save_pub(res$panel_tx, tx_prefix, width_mm = 203, height_mm = 152)
  message("Wrote: ", paste0(met_prefix, ".png"))
  message("Wrote: ", paste0(met_prefix, ".pdf"))
  message("Wrote: ", paste0(met_prefix, ".svg"))
  message("Wrote: ", paste0(tx_prefix, ".png"))
  message("Wrote: ", paste0(tx_prefix, ".pdf"))
  message("Wrote: ", paste0(tx_prefix, ".svg"))
}

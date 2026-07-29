#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

theme_slide_heatmap <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = base_size, face = "bold", colour = "#1F2937"),
      axis.text = element_text(size = base_size - 0.2, colour = "#1F2937"),
      axis.text.x = element_text(size = base_size - 0.2, face = "bold"),
      axis.text.y = element_text(size = base_size - 0.8),
      plot.title = element_text(size = base_size + 1.8, face = "bold", colour = "#111827"),
      plot.subtitle = element_text(size = base_size - 0.1, colour = "#4B5563"),
      strip.background = element_rect(fill = "#F3F4F6", colour = NA),
      strip.text.y.left = element_text(size = base_size - 0.1, face = "bold", angle = 0),
      legend.title = element_text(size = base_size - 0.1, face = "bold"),
      legend.text = element_text(size = base_size - 0.3),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.justification = "left",
      panel.grid = element_blank(),
      panel.spacing.y = grid::unit(2.0, "mm"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

save_pub_r <- function(plot, filename, width_mm = 380, height_mm = 165, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
}

wrap_label <- function(x, width = 20) {
  stringr::str_wrap(x, width = width)
}

root <- "/Users/mingfeichen/Manuscript"
infile <- file.path(root, "metabolite_adsorption_merged_long.csv")
out_dir <- file.path(getwd(), "outputs", "manual-20260601-a7", "presentations", "metabolite_adsorption_heatmap", "assets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_base <- file.path(out_dir, "metabolite_adsorption_heatmap_slide")

if (!file.exists(infile)) {
  stop("Input file not found: ", infile, call. = FALSE)
}

group_order <- c(
  "clay",
  "iron mineral",
  "Bacillus sediment",
  "Rhodano sediment"
)

group_labels <- c(
  "clay" = "Clay",
  "iron mineral" = "Iron mineral",
  "Bacillus sediment" = "Bacillus sediment",
  "Rhodano sediment" = "Rhodano sediment"
)

chem_order <- c(
  "Amino acids & peptides",
  "Nucleosides & bases",
  "Organic acids",
  "Polyamines & amines",
  "Vitamins & cofactors",
  "Aromatic / benzenoids",
  "Other / synthetic"
)

df <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    value = as.numeric(value),
    group = factor(group, levels = group_order),
    chem_group = factor(chem_group, levels = chem_order),
    group_label = factor(group_labels[as.character(group)], levels = unname(group_labels[group_order])),
    metabolite_label = wrap_label(metabolite, width = 20)
  )

key_summary <- df %>%
  group_by(key) %>%
  summarise(
    metabolite = first(metabolite),
    metabolite_label = first(metabolite_label),
    chem_group = first(chem_group),
    n_groups = n_distinct(group[value > 0]),
    max_value = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_groups >= 2, !is.na(metabolite_label), metabolite_label != "NA") %>%
  arrange(chem_group, desc(max_value), metabolite)

ordered_keys <- key_summary$key
label_map <- key_summary %>%
  distinct(key, metabolite_label)

plot_df <- df %>%
  filter(key %in% ordered_keys, !is.na(metabolite_label), metabolite_label != "NA") %>%
  select(-metabolite_label) %>%
  left_join(label_map, by = "key") %>%
  mutate(
    metabolite_label = factor(metabolite_label, levels = rev(unique(key_summary$metabolite_label))),
    group_label = factor(group_labels[as.character(group)], levels = unname(group_labels[group_order]))
  ) %>%
  filter(!is.na(metabolite_label))

p <- ggplot(plot_df, aes(x = group_label, y = metabolite_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(chem_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(
    colours = c("#F8FBFF", "#DCEAF6", "#A8CCE6", "#5C97C8", "#1F5D99"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    name = "Adsorption"
  ) +
  labs(
    title = "Adsorption profiles of retained metabolites across clay, iron mineral, and sediment samples",
    subtitle = "Only metabolites present in at least two groups are shown; values are percentages",
    x = "Source",
    y = NULL
  ) +
  theme_slide_heatmap(base_size = 9) +
  theme(
    axis.text.x = element_text(size = 8.3, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 6.1),
    strip.text.y.left = element_text(size = 8.8, face = "bold", angle = 0),
    legend.key.width = grid::unit(8, "mm"),
    legend.key.height = grid::unit(3.5, "mm")
  )

save_pub_r(p, out_base, width_mm = 380, height_mm = 170)

message("Wrote: ", out_base, ".png")
message("Wrote: ", out_base, ".pdf")
message("Wrote: ", out_base, ".svg")

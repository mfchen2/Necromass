#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(svglite)
  library(ragg)
})

theme_nature_contract <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.2, colour = "#222222"),
      axis.text.y = element_text(size = base_size - 0.6),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.4),
      strip.text.y.left = element_text(size = base_size - 0.1, face = "bold", angle = 0),
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      plot.title = element_text(size = base_size + 0.7, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.1, colour = "#4A4A4A"),
      panel.grid.major.x = element_line(colour = "#E9E9E9", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = grid::unit(2.2, "mm"),
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.margin = margin(6, 8, 8, 4),
      strip.placement = "outside"
    )
}

save_pub_r <- function(plot, filename, width_mm = 190, height_mm = 230, dpi = 600) {
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

wrap_label <- function(x, width = 18) {
  paste(strwrap(x, width = width), collapse = "\n")
}

root <- "/Users/mingfeichen/Manuscript"
infile <- file.path(root, "metabolite_adsorption_merged_long.csv")
out_base <- file.path(root, "metabolite_adsorption_merged_barplot")

if (!file.exists(infile)) {
  stop("Input file not found: ", infile, call. = FALSE)
}

df <- read.csv(infile, stringsAsFactors = FALSE) %>%
  mutate(
    value = as.numeric(value),
    group = factor(group, levels = c("clay", "iron mineral", "Bacillus sediment", "Rhodano sediment")),
    chem_group = factor(
      chem_group,
      levels = c(
        "Amino acids & peptides",
        "Nucleosides & bases",
        "Organic acids",
        "Polyamines & amines",
        "Vitamins & cofactors",
        "Aromatic / benzenoids",
        "Carbohydrates & sugars",
        "Polyols & glycans",
        "Other / synthetic"
      )
    )
  )

key_summary <- df %>%
  group_by(key) %>%
  summarise(
    metabolite = first(metabolite),
    chem_group = first(chem_group),
    n_groups = n_distinct(group[value > 0]),
    max_value = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_groups >= 2) %>%
  arrange(chem_group, desc(max_value), metabolite) %>%
  mutate(
    metabolite_label = vapply(metabolite, wrap_label, character(1), width = 18)
  )

ordered_keys <- key_summary$key
label_map <- setNames(key_summary$metabolite_label, key_summary$key)

plot_df <- df %>%
  filter(key %in% ordered_keys) %>%
  mutate(key = factor(key, levels = ordered_keys))

fill_cols <- c(
  "clay" = "#D58B4D",
  "iron mineral" = "#2A9D8F",
  "Bacillus sediment" = "#5B7996",
  "Rhodano sediment" = "#E76F51"
)

labeller_fun <- labeller(
  chem_group = label_wrap_gen(width = 16)
)

p <- ggplot(plot_df, aes(x = key, y = value, fill = group)) +
  geom_col(
    width = 0.72,
    position = position_dodge2(width = 0.78, preserve = "single", padding = 0.08),
    colour = "white",
    linewidth = 0.18
  ) +
  facet_grid(
    . ~ chem_group,
    scales = "free_x",
    space = "free_x",
    drop = TRUE,
    labeller = labeller_fun
  ) +
  scale_x_discrete(labels = label_map, expand = expansion(add = 0.3)) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.03)),
    labels = label_number(accuracy = 1)
  ) +
  scale_fill_manual(values = fill_cols, drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0.18))) +
  labs(
    title = NULL,
    subtitle = "Metabolites present in at least two groups; necromass values rescaled to percentage",
    x = "Metabolite",
    y = "Adsorption (%)",
    fill = NULL
  ) +
  theme_nature_contract() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.2, margin = margin(t = 2)),
    axis.text.y = element_text(size = 6.5),
    legend.margin = margin(t = 1, b = 1),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(5.8, "mm")
  )

save_pub_r(p, out_base, width_mm = 380, height_mm = 110)

message("Wrote: ", out_base, ".png")
message("Wrote: ", out_base, ".pdf")
message("Wrote: ", out_base, ".svg")

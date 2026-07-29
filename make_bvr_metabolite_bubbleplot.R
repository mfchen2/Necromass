#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
})

root <- "/Users/mingfeichen"
source_path <- file.path(root, "B_vs_R_ttest_targeted_metabolites_013026.csv")
out_dir <- file.path(
  root,
  "Manuscript",
  "outputs",
  "manual-20260608-a3",
  "presentations",
  "bvr_metabolite_bubbleplot",
  "assets"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_nature <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(face = "plain"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(size = rel(0.92)),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.42, "cm"),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey88", linewidth = 0.25),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = rel(0.9)),
      strip.placement = "outside"
    )
}

save_pub <- function(plot, basename, width = 14, height = 8, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      file.path(out_dir, paste0(basename, ".png")),
      plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      device = ragg::agg_png
    )
  } else {
    ggsave(file.path(out_dir, paste0(basename, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  ggsave(
    file.path(out_dir, paste0(basename, ".pdf")),
    plot,
    width = width,
    height = height,
    bg = "white",
    device = grDevices::cairo_pdf
  )
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(out_dir, paste0(basename, ".svg")), plot, width = width, height = height, device = svglite::svglite, bg = "white")
  }
}

clean_feature <- function(x) {
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\b2 4\\b", "2,4")
  x <- str_replace_all(x, "\\b3 4\\b", "3,4")
  x <- str_squish(x)
  x
}

contrast_to_label <- c(
  "B_0hr vs R_0hr" = "0h",
  "B_8hr vs R_24hr" = "mid",
  "B_24hr vs R_48hr" = "late",
  "B_Al vs R_Al" = "Al",
  "B_K vs R_K" = "Kanamycin",
  "B_P vs R_P" = "Phage"
)

superclass_order <- c(
  "Amino acids & derivatives",
  "Nucleic acids (bases/nucleosides/nucleotides)",
  "Central carbon & organic acids",
  "Lipids & osmolytes (quaternary amines)",
  "Polyamines & guanidino compounds",
  "Aromatic compounds & phenolics",
  "Cofactors (pterins & vitamins)",
  "Carbohydrates & amino sugars",
  "Methylation & sulfur salvage"
)

dat <- read_csv(source_path, show_col_types = FALSE) %>%
  mutate(
    contrast_label = recode(contrast, !!!contrast_to_label),
    sig = !is.na(padj) & padj < 0.05 & ttest_status == "ok",
    feature_label = clean_feature(feature),
    superclass = factor(superclass, levels = superclass_order)
  )

# Keep only recurrent significant features, then select a compact set for the plot.
feature_summary <- dat %>%
  filter(sig) %>%
  group_by(feature, feature_label, superclass) %>%
  summarise(
    n_sig = n_distinct(contrast),
    min_padj = min(padj, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig), min_padj, desc(max_abs_log2FC), superclass, feature_label)

selected_features <- feature_summary %>%
  filter(n_sig >= 4) %>%
  slice_head(n = 20)

if (nrow(selected_features) < 12) {
  selected_features <- feature_summary %>% slice_head(n = 20)
}

plot_df <- dat %>%
  semi_join(selected_features, by = c("feature", "feature_label", "superclass")) %>%
  mutate(
    feature_label = factor(
      feature_label,
      levels = selected_features %>%
        arrange(superclass, desc(n_sig), min_padj, desc(max_abs_log2FC), feature_label) %>%
        pull(feature_label) %>%
        unique()
    ),
    contrast_label = factor(contrast_label, levels = c("0h", "mid", "late", "Al", "Kanamycin", "Phage"))
  ) %>%
  arrange(superclass, feature_label, contrast_label)

max_abs <- max(abs(plot_df$log2FC), na.rm = TRUE)
max_abs <- max(max_abs, 8)

base <- ggplot(plot_df, aes(x = contrast_label, y = feature_label)) +
  geom_point(
    data = plot_df,
    aes(fill = log2FC, size = -log10(padj)),
    shape = 21,
    colour = "grey70",
    stroke = 0.15,
    alpha = 0.92
  ) +
  geom_point(
    data = filter(plot_df, sig),
    aes(fill = log2FC, size = -log10(padj)),
    shape = 21,
    colour = "black",
    stroke = 0.48,
    alpha = 1
  ) +
  facet_grid(superclass ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "#F7F7F7",
    high = "#C53030",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    oob = squish,
    name = expression(log[2]~FC~"(Bacillus vs Rhodanobacter)")
  ) +
  scale_size_continuous(
    range = c(1.8, 7.2),
    breaks = c(1, 2, 4, 6),
    name = expression(-log[10](padj))
  ) +
  labs(
    title = "Targeted metabolite differences between Bacillus and Rhodanobacter",
    subtitle = "Selected metabolites recur across at least four of the six contrasts; thick outlines indicate padj < 0.05.",
    x = NULL,
    y = NULL
  ) +
  theme_nature(base_size = 8.5) +
  theme(
    axis.text.x = element_text(size = 8.2, angle = 0, vjust = 0.5),
    axis.text.y = element_text(size = 7.8, lineheight = 0.95),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.1),
    legend.position = "right",
    legend.title = element_text(size = 8.2),
    legend.text = element_text(size = 7.7),
    legend.key.width = unit(0.45, "cm"),
    panel.spacing.y = unit(0.18, "lines"),
    plot.margin = margin(6, 10, 6, 6)
  )

write_csv(
  selected_features,
  file.path(out_dir, "bvr_metabolite_bubbleplot_selected_features.csv")
)
write_csv(
  plot_df,
  file.path(out_dir, "bvr_metabolite_bubbleplot_data.csv")
)

save_pub(base, "bvr_metabolite_bubbleplot", width = 14, height = 8, dpi = 600)

message("Wrote bubble plot outputs to: ", out_dir)

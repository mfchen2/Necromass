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

root <- "/Users/mingfeichen"
out_dir <- file.path(
  getwd(),
  "outputs",
  "manual-20260601-a6",
  "presentations",
  "metabolite_ttest_remake",
  "assets"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.05), margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(0.92), margin = margin(b = 6)),
      axis.title = element_text(face = "plain"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.border = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", size = rel(0.9))
    )
}

save_all <- function(plot, basename, width, height) {
  ggsave(file.path(out_dir, paste0(basename, ".png")), plot, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(basename, ".pdf")), plot, width = width, height = height)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(out_dir, paste0(basename, ".svg")), plot, width = width, height = height, device = svglite::svglite)
  }
}

normalize_feature <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  x <- str_replace_all(x, "_negative$|_positive$", "")
  x <- str_replace_all(x, "^malic_acid_peak_1$", "malic_acid")
  x
}

pretty_feature <- function(x) {
  x <- str_trim(as.character(x))
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x
}

wrap_feature <- function(x, width = 20) {
  vapply(str_wrap(x, width = width), identity, character(1))
}

treatment_display_timecourse <- function(phylogeny, treatment) {
  case_when(
    phylogeny == "B" & treatment == "8hr" ~ "mid",
    phylogeny == "B" & treatment == "24hr" ~ "late",
    phylogeny == "B" & treatment == "P" ~ "Phage",
    phylogeny == "B" & treatment == "K" ~ "Kana",
    phylogeny == "B" & treatment == "Al" ~ "Al",
    phylogeny == "R" & treatment == "24hr" ~ "mid",
    phylogeny == "R" & treatment == "48hr" ~ "late",
    phylogeny == "R" & treatment == "P" ~ "Phage",
    phylogeny == "R" & treatment == "K" ~ "Kana",
    phylogeny == "R" & treatment == "Al" ~ "Al",
    TRUE ~ NA_character_
  )
}

contrast_display <- function(x) {
  recode(
    x,
    `B_0hr vs R_0hr` = "0h",
    `B_8hr vs R_24hr` = "mid",
    `B_24hr vs R_48hr` = "late",
    `B_P vs R_P` = "Phage",
    `B_K vs R_K` = "Kana",
    `B_Al vs R_Al` = "Al"
  )
}

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

superclass_palette <- c(
  "Amino acids & derivatives" = "#F58518",
  "Nucleic acids (bases/nucleosides/nucleotides)" = "#4C78A8",
  "Central carbon & organic acids" = "#54A24B",
  "Lipids & osmolytes (quaternary amines)" = "#E45756",
  "Polyamines & guanidino compounds" = "#B279A2",
  "Aromatic compounds & phenolics" = "#72B7B2",
  "Cofactors (pterins & vitamins)" = "#FF9DA6",
  "Carbohydrates & amino sugars" = "#1F9E89",
  "Methylation & sulfur salvage" = "#9D9D9D"
)

timecourse_path <- file.path(root, "metabolite_log2FC_ttests_vs0hr.csv")
bvr_path <- file.path(root, "B_vs_R_ttest_targeted_metabolites_013026.csv")

timecourse <- read_csv(timecourse_path, show_col_types = FALSE) %>%
  mutate(
    key = normalize_feature(feature),
    phylogeny_label = recode(phylogeny, B = "Bacillus", R = "Rhodanobacter"),
    treatment_display = treatment_display_timecourse(phylogeny, treatment),
    direction = if_else(log2FC >= 0, "Enriched", "Depleted")
  )

bvr <- read_csv(bvr_path, show_col_types = FALSE) %>%
  rename_with(~ if_else(.x == "", "row_id", .x)) %>%
  mutate(
    key = normalize_feature(feature),
    contrast_display = contrast_display(contrast),
    direction = if_else(log2FC >= 0, "B enriched", "R enriched")
  )

superclass_map <- bvr %>%
  distinct(key, superclass) %>%
  filter(!is.na(superclass)) %>%
  mutate(superclass = factor(superclass, levels = superclass_order))

timecourse_sig <- timecourse %>%
  left_join(superclass_map, by = "key") %>%
  filter(ttest_status == "ok", !is.na(pvalue), pvalue < 0.05, !is.na(superclass)) %>%
  mutate(
    treatment_display = factor(treatment_display, levels = c("mid", "late", "Phage", "Kana", "Al")),
    direction = factor(direction, levels = c("Enriched", "Depleted")),
    phylogeny_label = factor(phylogeny_label, levels = c("Bacillus", "Rhodanobacter"))
  )

bar_summary <- timecourse_sig %>%
  count(phylogeny_label, direction, treatment_display, superclass, name = "n") %>%
  mutate(
    treatment_display = factor(treatment_display, levels = c("mid", "late", "Phage", "Kana", "Al")),
    direction = factor(direction, levels = c("Enriched", "Depleted")),
    phylogeny_label = factor(phylogeny_label, levels = c("Bacillus", "Rhodanobacter")),
    superclass = factor(superclass, levels = superclass_order)
  )

bar_plot <- ggplot(bar_summary, aes(x = treatment_display, y = n, fill = superclass)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.2) +
  facet_grid(direction ~ phylogeny_label, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = superclass_palette, drop = FALSE) +
  labs(
    title = "Significant targeted metabolites cluster into different chemical superclasses",
    subtitle = "Counts are based on t-tests versus the 0 h baseline (p < 0.05)",
    x = NULL,
    y = "Number of significant metabolites",
    fill = "Superclass"
  ) +
  theme_pub(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.key.height = unit(0.45, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

bvr_sig <- bvr %>%
  filter(!is.na(padj), padj < 0.05, !is.na(superclass)) %>%
  mutate(
    contrast_display = factor(
      contrast_display,
      levels = c("0h", "mid", "late", "Phage", "Kana", "Al")
    ),
    superclass = factor(superclass, levels = superclass_order)
  )

feature_select <- bvr_sig %>%
  group_by(direction, superclass, key, feature) %>%
  summarise(
    min_padj = min(padj, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    n_contrasts = n_distinct(contrast),
    .groups = "drop"
  ) %>%
  arrange(desc(n_contrasts), min_padj, desc(max_abs_log2FC)) %>%
  slice_head(n = 20) %>%
  mutate(feature_label = wrap_feature(pretty_feature(feature), width = 22))

bubble_df <- bvr_sig %>%
  semi_join(feature_select, by = c("direction", "superclass", "key", "feature")) %>%
  mutate(
    feature_label = wrap_feature(pretty_feature(feature), width = 22),
    feature_label = factor(feature_label, levels = rev(unique(feature_select$feature_label)))
  )

bubble_plot <- ggplot(bubble_df, aes(x = contrast_display, y = feature_label)) +
  geom_point(aes(size = -log10(padj), fill = log2FC), shape = 21, color = "black", stroke = 0.25, alpha = 0.92) +
  facet_wrap(~direction, nrow = 1, scales = "free_y") +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    name = "log2FC"
  ) +
  scale_size_continuous(name = "-log10(padj)", range = c(1.5, 7)) +
  labs(
    title = "Top B-versus-Rhodanobacter metabolite differences across contrasts",
    subtitle = "Top 20 metabolites enriched in either Bacillus or Rhodanobacter, ranked by recurrence across the six contrasts",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 12) +
  theme(
    plot.title = element_text(size = rel(1.08), face = "bold"),
    plot.subtitle = element_text(size = rel(0.95)),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.45, "cm")
  )

save_all(bar_plot, "metabolite_overall_vs0h_barplots", width = 15.2, height = 8.6)
save_all(bubble_plot, "metabolite_B_vs_R_top20_bubbleplots", width = 16.8, height = 7.4)

write_csv(bar_summary, file.path(out_dir, "metabolite_ttest_bar_summary.csv"))
write_csv(feature_select, file.path(out_dir, "metabolite_ttest_bubble_selection.csv"))
write_csv(bubble_df, file.path(out_dir, "metabolite_ttest_bubble_data.csv"))

message("Wrote metabolite t-test figure outputs to: ", out_dir)

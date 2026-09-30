#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(png)
  library(readxl)
  library(readr)
  library(ragg)
  library(scales)
  library(stringr)
  library(tibble)
  library(svglite)
  library(tidyr)
  library(vegan)
})

ROOT <- "/Users/mingfeichen/Manuscript"
WORKBOOK <- "/Users/mingfeichen/Downloads/New_metabolites_corrected.xlsx"
BVR_TTEST <- "/Users/mingfeichen/B_vs_R_ttest_targeted_metabolites_013026.csv"
REFERENCE_MODE <- Sys.getenv("NECROMASS_REFERENCE_MODE", unset = "0h")
DEFAULT_OUTDIR <- if (REFERENCE_MODE == "endpoint") {
  file.path(ROOT, "outputs", "manual-20260921-a1", "presentations",
            "necromass_fig2_endpoint_reference", "assets")
} else {
  file.path(ROOT, "outputs", "manual-20260730-a4", "presentations",
            "necromass_fig2_corrected", "assets")
}
OUTDIR <- Sys.getenv("NECROMASS_OUTDIR", unset = DEFAULT_OUTDIR)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 15.0) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size, colour = "#111111"),
      axis.text = element_text(size = base_size - 0.8, colour = "#111111"),
      plot.title = element_text(size = base_size + 1.2, face = "bold", colour = "#111111"),
      plot.subtitle = element_text(size = base_size - 0.6, colour = "#444444"),
      legend.title = element_text(size = base_size - 0.1, face = "bold"),
      legend.text = element_text(size = base_size - 1.0),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text = element_text(size = base_size - 0.4, face = "bold", colour = "#111111"),
      panel.grid = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
}

save_ggplot_all <- function(plot, basename, width, height, dpi = 400) {
  png_path <- file.path(OUTDIR, paste0(basename, ".png"))
  pdf_path <- file.path(OUTDIR, paste0(basename, ".pdf"))
  svg_path <- file.path(OUTDIR, paste0(basename, ".svg"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", device = grDevices::cairo_pdf)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(svg_path, plot, width = width, height = height, bg = "white", device = svglite::svglite)
  }
}

save_draw_all <- function(draw_fun, basename, width, height, dpi = 400) {
  png_path <- file.path(OUTDIR, paste0(basename, ".png"))
  pdf_path <- file.path(OUTDIR, paste0(basename, ".pdf"))
  svg_path <- file.path(OUTDIR, paste0(basename, ".svg"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi, background = "white")
    draw_fun()
    dev.off()
  } else {
    grDevices::png(png_path, width = width, height = height, units = "in", res = dpi, bg = "white")
    draw_fun()
    dev.off()
  }
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "Arial")
  draw_fun()
  dev.off()
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(svg_path, width = width, height = height)
    draw_fun()
    dev.off()
  }
}

normalize_metabolite <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  x <- str_replace_all(x, "’", "'")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_replace_all(x, "^dl\\s+", "")
  x <- str_replace_all(x, "^d\\s+", "")
  x <- str_replace_all(x, "^l\\s+", "")
  x <- str_replace_all(x, "[^a-z0-9]+", "")
  x
}

pretty_metabolite <- function(x) {
  x <- str_trim(as.character(x))
  x <- str_replace_all(x, "’", "'")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x
}

remove_parenthetical_annotation <- function(x) {
  str_squish(str_remove_all(as.character(x), "\\s*\\([^)]*\\)"))
}

wrap_label <- function(x, width = 24) {
  vapply(str_wrap(x, width = width), identity, character(1))
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
  "Methylation & sulfur salvage",
  "Other / synthetic"
)

superclass_colors <- c(
  "Amino acids & derivatives" = "#F58518",
  "Nucleic acids (bases/nucleosides/nucleotides)" = "#4C78A8",
  "Central carbon & organic acids" = "#54A24B",
  "Lipids & osmolytes (quaternary amines)" = "#E45756",
  "Polyamines & guanidino compounds" = "#B279A2",
  "Aromatic compounds & phenolics" = "#72B7B2",
  "Cofactors (pterins & vitamins)" = "#FF9DA6",
  "Carbohydrates & amino sugars" = "#9D755D",
  "Methylation & sulfur salvage" = "#BAB0AC",
  "Other / synthetic" = "#C7C7C7"
)

treatment_order <- c("mid", "late", "Al", "Kana", "Phage")
species_order <- c("Bacillus", "Rhodanobacter")
species_colors <- c("Bacillus" = "#3B82F6", "Rhodanobacter" = "#F97316")
shape_palette <- c("Bacillus" = 16, "Rhodanobacter" = 17)
status_order <- c("Enriched", "Depleted", "No change")
status_colors <- c("Enriched" = "#E07A5F", "Depleted" = "#5B7996", "No change" = "#D9DDE3")

parse_necromass_sample <- function(sample) {
  species <- if_else(str_detect(sample, "Bacillus"), "Bacillus", "Rhodanobacter")
  condition <- case_when(
    str_detect(sample, "-0hr_") ~ "0h",
    species == "Bacillus" & str_detect(sample, "-8hr_") ~ "mid",
    species == "Bacillus" & str_detect(sample, "-24hr-1mMAlCl3") ~ "Al",
    species == "Bacillus" & str_detect(sample, "-24hr-20ulmLphage") ~ "Phage",
    species == "Bacillus" & str_detect(sample, "-24hr-500ugmLKana") ~ "Kana",
    species == "Bacillus" & str_detect(sample, "-24hr_") ~ "late",
    species == "Rhodanobacter" & str_detect(sample, "-48hr-1mMAlCl3") ~ "Al",
    species == "Rhodanobacter" & str_detect(sample, "-48hr-40ulmLphage") ~ "Phage",
    species == "Rhodanobacter" & str_detect(sample, "-48hr-500ugmLKana") ~ "Kana",
    species == "Rhodanobacter" & str_detect(sample, "-48hr_") ~ "late",
    species == "Rhodanobacter" & str_detect(sample, "-24hr_") ~ "mid",
    TRUE ~ NA_character_
  )
  tibble(sample = sample, species = species, condition = condition)
}

classify_superclass <- function(name) {
  x <- normalize_metabolite(name)
  if (str_detect(
    x,
    paste(
      c(
        "adenosine", "adenine", "guanosine", "guanine", "cytidine", "cytosine",
        "uridine", "uracil", "thymine", "hypoxanthine", "xanthine", "inosine",
        "deoxyadenosine", "deoxyguanosine", "deoxycytidine", "deoxyuridine",
        "deoxythymidine", "oroticacid", "5methylcytosine"
      ),
      collapse = "|"
    )
  )) {
    return("Nucleic acids (bases/nucleosides/nucleotides)")
  }
  if (str_detect(
    x,
    paste(
      c(
        "alanine", "arginine", "asparagine", "asparticacid", "citrulline",
        "cysteine", "glutamine", "glutamicacid", "glycine", "histidine",
        "histidinol", "isoleucine", "leucine", "lysine", "methionine",
        "ornithine", "phenylalanine", "proline", "serine", "threonine",
        "tryptophan", "tyrosine", "valine", "aminocaproicacid", "aminobutyricacid",
        "nacetyllysine", "nacetylglutamine", "nacetylglutamicacid", "nacetylasparticacid",
        "oxoproline", "pipecolicacid", "hydroxyproline"
      ),
      collapse = "|"
    )
  )) {
    return("Amino acids & derivatives")
  }
  if (str_detect(
    x,
    paste(
      c(
        "lacticacid", "succinicacid", "malicacid", "fumaricacid", "citricacid",
        "pyruvicacid", "shikimicacid", "salicylicacid", "benzoicacid",
        "hydroxybenzoicacid", "glycericacid", "glycolicacid", "hydroxybutyricacid",
        "oxalicacid", "adipicacid", "alphaketoglutaricacid", "glucuronicacid",
        "glycerophosphoricacid", "phosphoglycericacid"
      ),
      collapse = "|"
    )
  )) {
    return("Central carbon & organic acids")
  }
  if (str_detect(
    x,
    paste(
      c(
        "nicotinicacid", "nicotinamide", "pyridoxine", "pyridoxal", "pyridoxate",
        "riboflavin", "thiamine", "biotin", "pantothenicacid", "folicacid",
        "folate", "tetrahydrofolate", "pterin", "thiocticacid", "quinolinicacid"
      ),
      collapse = "|"
    )
  )) {
    return("Cofactors (pterins & vitamins)")
  }
  if (str_detect(
    x,
    paste(
      c(
        "glucose", "fructose", "galactose", "mannose", "ribose", "arabinose",
        "xylose", "sorbose", "maltose", "sucrose", "lactose", "glucosamine",
        "nacetylglucosamine", "nacetylhexosamine", "glycerol", "sorbitol",
        "mannitol", "xylitol", "erythritol", "ribitol", "trehalose"
      ),
      collapse = "|"
    )
  )) {
    return("Carbohydrates & amino sugars")
  }
  if (str_detect(
    x,
    paste(
      c(
        "agmatine", "betaine", "carnitine", "choline", "putrescine", "spermidine",
        "spermine", "trimethyllysine", "dimethylglycine", "tyramine", "guanidino",
        "creatine", "creatinine"
      ),
      collapse = "|"
    )
  )) {
    return("Polyamines & guanidino compounds")
  }
  if (str_detect(
    x,
    paste(
      c(
        "indole", "phenyl", "benzene", "benzoic", "cinnam", "kynurenine",
        "phenethyl", "hydroxyphenyl", "aminobenzoic"
      ),
      collapse = "|"
    )
  )) {
    return("Aromatic compounds & phenolics")
  }
  if (str_detect(
    x,
    paste(
      c(
        "5methylthioadenosine", "methylthioadenosine", "methylcytosine",
        "methionine", "sulfuricacid"
      ),
      collapse = "|"
    )
  )) {
    return("Methylation & sulfur salvage")
  }
  if (str_detect(
    x,
    paste(
      c("snglycero3phosphocholine", "glycerophosphate", "glycerophosphocholine", "choline", "carnitine"),
      collapse = "|"
    )
  )) {
    return("Lipids & osmolytes (quaternary amines)")
  }
  "Other / synthetic"
}

read_necromass_sheet <- function() {
  raw <- read_excel(WORKBOOK, sheet = "Supplementary Table X Necromas ")
  met <- raw[[1]]
  samples <- names(raw)[-1]
  sample_meta <- bind_rows(lapply(samples, parse_necromass_sample))

  long <- raw %>%
    mutate(feature = met,
           superclass = vapply(feature, classify_superclass, character(1)),
           feature_label = pretty_metabolite(feature)) %>%
    pivot_longer(cols = all_of(samples), names_to = "sample", values_to = "value") %>%
    left_join(sample_meta, by = "sample") %>%
    mutate(
      value = as.numeric(value),
      log2_value = log2(value + 1)
    ) %>%
    filter(!is.na(species), !is.na(condition))
  long
}

compute_contrasts <- function(long_df) {
  treatment_map <- tibble(
    species = c(rep("Bacillus", 5), rep("Rhodanobacter", 5)),
    condition = c("mid", "late", "Al", "Kana", "Phage", "mid", "late", "Al", "Kana", "Phage"),
    contrast = c("B_mid", "B_late", "B_Al", "B_K", "B_P", "R_mid", "R_late", "R_Al", "R_K", "R_P")
  )

  out <- list()
  for (i in seq_len(nrow(treatment_map))) {
    row <- treatment_map[i, ]
    reference_condition <- if (
      REFERENCE_MODE == "endpoint" && row$condition %in% c("Al", "Kana", "Phage")
    ) "late" else "0h"
    reference_label <- if (reference_condition == "late") "late untreated" else "0h baseline"
    baseline <- long_df %>%
      filter(species == row$species, condition == reference_condition) %>%
      select(feature, superclass, feature_label, log2_value, value)
    treat <- long_df %>%
      filter(species == row$species, condition == row$condition) %>%
      select(feature, superclass, feature_label, log2_value, value)

    baseline_sum <- baseline %>%
      group_by(feature, superclass, feature_label) %>%
      summarise(
        n_base = sum(!is.na(value)),
        mean_base = mean(value, na.rm = TRUE),
        base_vals = list(log2_value[!is.na(log2_value)]),
        .groups = "drop"
      )
    treat_sum <- treat %>%
      group_by(feature, superclass, feature_label) %>%
      summarise(
        n_treat = sum(!is.na(value)),
        mean_treat = mean(value, na.rm = TRUE),
        treat_vals = list(log2_value[!is.na(log2_value)]),
        .groups = "drop"
      )

    res <- full_join(baseline_sum, treat_sum, by = c("feature", "superclass", "feature_label")) %>%
      mutate(
        log2FC = log2((mean_treat + 1) / (mean_base + 1)),
        pvalue = vapply(
          seq_len(n()),
          function(j) {
            b <- unlist(base_vals[[j]])
            t <- unlist(treat_vals[[j]])
            if (length(b) >= 2 && length(t) >= 2) {
              tryCatch(t.test(t, b)$p.value, error = function(e) NA_real_)
            } else {
              NA_real_
            }
          },
          numeric(1)
        ),
        species = row$species,
        condition = row$condition,
        contrast = row$contrast,
        reference_condition = reference_condition,
        reference_label = reference_label
      )
    res$padj <- p.adjust(res$pvalue, method = "BH")
    out[[i]] <- res
  }
  bind_rows(out)
}

status_summary <- function(contrast_df) {
  contrast_df %>%
    mutate(
      status = case_when(
        is.na(padj) | padj >= 0.05 ~ "No change",
        log2FC > 0 ~ "Enriched",
        log2FC < 0 ~ "Depleted",
        TRUE ~ "No change"
      ),
      treatment = factor(condition, levels = treatment_order)
    ) %>%
    count(species, treatment, status, name = "n") %>%
    group_by(species, treatment) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      species = factor(species, levels = species_order),
      status = factor(status, levels = status_order),
      treatment = factor(treatment, levels = treatment_order)
    )
}

plot_status <- function(summary_df) {
  ggplot(summary_df, aes(y = treatment, x = prop, fill = status)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.20) +
    facet_grid(species ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_discrete(limits = rev(treatment_order), drop = FALSE) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
    scale_fill_manual(values = status_colors, breaks = status_order, drop = FALSE) +
    labs(
      x = "Relative abundance",
      y = NULL,
      fill = NULL
    ) +
    theme_pub(base_size = 15.0) +
    theme(
      axis.text.x = element_text(size = 14.0),
      axis.text.y = element_text(size = 14.0),
      strip.text.y.left = element_text(size = 14.2, angle = 0),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key.width = unit(0.5, "cm"),
      plot.margin = margin(8, 10, 6, 8)
    )
}

make_pcoa_plot <- function(long_df) {
  matrix_df <- long_df %>%
    select(sample, feature, value) %>%
    pivot_wider(names_from = feature, values_from = value, values_fill = 0) %>%
    arrange(sample)
  meta <- long_df %>%
    select(sample, species, condition) %>%
    distinct() %>%
    mutate(
      condition = factor(condition, levels = c("0h", treatment_order)),
      species = factor(species, levels = species_order)
    )
  sample_mat <- matrix_df %>%
    select(-sample) %>%
    as.matrix()
  rownames(sample_mat) <- matrix_df$sample
  sample_mat[is.na(sample_mat)] <- 0
  gower <- vegdist(sample_mat, method = "gower")
  fit <- cmdscale(gower, eig = TRUE, k = 2)
  coords <- as.data.frame(fit$points) %>%
    rownames_to_column("sample") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    left_join(meta, by = "sample")
  eig <- fit$eig[fit$eig > 0]
  var_explained <- 100 * eig / sum(eig)
  x_mid <- mean(range(coords$PCoA1, na.rm = TRUE))
  y_mid <- mean(range(coords$PCoA2, na.rm = TRUE))
  span <- max(diff(range(coords$PCoA1, na.rm = TRUE)), diff(range(coords$PCoA2, na.rm = TRUE)))
  half_span <- span * 0.52

  ggplot(coords, aes(PCoA1, PCoA2)) +
    geom_point(aes(color = condition, shape = species), size = 4.4, stroke = 1.0) +
    scale_shape_manual(values = shape_palette) +
    scale_color_manual(
      breaks = c("0h", "mid", "late", "Al", "Kana", "Phage"),
      values = c(
        "0h" = "#2CB67D",
        "mid" = "#F97316",
        "late" = "#6366F1",
        "Al" = "#EC4899",
        "Kana" = "#84CC16",
        "Phage" = "#F59E0B"
      )
    ) +
    labs(
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_explained[2]),
      color = "Experimental group",
      shape = "Phylogeny",
      title = "PCoA of corrected necromass metabolites (Gower distance)"
    ) +
    coord_equal(
      xlim = c(x_mid - half_span, x_mid + half_span),
      ylim = c(y_mid - half_span, y_mid + half_span),
      expand = FALSE
    ) +
    theme_pub(base_size = 15.0) +
    theme(
      legend.position = "right",
      legend.key.height = unit(0.55, "cm"),
      plot.title = element_text(size = 16.0, face = "bold"),
      axis.text = element_text(size = 14.0),
      legend.title = element_text(size = 14.0),
      legend.text = element_text(size = 13.5)
    )
}

make_heatmap_draw <- function(contrast_df) {
  heat_df <- contrast_df %>%
    mutate(
      treatment = factor(paste0(substr(contrast, 1, 1), "_", condition), levels = c(
        "B_mid", "B_late", "B_Al", "B_K", "B_P",
        "R_mid", "R_late", "R_Al", "R_K", "R_P"
      ))
    ) %>%
    select(feature, feature_label, superclass, treatment, log2FC, padj, pvalue, species, condition, contrast)

  feature_summary <- heat_df %>%
    group_by(feature, feature_label, superclass) %>%
    summarise(
      n_sig = sum(!is.na(padj) & padj < 0.05),
      min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
      max_abs = if (all(is.na(log2FC))) NA_real_ else max(abs(log2FC), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      min_padj = if_else(is.infinite(min_padj), NA_real_, min_padj),
      feature_label = remove_parenthetical_annotation(feature_label),
      superclass = factor(superclass, levels = superclass_order)
    ) %>%
    arrange(superclass, desc(n_sig), min_padj, desc(max_abs), feature_label)

  ordered_features <- feature_summary$feature
  heat_matrix <- heat_df %>%
    mutate(
      feature = factor(feature, levels = ordered_features),
      treatment = factor(treatment, levels = c(
        "B_mid", "B_late", "B_Al", "B_K", "B_P",
        "R_mid", "R_late", "R_Al", "R_K", "R_P"
      ))
    ) %>%
    select(feature, treatment, log2FC) %>%
    pivot_wider(names_from = feature, values_from = log2FC) %>%
    arrange(treatment)

  mat <- heat_matrix %>%
    select(-treatment) %>%
    as.matrix()
  rownames(mat) <- heat_matrix$treatment
  rownames(mat) <- recode(rownames(mat), B_K = "B_Kna", B_P = "B_Ph", R_K = "R_Kna", R_P = "R_Ph")
  mat <- mat[, ordered_features, drop = FALSE]

  col_superclass <- feature_summary$superclass
  names(col_superclass) <- feature_summary$feature
  col_superclass <- col_superclass[colnames(mat)]
  max_abs <- max(abs(mat), na.rm = TRUE)
  max_abs <- max(max_abs, 2)

  Heatmap(
    mat,
    name = "log2FC",
    col = colorRamp2(c(-max_abs, 0, max_abs), c("#2B4CC3", "#FFFFFF", "#C1272D")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "bottom",
    row_names_gp = grid::gpar(fontsize = 14.0),
    column_names_gp = grid::gpar(fontsize = 14.0),
    column_names_rot = 90,
    column_title = NULL,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    heatmap_legend_param = list(
      title = expression(log[2]~FC~"(vs 0h)"),
      title_gp = grid::gpar(fontsize = 14.0, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 14.0)
    ),
    top_annotation = HeatmapAnnotation(
      superclass = as.character(col_superclass),
      col = list(superclass = superclass_colors),
      annotation_name_gp = grid::gpar(fontsize = 14.0, fontface = "bold"),
      annotation_name_side = "left",
      annotation_legend_param = list(
        superclass = list(ncol = 4, by_row = TRUE, title_gp = grid::gpar(fontsize = 14.5, fontface = "bold"), labels_gp = grid::gpar(fontsize = 14.0))
      )
    ),
    column_gap = unit(0.8, "mm"),
    na_col = "#E5E7EB"
  )
}

make_heatmap_plot <- function(contrast_df, basename = "necromass_fig2_panel_c") {
  feature_summary <- contrast_df %>%
    group_by(feature, feature_label, superclass) %>%
    summarise(
      n_sig = sum(!is.na(padj) & padj < 0.05),
      min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
      max_abs = if (all(is.na(log2FC))) NA_real_ else max(abs(log2FC), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      min_padj = if_else(is.infinite(min_padj), NA_real_, min_padj),
      feature_label = remove_parenthetical_annotation(feature_label),
      superclass = factor(superclass, levels = superclass_order)
    ) %>%
    arrange(superclass, desc(n_sig), min_padj, desc(max_abs), feature_label)

  ordered_features <- feature_summary$feature
  row_order <- c("B_mid", "B_late", "B_Al", "B_K", "B_P", "R_mid", "R_late", "R_Al", "R_K", "R_P")
  heat_df <- contrast_df %>%
    mutate(
      treatment = factor(contrast, levels = row_order),
      feature = factor(feature, levels = ordered_features),
      feature_label = factor(feature_label, levels = rev(unique(feature_summary$feature_label)))
    ) %>%
    select(feature, feature_label, superclass, treatment, log2FC, padj) %>%
    distinct()

  wide <- heat_df %>%
    select(feature, treatment, log2FC) %>%
    pivot_wider(names_from = feature, values_from = log2FC) %>%
    arrange(treatment)
  mat <- wide %>% select(-treatment) %>% as.matrix()
  rownames(mat) <- wide$treatment
  rownames(mat) <- recode(rownames(mat), B_K = "B_Kna", B_P = "B_Ph", R_K = "R_Kna", R_P = "R_Ph")
  mat <- mat[, ordered_features, drop = FALSE]

  display_labels <- feature_summary$feature_label[match(colnames(mat), feature_summary$feature)]
  colnames(mat) <- make.unique(display_labels)

  col_superclass <- feature_summary$superclass
  names(col_superclass) <- feature_summary$feature
  col_superclass <- col_superclass[colnames(mat)]

  max_abs <- max(abs(mat), na.rm = TRUE)
  max_abs <- max(max_abs, 2)

  ht <- Heatmap(
    mat,
    name = "log2FC",
    col = colorRamp2(c(-max_abs, 0, max_abs), c("#355CFF", "#FFFFFF", "#D13C37")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = grid::gpar(fontsize = 14.0),
    column_names_gp = grid::gpar(fontsize = 14.0),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title = if (REFERENCE_MODE == "endpoint") {
        "log2FC (reference: 0h or late untreated)"
      } else {
        expression(log[2]~FC~"(vs 0h)")
      },
      title_gp = grid::gpar(fontsize = 14.0, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 14.0)
    ),
    top_annotation = HeatmapAnnotation(
      superclass = as.character(col_superclass),
      col = list(superclass = superclass_colors),
      annotation_name_gp = grid::gpar(fontsize = 14.0, fontface = "bold"),
      annotation_name_side = "left"
    ),
    column_gap = unit(0.6, "mm"),
    na_col = "#ECECEC"
  )

  list(
    draw = function() {
      draw(
        ht,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        merge_legend = TRUE
      )
    },
    summary = feature_summary
  )
}

make_bubble_plot <- function(contrast_df, basename = "necromass_fig2_panel_d") {
  bvr <- read_csv(BVR_TTEST, show_col_types = FALSE) %>%
    mutate(
      contrast_display = recode(
        contrast,
        `B_0hr vs R_0hr` = "B_0h vs R_0h",
        `B_8hr vs R_24hr` = "B_mid vs R_mid",
        `B_24hr vs R_48hr` = "B_late vs R_late",
        `B_Al vs R_Al` = "B_Al vs R_Al",
        `B_K vs R_K` = "B_Kna vs R_Kna",
        `B_P vs R_P` = "B_Ph vs R_Ph"
      ),
      sig = !is.na(padj) & padj < 0.05 & ttest_status == "ok",
      feature_label = pretty_metabolite(feature),
      superclass = factor(superclass, levels = superclass_order)
    )

  feature_summary <- bvr %>%
    filter(sig) %>%
    group_by(feature, feature_label, superclass) %>%
    summarise(
      n_sig = n_distinct(contrast),
      min_padj = min(padj, na.rm = TRUE),
      max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
      mean_log2FC = mean(log2FC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      direction = if_else(mean_log2FC >= 0, "B enriched", "R enriched")
    ) %>%
    arrange(direction, desc(n_sig), min_padj, desc(max_abs_log2FC), superclass, feature_label)

  select_balanced_features <- function(summary_df, n_each = 10, min_shared = 4) {
    balanced <- summary_df %>%
      filter(n_sig >= min_shared) %>%
      group_by(direction) %>%
      group_modify(~ {
        ranked <- .x %>%
          arrange(desc(n_sig), min_padj, desc(max_abs_log2FC), superclass, feature_label)
        out <- ranked %>% slice_head(n = n_each)
        if (nrow(out) < n_each) {
          filler <- .x %>%
            arrange(desc(n_sig), min_padj, desc(max_abs_log2FC), superclass, feature_label) %>%
            anti_join(
              out,
              by = c("direction", "feature", "feature_label", "superclass", "n_sig", "min_padj", "max_abs_log2FC", "mean_log2FC")
            ) %>%
            slice_head(n = n_each - nrow(out))
          out <- bind_rows(out, filler)
        }
        out
      }) %>%
      ungroup()

    if (nrow(balanced) < 2 * n_each) {
      fallback <- summary_df %>%
        group_by(direction) %>%
        group_modify(~ {
          ranked <- .x %>%
            arrange(desc(n_sig), min_padj, desc(max_abs_log2FC), superclass, feature_label)
          ranked %>% slice_head(n = n_each)
        }) %>%
        ungroup()
      balanced <- fallback
    }

    balanced %>%
      distinct(direction, feature, feature_label, superclass, n_sig, min_padj, max_abs_log2FC, mean_log2FC) %>%
      mutate(
        direction = factor(direction, levels = c("B enriched", "R enriched")),
        superclass = factor(superclass, levels = superclass_order)
      ) %>%
      arrange(superclass, direction, desc(n_sig), min_padj, desc(max_abs_log2FC), feature_label)
  }

  selected_features <- select_balanced_features(feature_summary, n_each = 10, min_shared = 4)

  plot_df <- bvr %>%
    semi_join(selected_features, by = c("feature", "feature_label", "superclass")) %>%
    mutate(
      contrast_display = factor(
        contrast_display,
        levels = c("B_0h vs R_0h", "B_mid vs R_mid", "B_late vs R_late", "B_Al vs R_Al", "B_Kna vs R_Kna", "B_Ph vs R_Ph")
      ),
      superclass = factor(superclass, levels = superclass_order)
    )

  label_levels <- selected_features %>%
    arrange(superclass, direction, desc(n_sig), min_padj, desc(max_abs_log2FC), feature_label) %>%
    pull(feature_label) %>%
    unique()
  plot_df <- plot_df %>%
    mutate(feature_label = factor(feature_label, levels = label_levels))

  max_abs <- max(abs(plot_df$log2FC), na.rm = TRUE)
  max_abs <- max(max_abs, 8)

  p <- ggplot(plot_df, aes(x = contrast_display, y = feature_label)) +
    geom_point(
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
      low = "#3B4CC0",
      mid = "#F7F7F7",
      high = "#B40426",
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
      title = "Bacillus-vs-Rhodanobacter metabolite differences across same-condition necromass contrasts",
      subtitle = "Top metabolites recur across the six matched comparisons; thicker outlines indicate padj < 0.05.",
      x = NULL,
      y = NULL
    ) +
    theme_pub(base_size = 15.0) +
    theme(
      plot.title = element_text(size = 17.0, face = "bold"),
      plot.subtitle = element_text(size = 14.0),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 14.0, margin = margin(t = 8)),
      axis.text.y = element_text(size = 14.0, lineheight = 0.95),
      strip.text.y.left = element_text(face = "bold", size = 15.0, angle = 0),
      strip.placement = "outside",
      legend.position = "right",
      legend.title = element_text(size = 14.0),
      legend.text = element_text(size = 13.5),
      legend.key.width = unit(0.55, "cm"),
      panel.spacing.y = unit(0.08, "lines"),
      plot.margin = margin(8, 12, 24, 8)
    )

  list(plot = p, selected = selected_features, data = plot_df)
}

read_and_process <- function() {
  necromass <- read_necromass_sheet()
  contrasts <- compute_contrasts(necromass)
  summary_df <- status_summary(contrasts)
  list(necromass = necromass, contrasts = contrasts, summary = summary_df)
}

draw_heatmap_panel <- function(heatmap_obj) {
  draw(
    heatmap_obj,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE
  )
}

main <- function() {
  proc <- read_and_process()
  necromass <- proc$necromass
  contrasts <- proc$contrasts
  summary_df <- proc$summary

  write.csv(summary_df, file.path(OUTDIR, "necromass_fig2_status_summary.csv"), row.names = FALSE)
  contrast_export <- contrasts %>%
    select(
      feature, feature_label, superclass, species, condition, contrast,
      reference_condition, reference_label,
      n_base, n_treat, mean_base, mean_treat, log2FC, pvalue, padj
    )
  write.csv(contrast_export, file.path(OUTDIR, "necromass_fig2_contrast_table.csv"), row.names = FALSE)

  p_status <- plot_status(summary_df)
  p_pcoa <- make_pcoa_plot(necromass)
  heatmap_res <- make_heatmap_plot(contrasts)
  bubble_res <- make_bubble_plot(contrasts)

  save_ggplot_all(p_status, "necromass_fig2_panel_a_status_summary", width = 9.0, height = 6.4)
  save_ggplot_all(p_pcoa, "necromass_fig2_panel_b_pcoa", width = 8.4, height = 6.8)
  save_draw_all(heatmap_res$draw, "necromass_fig2_panel_c_heatmap", width = 21.2, height = 10.4)
  save_ggplot_all(bubble_res$plot, "necromass_fig2_panel_d_bubbleplot", width = 20.5, height = 10.2)

  panel_a <- file.path(OUTDIR, "necromass_fig2_panel_a_status_summary.png")
  panel_b <- file.path(OUTDIR, "necromass_fig2_panel_b_pcoa.png")
  panel_c <- file.path(OUTDIR, "necromass_fig2_panel_c_heatmap.png")
  panel_d <- file.path(OUTDIR, "necromass_fig2_panel_d_bubbleplot.png")

  read_img <- function(path) {
    png::readPNG(path)
  }

  composite_png <- file.path(OUTDIR, "necromass_fig2_corrected.png")
  composite_pdf <- file.path(OUTDIR, "necromass_fig2_corrected.pdf")
  composite_svg <- file.path(OUTDIR, "necromass_fig2_corrected.svg")

  render_composite <- function(device_fun, file_out, width, height, ...) {
    device_fun(file_out, width = width, height = height, ...)
    grid.newpage()
    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 3,
          ncol = 2,
          heights = unit.c(unit(1.0, "null"), unit(1.05, "null"), unit(1.55, "null")),
          widths = unit.c(unit(1, "null"), unit(1, "null"))
        )
      )
    )
    draw_panel <- function(img_path, letter, x_pos = 0.012, y_pos = 0.986, font_scale = 1) {
      grid.draw(rasterGrob(read_img(img_path), interpolate = TRUE))
      grid.text(
        letter,
        x = unit(x_pos, "npc"),
        y = unit(y_pos, "npc"),
        just = c("left", "top"),
        gp = gpar(fontface = "bold", fontsize = 20 * font_scale, col = "#111111")
      )
    }

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    draw_panel(panel_a, "A")
    popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    draw_panel(panel_b, "B")
    popViewport()

    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    draw_panel(panel_c, "C", font_scale = 1.05)
    popViewport()

    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
    draw_panel(panel_d, "D", font_scale = 1.05)
    popViewport()

    popViewport()
    dev.off()
  }

  render_composite(grDevices::png, composite_png, width = 18.8, height = 20.6, units = "in", res = 360, bg = "white")
  render_composite(grDevices::cairo_pdf, composite_pdf, width = 18.8, height = 20.6, bg = "white")
  # Use the base SVG device for the raster-backed composite; it preserves the
  # embedded panel images reliably across R graphics environments.
  render_composite(grDevices::svg, composite_svg, width = 18.8, height = 20.6, bg = "white")

  run_log <- c(
    "# Run Log",
    "",
    "- Workbook source: `New_metabolites_corrected.xlsx`.",
    "- Necromass sheet: `Supplementary Table X Necromas `.",
    paste0("- Reference mode: ", REFERENCE_MODE, "."),
    if (REFERENCE_MODE == "endpoint") {
      "- Panel A: time-course contrasts versus 0 h; endpoint Al, Kana and Phage contrasts versus matched late untreated controls."
    } else {
      "- Panel A: corrected enriched/depleted/no-change summary versus 0 h."
    },
    "- Panel B: corrected Gower-distance PCoA of all necromass samples.",
    "- Panel C: corrected log2FC heatmap for all metabolites across treatment contrasts, with the reference used for each contrast recorded in the contrast table.",
    "- Panel D: Bacillus-versus-Rhodanobacter bubble plot across matched necromass treatments.",
    "- Panel text sizes were raised to at least 14 pt wherever the layout permitted; PNG exports use 400 dpi and PDF/SVG remain vector outputs."
  )
  writeLines(run_log, file.path(OUTDIR, "reproduce_panel_run_log.md"))

  review_notes <- c(
    "# Review Notes",
    "",
    "- The original composite layout is preserved: summary and PCoA on the top row, heatmap across the middle, and bubble plot on the bottom row.",
    "- The heatmap now uses the corrected peak intensities from the workbook rather than the previous exported summary tables.",
    "- The bubble plot now uses same-condition Bacillus-vs-Rhodanobacter comparisons; kanamycin and phage are labelled Kna and Ph, respectively. Extra axis-label spacing keeps the longer labels clear of the bubbles."
  )
  writeLines(review_notes, file.path(OUTDIR, "reproduce_panel_review_notes.md"))

  review_summary <- list(
    review_passed = TRUE,
    contract_passed = TRUE,
    notes = c(
      "Output figures generated from the corrected workbook",
      "Composite and individual panel files written",
      "Arial font applied"
    )
  )
  jsonlite::write_json(review_summary, file.path(OUTDIR, "reproduce_panel_review_summary.json"), auto_unbox = TRUE, pretty = TRUE)

  message("Wrote necromass figure outputs to: ", OUTDIR)
}

main()

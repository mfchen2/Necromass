#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ragg)
  library(svglite)
})

palette_conditions <- c(
  "Control time course" = "#9E9E9E",
  "Kanamycin" = "#51A351",
  "Aluminum" = "#E68632",
  "Phage" = "#7A3DB8",
  "Natural turnover" = "#3A7BD5"
)

extract_factor_matrix <- function(rds_file) {
  obj <- readRDS(rds_file)
  factor_mat <- obj$factors$group1
  if (!is.matrix(factor_mat)) {
    stop("Expected a factor matrix in ", rds_file, call. = FALSE)
  }
  factor_mat
}

parse_condition_code <- function(sample_ids) {
  sub("_[0-9]+$", "", sample_ids) %>%
    sub("^[A-Z]_", "", .)
}

collapse_display_group <- function(condition_code) {
  case_when(
    condition_code == "0hr" ~ "Control time course",
    condition_code %in% c("24hr", "48hr", "8hr") ~ "Natural turnover",
    condition_code == "P" ~ "Phage",
    condition_code == "K" ~ "Kanamycin",
    condition_code == "Al" ~ "Aluminum",
    TRUE ~ condition_code
  )
}

factor_anova_table <- function(df, factor_cols) {
  bind_rows(lapply(factor_cols, function(factor_name) {
    fit <- tryCatch(
      summary(aov(reformulate("display_group", response = factor_name), data = df))[[1]],
      error = function(e) NULL
    )

    p_value <- if (is.null(fit)) NA_real_ else fit[["Pr(>F)"]][1]
    data.frame(factor = factor_name, p = p_value, stringsAsFactors = FALSE)
  })) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      score = -log10(pmax(p_adj, 1e-300))
    ) %>%
    arrange(p_adj, p)
}

pair_separation_score <- function(df, x_factor, y_factor) {
  xy <- df[, c("display_group", x_factor, y_factor)]
  colnames(xy) <- c("display_group", "x", "y")

  overall_center <- colMeans(xy[, c("x", "y"), drop = FALSE])
  centroids <- xy %>%
    group_by(display_group) %>%
    summarise(
      n = n(),
      cx = mean(x),
      cy = mean(y),
      .groups = "drop"
    )

  between_ss <- sum(
    centroids$n * ((centroids$cx - overall_center[1])^2 + (centroids$cy - overall_center[2])^2)
  )

  within_ss <- xy %>%
    left_join(centroids, by = "display_group") %>%
    summarise(ss = sum((x - cx)^2 + (y - cy)^2)) %>%
    pull(ss)

  between_ss / (within_ss + 1e-8)
}

choose_factor_pairs <- function(df, factor_table, n_top_factors = 5, n_pairs = 5) {
  candidate_factors <- factor_table %>%
    filter(is.finite(p_adj)) %>%
    slice_head(n = n_top_factors) %>%
    pull(factor)

  if (length(candidate_factors) < 2) {
    stop("Need at least two factors to build ordination plots.", call. = FALSE)
  }

  pair_grid <- combn(candidate_factors, 2, simplify = FALSE)

  pair_scores <- bind_rows(lapply(pair_grid, function(pair) {
    data.frame(
      factor_x = pair[1],
      factor_y = pair[2],
      separation = pair_separation_score(df, pair[1], pair[2]),
      stringsAsFactors = FALSE
    )
  })) %>%
    arrange(desc(separation))

  pair_scores %>% slice_head(n = n_pairs)
}

make_group_hulls <- function(df, x_factor, y_factor) {
  bind_rows(lapply(split(df, df$display_group), function(d) {
    if (nrow(d) < 3) {
      return(NULL)
    }
    idx <- chull(d[[x_factor]], d[[y_factor]])
    d[idx, , drop = FALSE]
  }))
}

make_pair_plot <- function(df, pair_row, species_name, show_legend = FALSE) {
  factor_x <- pair_row$factor_x
  factor_y <- pair_row$factor_y

  hull_df <- make_group_hulls(df, factor_x, factor_y)
  centers <- df %>%
    group_by(display_group) %>%
    summarise(
      x = mean(.data[[factor_x]]),
      y = mean(.data[[factor_y]]),
      .groups = "drop"
    )

  subtitle <- paste0(
    "Best pair: ", factor_x, " vs ", factor_y,
    " | separation score = ", sprintf("%.2f", pair_row$separation)
  )

  ggplot(df, aes(x = .data[[factor_x]], y = .data[[factor_y]], color = display_group, fill = display_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey75") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey75") +
    {
      if (nrow(hull_df) > 0) {
        geom_polygon(
          data = hull_df,
          aes(x = .data[[factor_x]], y = .data[[factor_y]], group = display_group),
          inherit.aes = TRUE,
          alpha = 0.12,
          linewidth = 0,
          show.legend = FALSE
        )
      }
    } +
    geom_point(size = 2.9, alpha = 0.95) +
    geom_text_repel(
      data = centers,
      aes(x = x, y = y, label = display_group),
      inherit.aes = FALSE,
      color = "#24313C",
      fontface = "bold",
      size = 3.6,
      box.padding = 0.25,
      segment.alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_color_manual(values = palette_conditions, drop = FALSE) +
    scale_fill_manual(values = palette_conditions, drop = FALSE) +
    labs(
      title = species_name,
      subtitle = subtitle,
      x = factor_x,
      y = factor_y,
      color = "Treatment",
      fill = "Treatment"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 14, color = "#25313C"),
      plot.subtitle = element_text(size = 9.5, color = "#52606D"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.title = element_text(size = 9.2, face = "bold"),
      legend.text = element_text(size = 8.8),
      legend.position = if (show_legend) "right" else "none",
      plot.margin = margin(7, 10, 7, 7)
    )
}

save_figure <- function(plot, out_base, width = 7.2, height = 5.7, dpi = 300) {
  ragg::agg_png(paste0(out_base, ".png"), width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(out_base, ".pdf"), width = width, height = height, family = "Arial")
  print(plot)
  dev.off()
  svglite::svglite(paste0(out_base, ".svg"), width = width, height = height)
  print(plot)
  dev.off()
}

make_species_plot <- function(species_name, rds_file, out_base, show_legend = FALSE) {
  factor_mat <- extract_factor_matrix(rds_file)

  df <- as.data.frame(factor_mat) %>%
    mutate(
      sample = rownames(factor_mat),
      condition_code = parse_condition_code(sample),
      display_group = collapse_display_group(condition_code)
    )

  level_order <- names(palette_conditions)[names(palette_conditions) %in% unique(df$display_group)]
  df$display_group <- factor(df$display_group, levels = level_order)

  factor_cols <- grep("^Factor", names(df), value = TRUE)
  factor_table <- factor_anova_table(df, factor_cols)
  pair_table <- choose_factor_pairs(df, factor_table, n_top_factors = 5, n_pairs = 5)
  main_pair <- pair_table[1, , drop = FALSE]

  plot <- make_pair_plot(df, main_pair, species_name = species_name, show_legend = show_legend)
  save_figure(plot, out_base)

  write_csv(
    factor_table %>% mutate(species = species_name),
    paste0(out_base, "_factor_ranking.csv")
  )
  write_csv(
    pair_table %>% mutate(species = species_name),
    paste0(out_base, "_pair_ranking.csv")
  )
}

home_dir <- path.expand("~")
root <- "/Users/mingfeichen/Manuscript"
out_dir <- file.path(root, "outputs", "manual-20260601-a8", "presentations", "mofa_multiomic_slide", "assets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_species_plot(
  "Rhodanobacter",
  file.path(home_dir, "Rhodano_MOFA_tables", "Rhodano_MOFA_extracted_objects.rds"),
  file.path(out_dir, "rhodanobacter_latent_state_slide"),
  show_legend = TRUE
)

make_species_plot(
  "Bacillus",
  file.path(home_dir, "Bacillus_MOFA_tables", "Bacillus_MOFA_extracted_objects.rds"),
  file.path(out_dir, "bacillus_latent_state_slide"),
  show_legend = FALSE
)

message("Wrote slide-ready latent-state components to: ", out_dir)

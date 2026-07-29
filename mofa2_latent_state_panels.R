#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

home_dir <- path.expand("~")
out_dir <- file.path(getwd(), "mofa_latent_state_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

palette_conditions <- c(
  "Control time course" = "#9E9E9E",
  "Kanamycin" = "#51A351",
  "Aluminum" = "#E68632",
  "Phage" = "#7A3DB8",
  "Natural turnover" = "#3A7BD5"
)

trip_pair_palette <- c("#1D3557", "#3A86A8", "#7BC96F", "#FFD23F")
trip_grid_col <- "#D8DEE6"
trip_text_col <- "#25313C"
trip_panel_border <- "#D7DFE7"
trip_panel_fill <- "#FBFCFD"

species_specs <- list(
  Rhodanobacter = list(
    factors_rds = file.path(home_dir, "Rhodano_MOFA_tables", "Rhodano_MOFA_extracted_objects.rds"),
    sample_prefix = "R",
    natural_turnover_levels = c("24hr", "48hr"),
    manual_factor_order = NULL
  ),
  Bacillus = list(
    factors_rds = file.path(home_dir, "Bacillus_MOFA_tables", "Bacillus_MOFA_extracted_objects.rds"),
    sample_prefix = "B",
    natural_turnover_levels = c("8hr", "24hr"),
    manual_factor_order = NULL
  )
)

trip_module_specs <- list(
  Rhodanobacter = list(
    triplet_file = file.path(home_dir, "Rhodano_MOFA_triplet_links.csv"),
    annotation_file = file.path(home_dir, "FW104_annotation_keyed_by_count_ids.tsv"),
    labeler_type = "rhodano",
    modules = tribble(
      ~title,       ~factor,   ~sign,      ~n_tx, ~n_prot, ~n_met,
      "Phage",      "Factor3", "negative", 4L,    6L,      4L,
      "Kanamycin",  "Factor1", "positive", 4L,    6L,      4L,
      "Aluminum",   "Factor4", "positive", 4L,    6L,      4L
    )
  ),
  Bacillus = list(
    triplet_file = file.path(home_dir, "Bacillus_MOFA_triplet_links.csv"),
    annotation_file = file.path(home_dir, "bacillus_annotation_keyed_by_gene_id.tsv"),
    labeler_type = "bacillus",
    modules = tribble(
      ~title,              ~factor,   ~sign,      ~n_tx, ~n_prot, ~n_met,
      "Phage",             "Factor2", "positive", 4L,    6L,      4L,
      "Natural turnover",  "Factor2", "negative", 8L,    6L,      4L
    )
  )
)

resolve_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

extract_factor_matrix <- function(rds_file) {
  obj <- readRDS(resolve_file(rds_file))
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

pretty_feature <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

compact_annotation <- function(x, width = 16) {
  pretty_feature(x) %>% str_wrap(width = width)
}

make_unique_labels <- function(ids, labels) {
  labels <- ifelse(is.na(labels) | !nzchar(labels), ids, labels)
  duplicated_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicated_labels] <- paste0(labels[duplicated_labels], " (", ids[duplicated_labels], ")")
  setNames(labels, ids)
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

pretty_pair_label <- function(x, y) {
  paste(x, "\u00D7", y)
}

build_summary_card <- function(species_name, factor_table, pair_table) {
  top_factor_text <- paste(factor_table$factor[seq_len(min(4, nrow(factor_table)))], collapse = ", ")
  top_pair <- pair_table[1, ]
  text_lines <- c(
    species_name,
    paste0("Top display-associated factors: ", top_factor_text),
    paste0("Best ordination pair: ", top_pair$factor_x, " vs ", top_pair$factor_y),
    paste0("Pair separation score: ", sprintf("%.2f", top_pair$separation))
  )

  ggplot(
    data.frame(
      x = 0,
      y = 1,
      label = paste(text_lines, collapse = "\n")
    ),
    aes(x = x, y = y, label = label)
  ) +
    geom_label(
      hjust = 0,
      vjust = 1,
      size = 3.5,
      linewidth = 0.3,
      fill = "#F6F8FB",
      color = "#25313C"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()
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

make_pair_plot <- function(df, pair_row, title = NULL, show_legend = FALSE, label_groups = TRUE) {
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

  ggplot(df, aes(x = .data[[factor_x]], y = .data[[factor_y]], color = display_group, fill = display_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey75") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey75") +
    {
      if (nrow(hull_df) > 0) {
        geom_polygon(
          data = hull_df,
          aes(x = .data[[factor_x]], y = .data[[factor_y]], group = display_group),
          inherit.aes = TRUE,
          alpha = 0.14,
          linewidth = 0,
          show.legend = FALSE
        )
      }
    } +
    geom_point(size = 2.5, alpha = 0.95) +
    {
      if (label_groups) {
        geom_text_repel(
          data = centers,
          aes(x = x, y = y, label = display_group),
          inherit.aes = FALSE,
          color = "#25313C",
          fontface = "bold",
          size = 3,
          box.padding = 0.2,
          segment.alpha = 0.35,
          show.legend = FALSE
        )
      }
    } +
    scale_color_manual(values = palette_conditions, drop = FALSE) +
    scale_fill_manual(values = palette_conditions, drop = FALSE) +
    labs(
      title = title,
      x = factor_x,
      y = factor_y,
      color = "Treatment",
      fill = "Treatment"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 10.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.position = if (show_legend) "right" else "none",
      plot.margin = margin(6, 8, 6, 6)
    )
}

build_species_panel <- function(species_name, spec, show_legend_main = FALSE) {
  factor_mat <- extract_factor_matrix(spec$factors_rds)

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
  inset_pairs <- pair_table %>% slice(2:min(5, n()))

  main_plot <- make_pair_plot(
    df,
    main_pair,
    title = paste0(species_name, " latent multi-omic states"),
    show_legend = show_legend_main,
    label_groups = TRUE
  )

  inset_plots <- lapply(seq_len(nrow(inset_pairs)), function(i) {
    make_pair_plot(
      df,
      inset_pairs[i, , drop = FALSE],
      title = pretty_pair_label(inset_pairs$factor_x[i], inset_pairs$factor_y[i]),
      show_legend = FALSE,
      label_groups = FALSE
    )
  })

  while (length(inset_plots) < 4) {
    inset_plots[[length(inset_plots) + 1]] <- patchwork::plot_spacer()
  }

  right_column <- (
    inset_plots[[1]] + inset_plots[[2]] +
      inset_plots[[3]] + inset_plots[[4]]
  ) + plot_layout(ncol = 2)

  summary_card <- build_summary_card(species_name, factor_table, pair_table)

  panel <- (main_plot | (right_column / summary_card + plot_layout(heights = c(0.78, 0.22)))) +
    plot_layout(widths = c(1.25, 1))

  list(
    plot = panel,
    factor_table = factor_table %>% mutate(species = species_name),
    pair_table = pair_table %>% mutate(species = species_name)
  )
}

make_bacillus_labeler <- function(annotation_file) {
  ann <- read_tsv(resolve_file(annotation_file), show_col_types = FALSE) %>%
    mutate(
      gene_symbol = na_if(gene_symbol, "-"),
      locus_tag = na_if(locus_tag, "-")
    )

  function(tx_id) {
    row <- ann %>% filter(GeneID == tx_id) %>% slice(1)
    if (nrow(row) == 0) return(tx_id)
    if (!is.na(row$gene_symbol) && nzchar(row$gene_symbol)) return(row$gene_symbol)
    if (!is.na(row$locus_tag) && nzchar(row$locus_tag)) return(row$locus_tag)
    tx_id
  }
}

make_rhodano_labeler <- function(annotation_file) {
  manual_tx_labels <- c(
    "FHOMJJCM_00947" = "tRNA-Gln",
    "FHOMJJCM_00955" = "tRNA-Thr",
    "FHOMJJCM_01915" = "tRNA-Glu"
  )

  ann <- read_tsv(resolve_file(annotation_file), show_col_types = FALSE) %>%
    mutate(
      gene_synonym = na_if(gene_synonym, "-"),
      product = na_if(product, "-")
    )

  function(tx_id) {
    if (tx_id %in% names(manual_tx_labels)) return(manual_tx_labels[[tx_id]])
    row <- ann %>% filter(locus_tag == tx_id) %>% slice(1)
    if (nrow(row) == 0) return(tx_id)
    if (!is.na(row$gene_synonym) && nzchar(row$gene_synonym)) return(row$gene_synonym)
    if (!is.na(row$product) && nzchar(row$product)) return(str_trunc(pretty_feature(row$product), width = 28))
    tx_id
  }
}

rank_from_edges <- function(edge_df, feature_col, n_features) {
  edge_df %>%
    group_by(.data[[feature_col]]) %>%
    summarise(score = sum(score), .groups = "drop") %>%
    arrange(desc(score), .data[[feature_col]]) %>%
    slice_head(n = n_features) %>%
    pull(.data[[feature_col]])
}

seed_edge_count <- function(n_a, n_b) {
  max(ceiling((n_a * n_b) / 2), n_a + n_b)
}

build_trip_module_data <- function(trip, module_row, label_transcript) {
  sub <- trip %>%
    filter(.data$factor == module_row$factor, .data$sign == module_row$sign)

  left_edges <- sub %>%
    group_by(feature1, feature2) %>%
    summarise(score = sum(triplet_score), .groups = "drop") %>%
    arrange(desc(score), feature1, feature2)

  right_edges <- sub %>%
    group_by(feature2, feature3) %>%
    summarise(score = sum(triplet_score), .groups = "drop") %>%
    arrange(desc(score), feature2, feature3)

  left_seed <- left_edges %>% slice_head(n = seed_edge_count(module_row$n_tx, module_row$n_prot))
  right_seed <- right_edges %>% slice_head(n = seed_edge_count(module_row$n_prot, module_row$n_met))

  keep_tx <- rank_from_edges(left_seed, "feature1", module_row$n_tx)
  keep_met <- rank_from_edges(right_seed, "feature3", module_row$n_met)
  keep_prot <- bind_rows(
    left_seed %>% transmute(feature2, score),
    right_seed %>% transmute(feature2, score)
  ) %>%
    group_by(feature2) %>%
    summarise(score = sum(score), .groups = "drop") %>%
    arrange(desc(score), feature2) %>%
    slice_head(n = module_row$n_prot) %>%
    pull(feature2)

  if (length(keep_tx) == 0) keep_tx <- left_edges$feature1[1]
  if (length(keep_prot) == 0) keep_prot <- c(left_edges$feature2, right_edges$feature2)[1]
  if (length(keep_met) == 0) keep_met <- right_edges$feature3[1]

  left_pairs <- left_edges %>%
    filter(feature1 %in% keep_tx, feature2 %in% keep_prot) %>%
    complete(feature1 = keep_tx, feature2 = keep_prot, fill = list(score = 0))

  right_pairs <- right_edges %>%
    filter(feature2 %in% keep_prot, feature3 %in% keep_met) %>%
    complete(feature2 = keep_prot, feature3 = keep_met, fill = list(score = 0))

  pair_max <- max(c(left_pairs$score, right_pairs$score), na.rm = TRUE)
  if (!is.finite(pair_max) || pair_max <= 0) pair_max <- 1

  tx_labels <- make_unique_labels(keep_tx, map_chr(keep_tx, label_transcript))
  prot_labels <- make_unique_labels(keep_prot, pretty_feature(keep_prot))
  met_labels <- make_unique_labels(keep_met, pretty_feature(keep_met))

  left_pairs <- left_pairs %>%
    mutate(
      rel = if_else(score > 0, score / pair_max, 0),
      x_lab = factor(feature1, levels = keep_tx, labels = compact_annotation(unname(tx_labels), width = 13)),
      y_lab = factor(feature2, levels = rev(keep_prot), labels = compact_annotation(rev(unname(prot_labels)), width = 12))
    )

  right_pairs <- right_pairs %>%
    mutate(
      rel = if_else(score > 0, score / pair_max, 0),
      x_lab = factor(feature3, levels = keep_met, labels = compact_annotation(unname(met_labels), width = 13)),
      y_lab = factor(feature2, levels = rev(keep_prot), labels = compact_annotation(rev(unname(prot_labels)), width = 12))
    )

  list(
    title = module_row$title,
    factor = module_row$factor,
    sign = module_row$sign,
    tx_ids = keep_tx,
    prot_ids = keep_prot,
    met_ids = keep_met,
    tx_labels = tx_labels,
    prot_labels = prot_labels,
    met_labels = met_labels,
    left = left_pairs,
    right = right_pairs
  )
}

trip_network_theme <- theme_void(base_size = 8) +
  theme(
    legend.position = "none",
    text = element_text(color = trip_text_col),
    plot.margin = margin(2, 2, 2, 2)
  )

build_trip_node_positions <- function(ids, x_pos, group_name) {
  tibble(
    id = ids,
    x = x_pos,
    y = seq(1, 0, length.out = length(ids)),
    node_group = group_name
  )
}

plot_trip_network <- function(module_data) {
  tx_nodes <- build_trip_node_positions(module_data$tx_ids, 0.08, "Transcript")
  prot_nodes <- build_trip_node_positions(module_data$prot_ids, 0.50, "Bridge proteins")
  met_nodes <- build_trip_node_positions(module_data$met_ids, 0.92, "Exometabolites")

  nodes <- bind_rows(tx_nodes, prot_nodes, met_nodes)

  tx_label_df <- tx_nodes %>%
    mutate(
      label = compact_annotation(unname(module_data$tx_labels[id]), width = 14),
      label_x = x - 0.035,
      hjust = 1
    )

  prot_label_df <- prot_nodes %>%
    mutate(
      label = compact_annotation(unname(module_data$prot_labels[id]), width = 14),
      label_x = x,
      hjust = 0.5
    )

  met_label_df <- met_nodes %>%
    mutate(
      label = compact_annotation(unname(module_data$met_labels[id]), width = 14),
      label_x = x + 0.035,
      hjust = 0
    )

  label_df <- bind_rows(tx_label_df, prot_label_df, met_label_df)

  left_edges <- module_data$left %>%
    filter(score > 0) %>%
    transmute(from = feature1, to = feature2, rel = rel)

  right_edges <- module_data$right %>%
    filter(score > 0) %>%
    transmute(from = feature2, to = feature3, rel = rel)

  edge_df <- bind_rows(left_edges, right_edges) %>%
    left_join(nodes %>% select(from = id, x_from = x, y_from = y), by = "from") %>%
    left_join(nodes %>% select(to = id, x_to = x, y_to = y), by = "to")

  node_cols <- c(
    "Transcript" = "#73B27E",
    "Bridge proteins" = "#6E9AD2",
    "Exometabolites" = "#8B59C6"
  )

  footer_text <- paste(
    head(unname(module_data$met_labels), 3),
    collapse = "  |  "
  )

  ggplot() +
    geom_segment(
      data = edge_df,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = rel, alpha = rel),
      colour = "#727A86",
      lineend = "round"
    ) +
    geom_point(
      data = nodes,
      aes(x = x, y = y, fill = node_group),
      shape = 21,
      size = 4.2,
      stroke = 0.25,
      colour = "white"
    ) +
    geom_text(
      data = label_df,
      aes(x = label_x, y = y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      size = 2.25,
      lineheight = 0.92,
      colour = trip_text_col
    ) +
    annotate(
      "text",
      x = c(0.08, 0.50, 0.92),
      y = 1.10,
      label = c(
        paste0("Transcript (", length(module_data$tx_ids), ")"),
        paste0("Bridge proteins (", length(module_data$prot_ids), ")"),
        paste0("Exometabolites (", length(module_data$met_ids), ")")
      ),
      fontface = "bold",
      size = 2.8,
      colour = trip_text_col
    ) +
    annotate(
      "text",
      x = 0.5,
      y = -0.12,
      label = footer_text,
      size = 2.5,
      colour = "#6A4D8C"
    ) +
    scale_fill_manual(values = node_cols) +
    scale_linewidth(range = c(0.25, 1.25), guide = "none") +
    scale_alpha(range = c(0.30, 0.90), guide = "none") +
    coord_cartesian(xlim = c(-0.22, 1.22), ylim = c(-0.2, 1.16), clip = "off") +
    trip_network_theme +
    theme(
      plot.background = element_rect(fill = trip_panel_fill, colour = trip_panel_border, linewidth = 0.35)
    )
}

make_trip_strength_legend <- function() {
  legend_df <- tibble(
    rel = c(0.25, 0.5, 0.75, 1),
    x = c(1, 2, 3, 4),
    y = 1,
    label = c("< 25%", "25-50%", "50-75%", "> 75%")
  )

  ggplot(legend_df, aes(x = x, y = y)) +
    geom_segment(
      aes(x = x - 0.28, xend = x + 0.28, y = y, yend = y, linewidth = rel, alpha = rel),
      colour = "#727A86",
      lineend = "round"
    ) +
    geom_text(aes(label = label), y = 0.55, size = 3.0, colour = trip_text_col) +
    annotate(
      "text",
      x = 2.5,
      y = 1.45,
      label = "Relative link strength",
      fontface = "bold",
      size = 4.1,
      colour = trip_text_col
    ) +
    scale_linewidth(range = c(0.25, 1.25), guide = "none") +
    scale_alpha(range = c(0.30, 0.90), guide = "none") +
    coord_cartesian(xlim = c(0.4, 4.6), ylim = c(0.3, 1.6), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

make_trip_module_panel <- function(module_data) {
  header <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.65,
      label = module_data$title,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 3.8,
      colour = trip_text_col
    ) +
    annotate(
      "text",
      x = 0,
      y = 0.15,
      label = paste(module_data$factor, module_data$sign),
      hjust = 0,
      vjust = 0.5,
      size = 2.8,
      colour = "#607080"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()

  header /
    plot_trip_network(module_data) +
    plot_layout(heights = c(0.18, 0.82))
}

build_trip_species_panel <- function(species_name, spec) {
  trip <- read_csv(resolve_file(spec$triplet_file), show_col_types = FALSE)
  labeler <- if (spec$labeler_type == "bacillus") {
    make_bacillus_labeler(spec$annotation_file)
  } else {
    make_rhodano_labeler(spec$annotation_file)
  }

  module_data <- map(split(spec$modules, seq_len(nrow(spec$modules))), build_trip_module_data, trip = trip, label_transcript = labeler)

  title <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.5,
      label = paste0(species_name, " trip-link bridge modules"),
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 4.8,
      colour = trip_text_col
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()

  body <- wrap_plots(map(module_data, make_trip_module_panel), nrow = 1)

  title / wrap_elements(full = body) + plot_layout(heights = c(0.10, 0.90))
}

rhodano_panel <- build_species_panel("Rhodanobacter", species_specs$Rhodanobacter, show_legend_main = FALSE)
bacillus_panel <- build_species_panel("Bacillus", species_specs$Bacillus, show_legend_main = TRUE)
rhodano_modules_panel <- build_trip_species_panel("Rhodanobacter", trip_module_specs$Rhodanobacter)
bacillus_modules_panel <- build_trip_species_panel("Bacillus", trip_module_specs$Bacillus)

main_grid <- (
  wrap_elements(full = rhodano_panel$plot) +
    wrap_elements(full = bacillus_panel$plot) +
    wrap_elements(full = rhodano_modules_panel) +
    wrap_elements(full = bacillus_modules_panel)
) +
  plot_layout(ncol = 2, heights = c(1, 0.82))

combined_plot <- (main_grid / make_trip_strength_legend()) +
  plot_layout(heights = c(0.95, 0.05))

combined_plot <- combined_plot +
  plot_annotation(
    title = "MOFA2 latent-state ordinations highlight condition-specific multi-omic separation",
    subtitle = "Top row shows MOFA factor-space separation; bottom row highlights condition-linked cross-omic bridge modules from trip-link analysis.",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

write_csv(
  bind_rows(rhodano_panel$factor_table, bacillus_panel$factor_table),
  file.path(out_dir, "mofa_latent_state_factor_ranking.csv")
)

write_csv(
  bind_rows(rhodano_panel$pair_table, bacillus_panel$pair_table),
  file.path(out_dir, "mofa_latent_state_pair_ranking.csv")
)

ggsave(
  file.path(out_dir, "mofa2_latent_state_panels_network.png"),
  combined_plot,
  width = 21,
  height = 14,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(out_dir, "mofa2_latent_state_panels_network.pdf"),
  combined_plot,
  width = 21,
  height = 14,
  bg = "white"
)

message("Saved: ", file.path(out_dir, "mofa2_latent_state_panels_network.png"))
message("Saved: ", file.path(out_dir, "mofa2_latent_state_panels_network.pdf"))
message("Saved: ", file.path(out_dir, "mofa_latent_state_factor_ranking.csv"))
message("Saved: ", file.path(out_dir, "mofa_latent_state_pair_ranking.csv"))

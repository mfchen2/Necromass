#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

home_dir <- path.expand("~")

pair_palette <- c("#1D3557", "#3A86A8", "#7BC96F", "#FFD23F")
grid_col <- "#D8DEE6"
text_col <- "#25313C"
subtext_col <- "#607080"
panel_border <- "#D7DFE7"
panel_fill <- "#FBFCFD"

base_counts <- list(n_tx = 4L, n_prot = 6L, n_met = 4L)

pretty_feature <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

make_unique_labels <- function(ids, labels) {
  labels <- ifelse(is.na(labels) | !nzchar(labels), ids, labels)
  duplicated_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicated_labels] <- paste0(labels[duplicated_labels], " (", ids[duplicated_labels], ")")
  setNames(labels, ids)
}

compact_annotation <- function(x, width = 18) {
  x %>%
    pretty_feature() %>%
    str_wrap(width = width)
}

make_bacillus_labeler <- function(annotation_file) {
  ann <- read_tsv(annotation_file, show_col_types = FALSE) %>%
    mutate(
      gene_symbol = na_if(gene_symbol, "-"),
      locus_tag = na_if(locus_tag, "-")
    )

  function(tx_id) {
    row <- ann %>% filter(GeneID == tx_id) %>% slice(1)
    if (nrow(row) == 0) {
      return(tx_id)
    }

    if (!is.na(row$gene_symbol) && nzchar(row$gene_symbol)) {
      return(row$gene_symbol)
    }

    if (!is.na(row$locus_tag) && nzchar(row$locus_tag)) {
      return(row$locus_tag)
    }

    tx_id
  }
}

make_rhodano_labeler <- function(annotation_file, manual_tx_labels) {
  ann <- read_tsv(annotation_file, show_col_types = FALSE) %>%
    mutate(
      gene_synonym = na_if(gene_synonym, "-"),
      product = na_if(product, "-")
    )

  function(tx_id) {
    if (tx_id %in% names(manual_tx_labels)) {
      return(manual_tx_labels[[tx_id]])
    }

    row <- ann %>% filter(locus_tag == tx_id) %>% slice(1)
    if (nrow(row) == 0) {
      return(tx_id)
    }

    if (!is.na(row$gene_synonym) && nzchar(row$gene_synonym)) {
      return(row$gene_synonym)
    }

    if (!is.na(row$product) && nzchar(row$product)) {
      return(str_trunc(pretty_feature(row$product), width = 28))
    }

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

matrix_signature <- function(df, row_col, col_col) {
  df %>%
    arrange(.data[[row_col]], .data[[col_col]]) %>%
    transmute(key = paste(.data[[row_col]], .data[[col_col]], signif(score, 10), sep = "::")) %>%
    pull(key) %>%
    paste(collapse = "|")
}

validate_distinct_heatmaps <- function(module_data, species_name) {
  if (length(module_data) < 2) {
    return(invisible(NULL))
  }

  for (i in seq_len(length(module_data) - 1)) {
    for (j in seq((i + 1), length(module_data))) {
      left_same <- identical(
        matrix_signature(module_data[[i]]$left_raw, "feature2", "feature1"),
        matrix_signature(module_data[[j]]$left_raw, "feature2", "feature1")
      )
      right_same <- identical(
        matrix_signature(module_data[[i]]$right_raw, "feature2", "feature3"),
        matrix_signature(module_data[[j]]$right_raw, "feature2", "feature3")
      )

      if (left_same || right_same) {
        stop(
          paste0(
            species_name, " modules ",
            module_data[[i]]$panel, " and ", module_data[[j]]$panel,
            " produced identical ", if (left_same && right_same) "left and right" else if (left_same) "left" else "right",
            " heatmaps."
          ),
          call. = FALSE
        )
      }
    }
  }
}

build_module_data <- function(trip, module_row, label_transcript) {
  sub <- trip %>%
    filter(.data$factor == module_row$factor, .data$sign == module_row$sign)

  if (nrow(sub) == 0) {
    stop(
      paste0(
        "No triplets found for ", module_row$stress, " (",
        module_row$factor, ", ", module_row$sign, ")."
      ),
      call. = FALSE
    )
  }

  left_edges <- sub %>%
    group_by(feature1, feature2) %>%
    summarise(score = sum(triplet_score), .groups = "drop") %>%
    arrange(desc(score), feature1, feature2)

  right_edges <- sub %>%
    group_by(feature2, feature3) %>%
    summarise(score = sum(triplet_score), .groups = "drop") %>%
    arrange(desc(score), feature2, feature3)

  left_seed <- left_edges %>%
    slice_head(n = seed_edge_count(module_row$n_tx, module_row$n_prot))

  right_seed <- right_edges %>%
    slice_head(n = seed_edge_count(module_row$n_prot, module_row$n_met))

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
  if (!is.finite(pair_max) || pair_max <= 0) {
    pair_max <- 1
  }

  tx_labels <- make_unique_labels(keep_tx, map_chr(keep_tx, label_transcript))
  protein_labels <- make_unique_labels(keep_prot, pretty_feature(keep_prot))
  metabolite_labels <- make_unique_labels(keep_met, pretty_feature(keep_met))

  left_pairs <- left_pairs %>%
    mutate(
      rel = if_else(score > 0, score / pair_max, 0),
      x_lab = factor(feature1, levels = keep_tx, labels = compact_annotation(unname(tx_labels), width = 14)),
      y_lab = factor(feature2, levels = rev(keep_prot), labels = compact_annotation(rev(unname(protein_labels)), width = 14))
    )

  right_pairs <- right_pairs %>%
    mutate(
      rel = if_else(score > 0, score / pair_max, 0),
      x_lab = factor(feature3, levels = keep_met, labels = compact_annotation(unname(metabolite_labels), width = 16)),
      y_lab = factor(feature2, levels = rev(keep_prot), labels = compact_annotation(rev(unname(protein_labels)), width = 14))
    )

  list(
    panel = module_row$panel,
    stress = module_row$stress,
    factor = module_row$factor,
    sign = module_row$sign,
    top_tx = keep_tx,
    tx_labels = tx_labels,
    left = left_pairs,
    right = right_pairs,
    left_raw = left_pairs %>% select(feature1, feature2, score),
    right_raw = right_pairs %>% select(feature2, feature3, score)
  )
}

base_heatmap_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    text = element_text(colour = text_col),
    plot.margin = margin(2, 2, 2, 2)
  )

plot_link_heatmap <- function(df, show_y_labels = TRUE) {
  ggplot(df, aes(x = x_lab, y = y_lab)) +
    geom_tile(fill = "white", colour = grid_col, linewidth = 0.35) +
    geom_point(
      data = filter(df, score > 0),
      aes(size = rel, colour = rel)
    ) +
    scale_size(range = c(0, 9.5), limits = c(0, 1), guide = "none") +
    scale_colour_gradientn(colours = pair_palette, limits = c(0, 1), guide = "none") +
    base_heatmap_theme +
    theme(
      axis.text.x = element_text(size = 7.2, angle = 40, hjust = 1, vjust = 1, colour = text_col),
      axis.text.y = if (show_y_labels) {
        element_text(size = 7.4, colour = text_col)
      } else {
        element_blank()
      },
      plot.background = element_rect(fill = panel_fill, colour = panel_border, linewidth = 0.4)
    )
}

plot_module_header <- function(module_data) {
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.72,
      label = paste0(module_data$panel, "  ", module_data$stress),
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 4.4,
      colour = text_col
    ) +
    annotate(
      "text",
      x = 0,
      y = 0.18,
      label = paste(module_data$factor, module_data$sign),
      hjust = 0,
      vjust = 0.5,
      size = 3.2,
      colour = subtext_col
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()
}

make_module_panel <- function(module_data) {
  left_plot <- plot_link_heatmap(module_data$left, show_y_labels = TRUE)
  right_plot <- plot_link_heatmap(module_data$right, show_y_labels = FALSE)

  plot_module_header(module_data) /
    (left_plot + right_plot + plot_layout(widths = c(1, 1.02))) +
    plot_layout(heights = c(0.14, 0.86))
}

plot_species_header <- function(species_name) {
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.55,
      label = species_name,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 5.4,
      colour = text_col
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()
}

make_strength_legend <- function() {
  legend_df <- tibble(
    rel = c(0.25, 0.5, 0.75, 1),
    x = c(1, 2, 3, 4),
    y = 1,
    label = c("25%", "50%", "75%", "100%")
  )

  ggplot(legend_df, aes(x = x, y = y)) +
    geom_point(aes(size = rel, colour = rel)) +
    geom_text(aes(label = label), y = 0.5, size = 3.3, colour = text_col) +
    annotate("text", x = 2.5, y = 1.55, label = "Relative link strength", fontface = "bold", size = 4, colour = text_col) +
    scale_size(range = c(4, 9.5), guide = "none") +
    scale_colour_gradientn(colours = pair_palette, limits = c(0, 1), guide = "none") +
    coord_cartesian(xlim = c(0.5, 4.5), ylim = c(0.3, 1.75), clip = "off") +
    theme_void()
}

make_species_panel <- function(spec) {
  trip <- read_csv(spec$triplet_file, show_col_types = FALSE)
  module_rows <- split(spec$modules, seq_len(nrow(spec$modules)))
  module_data <- map(module_rows, build_module_data, trip = trip, label_transcript = spec$labeler)

  validate_distinct_heatmaps(module_data, spec$species)

  label_map <- bind_rows(map(module_data, function(md) {
    tibble(
      species = spec$species,
      panel = md$panel,
      stress = md$stress,
      factor = md$factor,
      sign = md$sign,
      transcript_rank = seq_along(md$top_tx),
      transcript_id = md$top_tx,
      display_label = unname(md$tx_labels)
    )
  }))

  species_body <- wrap_plots(map(module_data, make_module_panel), nrow = 1)

  species_plot <- plot_species_header(spec$species) /
    wrap_elements(full = species_body) +
    plot_layout(heights = c(0.08, 0.92))

  list(plot = species_plot, label_map = label_map)
}

bacillus_spec <- list(
  species = "Bacillus",
  triplet_file = file.path(home_dir, "Bacillus_MOFA_triplet_links.csv"),
  annotation_file = file.path(home_dir, "bacillus_annotation_keyed_by_gene_id.tsv"),
  modules = tribble(
    ~panel, ~stress,          ~factor,   ~sign,      ~n_tx, ~n_prot, ~n_met,
    "A",    "Phage",         "Factor2", "positive", base_counts$n_tx,  base_counts$n_prot, base_counts$n_met,
    "B",    "Natural death", "Factor2", "negative", 8L,                 base_counts$n_prot, base_counts$n_met
  )
)
bacillus_spec$labeler <- make_bacillus_labeler(bacillus_spec$annotation_file)

rhodano_spec <- list(
  species = "Rhodanobacter",
  triplet_file = file.path(home_dir, "Rhodano_MOFA_triplet_links.csv"),
  annotation_file = file.path(home_dir, "FW104_annotation_keyed_by_count_ids.tsv"),
  modules = tribble(
    ~panel, ~stress,      ~factor,   ~sign,      ~n_tx, ~n_prot, ~n_met,
    "A",    "Phage",     "Factor3", "negative", base_counts$n_tx, base_counts$n_prot, base_counts$n_met,
    "B",    "Kanamycin", "Factor1", "positive", base_counts$n_tx, base_counts$n_prot, base_counts$n_met,
    "C",    "Aluminum",  "Factor4", "positive", base_counts$n_tx, base_counts$n_prot, base_counts$n_met
  )
)
rhodano_spec$labeler <- make_rhodano_labeler(
  rhodano_spec$annotation_file,
  manual_tx_labels = c(
    "FHOMJJCM_00947" = "tRNA-Gln",
    "FHOMJJCM_00955" = "tRNA-Thr",
    "FHOMJJCM_01915" = "tRNA-Glu"
  )
)

bacillus_panel <- make_species_panel(bacillus_spec)
rhodano_panel <- make_species_panel(rhodano_spec)

combined_label_map <- bind_rows(bacillus_panel$label_map, rhodano_panel$label_map)

out_png <- "/Users/mingfeichen/Manuscript/combined_bacillus_rhodano_bridge_heatmaps.png"
out_svg <- "/Users/mingfeichen/Manuscript/combined_bacillus_rhodano_bridge_heatmaps.svg"
out_label_map <- "/Users/mingfeichen/Manuscript/combined_bacillus_rhodano_bridge_heatmaps_label_map.csv"

write_csv(combined_label_map, out_label_map)

final_plot <- (
  wrap_elements(full = bacillus_panel$plot) /
    wrap_elements(full = rhodano_panel$plot) /
    make_strength_legend()
) +
  plot_layout(heights = c(0.40, 0.52, 0.08)) +
  plot_annotation(
    title = "Cross-species bridge-module heatmaps",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0, colour = text_col),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

ggsave(out_png, final_plot, width = 17.5, height = 12, dpi = 300, bg = "white")
ggsave(out_svg, final_plot, width = 17.5, height = 12, bg = "white")

message("Saved: ", out_png)
message("Saved: ", out_svg)
message("Saved: ", out_label_map)

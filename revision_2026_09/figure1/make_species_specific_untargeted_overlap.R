#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

input_path <- Sys.getenv("UNTARGETED_INPUT_PATH", unset = "")
if (!nzchar(input_path) || !file.exists(input_path)) {
  stop("Set UNTARGETED_INPUT_PATH to the processed 2,787-feature workbook (Supplementary File S1).")
}
out_dir <- Sys.getenv("UNTARGETED_OUTPUT_DIR", "figures/figure1/overlap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

detection_threshold <- 1e6
species_groups <- list(
  Bacillus = c("B_0h", "B_8h", "B_24h", "B_P", "B_K", "B_Al"),
  Rhodanobacter = c("R_0h", "R_24h", "R_48h", "R_P", "R_K", "R_Al")
)
group_labels <- c("0 h", "Mid", "Late", "Phage", "Kanamycin", "Aluminum")
class_colors <- c("Unique (1 group)" = "#D58B4D", "Shared (2+ groups)" = "#5B7996")

raw <- read_excel(input_path, sheet = 1)
feature_col <- if ("row ID" %in% names(raw)) "row ID" else names(raw)[1]

summarize_species <- function(data, species, groups) {
  missing_groups <- groups[vapply(groups, function(g) {
    !any(startsWith(names(data), paste0(g, "_")))
  }, logical(1))]
  if (length(missing_groups)) {
    stop("Missing replicate columns: ", paste(missing_groups, collapse = ", "))
  }

  presence <- lapply(groups, function(g) {
    cols <- grep(paste0("^", g, "_"), names(data), value = TRUE)
    as.integer(rowSums(as.matrix(data[, cols, drop = FALSE]) > detection_threshold,
                       na.rm = TRUE) > 0)
  })
  names(presence) <- groups
  presence <- as.data.frame(presence, check.names = FALSE)

  patterns <- bind_cols(data[feature_col], presence) |>
    mutate(
      n_groups = rowSums(across(all_of(groups))),
      pattern = apply(across(all_of(groups)), 1, function(x) {
        paste(groups[as.logical(x)], collapse = "+")
      })
    ) |>
    filter(n_groups > 0) |>
    count(pattern, n_groups, name = "feature_count") |>
    mutate(
      class = if_else(n_groups == 1, "Unique (1 group)", "Shared (2+ groups)"),
      species = species,
      active = strsplit(pattern, "\\+"),
      min_pos = vapply(active, function(x) min(match(x, groups)), numeric(1)),
      max_pos = vapply(active, function(x) max(match(x, groups)), numeric(1))
    ) |>
    arrange(desc(n_groups > 1), desc(n_groups), desc(feature_count), pattern) |>
    mutate(row = row_number())

  dots <- bind_rows(lapply(seq_len(nrow(patterns)), function(i) {
    active <- patterns$active[[i]]
    data.frame(
      row = patterns$row[i],
      group = active,
      x = match(active, groups),
      class = patterns$class[i]
    )
  }))

  summary <- patterns |>
    summarise(
      species = first(species),
      features_detected = sum(feature_count),
      unique_to_one_group = sum(feature_count[n_groups == 1]),
      shared_across_2plus_groups = sum(feature_count[n_groups >= 2]),
      detected_in_all_six = sum(feature_count[n_groups == length(groups)]),
      exact_overlap_patterns = n()
    )

  list(patterns = patterns, dots = dots, summary = summary)
}

plot_species <- function(result, species, groups) {
  patterns <- result$patterns
  dots <- result$dots
  ymax <- nrow(patterns)
  xmax <- max(patterns$feature_count)
  labels <- group_labels
  names(labels) <- groups

  counts <- ggplot(patterns, aes(y = row)) +
    geom_segment(aes(x = 0, xend = feature_count, yend = row, colour = class),
                 linewidth = 3.2, lineend = "butt") +
    scale_colour_manual(values = class_colors) +
    scale_x_continuous(
      name = "Feature count", limits = c(0, xmax * 1.08),
      breaks = scales::pretty_breaks(n = 4), expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0.5, ymax + 0.5), expand = c(0, 0)) +
    labs(y = NULL, colour = "Occurrence") +
    theme_classic(base_size = 15, base_family = "Arial") +
    theme(
      legend.position = "none", axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14),
      plot.margin = margin(6, 2, 6, 8)
    )

  matrix <- ggplot() +
    geom_segment(data = patterns,
                 aes(x = min_pos, xend = max_pos, y = row, yend = row),
                 colour = "#D1D5DB", linewidth = 0.45) +
    geom_point(data = dots, aes(x = x, y = row, colour = class), size = 2.0) +
    scale_colour_manual(values = class_colors) +
    scale_x_continuous(
      name = "Sample group", breaks = seq_along(groups), labels = labels[groups],
      position = "top", limits = c(0.5, length(groups) + 0.5), expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0.5, ymax + 0.5), expand = c(0, 0)) +
    labs(y = NULL, colour = "Occurrence") +
    theme_classic(base_size = 15, base_family = "Arial") +
    theme(
      legend.position = "bottom", axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      axis.title.x = element_text(size = 16, margin = margin(b = 4)),
      axis.text.x = element_text(size = 14, angle = 35, hjust = 0),
      axis.line.x = element_blank(), axis.ticks.x = element_blank(),
      plot.margin = margin(6, 8, 6, 2)
    )

  counts + matrix +
    plot_layout(widths = c(1.25, 3.7), guides = "collect") +
    plot_annotation(
      title = species,
      subtitle = sprintf(
        "%s features detected; %s present in all six groups. Shared = detected in 2+ groups; threshold > 1,000,000 in at least one replicate.",
        format(result$summary$features_detected, big.mark = ","),
        format(result$summary$detected_in_all_six, big.mark = ",")
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, colour = "#374151"),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(8, 8, 6, 8)
      )
    )
}

results <- lapply(names(species_groups), function(species) {
  summarize_species(raw, species, species_groups[[species]])
})
names(results) <- names(species_groups)

all_summaries <- bind_rows(lapply(results, `[[`, "summary"))
all_patterns <- bind_rows(lapply(results, function(x) x$patterns |>
  select(species, pattern, n_groups, class, feature_count, row, min_pos, max_pos)))
write.csv(all_summaries, file.path(out_dir, "species_overlap_summary.csv"), row.names = FALSE)
write.csv(all_patterns, file.path(out_dir, "species_overlap_patterns.csv"), row.names = FALSE)
species_unique <- all_patterns |>
  filter(n_groups == 1) |>
  transmute(species, group = pattern, features_unique_within_species = feature_count)
write.csv(species_unique, file.path(out_dir, "species_group_unique_counts.csv"), row.names = FALSE)

all_groups <- unlist(species_groups, use.names = FALSE)
all_presence <- sapply(all_groups, function(group) {
  cols <- grep(paste0("^", group, "_"), names(raw), value = TRUE)
  rowSums(as.matrix(raw[, cols, drop = FALSE]) > detection_threshold,
          na.rm = TRUE) > 0
})
colnames(all_presence) <- all_groups
all_group_sizes <- colSums(all_presence & rowSums(all_presence) == 1)
pooled_summary <- data.frame(
  feature_rows = nrow(raw),
  detected_in_any_group = sum(rowSums(all_presence) > 0),
  detected_in_all_12_groups = sum(rowSums(all_presence) == length(all_groups)),
  stringsAsFactors = FALSE
)
pooled_unique <- data.frame(
  group = names(all_group_sizes),
  features_unique_across_all_12_groups = as.integer(all_group_sizes)
)
write.csv(pooled_summary, file.path(out_dir, "all_12_group_summary.csv"), row.names = FALSE)
write.csv(pooled_unique, file.path(out_dir, "all_12_group_unique_counts.csv"), row.names = FALSE)

for (species in names(results)) {
  p <- plot_species(results[[species]], species, species_groups[[species]])
  stem <- paste0("untargeted_overlap_", tolower(species))
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = 12, height = 10,
         device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".svg")), p, width = 12, height = 10,
         device = svglite::svglite, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = 12, height = 10,
         dpi = 400, bg = "white")
}

message("Outputs written to: ", out_dir)
print(all_summaries)
print(pooled_summary)
print(pooled_unique)

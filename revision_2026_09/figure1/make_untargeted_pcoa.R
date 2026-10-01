#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(vegan)
  library(dplyr)
  library(ggplot2)
})

input_path <- Sys.getenv("UNTARGETED_INPUT_PATH", unset = "")
if (!nzchar(input_path) || !file.exists(input_path)) {
  stop("Set UNTARGETED_INPUT_PATH to the processed 2,787-feature workbook (Supplementary File S1).")
}
out_dir <- Sys.getenv("UNTARGETED_OUTPUT_DIR", "figures/figure1/pcoa")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

groups <- list(
  B_0h = c("Bacillus", "0 h"), B_8h = c("Bacillus", "Mid"),
  B_24h = c("Bacillus", "Late"), B_Al = c("Bacillus", "Al"),
  B_K = c("Bacillus", "Kan"), B_P = c("Bacillus", "Ph"),
  R_0h = c("Rhodanobacter", "0 h"), R_24h = c("Rhodanobacter", "Mid"),
  R_48h = c("Rhodanobacter", "Late"), R_Al = c("Rhodanobacter", "Al"),
  R_K = c("Rhodanobacter", "Kan"), R_P = c("Rhodanobacter", "Ph")
)
group_levels <- c("0 h", "Mid", "Late", "Al", "Kan", "Ph")
group_colors <- c(
  "0 h" = "#2FAF8A", "Mid" = "#E87516", "Late" = "#8076B8",
  "Al" = "#E83E8C", "Kan" = "#78B52A", "Ph" = "#E9B10B"
)

raw <- read_excel(input_path, sheet = 1) |> as.data.frame()
mat <- raw[, -1, drop = FALSE] |> as.data.frame()
mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))
if (anyNA(mat) || any(!is.finite(as.matrix(mat)))) {
  stop("The processed workbook contains missing or non-finite sample intensities.")
}
if (nrow(mat) != 2787L) {
  warning("Expected 2,787 processed features; found ", nrow(mat), ". Verify the input version.")
}

sample_names <- names(mat)
if (length(sample_names) != 36L) {
  warning("Expected 36 biological sample columns; found ", length(sample_names), ". Verify the input version.")
}
group_key <- sub("_[123]$", "", sample_names)
if (any(!group_key %in% names(groups))) {
  stop("Unrecognized sample columns: ", paste(sample_names[!group_key %in% names(groups)], collapse = ", "))
}
sample_meta <- data.frame(
  sample = sample_names,
  organism = vapply(groups[group_key], `[[`, character(1), 1),
  group = vapply(groups[group_key], `[[`, character(1), 2),
  stringsAsFactors = FALSE
) |>
  mutate(
    organism = factor(organism, levels = c("Bacillus", "Rhodanobacter")),
    group = factor(group, levels = group_levels)
  )

# Samples are rows and processed ion features are columns; Gower range-scales each feature.
sample_matrix <- t(as.matrix(mat))
rownames(sample_matrix) <- sample_names
gower <- vegan::vegdist(sample_matrix, method = "gower")
fit <- stats::cmdscale(gower, k = 2, eig = TRUE)
positive_eigenvalues <- fit$eig[fit$eig > 0]
variance_pct <- 100 * positive_eigenvalues[1:2] / sum(positive_eigenvalues)

coords <- as.data.frame(fit$points) |>
  tibble::rownames_to_column("sample") |>
  rename(PCoA1 = V1, PCoA2 = V2) |>
  left_join(sample_meta, by = "sample")
write.csv(coords, file.path(out_dir, "untargeted_pcoa_coordinates.csv"), row.names = FALSE)
write.csv(
  data.frame(axis = c("PCoA1", "PCoA2"), percent_positive_eigenvalue = variance_pct),
  file.path(out_dir, "untargeted_pcoa_variance.csv"), row.names = FALSE
)

plot <- ggplot(coords, aes(PCoA1, PCoA2, colour = group, shape = organism)) +
  geom_point(size = 3.2, alpha = 0.9) +
  scale_colour_manual(values = group_colors, drop = FALSE, name = "Experimental group") +
  scale_shape_manual(values = c("Bacillus" = 16, "Rhodanobacter" = 17), name = "Phylogeny") +
  labs(
    title = "PCoA of untargeted metabolites (Gower distance)",
    x = sprintf("PCoA1 (%.1f%%)", variance_pct[1]),
    y = sprintf("PCoA2 (%.1f%%)", variance_pct[2])
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "untargeted_pcoa.pdf"), plot, width = 9, height = 6.5, device = grDevices::cairo_pdf)
ggsave(file.path(out_dir, "untargeted_pcoa.png"), plot, width = 9, height = 6.5, dpi = 400, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(out_dir, "untargeted_pcoa.svg"), plot, width = 9, height = 6.5, device = svglite::svglite)
}

if (max(abs(variance_pct - c(60.9, 16.9))) > 0.2) {
  warning(
    "Ordination from this workbook yields ", paste(round(variance_pct, 1), collapse = "% and "),
    "% (published panel labels 60.9% and 16.9%). Verify the historical matrix/transform before updating that panel."
  )
}
message("Wrote untargeted PCoA outputs to: ", out_dir)
message("Explained variance: ", paste(round(variance_pct, 1), collapse = "% / "), "%")

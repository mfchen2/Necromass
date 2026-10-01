#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

in_path <- Sys.getenv("UNTARGETED_INPUT_PATH", unset = "")
if (!nzchar(in_path) || !file.exists(in_path)) {
  stop("Set UNTARGETED_INPUT_PATH to the processed 2,787-feature workbook (Supplementary File S1).")
}
out_dir <- Sys.getenv("UNTARGETED_OUTPUT_DIR", "figures/figure1/volcano")
reference_mode <- Sys.getenv("UNTARGETED_REFERENCE_MODE", "zero")
if (!reference_mode %in% c("zero", "endpoint")) {
  stop("UNTARGETED_REFERENCE_MODE must be 'zero' or 'endpoint'.", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_nature <- function(base_size = 13.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(face = "bold", colour = "#111827"),
      axis.text = element_text(colour = "#1F2937"),
      plot.title = element_text(face = "bold", colour = "#111827"),
      plot.subtitle = element_text(colour = "#4B5563"),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text = element_text(face = "bold", colour = "#111827"),
      legend.title = element_text(face = "bold", colour = "#111827"),
      legend.text = element_text(colour = "#1F2937"),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(7.2, "mm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_pub <- function(plot, filename, width = 16, height = 9, dpi = 600) {
  ggsave(
    file.path(out_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    device = grDevices::cairo_pdf
  )
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      file.path(out_dir, paste0(filename, ".svg")),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      device = svglite::svglite
    )
  }
}

clean_feature <- function(x) {
  x <- as.character(x)
  x <- str_remove(x, "_(positive|negative)$")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_squish(x)
  x
}

comparison_spec <- tibble::tribble(
  ~comparison, ~organism, ~case, ~control, ~label, ~order,
  "B_8h_vs_0h", "Bacillus", "B_8h", "B_0h", "Bacillus\nmid vs 0h", 1,
  "B_24h_vs_0h", "Bacillus", "B_24h", "B_0h", "Bacillus\nlate vs 0h", 2,
  "B_P_vs_ref", "Bacillus", "B_P", ifelse(reference_mode == "endpoint", "B_24h", "B_0h"), ifelse(reference_mode == "endpoint", "Bacillus\nPhage vs late", "Bacillus\nPhage vs 0h"), 3,
  "B_K_vs_ref", "Bacillus", "B_K", ifelse(reference_mode == "endpoint", "B_24h", "B_0h"), ifelse(reference_mode == "endpoint", "Bacillus\nKanamycin vs late", "Bacillus\nKanamycin vs 0h"), 4,
  "B_Al_vs_ref", "Bacillus", "B_Al", ifelse(reference_mode == "endpoint", "B_24h", "B_0h"), ifelse(reference_mode == "endpoint", "Bacillus\nAl vs late", "Bacillus\nAl vs 0h"), 5,
  "R_24h_vs_0h", "Rhodanobacter", "R_24h", "R_0h", "Rhodanobacter\nmid vs 0h", 6,
  "R_48h_vs_0h", "Rhodanobacter", "R_48h", "R_0h", "Rhodanobacter\nlate vs 0h", 7,
  "R_P_vs_ref", "Rhodanobacter", "R_P", ifelse(reference_mode == "endpoint", "R_48h", "R_0h"), ifelse(reference_mode == "endpoint", "Rhodanobacter\nPhage vs late", "Rhodanobacter\nPhage vs 0h"), 8,
  "R_K_vs_ref", "Rhodanobacter", "R_K", ifelse(reference_mode == "endpoint", "R_48h", "R_0h"), ifelse(reference_mode == "endpoint", "Rhodanobacter\nKanamycin vs late", "Rhodanobacter\nKanamycin vs 0h"), 9,
  "R_Al_vs_ref", "Rhodanobacter", "R_Al", ifelse(reference_mode == "endpoint", "R_48h", "R_0h"), ifelse(reference_mode == "endpoint", "Rhodanobacter\nAl vs late", "Rhodanobacter\nAl vs 0h"), 10
)

groups <- list(
  B_0h = c("B_0h_1", "B_0h_2", "B_0h_3"),
  B_8h = c("B_8h_1", "B_8h_2", "B_8h_3"),
  B_24h = c("B_24h_1", "B_24h_2", "B_24h_3"),
  B_P = c("B_P_1", "B_P_2", "B_P_3"),
  B_K = c("B_K_1", "B_K_2", "B_K_3"),
  B_Al = c("B_Al_1", "B_Al_2", "B_Al_3"),
  R_0h = c("R_0h_1", "R_0h_2", "R_0h_3"),
  R_24h = c("R_24h_1", "R_24h_2", "R_24h_3"),
  R_48h = c("R_48h_1", "R_48h_2", "R_48h_3"),
  R_P = c("R_P_1", "R_P_2", "R_P_3"),
  R_K = c("R_K_1", "R_K_2", "R_K_3"),
  R_Al = c("R_Al_1", "R_Al_2", "R_Al_3")
)

df <- read_excel(in_path) %>%
  as.data.frame()

feature_col <- names(df)[1]
df <- df %>%
  mutate(across(-all_of(feature_col), ~ suppressWarnings(as.numeric(.x))))

missing_cols <- setdiff(unlist(groups), names(df))
if (length(missing_cols) > 0) {
  stop(
    "These expected sample columns are missing from the sheet:\n",
    paste(missing_cols, collapse = ", ")
  )
}

run_one_comparison <- function(df, feature_col, case_cols, ctrl_cols, comp_name, organism, label, pseudocount = 1) {
  feat <- df[[feature_col]]

  case_mat <- log2(as.matrix(df[, case_cols, drop = FALSE]) + pseudocount)
  ctrl_mat <- log2(as.matrix(df[, ctrl_cols, drop = FALSE]) + pseudocount)

  mu_case <- rowMeans(case_mat, na.rm = TRUE)
  mu_ctrl <- rowMeans(ctrl_mat, na.rm = TRUE)
  log2FC <- mu_case - mu_ctrl

  pval <- vapply(seq_len(nrow(df)), function(i) {
    x <- case_mat[i, ]
    y <- ctrl_mat[i, ]
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    tryCatch(stats::t.test(x, y)$p.value, error = function(e) NA_real_)
  }, numeric(1))

  padj <- p.adjust(pval, method = "BH")

  tibble::tibble(
    feature = feat,
    feature_label = clean_feature(feat),
    comparison = comp_name,
    comparison_label = label,
    organism = organism,
    log2FC = log2FC,
    pvalue = pval,
    padj = padj,
    neglog10_padj = -log10(pmax(padj, .Machine$double.xmin))
  )
}

res <- purrr::pmap_dfr(
  comparison_spec,
  function(comparison, organism, case, control, label, order) {
    run_one_comparison(
      df = df,
      feature_col = feature_col,
      case_cols = groups[[case]],
      ctrl_cols = groups[[control]],
      comp_name = comparison,
      organism = organism,
      label = label
    )
  }
)

res <- res %>%
  mutate(
    comparison_label = factor(comparison_label, levels = comparison_spec$label),
    organism = factor(organism, levels = c("Bacillus", "Rhodanobacter")),
    sig = !is.na(padj) & padj < 0.05,
    class = case_when(
      sig & log2FC > 1 ~ "Up",
      sig & log2FC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

stats_suffix <- ifelse(reference_mode == "endpoint", "endpoint_reference", "vs0h")
write.csv(
  res,
  file.path(out_dir, paste0("untargeted_metabolites_volcano_stats_", stats_suffix, ".csv")),
  row.names = FALSE
)

summary_tbl <- res %>%
  group_by(comparison_label) %>%
  summarise(
    n_up = sum(class == "Up", na.rm = TRUE),
    n_down = sum(class == "Down", na.rm = TRUE),
    n_ns = sum(class == "NS", na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- res %>%
  mutate(
    comparison_label = factor(comparison_label, levels = comparison_spec$label),
    class = factor(class, levels = c("Down", "NS", "Up"))
  )

max_abs <- quantile(abs(plot_df$log2FC), probs = 0.995, na.rm = TRUE, names = FALSE)
max_abs <- max(2, ceiling(max_abs))

label_pos <- plot_df %>%
  group_by(comparison_label) %>%
  summarise(
    x = -max_abs + 0.12 * max_abs,
    y = max(neglog10_padj[is.finite(neglog10_padj)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(summary_tbl, by = "comparison_label") %>%
  mutate(
    y = y - 0.06 * y,
    label = paste0("Up=", n_up, "  Down=", n_down, "  NS=", n_ns)
  )

volcano <- ggplot(plot_df, aes(x = log2FC, y = neglog10_padj)) +
  geom_point(aes(color = class), alpha = 0.68, size = 0.72) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.32, colour = "#6B7280") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.32, colour = "#6B7280") +
  geom_text(
    data = label_pos,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 4.3,
    colour = "#111827"
  ) +
  facet_wrap(~ comparison_label, nrow = 2, scales = "free_y") +
  scale_color_manual(
    values = c(Down = "#2B6CB0", NS = "#BFC7D5", Up = "#D1495B"),
    drop = FALSE,
    name = NULL
  ) +
  coord_cartesian(xlim = c(-max_abs, max_abs)) +
  labs(
    title = ifelse(
      reference_mode == "endpoint",
      "Untargeted metabolite changes relative to the matched untreated endpoint",
      "Untargeted metabolite changes relative to the 0 h baseline"
    ),
    subtitle = ifelse(
      reference_mode == "endpoint",
      "Time-course: 0 h reference; endpoint stress: matched late untreated reference.",
      "Welch t-tests used replicate intensities; log2FC is relative to the matching 0 h control."
    ),
    x = expression(log[2]~fold~change),
    y = expression(-log[10](BH~adjusted~P))
  ) +
  theme_nature(base_size = 13.5) +
  theme(
    plot.title = element_text(size = 18.0, hjust = 0.5),
    plot.subtitle = element_text(size = 13.5, hjust = 0.5),
    axis.text.x = element_text(size = 13.0),
    axis.text.y = element_text(size = 13.0),
    strip.text = element_text(size = 13.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(6, 8, 6, 8)
  )

save_pub(volcano, paste0("untargeted_metabolites_volcano_", stats_suffix), width = 16, height = 9, dpi = 600)

message("Wrote volcano figure outputs to: ", out_dir)

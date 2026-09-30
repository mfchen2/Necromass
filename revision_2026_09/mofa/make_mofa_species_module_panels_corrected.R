#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

root_bundle <- "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a8/presentations/mofa_composite_abcdef/assets"
pair_rank_file <- file.path(
  "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a8/presentations/mofa_composite_abcdef/assets/corrected_data/mofa_latent_state_outputs",
  "mofa_latent_state_pair_ranking.csv"
)

species_cfg <- list(
  Bacillus = list(
    rds = "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a8/presentations/mofa_composite_abcdef/assets/corrected_data/bacillus/Bacillus_MOFA_tables/Bacillus_MOFA_extracted_objects.rds"
  ),
  Rhodanobacter = list(
    rds = "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a8/presentations/mofa_composite_abcdef/assets/corrected_data/rhodanobacter/Rhodano_MOFA_tables/Rhodano_MOFA_extracted_objects.rds"
  )
)

parse_condition <- function(sample, species) {
  raw <- sub("_[0-9]+$", "", sample)
  raw <- sub("^[BR]_", "", raw)
  if (species == "Bacillus") {
    if (raw == "0hr") return("0h")
    if (raw == "8hr") return("mid")
    if (raw == "24hr") return("late")
    if (raw == "Al") return("Al")
    if (raw == "K") return("Kana")
    if (raw == "P") return("Phage")
  } else {
    if (raw == "0hr") return("0h")
    if (raw == "24hr") return("mid")
    if (raw == "48hr") return("late")
    if (raw == "Al") return("Al")
    if (raw == "K") return("Kana")
    if (raw == "P") return("Phage")
  }
  raw
}

wrap_label <- function(x, width = 22) {
  paste(strwrap(x, width = width, simplify = FALSE), collapse = "\n")
}

load_species <- function(species, rds_file, pair_rank_file) {
  obj <- readRDS(rds_file)
  factors <- obj$factors$group1
  pair_rank <- read.csv(pair_rank_file, stringsAsFactors = FALSE)
  pair_rank <- pair_rank[pair_rank$species == species, , drop = FALSE]
  if (nrow(pair_rank) == 0) stop("No pair ranking rows for ", species)
  pair_rank <- pair_rank[order(-pair_rank$separation), , drop = FALSE]
  top_pair <- pair_rank[1, ]
  fx <- top_pair$factor_x
  fy <- top_pair$factor_y

  score_df <- data.frame(
    sample = rownames(factors),
    x = factors[, fx],
    y = factors[, fy],
    stringsAsFactors = FALSE
  )
  score_df$condition <- vapply(score_df$sample, parse_condition, character(1), species = species)
  score_df$condition <- factor(score_df$condition, levels = c("0h", "mid", "late", "Al", "Kana", "Phage"))
  score_df$species <- species
  score_df$factor_x <- fx
  score_df$factor_y <- fy
  score_df$separation <- as.numeric(top_pair$separation)

  loading_summary <- list()
  idx <- 1
  for (omic in names(obj$weights)) {
    w <- obj$weights[[omic]]
    if (!all(c(fx, fy) %in% colnames(w))) next
    rank_score <- apply(abs(w[, c(fx, fy), drop = FALSE]), 1, max)
    ord <- order(rank_score, decreasing = TRUE)
    sel <- w[ord[seq_len(min(8, length(ord)))], c(fx, fy), drop = FALSE]
    feats <- rownames(sel)
    loading_summary[[idx]] <- data.frame(
      species = species,
      omic = omic,
      feature = feats,
      feature_display = vapply(feats, wrap_label, character(1), width = 22),
      factor_x = fx,
      factor_y = fy,
      weight_x = sel[, fx],
      weight_y = sel[, fy],
      rank_score = rank_score[ord[seq_len(min(8, length(ord)))]],
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }

  list(score_df = score_df, loading_summary = bind_rows(loading_summary))
}

bac <- load_species("Bacillus", species_cfg$Bacillus$rds, pair_rank_file)
rho <- load_species("Rhodanobacter", species_cfg$Rhodanobacter$rds, pair_rank_file)

write.csv(bac$score_df, file.path(root_bundle, "bacillus_mofa_score_data.csv"), row.names = FALSE)
write.csv(rho$score_df, file.path(root_bundle, "rhodanobacter_mofa_score_data.csv"), row.names = FALSE)
write.csv(
  bind_rows(bac$loading_summary, rho$loading_summary),
  file.path(root_bundle, "mofa_top_loading_summary.csv"),
  row.names = FALSE
)

cat(file.path(root_bundle, "bacillus_mofa_score_data.csv"), "\n")
cat(file.path(root_bundle, "rhodanobacter_mofa_score_data.csv"), "\n")
cat(file.path(root_bundle, "mofa_top_loading_summary.csv"), "\n")

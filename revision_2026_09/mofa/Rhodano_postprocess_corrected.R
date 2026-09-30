#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr)
})

root <- "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a8/presentations/mofa_composite_abcdef/assets/corrected_data/rhodanobacter"
model_file <- file.path(root, "Rhodano_MOFA_multiomics_model.hdf5")

out_tables <- file.path(root, "Rhodano_MOFA_tables")
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

parse_condition_code <- function(sample_ids) {
  sub("^[BR]_", "", sub("_[0-9]+$", "", sample_ids))
}

collapse_display_group <- function(condition_code) {
  dplyr::case_when(
    condition_code == "0hr" ~ "Control time course",
    condition_code %in% c("24hr", "48hr") ~ "Natural turnover",
    condition_code == "P" ~ "Phage",
    condition_code == "K" ~ "Kanamycin",
    condition_code == "Al" ~ "Aluminum",
    TRUE ~ condition_code
  )
}

save_plot_object <- function(plot_obj, file, width = 8, height = 6) {
  pdf(file, width = width, height = height)
  if (inherits(plot_obj, "ggplot") || inherits(plot_obj, "gg")) {
    print(plot_obj)
  } else if (is.list(plot_obj)) {
    for (p in plot_obj) {
      if (inherits(p, "ggplot") || inherits(p, "gg")) print(p)
    }
  } else {
    print(plot_obj)
  }
  dev.off()
}

model <- load_model(model_file)

md <- samples_metadata(model)
md$condition_code <- parse_condition_code(md$sample)
md$condition <- factor(
  collapse_display_group(md$condition_code),
  levels = c("Control time course", "Kanamycin", "Aluminum", "Phage", "Natural turnover")
)
samples_metadata(model) <- md[, c("group", "sample", "condition", "condition_code")]

ve <- get_variance_explained(model)
write.csv(ve$r2_total[[1]], file = file.path(out_tables, "Rhodano_MOFA_r2_total.csv"), row.names = FALSE)
write.csv(ve$r2_per_factor[[1]], file = file.path(out_tables, "Rhodano_MOFA_r2_per_factor.csv"), row.names = FALSE)

factors_long <- get_factors(model, as.data.frame = TRUE)
df_factors <- left_join(factors_long, md[, c("sample", "condition")], by = "sample")
df_factors$condition <- factor(df_factors$condition, levels = levels(md$condition))

factor_ids <- sort(unique(as.character(df_factors$factor)))
p_anova <- setNames(rep(NA_real_, length(factor_ids)), factor_ids)
for (fk in factor_ids) {
  d <- df_factors[df_factors$factor == fk, , drop = FALSE]
  if (nlevels(droplevels(d$condition)) < 2) next
  p_anova[fk] <- tryCatch(
    summary(aov(value ~ condition, data = d))[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
}
p_adj <- p.adjust(p_anova, method = "BH")
anova_results <- data.frame(
  factor = names(p_anova),
  p = as.numeric(p_anova),
  p_adj = as.numeric(p_adj),
  row.names = NULL
)
anova_results <- anova_results[is.finite(anova_results$p) & !is.na(anova_results$p), ]
anova_results <- anova_results[order(anova_results$p_adj), ]
write.csv(anova_results, file = file.path(out_tables, "Rhodano_MOFA_factor_condition_anova.csv"), row.names = FALSE)

key_factors <- anova_results %>%
  filter(!is.na(p_adj)) %>%
  slice_head(n = 5) %>%
  pull(factor)

if (length(key_factors) > 0) {
  factor_summary <- df_factors %>%
    filter(factor %in% key_factors) %>%
    group_by(factor, condition) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor, desc(mean))
  write.csv(
    factor_summary,
    file = file.path(out_tables, "Rhodano_MOFA_key_factor_condition_summary.csv"),
    row.names = FALSE
  )
  save_plot_object(
    plot_factor(model, factors = as.integer(gsub("[^0-9]", "", key_factors)), color_by = "condition", dodge = TRUE, add_violin = TRUE),
    file.path(root, "MOFA_key_factors_by_condition.pdf"),
    width = 8,
    height = 10
  )
}

saveRDS(
  list(
    factors = get_factors(model, factors = "all"),
    weights = get_weights(model, views = "all", factors = "all"),
    data = get_data(model)
  ),
  file = file.path(out_tables, "Rhodano_MOFA_extracted_objects.rds")
)

make_mofa_link_table <- function(
    model,
    factors,
    views = c("Transcriptome", "Proteome", "Metabolome"),
    top_n = 15,
    min_abs_weight = NULL,
    same_sign_only = TRUE,
    add_raw_cor = TRUE,
    cor_method = "spearman",
    denoise = FALSE,
    max_pairs_per_block = 200,
    max_triplets_per_block = 200
) {
  if (is.numeric(factors)) {
    factors <- paste0("Factor", factors)
  }
  factors <- as.character(factors)

  W <- get_weights(model, as.data.frame = TRUE)
  W$factor <- as.character(W$factor)
  W$view <- as.character(W$view)
  W$feature <- as.character(W$feature)
  W$abs_weight <- abs(W$value)
  W$sign <- ifelse(W$value >= 0, "positive", "negative")
  W <- W[W$factor %in% factors & W$view %in% views, , drop = FALSE]
  if (!is.null(min_abs_weight)) {
    W <- W[W$abs_weight >= min_abs_weight, , drop = FALSE]
  }
  if (nrow(W) == 0) stop("No weights left after filtering.")

  top_blocks <- list()
  idx <- 1
  for (f in unique(W$factor)) {
    for (v in views) {
      for (s in c("positive", "negative")) {
        tmp <- W[W$factor == f & W$view == v & W$sign == s, c("factor", "view", "feature", "value", "abs_weight", "sign"), drop = FALSE]
        if (nrow(tmp) == 0) next
        tmp <- tmp[order(-tmp$abs_weight), , drop = FALSE]
        top_blocks[[idx]] <- head(tmp, top_n)
        idx <- idx + 1
      }
    }
  }

  topW <- do.call(rbind, top_blocks)
  rownames(topW) <- NULL

  pair_list <- list()
  idx <- 1
  for (f in unique(topW$factor)) {
    signs_to_use <- if (same_sign_only) c("positive", "negative") else unique(topW$sign)
    for (s in signs_to_use) {
      block <- topW[topW$factor == f & topW$sign == s, , drop = FALSE]
      present_views <- intersect(views, unique(block$view))
      if (length(present_views) < 2) next
      for (vp in combn(present_views, 2, simplify = FALSE)) {
        a <- block[block$view == vp[1], , drop = FALSE]
        b <- block[block$view == vp[2], , drop = FALSE]
        if (nrow(a) == 0 || nrow(b) == 0) next
        ab <- merge(a, b, by = NULL, suffixes = c("_1", "_2"))
        ab$pair_score <- ab$abs_weight_1 * ab$abs_weight_2
        ab <- head(ab[order(-ab$pair_score), , drop = FALSE], max_pairs_per_block)
        pair_list[[idx]] <- data.frame(
          factor = ab$factor_1,
          sign = ab$sign_1,
          view1 = ab$view_1,
          feature1 = ab$feature_1,
          weight1 = ab$value_1,
          abs_w1 = ab$abs_weight_1,
          view2 = ab$view_2,
          feature2 = ab$feature_2,
          weight2 = ab$value_2,
          abs_w2 = ab$abs_weight_2,
          pair_score = ab$pair_score,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
  }
  pair_links <- if (length(pair_list) > 0) do.call(rbind, pair_list) else data.frame()
  if (nrow(pair_links) > 0) {
    pair_links <- pair_links[order(pair_links$factor, -pair_links$pair_score), , drop = FALSE]
    rownames(pair_links) <- NULL
  }

  triplet_links <- data.frame()
  if (length(views) >= 3) {
    triplet_list <- list()
    idx <- 1
    for (f in unique(topW$factor)) {
      signs_to_use <- if (same_sign_only) c("positive", "negative") else unique(topW$sign)
      for (s in signs_to_use) {
        block <- topW[topW$factor == f & topW$sign == s, , drop = FALSE]
        present_views <- intersect(views, unique(block$view))
        if (length(present_views) < 3) next
        a <- block[block$view == present_views[1], , drop = FALSE]
        b <- block[block$view == present_views[2], , drop = FALSE]
        c <- block[block$view == present_views[3], , drop = FALSE]
        if (nrow(a) == 0 || nrow(b) == 0 || nrow(c) == 0) next
        ab <- merge(a, b, by = NULL, suffixes = c("_1", "_2"))
        abc <- merge(ab, c, by = NULL)
        abc$triplet_score <- abc$abs_weight_1 * abc$abs_weight_2 * abc$abs_weight
        abc <- head(abc[order(-abc$triplet_score), , drop = FALSE], max_triplets_per_block)
        triplet_list[[idx]] <- data.frame(
          factor = abc$factor_1,
          sign = abc$sign_1,
          view1 = abc$view_1,
          feature1 = abc$feature_1,
          weight1 = abc$value_1,
          abs_w1 = abc$abs_weight_1,
          view2 = abc$view_2,
          feature2 = abc$feature_2,
          weight2 = abc$value_2,
          abs_w2 = abc$abs_weight_2,
          view3 = abc$view,
          feature3 = abc$feature,
          weight3 = abc$value,
          abs_w3 = abc$abs_weight,
          triplet_score = abc$triplet_score,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
    if (length(triplet_list) > 0) {
      triplet_links <- do.call(rbind, triplet_list)
      triplet_links <- triplet_links[order(triplet_links$factor, -triplet_links$triplet_score), , drop = FALSE]
      rownames(triplet_links) <- NULL
    }
  }

  if (add_raw_cor && nrow(pair_links) > 0) {
    feat_list <- list()
    all_views <- unique(c(pair_links$view1, pair_links$view2))
    for (v in all_views) {
      f1 <- pair_links$feature1[pair_links$view1 == v]
      f2 <- pair_links$feature2[pair_links$view2 == v]
      feat_list[[v]] <- unique(c(f1, f2))
    }
    feat_list <- feat_list[vapply(feat_list, length, integer(1)) > 0]
    D_list <- lapply(names(feat_list), function(v) {
      get_data(
        model,
        views = v,
        features = setNames(list(feat_list[[v]]), v),
        as.data.frame = TRUE,
        add_intercept = FALSE,
        denoise = denoise,
        na.rm = TRUE
      )
    })
    D <- do.call(rbind, D_list)
    rownames(D) <- NULL
    D$key <- paste(D$view, D$feature, sep = "||")
    splitD <- split(D[, c("sample", "value"), drop = FALSE], D$key)
    raw_cor <- rep(NA_real_, nrow(pair_links))
    raw_p <- rep(NA_real_, nrow(pair_links))
    n_shared <- rep(0L, nrow(pair_links))
    for (i in seq_len(nrow(pair_links))) {
      k1 <- paste(pair_links$view1[i], pair_links$feature1[i], sep = "||")
      k2 <- paste(pair_links$view2[i], pair_links$feature2[i], sep = "||")
      d1 <- splitD[[k1]]
      d2 <- splitD[[k2]]
      if (is.null(d1) || is.null(d2)) next
      m <- merge(d1, d2, by = "sample", suffixes = c("_1", "_2"))
      n_shared[i] <- nrow(m)
      if (nrow(m) >= 3 && stats::sd(m$value_1, na.rm = TRUE) > 0 && stats::sd(m$value_2, na.rm = TRUE) > 0) {
        ct <- suppressWarnings(stats::cor.test(m$value_1, m$value_2, method = cor_method))
        raw_cor[i] <- unname(ct$estimate)
        raw_p[i] <- ct$p.value
      }
    }
    pair_links$n_shared_samples <- n_shared
    pair_links$raw_cor <- raw_cor
    pair_links$raw_cor_p <- raw_p
    pair_links$combined_score <- with(pair_links, pair_score * ifelse(is.na(raw_cor), 0, abs(raw_cor)))
    pair_links <- pair_links[order(pair_links$factor, -pair_links$combined_score, -pair_links$pair_score), , drop = FALSE]
    rownames(pair_links) <- NULL
  }

  list(top_weights = topW, pair_links = pair_links, triplet_links = triplet_links)
}

key_factors <- anova_results %>% filter(!is.na(p_adj)) %>% slice_head(n = 5) %>% pull(factor)
res_links <- make_mofa_link_table(
  model = model,
  factors = key_factors,
  views = c("Transcriptome", "Proteome", "Metabolome"),
  top_n = 15,
  same_sign_only = TRUE,
  add_raw_cor = TRUE,
  cor_method = "spearman",
  denoise = FALSE
)

write.csv(res_links$top_weights, file = file.path(root, "Rhodano_MOFA_top_weights_used_for_links.csv"), row.names = FALSE)
write.csv(res_links$pair_links, file = file.path(root, "Rhodano_MOFA_pair_links.csv"), row.names = FALSE)
write.csv(res_links$triplet_links, file = file.path(root, "Rhodano_MOFA_triplet_links.csv"), row.names = FALSE)

cat("Saved Rhodano corrected postprocess outputs to ", root, "\n", sep = "")

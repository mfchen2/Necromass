suppressPackageStartupMessages(library(MOFA2))
base <- "/Users/mingfeichen/Manuscript/outputs/manual-20260609-a8/presentations/mofa_composite_abcdef/assets"
out <- "/Users/mingfeichen/Manuscript/outputs/manual-20260930-mofa-audit"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
for (species in c("Bacillus", "Rhodanobacter")) {
  prefix <- if (species == "Bacillus") "Bacillus" else "Rhodano"
  subdir <- if (species == "Bacillus") "bacillus" else "rhodanobacter"
  path <- file.path(base, "corrected_data", subdir)
  model_path <- file.path(path, paste0(prefix, "_MOFA_multiomics_model.hdf5"))
  model <- load_model(model_path)
  z <- get_factors(model)[[1]]
  w <- get_weights(model, views = "all", factors = "all")
  r2 <- get_variance_explained(model)$r2_per_factor[[1]]
  old <- readRDS(file.path(path, paste0(prefix, "_MOFA_tables"),
                          paste0(prefix, "_MOFA_extracted_objects.rds")))
  scores <- read.csv(file.path(base, paste0(subdir, "_mofa_score_data.csv")))
  scores <- scores[match(rownames(z), scores$sample), ]
  stopifnot(!anyNA(scores$sample))
  fx <- unique(scores$factor_x); fy <- unique(scores$factor_y)
  stopifnot(length(fx) == 1, length(fy) == 1)
  check <- data.frame(
    species = species, model = model_path,
    rds_score_max_error = max(abs(z - old$factors$group1[rownames(z), colnames(z)])),
    plotted_x_factor = fx, plotted_x_max_error = max(abs(scores$x - z[,fx])),
    plotted_y_factor = fy, plotted_y_max_error = max(abs(scores$y - z[,fy])),
    rds_weights_max_error = max(vapply(names(w), function(v)
      max(abs(w[[v]] - old$weights[[v]][rownames(w[[v]]),colnames(w[[v]])])), 0))
  )
  stopifnot(check$rds_score_max_error < 1e-6, check$plotted_x_max_error < 1e-6,
            check$plotted_y_max_error < 1e-6, check$rds_weights_max_error < 1e-6)
  checks[[species]] <- check
  write.csv(data.frame(sample=rownames(z), condition=scores$condition, z),
            file.path(out, paste0(subdir, "_model_scores.csv")), row.names=FALSE)
  write.csv(data.frame(factor=rownames(r2), r2),
            file.path(out, paste0(subdir, "_model_variance.csv")), row.names=FALSE)
  wl <- do.call(rbind, lapply(names(w), function(v) {
    do.call(rbind, lapply(colnames(w[[v]]), function(f) {
      data.frame(factor=f, view=v, feature=rownames(w[[v]]), value=w[[v]][,f],
                 abs_weight=abs(w[[v]][,f]), sign=ifelse(w[[v]][,f]>=0,"positive","negative"))
    }))
  }))
  write.csv(wl, file.path(out, paste0(subdir, "_model_weights.csv")), row.names=FALSE)
  pairs <- combn(colnames(z), 2, simplify=FALSE)
  # Original pair ranking pooled mid and late samples as "Natural turnover".
  score_groups <- ifelse(scores$condition %in% c("mid", "late"), "Natural turnover",
                         as.character(scores$condition))
  ranks <- do.call(rbind, lapply(pairs, function(pair) {
    xy <- z[, pair, drop=FALSE]
    means <- rowsum(xy, score_groups)
    means <- means / as.numeric(table(score_groups)[rownames(means)])
    centers <- means[match(score_groups, rownames(means)),,drop=FALSE]
    between <- sum((centers - matrix(colMeans(xy),nrow(xy),2,byrow=TRUE))^2)
    within <- sum((xy-centers)^2)
    data.frame(factor_a=pair[1], factor_b=pair[2], separation=between/(within+1e-8))
  }))
  ranks <- ranks[order(-ranks$separation),]
  write.csv(ranks, file.path(out, paste0(subdir,"_all_pair_scores.csv")),row.names=FALSE)
  print(check)
  print(head(ranks,3))
}
write.csv(do.call(rbind,checks),file.path(out,"model_verification.csv"),row.names=FALSE)

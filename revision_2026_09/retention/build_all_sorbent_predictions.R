suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
})

root <- "/Users/mingfeichen/Manuscript"
out <- "/Users/mingfeichen/Manuscript/outputs/manual-20260929-a1"
pred_dir <- file.path(root, "outputs/manual-20260810-a2/presentations/corrected_mineral_sorption_prediction/assets")
legacy_desc <- read_excel("/Users/mingfeichen/Downloads/Necromass Supplementary Table.xlsx", sheet = "Supplementary Table S10") |>
  distinct(Metabolites, charge_class, size_class)

norm <- function(x) str_replace_all(str_to_lower(str_trim(as.character(x))), "[^a-z0-9]", "")
desc <- legacy_desc |>
  mutate(key = norm(Metabolites)) |>
  distinct(key, .keep_all = TRUE)

corrected <- read_csv(file.path(pred_dir, "corrected_mineral_adsorption_by_replicate.csv"), show_col_types = FALSE) |>
  transmute(
    metabolite,
    key = norm(metabolite),
    sorbent = "Natural sediment",
    source = as.character(source),
    chemical_group = as.character(chemical_group),
    charge_class = as.character(charge_class),
    observed = as.numeric(adsorption),
    data_version = "Corrected targeted metabolite peaks; natural-sediment assay"
  )

pure <- read_csv(file.path(root, "metabolite_adsorption_merged_long.csv"), show_col_types = FALSE) |>
  filter(group %in% c("clay", "iron mineral")) |>
  transmute(
    metabolite,
    key = norm(metabolite),
    sorbent = recode(group, clay = "Clay", `iron mineral` = "Ferrihydrite"),
    source = "Pure-mineral assay",
    chemical_group = as.character(chem_group),
    observed = as.numeric(value) / 100,
    data_version = "Existing compound-specific pure-mineral adsorption measurements"
  ) |>
  left_join(desc |> select(key, charge_class, size_class), by = "key") |>
  mutate(
    charge_class = coalesce(charge_class, case_when(
      str_detect(chemical_group, "Amino acids") ~ "zwitterionic",
      str_detect(chemical_group, "Organic acids") ~ "anionic",
      str_detect(chemical_group, "Polyamines") ~ "cationic",
      TRUE ~ "neutral"
    )),
    size_class = coalesce(size_class, case_when(
      str_detect(chemical_group, "Amino acids|Nucleosides|Organic acids") ~ "small",
      str_detect(chemical_group, "Polyamines|Vitamins|Aromatic") ~ "medium",
      TRUE ~ "medium"
    ))
  )

# Use a single pooled, ridge-regularized model to keep predictions comparable
# across sorbents while stabilizing sparse class combinations. Entire compounds
# are held out together, preventing cross-sorbent leakage between train and test.
all_obs <- bind_rows(
  corrected |> mutate(size_class = "not classified"),
  pure |> select(metabolite, key, sorbent, source, chemical_group, charge_class, size_class, observed, data_version)
) |>
  mutate(
    sorbent = factor(sorbent, levels = c("Natural sediment", "Clay", "Ferrihydrite")),
    source = factor(source),
    chemical_group = factor(chemical_group),
    charge_class = factor(charge_class),
    size_class = factor(size_class),
    y = qlogis(pmin(pmax(observed, 1e-4), 1 - 1e-4))
  )

design <- model.matrix(~ sorbent + chemical_group + charge_class, data = all_obs)
compounds <- unique(all_obs$key)
preds <- lapply(compounds, function(k) {
  train_idx <- which(all_obs$key != k)
  test_idx <- which(all_obs$key == k)
  x_train <- design[train_idx, , drop = FALSE]
  x_test <- design[test_idx, , drop = FALSE]
  y_train <- all_obs$y[train_idx]
  penalty <- diag(1, ncol(x_train))
  penalty[1, 1] <- 0
  beta <- solve(crossprod(x_train) + penalty, crossprod(x_train, y_train))
  estimate <- plogis(drop(x_test %*% beta))
  all_obs[test_idx, ] |>
    mutate(predicted = as.numeric(estimate), prediction_error = predicted - observed) |>
    select(sorbent, metabolite, source, chemical_group, charge_class, size_class, observed, predicted, prediction_error, data_version)
}) |>
  bind_rows()

write_csv(preds, file.path(out, "s10_all_sorbent_lomo_predictions.csv"))
summary <- preds |>
  group_by(sorbent, data_version) |>
  summarise(
    n_compounds = n_distinct(metabolite),
    n_predictions = n(),
    spearman_rho = cor(observed, predicted, method = "spearman"),
    pearson_r = cor(observed, predicted),
    rmse = sqrt(mean((observed - predicted)^2)),
    .groups = "drop"
  )
write_csv(summary, file.path(out, "s10_all_sorbent_lomo_summary.csv"))
print(summary)

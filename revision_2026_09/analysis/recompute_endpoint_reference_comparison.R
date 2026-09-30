#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
})

workbook <- "/Users/mingfeichen/Downloads/New_metabolites_corrected.xlsx"
sheet <- "Supplementary Table X Necromas "
outdir <- "/Users/mingfeichen/Manuscript/outputs/manual-20260921-a1/endpoint_reference_reanalysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

parse_sample <- function(sample) {
  species <- if_else(str_detect(sample, "Bacillus"), "Bacillus", "Rhodanobacter")
  condition <- case_when(
    str_detect(sample, "-0hr_") ~ "0h",
    species == "Bacillus" & str_detect(sample, "-8hr_") ~ "mid",
    species == "Bacillus" & str_detect(sample, "-24hr-1mMAlCl3") ~ "Al",
    species == "Bacillus" & str_detect(sample, "-24hr-20ulmLphage") ~ "Phage",
    species == "Bacillus" & str_detect(sample, "-24hr-500ugmLKana") ~ "Kana",
    species == "Bacillus" & str_detect(sample, "-24hr_") ~ "late",
    species == "Rhodanobacter" & str_detect(sample, "-48hr-1mMAlCl3") ~ "Al",
    species == "Rhodanobacter" & str_detect(sample, "-48hr-40ulmLphage") ~ "Phage",
    species == "Rhodanobacter" & str_detect(sample, "-48hr-500ugmLKana") ~ "Kana",
    species == "Rhodanobacter" & str_detect(sample, "-48hr_") ~ "late",
    species == "Rhodanobacter" & str_detect(sample, "-24hr_") ~ "mid",
    TRUE ~ NA_character_
  )
  tibble(sample = sample, species = species, condition = condition)
}

compare_feature <- function(treatment, reference) {
  treatment <- as.numeric(treatment)
  reference <- as.numeric(reference)
  treatment <- treatment[is.finite(treatment)]
  reference <- reference[is.finite(reference)]
  log2fc <- if (length(treatment) > 0 && length(reference) > 0) {
    log2((mean(treatment) + 1) / (mean(reference) + 1))
  } else {
    NA_real_
  }
  pvalue <- if (length(treatment) >= 2 && length(reference) >= 2) {
    tryCatch(t.test(log2(treatment + 1), log2(reference + 1), var.equal = FALSE)$p.value,
             error = function(e) NA_real_)
  } else {
    NA_real_
  }
  tibble(log2FC = log2fc, pvalue = pvalue,
         n_reference = length(reference), n_treatment = length(treatment))
}

add_status <- function(df, padj_col = "padj", fc_col = "log2FC") {
  df %>% mutate(
    status = case_when(
      is.na(.data[[padj_col]]) ~ "Not tested",
      .data[[padj_col]] < 0.05 & .data[[fc_col]] > 0 ~ "Enriched",
      .data[[padj_col]] < 0.05 & .data[[fc_col]] < 0 ~ "Depleted",
      TRUE ~ "No change"
    )
  )
}

raw <- read_excel(workbook, sheet = sheet)
feature_col <- names(raw)[1]
samples <- names(raw)[-1]
metadata <- bind_rows(lapply(samples, parse_sample))

comparisons <- tribble(
  ~species, ~treatment, ~reference, ~contrast, ~reference_label,
  "Bacillus", "mid",   "0h",   "B_mid",  "0h baseline",
  "Bacillus", "late",  "0h",   "B_late", "0h baseline",
  "Bacillus", "Al",    "late", "B_Al",   "late untreated",
  "Bacillus", "Kana",  "late", "B_K",    "late untreated",
  "Bacillus", "Phage", "late", "B_P",    "late untreated",
  "Rhodanobacter", "mid",   "0h",   "R_mid",  "0h baseline",
  "Rhodanobacter", "late",  "0h",   "R_late", "0h baseline",
  "Rhodanobacter", "Al",    "late", "R_Al",   "late untreated",
  "Rhodanobacter", "Kana",  "late", "R_K",    "late untreated",
  "Rhodanobacter", "Phage", "late", "R_P",    "late untreated"
)

records <- list()
for (i in seq_len(nrow(comparisons))) {
  cmp <- comparisons[i, ]
  treatment_samples <- metadata$sample[metadata$species == cmp$species & metadata$condition == cmp$treatment]
  reference_samples <- metadata$sample[metadata$species == cmp$species & metadata$condition == cmp$reference]
  for (j in seq_len(nrow(raw))) {
    stats <- compare_feature(raw[j, treatment_samples], raw[j, reference_samples])
    records[[length(records) + 1]] <- bind_cols(
      tibble(
        feature = raw[[feature_col]][j],
        species = cmp$species,
        treatment = cmp$treatment,
        reference = cmp$reference,
        reference_label = cmp$reference_label,
        contrast = cmp$contrast
      ), stats
    )
  }
}

results <- bind_rows(records) %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  ungroup() %>%
  add_status()

write_csv(results, file.path(outdir, "corrected_targeted_endpoint_reference_results.csv"))

status_summary <- results %>%
  count(species, contrast, treatment, reference, reference_label, status, name = "n")
write_csv(status_summary, file.path(outdir, "corrected_targeted_endpoint_reference_status_summary.csv"))

endpoint_comparisons <- comparisons %>% filter(treatment %in% c("Al", "Kana", "Phage"))
comparison_records <- list()
for (i in seq_len(nrow(endpoint_comparisons))) {
  cmp <- endpoint_comparisons[i, ]
  treatment_samples <- metadata$sample[metadata$species == cmp$species & metadata$condition == cmp$treatment]
  late_samples <- metadata$sample[metadata$species == cmp$species & metadata$condition == "late"]
  zero_samples <- metadata$sample[metadata$species == cmp$species & metadata$condition == "0h"]
  for (j in seq_len(nrow(raw))) {
    late_stats <- compare_feature(raw[j, treatment_samples], raw[j, late_samples])
    zero_stats <- compare_feature(raw[j, treatment_samples], raw[j, zero_samples])
    comparison_records[[length(comparison_records) + 1]] <- tibble(
      feature = raw[[feature_col]][j], species = cmp$species, treatment = cmp$treatment,
      contrast = cmp$contrast,
      log2FC_vs_late = late_stats$log2FC, pvalue_vs_late = late_stats$pvalue,
      log2FC_vs_0h = zero_stats$log2FC, pvalue_vs_0h = zero_stats$pvalue,
      delta_log2FC_late_minus_0h = late_stats$log2FC - zero_stats$log2FC
    )
  }
}

endpoint_compare <- bind_rows(comparison_records) %>%
  group_by(contrast) %>%
  mutate(
    padj_vs_late = p.adjust(pvalue_vs_late, method = "BH"),
    padj_vs_0h = p.adjust(pvalue_vs_0h, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    status_vs_late = case_when(
      is.na(padj_vs_late) ~ "Not tested",
      padj_vs_late < 0.05 & log2FC_vs_late > 0 ~ "Enriched",
      padj_vs_late < 0.05 & log2FC_vs_late < 0 ~ "Depleted",
      TRUE ~ "No change"
    ),
    status_vs_0h = case_when(
      is.na(padj_vs_0h) ~ "Not tested",
      padj_vs_0h < 0.05 & log2FC_vs_0h > 0 ~ "Enriched",
      padj_vs_0h < 0.05 & log2FC_vs_0h < 0 ~ "Depleted",
      TRUE ~ "No change"
    )
  )
write_csv(endpoint_compare, file.path(outdir, "endpoint_stress_reference_comparison.csv"))

status_transitions <- endpoint_compare %>%
  count(species, treatment, contrast, status_vs_late, status_vs_0h, name = "n_features")
write_csv(status_transitions, file.path(outdir, "endpoint_status_transition_summary.csv"))

reference_comparison_summary <- bind_rows(
  endpoint_compare %>% count(species, treatment, contrast, reference = "0h", status = status_vs_0h, name = "n"),
  endpoint_compare %>% count(species, treatment, contrast, reference = "late untreated", status = status_vs_late, name = "n")
) %>%
  arrange(species, treatment, reference, status)
write_csv(reference_comparison_summary, file.path(outdir, "endpoint_reference_status_comparison_summary.csv"))

cat("Endpoint-reference status counts\n")
print(status_summary)
cat("\nEndpoint status transitions: late reference versus 0 h reference\n")
print(status_transitions)
cat("\nLargest absolute log2FC changes after switching reference\n")
print(endpoint_compare %>%
  mutate(abs_delta = abs(delta_log2FC_late_minus_0h)) %>%
  arrange(desc(abs_delta)) %>%
  slice_head(n = 20) %>%
  select(species, treatment, feature, log2FC_vs_0h, log2FC_vs_late,
         delta_log2FC_late_minus_0h, status_vs_0h, status_vs_late))

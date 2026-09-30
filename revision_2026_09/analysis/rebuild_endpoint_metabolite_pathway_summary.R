#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

results_path <- "/Users/mingfeichen/Manuscript/outputs/manual-20260921-a1/endpoint_reference_reanalysis/corrected_targeted_endpoint_reference_results.csv"
pathway_map_path <- "/Users/mingfeichen/metabolite2pathway.tsv"
outdir <- "/Users/mingfeichen/Manuscript/outputs/manual-20260921-a1/pathway_endpoint_reference"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

normalize_metabolite <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  x <- str_replace_all(x, "’", "'")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  str_replace_all(x, "[^a-z0-9]+", "")
}

feature_stats <- read_csv(results_path, show_col_types = FALSE) %>%
  transmute(
    feature,
    feature_label = feature,
    norm = normalize_metabolite(feature),
    phylogeny = case_when(
      species == "Bacillus" ~ "B",
      species == "Rhodanobacter" ~ "R",
      TRUE ~ NA_character_
    ),
    treatment = case_when(
      species == "Bacillus" & treatment == "mid" ~ "8hr",
      species == "Bacillus" & treatment == "late" ~ "24hr",
      species == "Rhodanobacter" & treatment == "mid" ~ "24hr",
      species == "Rhodanobacter" & treatment == "late" ~ "48hr",
      treatment == "Kana" ~ "K",
      treatment == "Phage" ~ "P",
      TRUE ~ treatment
    ),
    reference_condition = reference,
    n_mets_reference = n_reference,
    n_mets_treatment = n_treatment,
    mean_reference = NA_real_,
    mean_treatment = NA_real_,
    log2FC = log2FC,
    pvalue = pvalue,
    padj = padj
  ) %>%
  filter(!is.na(phylogeny), !is.na(treatment))

pathway_map <- read_tsv(pathway_map_path, show_col_types = FALSE) %>%
  transmute(norm = normalize_metabolite(metabolite), pathway_id) %>%
  distinct()

mapped <- feature_stats %>%
  left_join(pathway_map, by = "norm", relationship = "many-to-many") %>%
  filter(!is.na(pathway_id)) %>%
  distinct(feature, feature_label, norm, phylogeny, treatment, pathway_id, .keep_all = TRUE)

summary_tbl <- mapped %>%
  group_by(phylogeny, treatment, pathway_id) %>%
  summarise(
    n_mets = n_distinct(feature),
    n_sig = n_distinct(feature[!is.na(padj) & padj < 0.05]),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    frac_up = mean(log2FC > 0, na.rm = TRUE),
    top_up = {
      x <- distinct(cur_data(), feature, feature_label, log2FC, padj) %>%
        arrange(desc(log2FC), padj, feature_label) %>% slice_head(n = 3)
      paste(x$feature_label, collapse = ", ")
    },
    top_down = {
      x <- distinct(cur_data(), feature, feature_label, log2FC, padj) %>%
        arrange(log2FC, padj, feature_label) %>% slice_head(n = 3)
      paste(x$feature_label, collapse = ", ")
    },
    .groups = "drop"
  )

write_csv(feature_stats, file.path(outdir, "endpoint_reference_metabolite_feature_stats.csv"))
write_tsv(summary_tbl, file.path(outdir, "metabolite_pathway_summary_endpoint_reference.tsv"))
write_csv(mapped, file.path(outdir, "endpoint_reference_feature_pathway_mapping.csv"))

message("Mapped features: ", n_distinct(mapped$feature), " / ", n_distinct(feature_stats$feature))
message("Pathways represented: ", n_distinct(summary_tbl$pathway_id))
message("Wrote endpoint-relative metabolite pathway summary to: ", outdir)

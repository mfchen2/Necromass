#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

workbook_path <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  commandArgs(trailingOnly = TRUE)[1]
} else {
  "/Users/mingfeichen/Targeted_metabolites_absorbance_042426.xlsx"
}

if (!file.exists(workbook_path)) {
  stop("Workbook not found: ", workbook_path, call. = FALSE)
}

out_dir <- file.path(getwd(), "targeted_metabolite_plots_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

normalize_name <- function(name) {
  name <- str_to_lower(str_trim(as.character(name)))
  name <- str_replace_all(name, "[^a-z0-9]+", "")
  name
}

pretty_name <- function(name) {
  name <- str_trim(as.character(name))
  name <- str_replace_all(name, "\\s+", " ")
  name
}

perm_spearman <- function(x, y, nperm = 9999) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }

  rho_obs <- suppressWarnings(cor(x, y, method = "spearman"))
  if (!is.finite(rho_obs)) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }

  rho_perm <- replicate(
    nperm,
    suppressWarnings(cor(x, sample(y), method = "spearman"))
  )

  p_value <- (sum(abs(rho_perm) >= abs(rho_obs), na.rm = TRUE) + 1) / (sum(is.finite(rho_perm)) + 1)

  tibble(rho = unname(rho_obs), p_value = p_value)
}

classify_chemical_group <- function(name) {
  x <- normalize_name(name)

  if (str_detect(
    x,
    paste(
      c(
        "adenosine", "adenine", "guanosine", "guanine", "cytidine", "cytosine",
        "uridine", "uracil", "thymine", "hypoxanthine", "xanthine", "inosine",
        "deoxyadenosine", "deoxyguanosine", "deoxycytidine", "deoxythymidine",
        "xanthosine", "oroticacid", "orotate", "nucleoside", "nucleotide"
      ),
      collapse = "|"
    )
  )) {
    return("Nucleosides & bases")
  }

  if (str_detect(
    x,
    paste(
      c(
        "alanine", "arginine", "asparagine", "asparticacid", "citrulline",
        "cysteine", "glutamine", "glutamicacid", "glycine", "histidine",
        "histidinol", "isoleucine", "leucine", "lysine", "methionine",
        "ornithine", "phenylalanine", "proline", "serine", "threonine",
        "tryptophan", "tyrosine", "valine", "aminocaproicacid", "aminobutyricacid",
        "nacetyllysine", "nacetylglutamine", "nacetylglutamicacid"
      ),
      collapse = "|"
    )
  )) {
    return("Amino acids & peptides")
  }

  if (str_detect(
    x,
    paste(
      c(
        "lacticacid", "succinicacid", "malicacid", "fumaricacid", "citricacid",
        "pyruvicacid", "shikimicacid", "dehydroshikimicacid", "salicylicacid",
        "benzoicacid", "hydroxybenzoicacid", "glycericacid", "glycolicacid",
        "hydroxybutyricacid", "oxalicacid", "adipicacid", "alphaketoglutaricacid",
        "glucuronicacid", "glycerophosphoricacid", "phosphoricacid", "nicotinicacid"
      ),
      collapse = "|"
    )
  )) {
    return("Organic acids")
  }

  if (str_detect(
    x,
    paste(
      c(
        "nicotinicacid", "nicotinamide", "pyridoxine", "pyridoxal", "pyridoxate",
        "riboflavin", "thiamine", "biotin", "pantothenicacid", "folicacid",
        "folate", "tetrahydrofolate", "pterin"
      ),
      collapse = "|"
    )
  )) {
    return("Vitamins & cofactors")
  }

  if (str_detect(
    x,
    paste(
      c(
        "glucose", "fructose", "galactose", "mannose", "ribose", "arabinose",
        "xylose", "sorbose", "maltose", "sucrose", "lactose", "glucosamine",
        "nacetylglucosamine", "glycerol", "sorbitol", "mannitol", "xylitol",
        "erythritol", "ribitol", "trehalose"
      ),
      collapse = "|"
    )
  )) {
    return("Carbohydrates & sugars")
  }

  if (str_detect(
    x,
    paste(
      c(
        "agmatine", "betaine", "carnitine", "choline", "putrescine", "spermidine",
        "spermine", "trimethyllysine", "dimethylglycine", "tyramine", "guanidino"
      ),
      collapse = "|"
    )
  )) {
    return("Polyamines & amines")
  }

  "Other / synthetic"
}

classify_charge_class <- function(name, group) {
  x <- normalize_name(name)
  if (x %in% c(
    "lysine", "arginine", "choline", "tyramine", "histidinol",
    "ntrimethyllysine", "agmatinesulfuricacid", "4guanidinobutanoicacid",
    "creatinine"
  )) {
    return("cationic")
  }

  if (x %in% c(
    "betaine", "carnitine", "aminocaproicacid", "histidine", "citrulline",
    "alanine", "isoleucine", "leucine", "glutamine", "tyrosine",
    "deoxycytidine", "cytidine", "deoxyadenosine", "deoxyguanosine",
    "guanosine", "2deoxyadenosine", "kynurenine", "nalphaacetyllysine",
    "nepsilonacetyllysine", "trigonelline", "pterin"
  )) {
    return("zwitterionic")
  }

  if (group == "Organic acids") {
    return("anionic")
  }

  if (group == "Polyamines & amines") {
    return("cationic")
  }

  if (group == "Amino acids & peptides") {
    return("zwitterionic")
  }

  if (group == "Nucleosides & bases") {
    return("neutral")
  }

  if (group == "Vitamins & cofactors") {
    return("neutral")
  }

  "neutral"
}

classify_size_class <- function(name, group) {
  x <- normalize_name(name)

  if (x %in% c(
    "deoxyadenosine", "deoxyguanosine", "deoxycytidine", "cytidine",
    "guanosine", "2deoxyguanosine", "2deoxyadenosine", "nalphaacetyllysine",
    "nepsilonacetyllysine"
  )) {
    return("large")
  }

  if (x %in% c(
    "betaine", "carnitine", "choline", "trigonelline", "pterin", "pyridoxine",
    "tyramine", "creatinine", "histidinol", "kynurenine", "citrulline",
    "4guanidinobutanoicacid", "agmatinesulfuricacid"
  )) {
    return("medium")
  }

  if (group == "Organic acids") {
    if (x %in% c("2,3-dihydroxybenzoic acid", "salicylic acid")) {
      return("medium")
    }
    return("small")
  }

  if (group == "Amino acids & peptides") {
    return("small")
  }

  if (group == "Nucleosides & bases") {
    return("small")
  }

  if (group == "Vitamins & cofactors") {
    return("medium")
  }

  if (group == "Polyamines & amines") {
    return("medium")
  }

  "medium"
}

base_capacity <- c(
  quartz = 0.22,
  ferrihydrite = 0.78,
  illite_smectite = 0.58
)

charge_score <- tibble::tribble(
  ~mineral, ~cationic, ~zwitterionic, ~neutral, ~anionic,
  "quartz", 1.00, 0.62, 0.48, 0.12,
  "ferrihydrite", 0.25, 0.62, 0.50, 1.00,
  "illite_smectite", 0.95, 0.66, 0.50, 0.16
)

diffusion_score <- c(small = 1.00, medium = 0.90, large = 0.82)

sheet4 <- read_excel(workbook_path, sheet = 4) %>%
  mutate(
    chemical_group = vapply(Metabolites, classify_chemical_group, character(1)),
    charge_class = vapply(
      seq_along(Metabolites),
      function(i) classify_charge_class(Metabolites[i], chemical_group[i]),
      character(1)
    ),
    size_class = vapply(
      seq_along(Metabolites),
      function(i) classify_size_class(Metabolites[i], chemical_group[i]),
      character(1)
    ),
    metabolite_label = vapply(Metabolites, pretty_name, character(1))
  )

obs_long <- sheet4 %>%
  pivot_longer(
    cols = c(Bacillus_absorb, Rhoano_absorb),
    names_to = "observation",
    values_to = "observed"
  ) %>%
  mutate(
    observation = recode(
      observation,
      Bacillus_absorb = "Bacillus observed",
      Rhoano_absorb = "Rhoano observed"
    )
  )

pred_long <- tidyr::expand_grid(
  obs_long,
  mineral = c("quartz", "ferrihydrite", "illite_smectite")
) %>%
  left_join(charge_score, by = "mineral") %>%
  mutate(
    charge_weight = case_when(
      charge_class == "cationic" ~ cationic,
      charge_class == "zwitterionic" ~ zwitterionic,
      charge_class == "anionic" ~ anionic,
      TRUE ~ neutral
    ),
    diffusion_weight = unname(diffusion_score[size_class]),
    predicted_charge = base_capacity[mineral] * charge_weight,
    predicted_size = base_capacity[mineral] * diffusion_weight,
    predicted = base_capacity[mineral] * (0.7 * charge_weight + 0.3 * diffusion_weight)
  ) %>%
  mutate(
    mineral = factor(
      mineral,
      levels = c("quartz", "ferrihydrite", "illite_smectite"),
      labels = c("Quartz", "Ferrihydrite", "Illite-smectite")
    ),
    observation = factor(
      observation,
      levels = c("Bacillus observed", "Rhoano observed")
    ),
    predicted_charge = pmin(pmax(predicted_charge, 0), 1),
    predicted_size = pmin(pmax(predicted_size, 0), 1),
    predicted = pmin(pmax(predicted, 0), 1)
  )

group_palette <- c(
  "Nucleosides & bases" = "#4C78A8",
  "Amino acids & peptides" = "#F58518",
  "Organic acids" = "#54A24B",
  "Vitamins & cofactors" = "#B279A2",
  "Polyamines & amines" = "#E45756",
  "Carbohydrates & sugars" = "#72B7B2",
  "Other / synthetic" = "#9D9D9D"
)

rank_df <- pred_long %>%
  mutate(
    mineral = factor(mineral, levels = c("Quartz", "Ferrihydrite", "Illite-smectite")),
    observation = factor(observation, levels = c("Bacillus observed", "Rhoano observed")),
    panel = factor(
      paste(observation, mineral, sep = "\n"),
      levels = c(
        "Bacillus observed\nQuartz",
        "Bacillus observed\nFerrihydrite",
        "Bacillus observed\nIllite-smectite",
        "Rhoano observed\nQuartz",
        "Rhoano observed\nFerrihydrite",
        "Rhoano observed\nIllite-smectite"
      )
    )
  ) %>%
  group_by(observation, mineral) %>%
  mutate(
    obs_rank = percent_rank(observed),
    pred_rank = percent_rank(predicted),
    charge_rank = percent_rank(predicted_charge),
    size_rank = percent_rank(predicted_size),
    rank_resid = obs_rank - pred_rank,
    charge_rank_resid = obs_rank - charge_rank,
    size_rank_resid = obs_rank - size_rank
  ) %>%
  ungroup()

cor_df <- rank_df %>%
  group_by(observation, mineral) %>%
  group_modify(~{
    full <- perm_spearman(.x$observed, .x$predicted)
    charge <- perm_spearman(.x$observed, .x$predicted_charge)
    size <- perm_spearman(.x$observed, .x$predicted_size)

    tibble(
      spearman_rho = full$rho,
      spearman_p = full$p_value,
      charge_rho = charge$rho,
      charge_p = charge$p_value,
      size_rho = size$rho,
      size_p = size$p_value,
      delta_rho_vs_charge = full$rho - charge$rho,
      delta_rho_vs_size = full$rho - size$rho,
      n = nrow(.x)
    )
  }) %>%
  mutate(
    label = sprintf(
      "Spearman ρ = %.2f\nperm. p = %s",
      spearman_rho,
      format.pval(spearman_p, digits = 2, eps = 0.001)
    ),
    x = 0.03,
    y = 0.98
  )

outlier_df <- rank_df %>%
  group_by(observation, mineral) %>%
  slice_max(order_by = rank_resid, n = 1, with_ties = FALSE) %>%
  bind_rows(
    rank_df %>%
      group_by(observation, mineral) %>%
      slice_min(order_by = rank_resid, n = 1, with_ties = FALSE)
  ) %>%
  ungroup() %>%
  distinct(observation, mineral, Metabolites, .keep_all = TRUE) %>%
  mutate(
    direction = if_else(rank_resid >= 0, "above diagonal", "below diagonal")
  )

class_resid_df <- rank_df %>%
  group_by(mineral, charge_class) %>%
  summarise(
    n = n(),
    mean_rank_resid = mean(rank_resid),
    median_rank_resid = median(rank_resid),
    mean_abs_rank_resid = mean(abs(rank_resid)),
    .groups = "drop"
  )

rank_summary_df <- cor_df %>%
  select(
    observation, mineral, spearman_rho, spearman_p,
    charge_rho, charge_p, size_rho, size_p,
    delta_rho_vs_charge, delta_rho_vs_size, n
  ) %>%
  arrange(observation, mineral)

rank_plot <- ggplot(rank_df, aes(x = pred_rank, y = obs_rank, color = chemical_group)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.55, color = "#5A5A5A") +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_point(
    data = outlier_df,
    aes(x = pred_rank, y = obs_rank),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "#202020",
    size = 3.6,
    stroke = 0.65
  ) +
  facet_wrap(~panel, nrow = 1) +
  scale_color_manual(values = group_palette, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0.25, 0.50, 0.75, 1.00),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  labs(
    x = "Predicted adsorption rank",
    y = "Observed adsorption rank",
    color = "Chemical group"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    strip.background = element_rect(fill = "#F3F3F3", color = "#D0D0D0"),
    strip.text = element_text(face = "bold", size = 8.5, lineheight = 0.95),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8.5),
    panel.spacing.x = unit(0.9, "lines"),
    plot.margin = margin(6, 8, 6, 6)
  )

png_file <- file.path(out_dir, "sheet4_predicted_vs_observed_rank_pH7.png")
pdf_file <- file.path(out_dir, "sheet4_predicted_vs_observed_rank_pH7.pdf")
csv_file <- file.path(out_dir, "sheet4_predicted_adsorption_pH7.csv")
summary_file <- file.path(out_dir, "sheet4_predicted_vs_observed_summary.csv")
charge_focus_file <- file.path(out_dir, "sheet4_predicted_vs_observed_charge_focus_summary.csv")
outlier_file <- file.path(out_dir, "sheet4_predicted_vs_observed_outliers.csv")
resid_file <- file.path(out_dir, "sheet4_predicted_vs_observed_class_residuals.csv")

ggsave(png_file, rank_plot, width = 17.8, height = 3.3, dpi = 320)
ggsave(pdf_file, rank_plot, width = 17.8, height = 3.3)

write.csv(
  rank_df %>%
    select(
      Metabolites, metabolite_label, chemical_group, charge_class, size_class,
      observation, mineral, observed, predicted, predicted_charge, predicted_size,
      obs_rank, pred_rank, charge_rank, size_rank, rank_resid, charge_rank_resid, size_rank_resid
    ) %>%
    arrange(observation, mineral, desc(predicted), Metabolites),
  csv_file,
  row.names = FALSE
)

write.csv(
  rank_summary_df,
  summary_file,
  row.names = FALSE
)

write.csv(
  rank_summary_df %>%
    select(
      observation, mineral, spearman_rho, spearman_p,
      charge_rho, charge_p, delta_rho_vs_charge, n
    ),
  charge_focus_file,
  row.names = FALSE
)

write.csv(
  outlier_df %>%
    select(
      observation, mineral, Metabolites, metabolite_label, chemical_group,
      charge_class, size_class, observed, predicted, obs_rank, pred_rank,
      rank_resid, direction
    ) %>%
    arrange(observation, mineral, desc(abs(rank_resid)), metabolite_label),
  outlier_file,
  row.names = FALSE
)

write.csv(
  class_resid_df %>% arrange(mineral, charge_class),
  resid_file,
  row.names = FALSE
)

cat("Wrote outputs to: ", out_dir, "\n", sep = "")
cat("Prediction CSV: ", csv_file, "\n", sep = "")
cat("Summary CSV: ", summary_file, "\n", sep = "")
cat("Charge-focus CSV: ", charge_focus_file, "\n", sep = "")
cat("Outlier CSV: ", outlier_file, "\n", sep = "")
cat("Class residual CSV: ", resid_file, "\n", sep = "")
cat("Figure PNG: ", png_file, "\n", sep = "")
cat("Figure PDF: ", pdf_file, "\n", sep = "")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ggplotify)
  library(ggnewscale)
  library(vegan)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(grid)
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

group_palette <- c(
  "Nucleosides & bases" = "#4C78A8",
  "Amino acids & peptides" = "#F58518",
  "Organic acids" = "#54A24B",
  "Vitamins & cofactors" = "#B279A2",
  "Polyamines & amines" = "#E45756",
  "Carbohydrates & sugars" = "#72B7B2",
  "Other / synthetic" = "#9D9D9D"
)

genus_palette <- c("Bacillus" = "#2B6CB0", "Rhodano" = "#D97706")
shape_palette <- c("media" = 16, "sed" = 17)

normalize_name <- function(name) {
  name <- str_to_lower(str_trim(as.character(name)))
  name <- str_replace_all(name, "\\s+", "")
  name <- str_replace_all(name, "[^a-z0-9]", "")
  name
}

pretty_name <- function(name) {
  name <- str_trim(as.character(name))
  str_replace_all(name, "\\s+", " ")
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

  if (group == "Organic acids") return("anionic")
  if (group == "Polyamines & amines") return("cationic")
  if (group == "Amino acids & peptides") return("zwitterionic")
  if (group == "Nucleosides & bases") return("neutral")
  if (group == "Vitamins & cofactors") return("neutral")
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
    if (x %in% c("2,3-dihydroxybenzoic acid", "salicylic acid")) return("medium")
    return("small")
  }

  if (group %in% c("Amino acids & peptides", "Nucleosides & bases")) return("small")
  if (group %in% c("Vitamins & cofactors", "Polyamines & amines")) return("medium")
  "medium"
}

make_pcoa_plot <- function(sheet2) {
  sample_matrix <- as.data.frame(sheet2)
  sample_matrix[[1]] <- NULL
  sample_matrix <- as.data.frame(lapply(sample_matrix, function(x) replace(x, is.na(x), 0)))

  sample_meta <- tibble(sample = colnames(sample_matrix)) %>%
    mutate(
      genus = if_else(str_starts(sample, "Bacillus"), "Bacillus", "Rhodano"),
      habitat = if_else(str_detect(sample, "media"), "media", "sed"),
      group = paste(genus, habitat, sep = "_")
    ) %>%
    mutate(
      genus = factor(genus, levels = c("Bacillus", "Rhodano")),
      habitat = factor(habitat, levels = c("media", "sed")),
      group = factor(group, levels = c("Bacillus_media", "Bacillus_sed", "Rhodano_media", "Rhodano_sed"))
    )

  rownames(sample_matrix) <- rownames(sheet2)
  sample_matrix <- t(as.matrix(sample_matrix))

  gower_dist <- vegdist(sample_matrix, method = "gower")
  pcoa_fit <- cmdscale(gower_dist, k = 2, eig = TRUE)
  pcoa_coords <- as.data.frame(pcoa_fit$points) %>%
    rownames_to_column("sample") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    left_join(sample_meta, by = "sample")

  pos_eig <- pcoa_fit$eig[pcoa_fit$eig > 0]
  var_explained <- 100 * pos_eig / sum(pos_eig)

  ggplot(pcoa_coords, aes(PCoA1, PCoA2, color = genus, shape = habitat)) +
    geom_point(size = 4.1, stroke = 1.0) +
    geom_text_repel(aes(label = sample), size = 3.0, show.legend = FALSE, max.overlaps = Inf) +
    scale_color_manual(values = genus_palette) +
    scale_shape_manual(values = shape_palette) +
    coord_equal() +
    labs(
      title = NULL,
      subtitle = NULL,
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_explained[2]),
      color = "Genus",
      shape = "Habitat"
    ) +
    theme_classic(base_size = 10.5) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right",
      axis.title = element_text(face = "bold", size = 10.5),
      axis.text = element_text(size = 9)
    )
}

make_heatmap_panel <- function(sheet3) {
  sheet3_long <- sheet3 %>%
    mutate(
      chemical_group = vapply(Metabolites, classify_chemical_group, character(1)),
      Bacillus_absorb = pmax(Bacillus_absorb, 0),
      Rhoano_absorb = pmax(Rhoano_absorb, 0),
      max_absorb = pmax(Bacillus_absorb, Rhoano_absorb)
    ) %>%
    arrange(
      factor(chemical_group, levels = names(group_palette)),
      desc(max_absorb),
      Metabolites
    ) %>%
    mutate(metabolite_label = pretty_name(Metabolites))

  annot_df <- sheet3_long %>%
    select(metabolite_label, chemical_group)

  heat_long <- sheet3_long %>%
    select(metabolite_label, Bacillus_absorb, Rhoano_absorb) %>%
    pivot_longer(
      cols = c(Bacillus_absorb, Rhoano_absorb),
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = factor(source, levels = c("Bacillus_absorb", "Rhoano_absorb")),
      source = recode(source, Bacillus_absorb = "Bacillus", Rhoano_absorb = "Rhodano"),
      metabolite_label = factor(metabolite_label, levels = sheet3_long$metabolite_label)
    )

  top_bar <- ggplot(annot_df, aes(x = factor(metabolite_label, levels = sheet3_long$metabolite_label), y = 1, fill = chemical_group)) +
    geom_tile(height = 1) +
    scale_fill_manual(values = group_palette, guide = "none") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void(base_size = 8) +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )

  heat_plot <- ggplot(heat_long, aes(x = metabolite_label, y = source, fill = value)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradientn(
      colors = c("#FFFFFF", "#BFD7EA", "#08306B"),
      limits = c(0, max(heat_long$value, na.rm = TRUE)),
      name = "Absorbance"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = c("Bacillus", "Rhodano")) +
    labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
    theme_classic(base_size = 8.5) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2.7),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 7.5),
      plot.margin = margin(0, 2, 0, 2)
    )

  (top_bar / heat_plot) +
    plot_layout(heights = c(0.10, 1))
}

make_sheet3_scatter <- function(sheet4) {
  sheet4_long <- sheet4 %>%
    mutate(
      chemical_group = vapply(Metabolites, classify_chemical_group, character(1)),
      metabolite_label = pretty_name(Metabolites)
    )

  scatter_cor <- cor.test(
    sheet4_long$Bacillus_absorb,
    sheet4_long$Rhoano_absorb,
    method = "pearson"
  )

  scatter_label <- sprintf(
    "Pearson r = %.2f\np = %.2e",
    unname(scatter_cor$estimate),
    scatter_cor$p.value
  )

  ggplot(
    sheet4_long,
    aes(x = Bacillus_absorb, y = Rhoano_absorb, color = chemical_group)
  ) +
    geom_point(size = 3.2, alpha = 0.9) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      level = 0.95,
      color = "#1A1A1A",
      fill = "#A9A9A9",
      linewidth = 0.8,
      inherit.aes = FALSE,
      aes(x = Bacillus_absorb, y = Rhoano_absorb),
      data = sheet4_long
    ) +
    annotate(
      "label",
      x = max(sheet4_long$Bacillus_absorb, na.rm = TRUE) * 0.03,
      y = max(sheet4_long$Rhoano_absorb, na.rm = TRUE) * 0.98,
      label = scatter_label,
      hjust = 0,
      vjust = 1,
      size = 3.4,
      linewidth = 0.25,
      fill = "white",
      color = "#1A1A1A"
    ) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    coord_equal(
      xlim = c(0, max(sheet4_long$Bacillus_absorb, na.rm = TRUE) * 1.05),
      ylim = c(0, max(sheet4_long$Rhoano_absorb, na.rm = TRUE) * 1.05),
      expand = TRUE
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Bacillus absorbance",
      y = "Rhoano absorbance",
      color = "Chemical group"
    ) +
    theme_classic(base_size = 10.5) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "none",
      axis.title = element_text(face = "bold", size = 10.5),
      axis.text = element_text(size = 9)
    )
}

make_predicted_scatter <- function(sheet4) {
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

  sheet4 <- sheet4 %>%
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

  pred_long <- tidyr::expand_grid(obs_long, mineral = c("quartz", "ferrihydrite", "illite_smectite")) %>%
    left_join(charge_score, by = "mineral") %>%
    mutate(
      charge_weight = case_when(
        charge_class == "cationic" ~ cationic,
        charge_class == "zwitterionic" ~ zwitterionic,
        charge_class == "anionic" ~ anionic,
        TRUE ~ neutral
      ),
      diffusion_weight = unname(diffusion_score[size_class]),
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
      predicted = pmin(pmax(predicted, 0), 1)
    )

  cor_df <- pred_long %>%
    group_by(observation, mineral) %>%
    group_modify(~{
      test <- cor.test(.x$observed, .x$predicted, method = "pearson")
      tibble(
        r = unname(test$estimate),
        p_value = test$p.value,
        r2 = unname(test$estimate)^2,
        n = nrow(.x)
      )
    }) %>%
    ungroup() %>%
    mutate(
      label = sprintf(
        "Pearson r = %.2f\np = %s\nR^2 = %.2f\nn = %d",
        r,
        format.pval(p_value, digits = 2, eps = 0.001),
        r2,
        n
      ),
      x = 0.05,
      y = 0.95
    )

  ggplot(pred_long, aes(x = predicted, y = observed, color = chemical_group)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.55, color = "#5A5A5A") +
    geom_point(size = 2.4, alpha = 0.9) +
    geom_text(
      data = cor_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      lineheight = 0.95,
      color = "#202020"
    ) +
    facet_grid(
      observation ~ mineral,
      labeller = labeller(
        observation = c(
          "Bacillus observed" = "Bacillus absorbance",
          "Rhoano observed" = "Rhoano absorbance"
        )
      )
    ) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Predicted adsorption",
      y = "Observed adsorption",
      color = "Chemical group"
    ) +
    theme_classic(base_size = 10.2) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      strip.background = element_rect(fill = "#F3F3F3", color = "#D0D0D0"),
      strip.text = element_text(face = "bold", size = 9.2),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 9.2),
      legend.text = element_text(size = 8.7),
      axis.title = element_text(face = "bold", size = 10.2),
      axis.text = element_text(size = 9)
    )
}

sheet2 <- read_excel(workbook_path, sheet = 2)
sheet3 <- read_excel(workbook_path, sheet = 3)
sheet4 <- read_excel(workbook_path, sheet = 4)

pcoa_plot <- make_pcoa_plot(sheet2)
heatmap_panel <- wrap_elements(full = make_heatmap_panel(sheet3))
scatter_plot <- make_sheet3_scatter(sheet4)
predicted_plot <- make_predicted_scatter(sheet4)

composite <- (pcoa_plot | scatter_plot) / heatmap_panel / predicted_plot +
  plot_layout(heights = c(0.90, 0.70, 1.60)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.99)
  )

png_file <- file.path(out_dir, "targeted_metabolite_composite_AD.png")
pdf_file <- file.path(out_dir, "targeted_metabolite_composite_AD.pdf")

ggsave(png_file, composite, width = 14.0, height = 12.8, dpi = 320)
ggsave(pdf_file, composite, width = 14.0, height = 12.8)

cat("Wrote outputs to: ", out_dir, "\n", sep = "")
cat("Composite PNG: ", png_file, "\n", sep = "")
cat("Composite PDF: ", pdf_file, "\n", sep = "")

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
  library(vegan)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
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
        "hydroxybutyricacid", "oxalicacid", "adipicacid", "alpha-ketoglutaricacid",
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
        "folate", "tetrahydrofolate"
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

sheet2 <- read_excel(workbook_path, sheet = 2)
sheet3 <- read_excel(workbook_path, sheet = 3)
sheet4 <- read_excel(workbook_path, sheet = 4)

sample_names <- sheet2[[1]]
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
    group = factor(group, levels = c(
      "Bacillus_media", "Bacillus_sed", "Rhodano_media", "Rhodano_sed"
    ))
  )

rownames(sample_matrix) <- rownames(sheet2)
sample_matrix <- as.matrix(sample_matrix)
sample_matrix <- t(sample_matrix)

gower_dist <- vegdist(sample_matrix, method = "gower")
pcoa_fit <- cmdscale(gower_dist, k = 2, eig = TRUE)
pcoa_coords <- as.data.frame(pcoa_fit$points) %>%
  rownames_to_column("sample") %>%
  rename(PCoA1 = V1, PCoA2 = V2) %>%
  left_join(sample_meta, by = "sample")

pos_eig <- pcoa_fit$eig[pcoa_fit$eig > 0]
var_explained <- 100 * pos_eig / sum(pos_eig)

pcoa_permanova <- adonis2(gower_dist ~ genus * habitat, data = sample_meta, permutations = 999, by = "terms") %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(!term %in% c("Residual", "Total")) %>%
  rename(p_value = `Pr(>F)`)

group_permanova <- adonis2(gower_dist ~ group, data = sample_meta, permutations = 999, by = "terms") %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(!term %in% c("Residual", "Total")) %>%
  rename(p_value = `Pr(>F)`)

pairwise_permanova <- function(dist_obj, meta, group_col = "group") {
  grp <- meta[[group_col]]
  pairs <- combn(levels(grp), 2, simplify = FALSE)

  bind_rows(lapply(pairs, function(pair) {
    keep <- grp %in% pair
    sub_meta <- meta[keep, , drop = FALSE]
    sub_meta[[group_col]] <- droplevels(sub_meta[[group_col]])
    sub_dist <- as.dist(as.matrix(dist_obj)[keep, keep])
    fit <- adonis2(sub_dist ~ group, data = data.frame(group = sub_meta[[group_col]]), permutations = 999)
    tibble(
      group1 = pair[1],
      group2 = pair[2],
      Df = fit$Df[1],
      SumOfSqs = fit$SumOfSqs[1],
      R2 = fit$R2[1],
      F = fit$F[1],
      p = fit$`Pr(>F)`[1]
    )
  })) %>%
    mutate(p_adj = p.adjust(p, method = "BH")) %>%
    arrange(p_adj, p)
}

pairwise_group_permanova <- pairwise_permanova(gower_dist, sample_meta)

pcoa_plot <- ggplot(pcoa_coords, aes(PCoA1, PCoA2, color = genus, shape = habitat)) +
  geom_point(size = 4.2, stroke = 1.1) +
  geom_text_repel(aes(label = sample), size = 3.2, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = genus_palette) +
  scale_shape_manual(values = shape_palette) +
  coord_cartesian(
    xlim = range(pcoa_coords$PCoA1) + c(-1, 1) * diff(range(pcoa_coords$PCoA1)) * 0.08,
    ylim = range(pcoa_coords$PCoA2) + c(-1, 1) * diff(range(pcoa_coords$PCoA2)) * 0.08,
    expand = FALSE
  ) +
  labs(
    title = "PCoA of Sheet 2 using Gower distance",
    x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
    y = sprintf("PCoA 2 (%.1f%%)", var_explained[2]),
    color = "Genus",
    shape = "Habitat"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    aspect.ratio = 1
  )

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

heat_matrix <- sheet3_long %>%
  select(Metabolites, Bacillus_absorb, Rhoano_absorb) %>%
  as.data.frame()
rownames(heat_matrix) <- sheet3_long$metabolite_label
heat_matrix <- as.matrix(heat_matrix[, c("Bacillus_absorb", "Rhoano_absorb")])
heat_matrix <- t(heat_matrix)
rownames(heat_matrix) <- c("Bacillus", "Rhodano")
colnames(heat_matrix) <- sheet3_long$metabolite_label

heat_group <- factor(sheet3_long$chemical_group, levels = names(group_palette))
names(heat_group) <- sheet3_long$metabolite_label

heat_col_fun <- colorRamp2(
  c(0, max(heat_matrix, na.rm = TRUE) * 0.5, max(heat_matrix, na.rm = TRUE)),
  c("#FFFFFF", "#BFD7EA", "#08306B")
)

heat_annotation <- HeatmapAnnotation(
  `Chemical group` = heat_group,
  col = list(`Chemical group` = group_palette),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_legend_param = list(title = "Chemical group")
)

heatmap_obj <- Heatmap(
  heat_matrix,
  name = "Absorbance",
  col = heat_col_fun,
  top_annotation = heat_annotation,
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 90,
  column_names_centered = FALSE,
  column_names_gp = grid::gpar(fontsize = 8.5),
  column_names_max_height = unit(65, "mm"),
  column_split = heat_group,
  column_gap = unit(1.2, "mm"),
  heatmap_height = unit(80, "mm"),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_title = "Sheet 3 absorbance after setting negatives to 0",
  row_title = "Sediment source",
  heatmap_legend_param = list(title = "Absorbance"),
  border = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 0.25)
)

absorbance_legend <- Legend(
  title = "Absorbance",
  col_fun = heat_col_fun,
  direction = "horizontal",
  legend_width = unit(32, "mm"),
  labels_gp = grid::gpar(fontsize = 8),
  title_gp = grid::gpar(fontsize = 9, fontface = "bold")
)

group_legend <- Legend(
  title = "Chemical group",
  at = names(group_palette),
  labels = names(group_palette),
  legend_gp = grid::gpar(fill = group_palette),
  direction = "horizontal",
  nrow = 1,
  by_row = TRUE,
  labels_gp = grid::gpar(fontsize = 8),
  title_gp = grid::gpar(fontsize = 9, fontface = "bold")
)

legend_row <- packLegend(
  absorbance_legend,
  group_legend,
  direction = "horizontal",
  gap = unit(6, "mm")
)

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

scatter_plot <- ggplot(
  sheet4_long,
  aes(x = Bacillus_absorb, y = Rhoano_absorb, color = chemical_group)
) +
  geom_point(size = 3.4, alpha = 0.9) +
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
  geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.6, color = "#666666") +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.6, color = "#666666") +
  annotate(
    "text",
    x = max(sheet4_long$Bacillus_absorb, na.rm = TRUE),
    y = 0.055,
    label = "5% threshold",
    hjust = 1,
    vjust = -0.2,
    size = 3.5,
    color = "#555555"
  ) +
  annotate(
    "label",
    x = max(sheet4_long$Bacillus_absorb, na.rm = TRUE) * 0.03,
    y = max(sheet4_long$Rhoano_absorb, na.rm = TRUE) * 0.98,
    label = scatter_label,
    hjust = 0,
    vjust = 1,
    size = 3.6,
    label.size = 0.25,
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
    title = "Sheet 4 scatter plot of metabolite absorbance",
    x = "Bacillus absorbance",
    y = "Rhoano absorbance",
    color = "Chemical group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

pcoa_png <- file.path(out_dir, "sheet2_pcoa_gower.png")
pcoa_pdf <- file.path(out_dir, "sheet2_pcoa_gower.pdf")
heat_png <- file.path(out_dir, "sheet3_absorbance_heatmap.png")
heat_pdf <- file.path(out_dir, "sheet3_absorbance_heatmap.pdf")
scatter_png <- file.path(out_dir, "sheet4_scatter_groups.png")
scatter_pdf <- file.path(out_dir, "sheet4_scatter_groups.pdf")

ggsave(pcoa_png, pcoa_plot, width = 7.2, height = 7.2, dpi = 320)
ggsave(pcoa_pdf, pcoa_plot, width = 7.2, height = 7.2)
png(heat_png, width = 8800, height = 2400, res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(
  layout = grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit.c(grid::unit(1, "null"), ComplexHeatmap:::height(legend_row) + grid::unit(3, "mm"))
  )
))
grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heatmap_obj, newpage = FALSE)
grid::popViewport()
grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.draw(legend_row)
grid::popViewport()
grid::popViewport()
dev.off()
pdf(heat_pdf, width = 30, height = 8.5)
grid::grid.newpage()
grid::pushViewport(grid::viewport(
  layout = grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit.c(grid::unit(1, "null"), ComplexHeatmap:::height(legend_row) + grid::unit(3, "mm"))
  )
))
grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heatmap_obj, newpage = FALSE)
grid::popViewport()
grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
grid::grid.draw(legend_row)
grid::popViewport()
grid::popViewport()
dev.off()
ggsave(scatter_png, scatter_plot, width = 7.4, height = 6.1, dpi = 320)
ggsave(scatter_pdf, scatter_plot, width = 7.4, height = 6.1)

write.csv(pcoa_permanova, file.path(out_dir, "sheet2_permanova_main_effects.csv"), row.names = FALSE)
write.csv(group_permanova, file.path(out_dir, "sheet2_permanova_group_model.csv"), row.names = FALSE)
write.csv(pairwise_group_permanova, file.path(out_dir, "sheet2_permanova_pairwise_groups.csv"), row.names = FALSE)
write.csv(pcoa_coords, file.path(out_dir, "sheet2_pcoa_coordinates.csv"), row.names = FALSE)

group_counts <- sheet4_long %>%
  count(chemical_group, sort = TRUE)
write.csv(group_counts, file.path(out_dir, "chemical_group_counts_sheet4.csv"), row.names = FALSE)

cat("Wrote outputs to: ", out_dir, "\n", sep = "")
cat("\nPERMANOVA main effects:\n")
print(pcoa_permanova)
cat("\nPERMANOVA group model:\n")
print(group_permanova)
cat("\nPairwise PERMANOVA:\n")
print(pairwise_group_permanova)

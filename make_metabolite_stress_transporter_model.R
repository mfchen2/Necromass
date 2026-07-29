#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

theme_nature <- function(base_size = 7.8, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(face = "bold", size = base_size, colour = "#111827"),
      axis.text = element_text(size = base_size - 0.35, colour = "#1F2937"),
      plot.title = element_text(face = "bold", size = base_size + 1.1, colour = "#111827"),
      plot.subtitle = element_text(size = base_size - 0.15, colour = "#4B5563"),
      strip.background = element_rect(fill = "#F3F4F6", colour = "#D1D5DB"),
      strip.text = element_text(face = "bold", size = base_size - 0.05, colour = "#111827"),
      legend.title = element_text(face = "bold", size = base_size - 0.1, colour = "#111827"),
      legend.text = element_text(size = base_size - 0.2, colour = "#1F2937"),
      legend.key.height = unit(3.4, "mm"),
      legend.key.width = unit(8, "mm"),
      plot.margin = margin(4, 6, 4, 6)
    )
}

save_pub <- function(plot, filename, width_mm = 183, height_mm = 163, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot)
  dev.off()
}

clean_label <- function(x) {
  x %>%
    str_replace_all("\n", " ") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

out_dir <- "/Users/mingfeichen/Manuscript/outputs/manual-20260608-a1/presentations/metabolite_stress_transporter_model/assets"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

met_file <- "/Users/mingfeichen/Manuscript/outputs/manual-20260601-a6/presentations/metabolite_ttest_remake/assets/metabolite_ttest_bubble_data.csv"
tx_file <- "/Users/mingfeichen/Manuscript/outputs/manual-20260603-a4/presentations/transcriptome_pathway_modules/assets/transcriptome_pathway_module_data.csv"

palette_contract <- c(
  phage = "#8C6BB1",
  stress = "#D58B4D",
  up = "#D1495B",
  down = "#3B73B9",
  neutral = "#EEF2F7",
  dark = "#111827",
  mid = "#4B5563",
  light = "#CBD5E1"
)

theme_set(theme_nature())

# -------------------------
# Panel a: mechanism sketch
# -------------------------
schem_nodes <- tibble::tribble(
  ~panel, ~x, ~y, ~label, ~fill, ~text_colour,
  "phage", 2.7, 2.65, "Phage\nlysis", palette_contract["phage"], "white",
  "phage", 2.7, 1.60, "Immediate\nrelease", "#F2EEF8", palette_contract["dark"],
  "phage", 2.7, 0.50, "Short\npulse", "#FFFFFF", palette_contract["dark"],
  "stress", 7.3, 2.65, "Al/Kana\nstress", palette_contract["stress"], "white",
  "stress", 7.3, 1.60, "Membrane\nleak", "#FDF5EC", palette_contract["dark"],
  "stress", 7.3, 0.50, "Sustained\nrelease", "#FFFFFF", palette_contract["dark"]
)

schem_boxes <- tibble::tribble(
  ~xmin, ~xmax, ~ymin, ~ymax, ~fill,
  1.3, 4.1, 2.10, 3.10, palette_contract["phage"],
  1.3, 4.1, 1.05, 1.95, "#F5F0FB",
  1.3, 4.1, -0.05, 0.95, "#FFFFFF",
  5.9, 8.7, 2.10, 3.10, palette_contract["stress"],
  5.9, 8.7, 1.05, 1.95, "#FFF2E4",
  5.9, 8.7, -0.05, 0.95, "#FFFFFF"
)

curve_df <- tibble::tribble(
  ~x, ~xend, ~y, ~yend, ~colour,
  2.7, 2.7, 2.35, 2.05, palette_contract["phage"],
  2.7, 2.7, 2.35, 2.05, palette_contract["stress"]
)

spark_phage <- tibble(
  x = c(1.45, 1.85, 2.1, 2.35, 2.55, 2.65, 2.75, 2.95, 3.2, 3.55, 3.9),
  y = c(0.30, 0.32, 0.35, 0.45, 1.65, 0.48, 0.38, 0.34, 0.32, 0.30, 0.29)
)
spark_stress <- tibble(
  x = c(6.15, 6.45, 6.75, 7.05, 7.35, 7.65, 7.95, 8.25, 8.55),
  y = c(0.30, 0.38, 0.52, 0.60, 0.62, 0.58, 0.52, 0.42, 0.34)
)

panel_a <- ggplot() +
  geom_rect(data = schem_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            colour = NA) +
  scale_fill_identity() +
  geom_curve(
    aes(x = 2.7, y = 2.35, xend = 2.7, yend = 2.05),
    curvature = 0.0, linewidth = 0.55, arrow = arrow(length = unit(1.6, "mm")),
    colour = palette_contract["phage"]
  ) +
  geom_curve(
    aes(x = 7.3, y = 2.35, xend = 7.3, yend = 2.05),
    curvature = 0.0, linewidth = 0.55, arrow = arrow(length = unit(1.6, "mm")),
    colour = palette_contract["stress"]
  ) +
  geom_line(data = spark_phage, aes(x = x, y = y), linewidth = 0.8, colour = palette_contract["phage"]) +
  geom_line(data = spark_stress, aes(x = x, y = y), linewidth = 0.8, colour = palette_contract["stress"]) +
  geom_point(data = spark_phage, aes(x = x, y = y), size = 1.2, colour = palette_contract["phage"]) +
  geom_point(data = spark_stress, aes(x = x, y = y), size = 1.2, colour = palette_contract["stress"]) +
  geom_label(
    data = dplyr::filter(schem_nodes, y == 2.7),
    aes(x = x, y = y, label = label, fill = fill, colour = text_colour),
    linewidth = 0,
    size = 3.05,
    fontface = "bold",
    label.padding = unit(0.18, "lines"),
    lineheight = 0.95
  ) +
  geom_label(
    data = dplyr::filter(schem_nodes, y != 2.7),
    aes(x = x, y = y, label = label, fill = fill, colour = text_colour),
    linewidth = 0,
    size = 3.05,
    fontface = "plain",
    label.padding = unit(0.18, "lines"),
    lineheight = 0.95
  ) +
  scale_colour_identity() +
  annotate("segment", x = 1.0, xend = 4.4, y = 0.0, yend = 0.0, linewidth = 0.35, colour = palette_contract["light"]) +
  annotate("segment", x = 5.6, xend = 8.9, y = 0.0, yend = 0.0, linewidth = 0.35, colour = palette_contract["light"]) +
  coord_cartesian(xlim = c(0.5, 9.5), ylim = c(-0.15, 3.45), clip = "off") +
  theme_void(base_size = 7.8) +
  theme(
    plot.margin = margin(4, 8, 2, 8)
  )

# -------------------------
# Panel b: metabolite matrix
# -------------------------
met_catalog <- tibble::tribble(
  ~feature_label, ~superclass, ~order,
  "Betaine", "Lipids & osmolytes", 1,
  "Carnitine", "Lipids & osmolytes", 2,
  "sn-Glycero-3-Phosphocholine", "Lipids & osmolytes", 3,
  "N-trimethyllysine", "Amino acids & derivatives", 4,
  "N-alpha-acetyl-lysine", "Amino acids & derivatives", 5,
  "N-Acetyl-Glutamine", "Amino acids & derivatives", 6,
  "4-guanidinobutanoic acid", "Polyamines & guanidino compounds", 7,
  "Adenosine", "Nucleic acids", 8,
  "Guanosine", "Nucleic acids", 9,
  "Cytosine", "Nucleic acids", 10,
  "Xanthine", "Nucleic acids", 11,
  "Orotic acid", "Nucleic acids", 12,
  "5-methylcytosine", "Nucleic acids", 13,
  "3-methyladenine", "Nucleic acids", 14,
  "pterin", "Cofactors & vitamins", 15,
  "2 4-dihydroxypteridine", "Cofactors & vitamins", 16,
  "2-hydroxybutyric acid", "Central carbon & organic acids", 17,
  "glycerol 2-phosphoric acid", "Central carbon & organic acids", 18,
  "alpha-hydroxyisobutyric acid", "Central carbon & organic acids", 19,
  "4-hydroxybenzoic acid", "Aromatic compounds & phenolics", 20
) %>%
  mutate(feature_label = clean_label(feature_label))

met_raw <- read_csv(met_file, show_col_types = FALSE) %>%
  mutate(
    contrast_display = factor(contrast_display, levels = c("Al", "Kana")),
    feature_label = clean_label(feature_label),
    superclass = case_when(
      superclass == "Lipids & osmolytes (quaternary amines)" ~ "Lipids & osmolytes",
      superclass == "Amino acids & derivatives" ~ "Amino acids & derivatives",
      superclass == "Nucleic acids (bases/nucleosides/nucleotides)" ~ "Nucleic acids",
      superclass == "Cofactors (pterins & vitamins)" ~ "Cofactors & vitamins",
      superclass == "Central carbon & organic acids" ~ "Central carbon & organic acids",
      superclass == "Polyamines & guanidino compounds" ~ "Polyamines & guanidino compounds",
      superclass == "Aromatic compounds & phenolics" ~ "Aromatic compounds & phenolics",
      TRUE ~ superclass
    )
  ) %>%
  filter(contrast_display %in% c("Al", "Kana"), padj < 0.05, feature_label %in% met_catalog$feature_label)

met_sel <- met_raw %>%
  left_join(met_catalog, by = "feature_label") %>%
  group_by(feature_label, contrast_display) %>%
  arrange(padj, desc(abs(log2FC)), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    feature_label = factor(feature_label, levels = rev(met_catalog$feature_label[order(met_catalog$order)])),
    superclass = factor(met_catalog$superclass[match(as.character(feature_label), met_catalog$feature_label)],
                        levels = c(
                          "Lipids & osmolytes",
                          "Amino acids & derivatives",
                          "Polyamines & guanidino compounds",
                          "Nucleic acids",
                          "Cofactors & vitamins",
                          "Central carbon & organic acids",
                          "Aromatic compounds & phenolics"
                        ))
  )

met_base <- met_catalog %>%
  mutate(
    feature_label = factor(feature_label, levels = rev(met_catalog$feature_label[order(met_catalog$order)])),
    superclass = factor(
      superclass,
      levels = c(
        "Lipids & osmolytes",
        "Amino acids & derivatives",
        "Polyamines & guanidino compounds",
        "Nucleic acids",
        "Cofactors & vitamins",
        "Central carbon & organic acids",
        "Aromatic compounds & phenolics"
      )
    )
  )

met_plot <- tidyr::expand_grid(
  met_base,
  contrast_display = factor(c("Al", "Kana"), levels = c("Al", "Kana"))
) %>%
  left_join(
    met_sel %>%
      select(feature_label, contrast_display, log2FC, padj),
    by = c("feature_label", "contrast_display")
  )

panel_b <- ggplot(met_plot, aes(x = contrast_display, y = feature_label)) +
  geom_tile(fill = "#FAFBFD", colour = "white", linewidth = 0.22) +
  geom_point(
    data = met_plot %>% filter(!is.na(log2FC)),
    aes(fill = log2FC, size = -log10(padj)),
    shape = 21, colour = "#1F2937", stroke = 0.22, alpha = 0.96
  ) +
  facet_grid(superclass ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_fill_gradient2(
    low = "#3B73B9",
    mid = "white",
    high = "#D1495B",
    midpoint = 0,
    name = "log2FC\n(B vs R)"
  ) +
  scale_size_continuous(range = c(1.8, 5.0), name = "-log10(padj)") +
  labs(
    title = "Al/Kana-enriched exometabolites",
    x = NULL,
    y = "Representative enriched metabolites"
  ) +
  theme_nature(base_size = 7.3) +
  theme(
    plot.title = element_text(size = 8.9, face = "bold", hjust = 0),
    axis.text.x = element_text(size = 7.2),
    axis.text.y = element_text(size = 6.9),
    strip.text.y.left = element_text(angle = 0, size = 7.3),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(colour = "#EEF2F7", linewidth = 0.22),
    panel.spacing.y = unit(0.10, "in"),
    plot.margin = margin(3, 8, 3, 8)
  )

# -------------------------
# Panel c: transporter transcriptome
# -------------------------
tx_raw <- read_csv(tx_file, show_col_types = FALSE) %>%
  mutate(
    treatment_core = str_remove(treatment_display, "^[BR] "),
    treatment_core = factor(treatment_core, levels = c("mid", "late", "Al", "Kana", "Phage")),
    organism = factor(organism, levels = c("Bacillus", "Rhodanobacter")),
    pathway_label = paste0(pathway_name, " (", pathway_id, ")")
  )

tx_catalog <- tibble::tribble(
  ~organism, ~pathway_id, ~pathway_name, ~pathway_class, ~order,
  "Bacillus", "map02010", "ABC transporters", "Transport / signaling", 1,
  "Bacillus", "map01110", "Biosynthesis of secondary metabolites", "Secondary metabolism / transport", 2,
  "Rhodanobacter", "map02020", "Two-component system", "Transport / signaling", 1,
  "Rhodanobacter", "map02030", "Bacterial chemotaxis", "Transport / signaling", 2,
  "Rhodanobacter", "map02040", "Flagellar assembly", "Transport / signaling", 3,
  "Rhodanobacter", "map02025", "Biofilm formation", "Transport / signaling", 4,
  "Rhodanobacter", "map01110", "Biosynthesis of secondary metabolites", "Secondary metabolism / transport", 5
)
tx_catalog_meta <- tx_catalog %>%
  select(organism, pathway_id, pathway_class, order)

tx_levels <- c(
  "ABC transporters (map02010)",
  "Biosynthesis of secondary metabolites (map01110)",
  "Two-component system (map02020)",
  "Bacterial chemotaxis (map02030)",
  "Flagellar assembly (map02040)",
  "Biofilm formation (map02025)"
)

tx_sel <- tx_raw %>%
  filter(str_detect(pathway_id, "^map")) %>%
  semi_join(tx_catalog_meta, by = c("organism", "pathway_id")) %>%
  left_join(tx_catalog_meta, by = c("organism", "pathway_id")) %>%
  mutate(
    pathway_label = paste0(pathway_name, " (", pathway_id, ")"),
    pathway_label = factor(pathway_label, levels = rev(tx_levels))
  ) %>%
  group_by(organism, pathway_id, treatment_core) %>%
  arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

tx_base <- tx_catalog %>%
  semi_join(tx_sel %>% distinct(organism, pathway_id), by = c("organism", "pathway_id")) %>%
  mutate(
    pathway_label = factor(
      paste0(pathway_name, " (", pathway_id, ")"),
      levels = rev(tx_levels)
    )
  ) %>%
  select(organism, pathway_label, pathway_class, order)

tx_plot <- tidyr::expand_grid(
  tx_base,
  treatment_core = factor(c("mid", "late", "Al", "Kana", "Phage"), levels = c("mid", "late", "Al", "Kana", "Phage"))
) %>%
  left_join(
    tx_sel %>%
      select(organism, pathway_label, treatment_core, NES, padj),
    by = c("organism", "pathway_label", "treatment_core")
  )

panel_c <- ggplot(tx_plot, aes(x = treatment_core, y = pathway_label)) +
  geom_tile(fill = "#FAFBFD", colour = "white", linewidth = 0.22) +
  geom_point(
    aes(fill = NES, size = -log10(padj)),
    shape = 21, colour = "#1F2937", stroke = 0.22, alpha = 0.96
  ) +
  facet_wrap(~ organism, ncol = 1, scales = "free_y") +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_fill_gradient2(
    low = "#3B73B9",
    mid = "white",
    high = "#D1495B",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(range = c(1.8, 5.0), name = "-log10(padj)") +
  labs(
    title = "Transport and sensing pathways",
    x = NULL,
    y = "Selected pathways"
  ) +
  theme_nature(base_size = 7.3) +
  theme(
    plot.title = element_text(size = 8.9, face = "bold", hjust = 0),
    axis.text.x = element_text(size = 7.2),
    axis.text.y = element_text(size = 6.9),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(colour = "#EEF2F7", linewidth = 0.22),
    panel.spacing.y = unit(0.12, "in"),
    plot.margin = margin(3, 8, 3, 8)
  )

write_csv(met_plot, file.path(out_dir, "metabolite_stress_transporter_model_metabolites.csv"))
write_csv(tx_sel, file.path(out_dir, "metabolite_stress_transporter_model_transcriptome.csv"))

figure_pair <- (panel_b | panel_c) +
  plot_layout(widths = c(1.2, 0.85), guides = "keep") +
  plot_annotation(
    title = "Al/Kana stress couples exometabolite accumulation with transport remodeling",
    subtitle = "The metabolite panel highlights the released pool; the transcriptome panel shows matching transport and sensing responses.",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11.2, colour = "#111827", hjust = 0.5),
      plot.subtitle = element_text(size = 8.1, colour = "#4B5563", hjust = 0.5)
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 9))

save_pub(panel_b, file.path(out_dir, "metabolite_stress_metabolites"), width_mm = 203, height_mm = 152)
save_pub(panel_c, file.path(out_dir, "metabolite_stress_transporters"), width_mm = 203, height_mm = 152)
save_pub(figure_pair, file.path(out_dir, "metabolite_stress_linked_pair"), width_mm = 406, height_mm = 165)

message("Wrote: ", file.path(out_dir, "metabolite_stress_metabolites.png"))
message("Wrote: ", file.path(out_dir, "metabolite_stress_metabolites.pdf"))
message("Wrote: ", file.path(out_dir, "metabolite_stress_metabolites.svg"))
message("Wrote: ", file.path(out_dir, "metabolite_stress_transporters.png"))
message("Wrote: ", file.path(out_dir, "metabolite_stress_transporters.pdf"))
message("Wrote: ", file.path(out_dir, "metabolite_stress_transporters.svg"))
message("Wrote: ", file.path(out_dir, "metabolite_stress_linked_pair.png"))
message("Wrote: ", file.path(out_dir, "metabolite_stress_linked_pair.pdf"))
message("Wrote: ", file.path(out_dir, "metabolite_stress_linked_pair.svg"))

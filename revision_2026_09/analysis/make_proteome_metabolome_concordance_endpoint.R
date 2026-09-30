#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

root <- "/Users/mingfeichen/Manuscript"
out <- file.path(root, "outputs/manual-20260921-a1/presentations/proteome_metabolome_concordance_endpoint_matched/assets")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

met_path <- file.path(root, "outputs/manual-20260921-a1/presentations/pathway_metabolite_relationships_endpoint_matched/assets/figure_a_scatter_points.csv")
b_prot_path <- file.path(root, "outputs/manual-20260824-a1/presentations/proteome_mineral_support/assets/Bacillus_Proteome_fgsea_Results.csv")
r_prot_path <- file.path(root, "outputs/manual-20260824-a1/presentations/proteome_mineral_support/assets/Rhodano_Proteome_fgsea_Results.csv")

met <- read_csv(met_path, show_col_types = FALSE) %>%
  transmute(
    organism,
    treatment_raw = as.character(treatment),
    pathway_id,
    met_score = as.numeric(met_score),
    n_sig = as.numeric(n_sig)
  )

read_proteome <- function(path, organism) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(str_detect(as.character(pathway), "^map")) %>%
    transmute(
      organism = organism,
      treatment_raw = as.character(Condition),
      pathway_id = as.character(pathway),
      protein_score = as.numeric(NES),
      protein_padj = as.numeric(padj)
    ) %>%
    distinct(organism, treatment_raw, pathway_id, .keep_all = TRUE)
}

prot <- bind_rows(
  read_proteome(b_prot_path, "Bacillus"),
  read_proteome(r_prot_path, "Rhodanobacter")
)

display_treatment <- function(organism, treatment) {
  case_when(
    organism == "Bacillus" & treatment == "8hr" ~ "mid",
    organism == "Bacillus" & treatment == "24hr" ~ "late",
    organism == "Rhodanobacter" & treatment == "24hr" ~ "mid",
    organism == "Rhodanobacter" & treatment == "48hr" ~ "late",
    TRUE ~ treatment
  )
}

joined <- prot %>%
  inner_join(met, by = c("organism", "treatment_raw", "pathway_id")) %>%
  mutate(
    treatment = display_treatment(organism, treatment_raw),
    relationship = case_when(
      is.na(protein_score) | is.na(met_score) ~ "missing",
      protein_score == 0 | met_score == 0 ~ "weak/zero",
      sign(protein_score) == sign(met_score) ~ "same direction",
      TRUE ~ "opposite direction"
    ),
    # Retain pathways supported by either omic layer; the significance fields
    # remain available in the exported table for stricter downstream filtering.
    informative = protein_padj <= 0.05 | n_sig >= 1
  )

write_csv(joined, file.path(out, "proteome_metabolome_endpoint_scatter_data.csv"))

plot_df <- joined %>% filter(informative)
organism_levels <- c("Bacillus", "Rhodanobacter")
treatment_levels <- c("mid", "late", "Al", "K", "P")
plot_df <- plot_df %>%
  mutate(
    organism = factor(organism, levels = organism_levels),
    treatment = factor(treatment, levels = treatment_levels),
    relationship = factor(relationship, levels = c("same direction", "opposite direction", "weak/zero"))
  )

relationship_cols <- c("same direction" = "#D95F6A", "opposite direction" = "#3F78B5", "weak/zero" = "#9CA3AF")

panel_a <- ggplot(plot_df, aes(x = protein_score, y = met_score, color = relationship, size = pmax(n_sig, 1))) +
  geom_hline(yintercept = 0, color = "grey88", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey88", linewidth = 0.35) +
  geom_point(alpha = 0.82, stroke = 0.2) +
  facet_grid(organism ~ treatment, drop = FALSE, scales = "free") +
  scale_color_manual(values = relationship_cols, drop = FALSE) +
  scale_size_continuous(range = c(2.0, 7.0), name = "n significant\nmetabolites") +
  labs(
    title = "Proteome–exometabolome pathway concordance",
    subtitle = "Available proteome FGSEA contrasts joined to the endpoint-matched exometabolome",
    x = "Proteome pathway enrichment (NES)",
    y = "Observed exometabolite pathway score",
    color = "Relationship"
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.spacing = unit(0.65, "lines")
  )

examples <- plot_df %>%
  mutate(rank_score = abs(protein_score) + abs(met_score)) %>%
  filter(relationship %in% c("same direction", "opposite direction")) %>%
  group_by(organism, relationship) %>%
  arrange(desc(rank_score), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(label = paste0(pathway_id, " ", organism, " / ", treatment))

write_csv(examples, file.path(out, "proteome_metabolome_endpoint_representative_examples.csv"))

panel_b <- ggplot(examples, aes(x = met_score, y = reorder(label, met_score), xend = protein_score, yend = reorder(label, met_score), color = relationship)) +
  geom_segment(linewidth = 0.8, alpha = 0.85) +
  geom_point(aes(x = met_score), shape = 21, fill = "white", size = 3.0, stroke = 0.7) +
  geom_point(aes(x = protein_score), shape = 22, fill = "white", size = 3.0, stroke = 0.7) +
  facet_wrap(~organism, ncol = 2, scales = "free_y") +
  scale_color_manual(values = relationship_cols) +
  labs(
    title = "Representative proteome–exometabolome pathway relationships",
    x = "Pathway score (metabolome left; proteome right)",
    y = NULL,
    color = "Relationship"
  ) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

fig <- panel_a / panel_b + plot_layout(heights = c(1.35, 1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

save_figure <- function(plot, base, width = 14, height = 13) {
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(base, ".svg"), plot, width = width, height = height, device = svglite::svglite)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 400, device = ragg::agg_png)
}

save_figure(fig, file.path(out, "proteome_metabolome_concordance"))

run_log <- c(
  "Proteome–exometabolome concordance using endpoint-matched metabolome scores",
  paste0("Joined informative pathway observations: ", nrow(plot_df)),
  paste0("Available proteome pathway exports: Bacillus conditions = ", paste(sort(unique(prot$treatment_raw[prot$organism == "Bacillus"])), collapse = ", "), "; Rhodanobacter conditions = ", paste(sort(unique(prot$treatment_raw[prot$organism == "Rhodanobacter"])), collapse = ", ")),
  "Bacillus K is left blank because no Bacillus K proteome FGSEA export was available.",
  "The archived proteome FGSEA exports do not include explicit reference metadata; this limitation is retained in the interpretation."
)
writeLines(run_log, file.path(out, "run_log.txt"))

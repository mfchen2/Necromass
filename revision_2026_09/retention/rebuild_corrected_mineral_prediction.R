#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(readr); library(lme4); library(ragg); library(svglite)
})

input <- "/Users/mingfeichen/Downloads/New_metabolites_corrected.xlsx"
out <- "/Users/mingfeichen/Manuscript/outputs/manual-20260810-a2/presentations/corrected_mineral_sorption_prediction/assets"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

norm <- function(x) str_replace_all(str_to_lower(str_trim(as.character(x))), "[^a-z0-9]", "")
group <- function(x) {
  z <- norm(x)
  if (str_detect(z, "adenosine|adenine|guanosine|guanine|cytidine|cytosine|uridine|uracil|thymine|hypoxanthine|xanthine|inosine|deoxy|orotic")) return("Nucleosides & bases")
  if (str_detect(z, "alanine|arginine|asparagine|aspartic|citrulline|cysteine|glutamine|glutamic|glycine|histidine|histidinol|isoleucine|leucine|lysine|methionine|ornithine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine|aminocaproic|aminobutyric|acetyllysine|acetylglutamine")) return("Amino acids & peptides")
  if (str_detect(z, "lactic|succinic|malic|fumaric|citric|pyruvic|shikimic|salicylic|benzoic|hydroxybenzoic|glyceric|glycolic|hydroxybutyric|oxalic|adipic|ketoglutaric|glucuronic|glycerophosphoric|phosphoric|nicotinic")) return("Organic acids")
  if (str_detect(z, "nicotinamide|pyridox|riboflavin|thiamine|biotin|pantothenic|folic|folate|tetrahydrofolate|pterin")) return("Vitamins & cofactors")
  if (str_detect(z, "glucose|fructose|galactose|mannose|ribose|arabinose|xylose|maltose|sucrose|lactose|glucosamine|glycerol|sorbitol|mannitol|trehalose")) return("Carbohydrates & sugars")
  if (str_detect(z, "agmatine|betaine|carnitine|choline|putrescine|spermidine|spermine|trimethyllysine|dimethylglycine|tyramine|guanidino")) return("Polyamines & amines")
  "Other / synthetic"
}

charge <- function(x, g) {
  z <- norm(x)
  if (z %in% c("lysine","arginine","choline","tyramine","histidinol","ntrimethyllysine","4guanidinobutanoicacid","creatinine") || g == "Polyamines & amines") return("cationic")
  if (g == "Organic acids") return("anionic")
  if (g == "Amino acids & peptides") return("zwitterionic")
  if (z %in% c("betaine","carnitine","histidine","citrulline","alanine","isoleucine","leucine","glutamine","tyrosine","cytidine","guanosine","pterin","nalphaacetyllysine")) return("zwitterionic")
  "neutral"
}

mineral <- read_excel(input, sheet = "Supplememntary Table Y Mineral")
names(mineral)[1] <- "metabolite"
long <- mineral %>%
  pivot_longer(-metabolite, names_to = "sample", values_to = "peak") %>%
  mutate(source = if_else(str_detect(sample, "Bacillus"), "Bacillus", "Rhodanobacter"),
         matrix = if_else(str_detect(sample, "-NA_"), "NA", "Sed"),
         replicate = str_extract(sample, "_[123]$"))

wide <- long %>%
  select(metabolite, source, matrix, replicate, peak) %>%
  pivot_wider(names_from = matrix, values_from = peak) %>%
  rename(NA_peak = `NA`, Sed_peak = Sed) %>%
  group_by(metabolite, source) %>%
  summarise(NA_peak = mean(NA_peak, na.rm = TRUE), Sed_peak = mean(Sed_peak, na.rm = TRUE),
            n_NA = sum(is.finite(NA_peak)), n_Sed = sum(is.finite(Sed_peak)), .groups = "drop") %>%
  mutate(raw_partitioning = (NA_peak - Sed_peak) / NA_peak,
         adsorption = pmin(pmax(raw_partitioning, 0), 1),
         chemical_group = vapply(metabolite, group, character(1)),
         charge_class = vapply(seq_along(metabolite), function(i) charge(metabolite[i], chemical_group[i]), character(1)))

audit <- wide %>% summarise(n_metabolites = n_distinct(metabolite), n_rows = n(), missing_NA = sum(is.na(NA_peak)), missing_Sed = sum(is.na(Sed_peak)), negative_raw = sum(raw_partitioning < 0, na.rm=TRUE), above_one_raw = sum(raw_partitioning > 1, na.rm=TRUE))
write_csv(wide, file.path(out, "corrected_mineral_adsorption_by_replicate.csv"))
write_csv(audit, file.path(out, "corrected_mineral_data_audit.csv"))

summary <- wide %>% group_by(source, chemical_group, charge_class) %>% summarise(n=n(), mean_adsorption=mean(adsorption, na.rm=TRUE), sd_adsorption=sd(adsorption, na.rm=TRUE), .groups="drop")
write_csv(summary, file.path(out, "corrected_mineral_adsorption_summary.csv"))

model <- lmer(qlogis(pmin(pmax(adsorption, 1e-4), 1-1e-4)) ~ charge_class + source + (1|metabolite), data=wide)
ci <- suppressMessages(confint(model, parm="beta_", method="Wald"))
coef <- data.frame(term=names(fixef(model)), estimate=unname(fixef(model)), conf_low=ci[,1], conf_high=ci[,2])
write_csv(coef, file.path(out, "corrected_charge_class_model_coefficients.csv"))
capture.output(summary(model), file=file.path(out, "corrected_charge_class_model_summary.txt"))

metabs <- unique(wide$metabolite)
cv <- lapply(metabs, function(m) {
  train <- filter(wide, metabolite != m); test <- filter(wide, metabolite == m)
  fit <- lmer(qlogis(pmin(pmax(adsorption, 1e-4), 1-1e-4)) ~ charge_class + source + (1|metabolite), data=train)
  pred <- plogis(predict(fit, newdata=test, allow.new.levels=TRUE))
  data.frame(metabolite=m, source=test$source, observed=test$adsorption, predicted=pred)
}) %>% bind_rows()
cv_summary <- data.frame(n_metabolites=length(metabs), n_rows=nrow(cv), spearman_rho=cor(cv$observed, cv$predicted, method="spearman"), pearson_r=cor(cv$observed, cv$predicted), rmse=sqrt(mean((cv$observed-cv$predicted)^2)))
write_csv(cv, file.path(out, "corrected_leave_one_metabolite_out_predictions.csv"))
write_csv(cv_summary, file.path(out, "corrected_leave_one_metabolite_out_summary.csv"))

pred_plot <- ggplot(cv, aes(predicted, observed, color=source)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey40") + geom_point(size=2.2, alpha=.8) +
  facet_wrap(~source, nrow=1) + scale_x_continuous(limits=c(0,1), labels=scales::label_percent()) + scale_y_continuous(limits=c(0,1), labels=scales::label_percent()) +
  labs(title="Descriptor-based adsorption predictions partially recover corrected observations", subtitle=sprintf("Leave-one-metabolite-out validation; Spearman ρ = %.2f; RMSE = %.2f", cv_summary$spearman_rho, cv_summary$rmse), x="Predicted adsorption", y="Observed adsorption", color="Source") +
  theme_classic(base_size=13) + theme(plot.title=element_text(size=15, face="bold"), plot.subtitle=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=13), strip.text=element_text(size=13, face="bold"), legend.position="bottom", legend.text=element_text(size=12), legend.title=element_text(size=12))
ggsave(file.path(out,"corrected_prediction_vs_observation.pdf"), pred_plot, width=10, height=5.8, device=cairo_pdf)
ggsave(file.path(out,"corrected_prediction_vs_observation.png"), pred_plot, width=10, height=5.8, dpi=400, device=ragg::agg_png)
ggsave(file.path(out,"corrected_prediction_vs_observation.svg"), pred_plot, width=10, height=5.8, device=svglite::svglite)

p <- ggplot(wide, aes(charge_class, adsorption, color=charge_class)) +
  geom_jitter(width=.12, height=0, alpha=.45, size=1.7) +
  stat_summary(fun=mean, geom="point", color="black", size=3) +
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", color="black", width=.18) +
  facet_wrap(~source, nrow=1) + scale_y_continuous(limits=c(0,1), labels=scales::label_percent()) +
  labs(title="Corrected adsorption is associated with metabolite charge class", subtitle="Adsorption calculated from matched corrected peak heights: (NA - sediment) / NA", x="Operational charge class", y="Estimated adsorption") +
  theme_classic(base_size=13) + theme(plot.title=element_text(size=15, face="bold"), plot.subtitle=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=13), strip.text=element_text(size=13, face="bold"), legend.position="none")
ggsave(file.path(out,"corrected_charge_class_adsorption.pdf"), p, width=11, height=5.8, device=cairo_pdf)
ggsave(file.path(out,"corrected_charge_class_adsorption.png"), p, width=11, height=5.8, dpi=400, device=ragg::agg_png)
ggsave(file.path(out,"corrected_charge_class_adsorption.svg"), p, width=11, height=5.8, device=svglite::svglite)

cat("Wrote corrected analysis to ", out, "\n", sep="")
print(audit)
print(coef)

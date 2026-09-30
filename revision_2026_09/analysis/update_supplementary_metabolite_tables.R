#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(stringr); library(openxlsx)})

source_xlsx <- "/Users/mingfeichen/Downloads/New_metabolites_corrected.xlsx"
template_xlsx <- "/Users/mingfeichen/Downloads/Necromass Supplementary Table (1).xlsx"
out_xlsx <- "/Users/mingfeichen/Manuscript/outputs/manual-20260818-a1/Necromass_Supplementary_Table_updated.xlsx"
dir.create(dirname(out_xlsx), recursive=TRUE, showWarnings=FALSE)

pretty <- function(x) str_squish(as.character(x))
safe_t <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  tryCatch(t.test(x, y, var.equal=FALSE)$p.value, error=function(e) NA_real_)
}
safe_tstat <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  tryCatch(unname(t.test(x, y, var.equal=FALSE)$statistic), error=function(e) NA_real_)
}
mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm=TRUE)

raw <- read_excel(source_xlsx, sheet="Supplementary Table X Necromas ")
names(raw)[1] <- "feature"
raw$feature <- pretty(raw$feature)
sample_cols <- names(raw)[-1]

meta <- tibble(sample=sample_cols) %>% mutate(
  species=case_when(str_detect(sample,"Bacillus") ~ "B", str_detect(sample,"Rhodano") ~ "R", TRUE ~ NA_character_),
  treatment=case_when(
    str_detect(sample,"0hr") ~ "0hr",
    str_detect(sample,"24hr-1mMAlCl3|48hr-1mMAlCl3") ~ "Al",
    str_detect(sample,"phage") ~ "Phage",
    str_detect(sample,"Kana") ~ "Kanamycin",
    str_detect(sample,"8hr") & !str_detect(sample,"48hr") ~ "8hr",
    str_detect(sample,"48hr") ~ "48hr",
    str_detect(sample,"24hr") ~ "24hr",
    TRUE ~ NA_character_),
  condition=case_when(
    treatment %in% c("0hr","8hr","24hr","48hr") ~ treatment,
    treatment=="Al" ~ "Al",
    treatment=="Phage" ~ "Phage",
    treatment=="Kanamycin" ~ "Kanamycin",
    TRUE ~ NA_character_))

long <- raw %>% pivot_longer(-feature, names_to="sample", values_to="intensity") %>% left_join(meta, by="sample")
long$intensity <- as.numeric(long$intensity)

superclass <- function(x) {
  z <- str_to_lower(str_replace_all(x,"[^a-z0-9]", ""))
  case_when(
    str_detect(z,"aden|guan|cytos|cytid|urac|urid|xanth|hypox|inos|deoxy|orotic|thym") ~ "Nucleic acids (bases/nucleosides/nucleotides)",
    str_detect(z,"alan|argin|aspart|glutam|glycin|histid|isoleuc|leucin|lysine|methion|phenylalan|prolin|serin|threon|tryptoph|tyros|valin|amino") ~ "Amino acids & derivatives",
    str_detect(z,"lactic|succin|malic|fumar|citric|pyruv|shikim|benzo|hydroxy|glycer|phosph|organic") ~ "Central carbon & organic acids",
    str_detect(z,"pyridox|pterin|nicotin|riboflav|thiamin|vitamin|folat") ~ "Cofactors (pterins & vitamins)",
    str_detect(z,"betaine|carnitin|cholin|guanid|trimethyl|polyamin|spermid|putresc") ~ "Lipids & osmolytes (quaternary amines)",
    str_detect(z,"glucose|fructose|sugar|glycerol|mannitol|sorbitol") ~ "Carbohydrates & amino sugars",
    TRUE ~ "Other / synthetic")
}

feature_info <- long %>% distinct(feature) %>% mutate(superclass=vapply(feature, superclass, character(1)))

within <- long %>% filter(!is.na(species), !is.na(condition)) %>%
  mutate(contrast=case_when(
    species=="B" & condition=="8hr" ~ "B_mid",
    species=="B" & condition=="24hr" ~ "B_late",
    species=="R" & condition=="24hr" ~ "R_mid",
    species=="R" & condition=="48hr" ~ "R_late",
    species=="B" & condition=="Al" ~ "B_Al",
    species=="R" & condition=="Al" ~ "R_Al",
    species=="B" & condition=="Phage" ~ "B_P",
    species=="R" & condition=="Phage" ~ "R_P",
    species=="B" & condition=="Kanamycin" ~ "B_K",
    species=="R" & condition=="Kanamycin" ~ "R_K",
    TRUE ~ NA_character_)) %>% filter(!is.na(contrast))

base_for <- function(sp) filter(long, species==sp, condition=="0hr")
calc_within <- function(f, sp, tr, label) {
  b <- long %>% filter(feature==f, species==sp, condition=="0hr") %>% pull(intensity)
  y <- long %>% filter(feature==f, species==sp, condition==tr) %>% pull(intensity)
  mb <- mean_or_na(b); my <- mean_or_na(y); pc <- 1
  tibble(feature=f, phylogeny=sp, treatment=label, mean_0hr=mb, mean_treat=my,
         log2FC=log2((my+pc)/(mb+pc)), tstat=safe_tstat(y,b), pvalue=safe_t(y,b), n_0hr=sum(is.finite(b)), n_treat=sum(is.finite(y)), pseudocount=pc)
}
tr_map <- tribble(~sp,~tr,~label,"B","8hr","mid","B","24hr","late","B","Al","Al","B","Phage","Phage","B","Kanamycin","Kanamycin","R","24hr","mid","R","48hr","late","R","Al","Al","R","Phage","Phage","R","Kanamycin","Kanamycin")
s3 <- crossing(feature=unique(long$feature), tr_map) %>% rowwise() %>% do(calc_within(.$feature,.$sp,.$tr,.$label)) %>% ungroup() %>% mutate(ttest_status=case_when(is.na(pvalue)~"not_tested", pvalue<0.05 & log2FC>0~"up", pvalue<0.05 & log2FC<0~"down", TRUE~"not_significant")) %>% left_join(feature_info,by="feature") %>% select(feature,phylogeny,treatment,mean_0hr,mean_treat,log2FC,tstat,pvalue,n_0hr,n_treat,pseudocount,ttest_status,superclass)

pair_map <- tribble(~label,~sp,~tr,"0h","B","0hr","mid","B","8hr","late","B","24hr","Al","B","Al","K","B","Kanamycin","P","B","Phage")
calc_pair <- function(f,label) {
  z <- pair_map %>% filter(label==!!label)
  a <- long %>% filter(feature==f, species=="B", condition==z$tr) %>% pull(intensity)
  rtr <- case_when(label=="0h" ~ "0hr", label=="mid" ~ "24hr", label=="late" ~ "48hr", label=="Al" ~ "Al", label=="K" ~ "Kanamycin", TRUE ~ "Phage")
  b <- long %>% filter(feature==f, species=="R", condition==rtr) %>% pull(intensity)
  ma<-mean_or_na(a); mb<-mean_or_na(b); pc<-1
  tibble(contrast=paste0("B_",label," vs R_",label),A=paste0("B_",label),B=paste0("R_",label),feature=f,mean_A=ma,mean_B=mb,log2FC=log2((ma+pc)/(mb+pc)),tstat=safe_tstat(a,b),pvalue=safe_t(a,b),n_A=sum(is.finite(a)),n_B=sum(is.finite(b)),pseudocount=pc)
}
s4 <- crossing(feature=unique(long$feature), label=pair_map$label) %>% rowwise() %>% do(calc_pair(.$feature,.$label)) %>% ungroup() %>% mutate(padj=p.adjust(pvalue,"BH"),direction=case_when(padj<0.05 & log2FC>0~"Bacillus enriched",padj<0.05 & log2FC<0~"Rhodanobacter enriched",TRUE~"not significant"),superclass=vapply(feature, superclass, character(1))) %>% select(contrast,A,B,feature,mean_A,mean_B,log2FC,tstat,pvalue,n_A,n_B,pseudocount,padj,direction,superclass)

mineral <- read_excel(source_xlsx, sheet="Supplememntary Table Y Mineral"); names(mineral)[1] <- "Metabolites"; mineral$Metabolites <- pretty(mineral$Metabolites)
min_long <- mineral %>% pivot_longer(-Metabolites,names_to="sample",values_to="peak") %>% mutate(species=if_else(str_detect(sample,"Bacillus"),"Bacillus","Rhodanobacter"),matrix=if_else(str_detect(sample,"-NA_"),"NA","Sed"))
s9 <- min_long %>% group_by(Metabolites,species,matrix) %>% summarise(mean_peak=mean(peak,na.rm=TRUE),.groups="drop") %>% pivot_wider(names_from=matrix,values_from=mean_peak) %>% mutate(adsorption=(`NA`-Sed)/`NA`) %>% select(Metabolites,species,adsorption) %>% pivot_wider(names_from=species,values_from=adsorption) %>% transmute(Metabolites,Bacillus_absorb=Bacillus,Rhoano_absorb=Rhodanobacter)

wb <- loadWorkbook(template_xlsx)
replace_sheet <- function(wb, name, title, dat) {
  sh <- wb[[name]]; deleteData(wb, sheet=name, cols=1:30, rows=1:max(1000, nrow(dat)+5), gridExpand=TRUE)
  writeData(wb, name, title, startRow=1, startCol=1); writeData(wb, name, dat, startRow=2, startCol=1, colNames=TRUE)
  addStyle(wb, name, createStyle(fontName="Arial", fontSize=10), rows=1:(nrow(dat)+2), cols=1:ncol(dat), gridExpand=TRUE, stack=TRUE)
  addStyle(wb, name, createStyle(textDecoration="bold", fontName="Arial"), rows=1:2, cols=1:ncol(dat), gridExpand=TRUE, stack=TRUE)
}
replace_sheet(wb,"Supplementary Table S3","Supplementary Table S3. Corrected targeted metabolite differential abundance results for within-species contrasts.",s3)
replace_sheet(wb,"Supplementary Table S4","Supplementary Table S4. Corrected targeted metabolite results for matched Bacillus versus Rhodanobacter contrasts.",s4)
replace_sheet(wb,"Supplementary Table S9","Supplementary Table S9. Corrected targeted metabolite partitioning after sediment exposure.",s9)
saveWorkbook(wb,out_xlsx,overwrite=TRUE)
write.csv(s3,sub(".xlsx$","_S3.csv",out_xlsx),row.names=FALSE); write.csv(s4,sub(".xlsx$","_S4.csv",out_xlsx),row.names=FALSE); write.csv(s9,sub(".xlsx$","_S9.csv",out_xlsx),row.names=FALSE)
cat("Wrote ",out_xlsx,"\n",sep="")

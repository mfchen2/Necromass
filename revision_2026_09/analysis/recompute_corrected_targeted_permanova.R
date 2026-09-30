#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(stringr); library(tibble); library(vegan); library(openxlsx)})
input <- "/Users/mingfeichen/Downloads/New_metabolites_corrected.xlsx"
out_xlsx <- "/Users/mingfeichen/Manuscript/outputs/manual-20260818-a1/Necromass_Supplementary_Table_updated.xlsx"
out_dir <- "/Users/mingfeichen/Manuscript/outputs/manual-20260818-a1"
d <- read_excel(input, sheet="Supplementary Table X Necromas "); names(d)[1] <- "feature"
mat <- d %>% select(-feature) %>% as.data.frame(); rownames(mat) <- d$feature; mat <- as.matrix(mat); mat[is.na(mat)] <- 0; mat <- t(mat)
meta <- tibble(sample=rownames(mat)) %>% mutate(
  BR=if_else(str_detect(sample,"Bacillus"),"Bacillus","Rhodanobacter"),
  Exp_group=case_when(str_detect(sample,"0hr")~"0h", str_detect(sample,"24hr-1mMAlCl3|48hr-1mMAlCl3")~"Al", str_detect(sample,"phage")~"P", str_detect(sample,"Kana")~"K", BR=="Bacillus" & str_detect(sample,"8hr")~"mid", BR=="Bacillus" & str_detect(sample,"24hr")~"late", BR=="Rhodanobacter" & str_detect(sample,"24hr")~"mid", BR=="Rhodanobacter" & str_detect(sample,"48hr")~"late", TRUE~NA_character_),
  BR=factor(BR,levels=c("Bacillus","Rhodanobacter")), Exp_group=factor(Exp_group,levels=c("0h","mid","late","Al","K","P")))
stopifnot(!anyNA(meta$Exp_group), nrow(meta)==36)
dist <- vegdist(mat, method="gower")
fit <- adonis2(dist ~ BR + Exp_group, data=meta, permutations=999, by="terms")
permanova <- as.data.frame(fit) %>% rownames_to_column("term") %>% filter(term %in% c("BR","Exp_group","Residual","Total")) %>% rename(p_value=`Pr(>F)`)
write.csv(permanova,file.path(out_dir,"corrected_targeted_permanova.csv"),row.names=FALSE)
pcoa <- cmdscale(dist,k=2,eig=TRUE)
xy <- as.data.frame(pcoa$points) %>% rownames_to_column("sample") %>% rename(PCoA1=V1,PCoA2=V2) %>% left_join(meta,by="sample")
cent <- xy %>% group_by(Exp_group,BR) %>% summarise(cx=mean(PCoA1),cy=mean(PCoA2),.groups="drop")
overall <- xy %>% group_by(Exp_group) %>% summarise(cx=mean(PCoA1),cy=mean(PCoA2),.groups="drop")
base <- overall %>% filter(Exp_group=="0h")
dist_summary <- bind_rows(lapply(levels(meta$Exp_group), function(g) { z <- cent %>% filter(as.character(Exp_group)==g); q <- overall %>% filter(as.character(Exp_group)==g); tibble(Exp_group=g, Exp=as.numeric(sqrt((q$cx-base$cx)^2+(q$cy-base$cy)^2)), BvsR=as.numeric(sqrt((z$cx[z$BR=="Bacillus"]-z$cx[z$BR=="Rhodanobacter"])^2+(z$cy[z$BR=="Bacillus"]-z$cy[z$BR=="Rhodanobacter"])^2))) }))
write.csv(xy,file.path(out_dir,"corrected_targeted_pcoa_coordinates.csv"),row.names=FALSE); write.csv(dist_summary,file.path(out_dir,"corrected_targeted_centroid_distances.csv"),row.names=FALSE)
wb <- loadWorkbook(out_xlsx)
deleteData(wb,"Supplementary Table S1",cols=1:10,rows=1:40,gridExpand=TRUE); writeData(wb,"Supplementary Table S1","Supplementary Table S1. PERMANOVA results for the corrected targeted metabolite profiles.",startRow=1); writeData(wb,"Supplementary Table S1",permanova,startRow=2,colNames=TRUE)
deleteData(wb,"Supplementary Table S2",cols=1:6,rows=1:30,gridExpand=TRUE); writeData(wb,"Supplementary Table S2","Supplementary Table S2. Centroid-distance summary for corrected targeted metabolite profiles.",startRow=1); writeData(wb,"Supplementary Table S2",dist_summary,startRow=2,colNames=TRUE)
saveWorkbook(wb,out_xlsx,overwrite=TRUE); print(permanova); print(dist_summary)

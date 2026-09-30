#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readxl); library(readr); library(openxlsx)})
template <- "/Users/mingfeichen/Downloads/Necromass Supplementary Table (1).xlsx"
out <- "/Users/mingfeichen/Manuscript/outputs/manual-20260818-a1/Necromass_Supplementary_Table_updated.xlsx"
wb <- loadWorkbook(out)
orig1 <- read_excel(template, sheet="Supplementary Table S1", col_names=FALSE)
orig2 <- read_excel(template, sheet="Supplementary Table S2", col_names=FALSE)
perm <- read_csv("/Users/mingfeichen/Manuscript/outputs/manual-20260818-a1/corrected_targeted_permanova.csv", show_col_types=FALSE)
cent <- read_csv("/Users/mingfeichen/Manuscript/outputs/manual-20260818-a1/corrected_targeted_centroid_distances.csv", show_col_types=FALSE)
deleteData(wb,"Supplementary Table S1",cols=1:10,rows=1:40,gridExpand=TRUE)
writeData(wb,"Supplementary Table S1",as.data.frame(orig1[1:10,]),startRow=1,colNames=FALSE)
writeData(wb,"Supplementary Table S1","Targeted (corrected)",startRow=11,colNames=FALSE)
writeData(wb,"Supplementary Table S1",perm,startRow=12,colNames=TRUE)
deleteData(wb,"Supplementary Table S2",cols=1:6,rows=1:30,gridExpand=TRUE)
writeData(wb,"Supplementary Table S2",as.data.frame(orig2[1:10,]),startRow=1,colNames=FALSE)
writeData(wb,"Supplementary Table S2","Targeted (corrected)",startRow=11,colNames=FALSE)
writeData(wb,"Supplementary Table S2",cent,startRow=12,colNames=TRUE)
saveWorkbook(wb,out,overwrite=TRUE)
cat("Restored untargeted S1/S2 blocks and retained corrected targeted blocks.\n")

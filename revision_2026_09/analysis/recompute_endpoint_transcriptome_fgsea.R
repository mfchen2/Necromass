#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

out_dir <- Sys.getenv(
  "TRANSCRIPTOME_FGSEA_OUTPUT_DIR",
  unset = "/Users/mingfeichen/Manuscript/outputs/manual-20260921-a1/transcriptome_fgsea_endpoint_reference"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_bacillus_counts <- function() {
  x <- read.csv("/Users/mingfeichen/Bacillus_RNAseq/counts/gene_counts.csv",
                check.names = FALSE, stringsAsFactors = FALSE)
  ids <- x$Geneid
  mat <- as.matrix(x[, setdiff(names(x), "Geneid"), drop = FALSE])
  colnames(mat) <- colnames(mat) |>
    str_replace_all("_0hr_", "_0h_") |>
    str_replace_all("_8hr_", "_8h_") |>
    str_replace_all("_24hr_", "_24h_")
  rownames(mat) <- ids
  storage.mode(mat) <- "integer"
  meta <- read.csv("/Users/mingfeichen/Bacillus_RNAseq/data/meta/sample_info.csv",
                   check.names = FALSE, stringsAsFactors = FALSE)
  list(counts = mat, meta = meta,
       gene_sets = readRDS("/Users/mingfeichen/kegg_pathway_gene_sets_counts_ids.rds"),
       organism = "Bacillus", phylogeny = "B")
}

read_rhodano_counts <- function() {
  x <- read.delim("/Users/mingfeichen/Rhodano_mRNAseq/counts/gene_counts_prokka.txt",
                  comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  sample_cols <- if ("Length" %in% names(x)) names(x)[7:ncol(x)] else names(x)[2:ncol(x)]
  mat <- as.matrix(x[, sample_cols, drop = FALSE])
  colnames(mat) <- basename(colnames(mat))
  colnames(mat) <- str_remove(colnames(mat), "[.]sorted[.]bam$")
  colnames(mat) <- str_remove(colnames(mat), "[.]bam$")
  rownames(mat) <- x$Geneid
  storage.mode(mat) <- "integer"
  meta <- read.csv("/Users/mingfeichen/Rhodano_mRNAseq/data/meta/sample_info.csv",
                   check.names = FALSE, stringsAsFactors = FALSE)
  list(counts = mat, meta = meta,
       gene_sets = readRDS("/Users/mingfeichen/FW104_kegg_pathway_gene_sets_counts_ids.rds"),
       organism = "Rhodanobacter", phylogeny = "R")
}

run_species <- function(spec) {
  meta <- spec$meta[match(colnames(spec$counts), spec$meta$sample), , drop = FALSE]
  stopifnot(all(meta$sample == colnames(spec$counts)))
  rownames(meta) <- meta$sample
  meta$condition <- factor(meta$condition)

  dds <- DESeqDataSetFromMatrix(
    countData = spec$counts,
    colData = meta,
    design = ~ condition
  )
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)

  contrasts <- if (spec$organism == "Bacillus") {
    c("8h", "24h", "Al", "Kanamycin", "P")
  } else {
    c("24h", "48h", "Al", "Kanamycin", "P")
  }
  refs <- ifelse(contrasts %in% c("Al", "Kanamycin", "P"),
                 ifelse(spec$organism == "Bacillus", "24h", "48h"), "0h")

  all_fgsea <- list()
  all_de <- list()
  for (i in seq_along(contrasts)) {
    ct <- contrasts[[i]]
    ref <- refs[[i]]
    label_ct <- if (ct == "Kanamycin") "K" else if (ct == "P") "P" else ct
    contrast_name <- paste0(label_ct, "_vs_", ref)
    res <- results(dds, contrast = c("condition", ct, ref))
    de <- as.data.frame(res) %>%
      tibble::rownames_to_column("GeneID") %>%
      mutate(organism = spec$organism, phylogeny = spec$phylogeny,
             contrast = contrast_name, reference = ref)
    write_csv(de, file.path(out_dir, paste0(spec$phylogeny, "_DE_", contrast_name, ".csv")))
    all_de[[contrast_name]] <- de

    ranks <- de$stat
    names(ranks) <- de$GeneID
    ranks <- ranks[is.finite(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    sets <- lapply(spec$gene_sets, intersect, y = names(ranks))
    sets <- sets[lengths(sets) >= 5]
    fg <- fgsea(pathways = sets, stats = ranks, minSize = 5, maxSize = 500)
    fg <- as.data.frame(fg)
    names(fg)[names(fg) == "pathway"] <- "set_id"
    fg <- fg %>%
      transmute(contrast = contrast_name, set_id, padj = padj, NES = NES,
                pval = pval, ES = ES, size = size,
                leadingEdge = vapply(leadingEdge, paste, collapse = "|", FUN.VALUE = character(1)),
                organism = spec$organism, phylogeny = spec$phylogeny, reference = ref)
    write_csv(fg, file.path(out_dir, paste0(spec$phylogeny, "_FGSEA_", contrast_name, ".csv")))
    all_fgsea[[contrast_name]] <- fg
  }
  bind_rows(all_fgsea)
}

bacillus <- run_species(read_bacillus_counts())
rhodano <- run_species(read_rhodano_counts())
combined <- bind_rows(bacillus, rhodano)
write_csv(combined, file.path(out_dir, "transcriptome_fgsea_endpoint_reference_all.csv"))
writeLines(c(
  "Transcriptome FGSEA recomputation with endpoint-matched references",
  paste0("Output directory: ", normalizePath(out_dir, winslash = "/")),
  "Bacillus: 8h and 24h versus 0h; Al, K and P versus 24h.",
  "Rhodanobacter: 24h and 48h versus 0h; Al, K and P versus 48h.",
  "Gene-level statistics were generated with DESeq2 from the archived count matrices.",
  "Pathway sets were the archived organism-specific KEGG pathway RDS files."
), file.path(out_dir, "README.txt"))

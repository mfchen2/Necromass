#!/usr/bin/env bash
set -euo pipefail

REF_DIR="data/ref"
BAM_DIR="bam_refseq"
COUNTS_DIR="counts"
GFF_FILE="${REF_DIR}/FW104_10B01.gff"   # or genes.gtf

mkdir -p "${COUNTS_DIR}"

# Adjust -s according to strandedness: 0 (unstranded), 1 (forward), 2 (reverse)
STRAND=2

featureCounts \
  -T 8 \
  -p \
  -s "${STRAND}" \
  -t gene \
  -g locus_tag \
  -a "${GFF_FILE}" \
  -o "${COUNTS_DIR}/gene_counts_refseq.txt" \
  "${BAM_DIR}"/*.sorted.bam

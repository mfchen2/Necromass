#!/usr/bin/env bash
set -euo pipefail

RAW_DIR="data/raw"      # contains subfolders R_0h_1, R_0h_2, ...
TRIM_DIR="trimmed"
QC_RAW_DIR="qc/raw_qc"
QC_TRIM_DIR="qc/trimmed_qc"

mkdir -p "${TRIM_DIR}" "${QC_RAW_DIR}" "${QC_TRIM_DIR}"

# Clean old lists if rerunning
rm -f samples.txt raw_R1_files.txt

######################################
# 1) Find all R1 files recursively
#    pattern: *_1.fq.gz (e.g. R_0h_1_1.fq.gz)
######################################
find "${RAW_DIR}" -type f -name '*_1.fq.gz' | sort > raw_R1_files.txt

if [ ! -s raw_R1_files.txt ]; then
  echo "ERROR: No R1 FASTQ files matching '*_1.fq.gz' found under ${RAW_DIR}" >&2
  exit 1
fi

# Start a fresh samples file
: > samples.txt

######################################
# 2) Raw QC with FastQC
######################################
while read -r R1; do
  # R1 full path, e.g. data/raw/R_0h_1/R_0h_1_1.fq.gz
  R2="${R1/_1.fq.gz/_2.fq.gz}"

  if [ ! -f "${R2}" ]; then
    echo "WARNING: R2 not found for ${R1} (expected ${R2}) – skipping." >&2
    continue
  fi

  # SAMPLE is basename without the trailing _1.fq.gz
  # e.g. R_0h_1_1.fq.gz -> R_0h_1
  SAMPLE="$(basename "${R1}" _1.fq.gz)"

  echo "${SAMPLE}" >> samples.txt

  fastqc -o "${QC_RAW_DIR}" "${R1}" "${R2}"
done < raw_R1_files.txt

# Deduplicate sample names (just in case)
sort -u samples.txt -o samples.txt

######################################
# 3) Trimming with fastp
######################################
while read -r R1; do
  R2="${R1/_1.fq.gz/_2.fq.gz}"
  SAMPLE="$(basename "${R1}" _1.fq.gz)"

  if [ ! -f "${R2}" ]; then
    echo "WARNING: R2 not found for ${R1} (expected ${R2}) – skipping trimming." >&2
    continue
  fi

  fastp \
    -i "${R1}" \
    -I "${R2}" \
    -o "${TRIM_DIR}/${SAMPLE}_trimmed_1.fq.gz" \
    -O "${TRIM_DIR}/${SAMPLE}_trimmed_2.fq.gz" \
    -h "${QC_TRIM_DIR}/${SAMPLE}_fastp.html" \
    -j "${QC_TRIM_DIR}/${SAMPLE}_fastp.json"

done < raw_R1_files.txt

######################################
# 4) FastQC on trimmed reads
######################################
while read -r SAMPLE; do
  fastqc -o "${QC_TRIM_DIR}" \
    "${TRIM_DIR}/${SAMPLE}_trimmed_1.fq.gz" \
    "${TRIM_DIR}/${SAMPLE}_trimmed_2.fq.gz"
done < samples.txt

######################################
# 5) MultiQC
######################################
multiqc "${QC_RAW_DIR}" "${QC_TRIM_DIR}" -o qc

#!/usr/bin/env bash
set -euo pipefail

REF_DIR="data/ref"
TRIM_DIR="trimmed"
BAM_DIR="bam_refseq"

GENOME_FA="${REF_DIR}/FW104_10B01_refseq.fna"
BT2_PREFIX="${REF_DIR}/genome_bt2"

mkdir -p "${BAM_DIR}"

######################################
# 1) Build Bowtie2 index (once)
######################################
if [ ! -e "${BT2_PREFIX}.1.bt2" ]; then
  echo "[INFO] Building Bowtie2 index for ${GENOME_FA} ..."
  bowtie2-build "${GENOME_FA}" "${BT2_PREFIX}"
else
  echo "[INFO] Bowtie2 index already exists at prefix ${BT2_PREFIX}"
fi

######################################
# 2) Align each sample, but SKIP finished ones
######################################
while read -r SAMPLE; do
  R1="${TRIM_DIR}/${SAMPLE}_trimmed_1.fq.gz"
  R2="${TRIM_DIR}/${SAMPLE}_trimmed_2.fq.gz"

  SAM_OUT="${BAM_DIR}/${SAMPLE}.sam"
  BAM_SORTED="${BAM_DIR}/${SAMPLE}.sorted.bam"
  BAM_INDEX="${BAM_SORTED}.bai"

  # --- Skip logic: if sorted BAM AND index exist, we assume sample is done
  if [ -f "${BAM_SORTED}" ] && [ -f "${BAM_INDEX}" ]; then
    echo "[SKIP] ${SAMPLE}: found ${BAM_SORTED} and ${BAM_INDEX}, assuming complete."
    continue
  fi

  # Basic sanity check
  if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "[WARN] Trimmed FASTQs missing for ${SAMPLE} (${R1} or ${R2}). Skipping." >&2
    continue
  fi

  echo "[RUN ] ${SAMPLE}: aligning and creating sorted BAM ..."

  # Clean up any old partial SAM from previous failed run
  rm -f "${SAM_OUT}"

  # Align
  bowtie2 \
    -x "${BT2_PREFIX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    -p 8 \
    -S "${SAM_OUT}"

  # Convert + sort
  samtools view -bS "${SAM_OUT}" \
    | samtools sort -o "${BAM_SORTED}" -

  # Index
  samtools index "${BAM_SORTED}"

  # Remove SAM to save space
  rm -f "${SAM_OUT}"

  echo "[DONE] ${SAMPLE}"
done < samples.txt

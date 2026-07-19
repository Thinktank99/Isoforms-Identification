#!/usr/bin/env bash
#SBATCH -J star_cdna_test
#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

set -euo pipefail

source ../scripts/env_conda.sh
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

REF=../ref
STAR_DIR="${REF}/star_index"

sample="ERR9871746"
CDNA_FASTQ="../results/06_align_count/${sample}/${sample}.cdna.fastq.gz"
OUT_DIR="../results/07_mapping_test_cdna"
mkdir -p "$OUT_DIR"

echo "Sample: $sample"
ls -lh "$CDNA_FASTQ"

STAR \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --genomeDir "$STAR_DIR" \
  --readFilesIn "$CDNA_FASTQ" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_DIR}/${sample}_" \
  --outSAMtype None

echo "==== Mapping results ===="
cat "${OUT_DIR}/${sample}_Log.final.out"

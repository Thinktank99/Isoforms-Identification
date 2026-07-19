#!/usr/bin/env bash
#SBATCH -J star_test_all
#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o ../logs/%x_%A_%a.out
#SBATCH -e ../logs/%x_%A_%a.err
#SBATCH --array=1-4

set -euo pipefail

source ../scripts/env_conda.sh

unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

REF=../ref
STAR_DIR="${REF}/star_index"

OUT_DIR=../results/07A_mapping_test
mkdir -p "$OUT_DIR"

SAMPLES=../scripts/samples.tsv

line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES")
sample=$(echo "$line" | cut -f1)
r1=$(echo "$line" | cut -f2)
r2=$(echo "$line" | cut -f3)

echo "Sample: $sample"
echo "R1: $r1"
echo "R2: $r2"
ls -lh "$r1" "$r2"

STAR \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --genomeDir "$STAR_DIR" \
  --readFilesIn "$r1" "$r2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_DIR}/${sample}_" \
  --outSAMtype None

echo "DONE: $sample"
echo "Log file:"
ls -lh "${OUT_DIR}/${sample}_Log.final.out"

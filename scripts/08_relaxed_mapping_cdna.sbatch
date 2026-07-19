#!/usr/bin/env bash
#SBATCH -J star_relaxed
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

SAMPLES=../scripts/samples.tsv
OUT_DIR=../results/08_relaxed_mapping
mkdir -p "$OUT_DIR"

line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES")
sample=$(echo "$line" | cut -f1)

CDNA_FASTQ="../results/06_align_count/${sample}/${sample}.cdna.fastq.gz"

echo "Sample: $sample"
echo "cDNA FASTQ: $CDNA_FASTQ"
ls -lh "$CDNA_FASTQ"

STAR \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --genomeDir "$STAR_DIR" \
  --readFilesIn "$CDNA_FASTQ" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_DIR}/${sample}_" \
  --outSAMtype None \
  --outFilterMultimapNmax 50 \
  --winAnchorMultimapNmax 100 \
  --outFilterMatchNmin 25 \
  --outFilterMatchNminOverLread 0.20 \
  --outFilterScoreMinOverLread 0.20 \
  --outFilterMismatchNoverLmax 0.10 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000

echo "==== Relaxed mapping results for $sample ===="
cat "${OUT_DIR}/${sample}_Log.final.out"

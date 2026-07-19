#!/usr/bin/env bash
#SBATCH -J fastqtosam
#SBATCH -p medium
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 02:30:00
#SBATCH -o ../logs/%x_%A_%a.out
#SBATCH -e ../logs/%x_%A_%a.err
#SBATCH --array=1-4

set -euo pipefail
source ../scripts/env_conda.sh

echo "Node: $(hostname)"
echo "Env:  ${CONDA_DEFAULT_ENV:-NA}"

SAMPLES=../scripts/samples.tsv
TRIM=../results/01_trim
OUT=../results/02_picard
mkdir -p "$OUT"

# Big temp dir (avoid /tmp on compute node)
export TMPDIR="$(readlink -f ../tmp)"
mkdir -p "$TMPDIR"
echo "TMPDIR=$TMPDIR"
df -h "$TMPDIR" || true

line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES")
sample=$(echo "$line" | cut -f1)

r1="${TRIM}/${sample}_R1.trim.fastq.gz"
r2="${TRIM}/${sample}_R2.trim.fastq.gz"
ubam="${OUT}/${sample}.unmapped.bam"

ls -lh "$r1" "$r2"

if command -v picard >/dev/null 2>&1; then
  PICARD_CMD="picard"
elif [ -f "$CONDA_PREFIX"/share/picard-*/picard.jar ]; then
  PICARD_JAR_PATH=$(ls -1 "$CONDA_PREFIX"/share/picard-*/picard.jar | head -n 1)
  PICARD_CMD="java -Xmx20g -Djava.io.tmpdir=$TMPDIR -jar $PICARD_JAR_PATH"
else
  echo "ERROR: Picard not found."
  exit 1
fi

echo "Using Picard as: $PICARD_CMD"

$PICARD_CMD FastqToSam \
  F1="$r1" \
  F2="$r2" \
  O="$ubam" \
  SM="$sample" \
  RG="RG_${sample}" \
  PL="ILLUMINA" \
  TMP_DIR="$TMPDIR"

ls -lh "$ubam"

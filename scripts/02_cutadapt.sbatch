#!/usr/bin/env bash
#SBATCH -J cutadapt
#SBATCH -p medium
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -o ../logs/%x_%A_%a.out
#SBATCH -e ../logs/%x_%A_%a.err
#SBATCH --array=1-4

set -euo pipefail
source ../scripts/env_conda.sh

echo "Node: $(hostname)"
echo "Env:  ${CONDA_DEFAULT_ENV:-NA}"
which cutadapt
cutadapt --version

SAMPLES=../scripts/samples.tsv
OUT_DIR=../results/01_trim
mkdir -p "$OUT_DIR"

line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES")
sample=$(echo "$line" | cut -f1)
r1=$(echo "$line" | cut -f2)
r2=$(echo "$line" | cut -f3)

out_r1="${OUT_DIR}/${sample}_R1.trim.fastq.gz"
out_r2="${OUT_DIR}/${sample}_R2.trim.fastq.gz"

# Common Illumina adapter candidates (R1)
A1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"   # TruSeq / universal
A2="CTGTCTCTTATACACATCT"                 # Nextera transposase (short)
A3="AGATCGGAAGAG"                         # partial Illumina (very common)

# Common Illumina adapter candidates (R2)
B1="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"   # TruSeq
B2="CTGTCTCTTATACACATCT"                 # Nextera
B3="AGATCGGAAGAG"                         # partial

# poly-A / poly-T tails (sometimes seen)
POLYA="AAAAAAAAAAAAAAAAAAAA"
POLYT="TTTTTTTTTTTTTTTTTTTT"

cutadapt \
  -j ${SLURM_CPUS_PER_TASK} \
  -q 20,20 --minimum-length 25 \
  -a "$A1" -a "$A2" -a "$A3" -a "$POLYA" \
  -A "$B1" -A "$B2" -A "$B3" -A "$POLYT" \
  --times 2 \
  --trim-n \
  -o "$out_r1" -p "$out_r2" \
  "$r1" "$r2"

echo "Done: $sample"
ls -lh "$out_r1" "$out_r2"

#!/usr/bin/env bash
#SBATCH -J fastqc
#SBATCH -p medium          # e.g. standard | medium | long | short (use your cluster's actual name)
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 01:30:00
#SBATCH -o ../logs/%x_%A_%a.out
#SBATCH -e ../logs/%x_%A_%a.err
#SBATCH --array=1-4

set -euo pipefail

# 1) Activate your environment inside the job
source ../scripts/env_conda.sh

# 2) Basic sanity check (prints paths into log)
echo "Node: $(hostname)"
echo "Env:  ${CONDA_DEFAULT_ENV:-NA}"
which fastqc
fastqc --version

# 3) Inputs from samples.tsv
SAMPLES=../scripts/samples.tsv
OUT=../results/00_qc
mkdir -p "$OUT"

# Read the SLURM array line (+1 to skip header)
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES")

sample=$(echo "$line" | cut -f1)
r1=$(echo "$line" | cut -f2)
r2=$(echo "$line" | cut -f3)

echo "Sample: $sample"
echo "R1: $r1"
echo "R2: $r2"

# 4) Run FastQC
fastqc -t ${SLURM_CPUS_PER_TASK} -o "$OUT" "$r1" "$r2"

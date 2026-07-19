#!/usr/bin/env bash
#SBATCH -J qc_sct
#SBATCH -p medium
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -o ../logs/%x_%A_%a.out
#SBATCH -e ../logs/%x_%A_%a.err
#SBATCH --array=1-4

set -euo pipefail

source ../scripts/env_conda.sh
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

INPUT_DIR=../results/07_seurat_raw
OUTPUT_DIR=../results/08_seurat_qc
mkdir -p "$OUTPUT_DIR"

FILES=(${INPUT_DIR}/*.seurat.raw.rds)

input="${FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
sample=$(basename "$input" .seurat.raw.rds)
output="${OUTPUT_DIR}/${sample}.qc_sct.rds"

Rscript ../scripts/08_qc_doubletfinder_sct.R "$input" "$output"

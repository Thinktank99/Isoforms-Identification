#!/usr/bin/env bash
#SBATCH -J build_seurat
#SBATCH -p medium
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

set -euo pipefail

source ../scripts/env_conda.sh
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

INPUT_DIR=../results/06_align_count
OUTPUT_DIR=../results/07_seurat_raw

mkdir -p "$OUTPUT_DIR"

Rscript ../scripts/06_build_seurat_objects.R "$INPUT_DIR" "$OUTPUT_DIR"

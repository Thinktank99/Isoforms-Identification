#!/usr/bin/env bash
#SBATCH -J markers
#SBATCH -p medium
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

set -euo pipefail

source ../scripts/env_conda.sh
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

INPUT_RDS=../results/09_integrated/integrated_seurat.rds
OUTPUT_TSV=../results/10_markers/all_markers.tsv

mkdir -p ../results/10_markers

Rscript ../scripts/10_find_markers.R "$INPUT_RDS" "$OUTPUT_TSV"

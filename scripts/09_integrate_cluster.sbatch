#!/usr/bin/env bash
#SBATCH -J integrate
#SBATCH -p medium
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 18:00:00
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

set -euo pipefail

source ../scripts/env_conda.sh
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

INPUT_DIR=../results/08_seurat_qc
OUTPUT_RDS=../results/09_integrated/integrated_seurat.rds

mkdir -p ../results/09_integrated

Rscript ../scripts/09_integrate_cluster.R "$INPUT_DIR" "$OUTPUT_RDS"

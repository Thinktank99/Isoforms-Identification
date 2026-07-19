#!/usr/bin/env bash
#SBATCH -J buildref
#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

set -euo pipefail
source ../scripts/env_conda.sh

unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

export TMPDIR="$(readlink -f ../tmp)"
mkdir -p "$TMPDIR"

REF=../ref
GENOME="${REF}/genome.fna"
GTF="${REF}/genomic.gtf"

DROPSEQ_DIR="${REF}/dropseq"
STAR_DIR="${REF}/star_index"
mkdir -p "$DROPSEQ_DIR" "$STAR_DIR"

DICT="${REF}/genome.dict"
REFLAT="${DROPSEQ_DIR}/genes.refFlat"
REDUCED="${DROPSEQ_DIR}/genes.reduced.gtf"

echo "GENOME=$GENOME"
echo "GTF=$GTF"
ls -lh "$GENOME" "$GTF"
which picard STAR ConvertToRefFlat ReduceGtf CreateIntervalsFiles
java -version

# 1) dict
if [ -s "$DICT" ]; then
  echo "[SKIP] genome.dict exists"
else
  picard CreateSequenceDictionary R="$GENOME" O="$DICT" TMP_DIR="$TMPDIR"
fi

# 2) refFlat
if [ -s "$REFLAT" ]; then
  echo "[SKIP] genes.refFlat exists"
else
  ConvertToRefFlat ANNOTATIONS_FILE="$GTF" SEQUENCE_DICTIONARY="$DICT" OUTPUT="$REFLAT"
fi

# 3) ReduceGtf (your param style)
if [ -s "$REDUCED" ]; then
  echo "[SKIP] genes.reduced.gtf exists"
else
  ReduceGtf SEQUENCE_DICTIONARY="$DICT" GTF="$GTF" OUTPUT="$REDUCED"
fi

# 4) intervals (OUTPUT is a DIRECTORY)
if [ -s "${DROPSEQ_DIR}/genes.intervals" ]; then
  echo "[SKIP] genes.intervals exists"
else
  CreateIntervalsFiles SEQUENCE_DICTIONARY="$DICT" REDUCED_GTF="$REDUCED" PREFIX=my OUTPUT="$DROPSEQ_DIR"
fi

# 5) STAR index (skip if complete unless forced)
SJDB_OVERHANG=99
need_star=0
for f in Genome SA SAindex genomeParameters.txt; do
  [ -s "${STAR_DIR}/${f}" ] || need_star=1
done
[ "${FORCE_STAR:-0}" -eq 1 ] && need_star=1

if [ "$need_star" -eq 0 ]; then
  echo "[SKIP] STAR index looks complete"
else
  rm -rf "$STAR_DIR"
  mkdir -p "$STAR_DIR"
  STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
    --runMode genomeGenerate \
    --genomeDir "$STAR_DIR" \
    --genomeFastaFiles "$GENOME" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang "$SJDB_OVERHANG" \
    --limitGenomeGenerateRAM 75000000000
fi

echo "DONE:"
ls -lh "$DICT"
ls -lh "$DROPSEQ_DIR" | head -n 50

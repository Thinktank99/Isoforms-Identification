#!/usr/bin/env bash
#SBATCH -J tag_xm_xc
#SBATCH -p medium
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o ../logs/%x_%A_%a.out
#SBATCH -e ../logs/%x_%A_%a.err
#SBATCH --array=1-4

set -uo pipefail

# Activate env inside job
source ../scripts/env_conda.sh

# Avoid cluster-injected JVM options
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

# Big temp (avoid /tmp space issues)
export TMPDIR="$(readlink -f ../tmp)"
mkdir -p "$TMPDIR"
echo "TMPDIR=$TMPDIR"
df -h "$TMPDIR" || true

SAMPLES=../scripts/samples.tsv
RAW_DIR=../data/fastq
OUT_DIR=../results/05_tagged
mkdir -p "$OUT_DIR"

# Get sample for this array task (+1 skips header)
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES")
sample=$(echo "$line" | cut -f1)
r1=$(echo "$line" | cut -f2)
r2=$(echo "$line" | cut -f3)

echo "Sample: $sample"
echo "R1: $r1"
echo "R2: $r2"
ls -lh "$r1" "$r2"

ubam="${OUT_DIR}/${sample}.unmapped.bam"
bam_xm="${OUT_DIR}/${sample}.xm.bam"
bam_xc="${OUT_DIR}/${sample}.xm_xc.bam"
bam_filt="${OUT_DIR}/${sample}.xm_xc.filtered.bam"

########################################
# 1) FASTQ -> uBAM
########################################

if [ ! -s "$ubam" ]; then
  echo "[RUN] picard FastqToSam"
  picard FastqToSam \
    F1="$r1" \
    F2="$r2" \
    O="$ubam" \
    SM="$sample" \
    RG="RG_${sample}" \
    PL="ILLUMINA" \
    TMP_DIR="$TMPDIR"
else
  echo "[SKIP] uBAM exists: $ubam"
fi

########################################
# 2) Tag UMI (XM) from READ2 bases 1-10
########################################

if [ ! -s "$bam_xm" ]; then
  echo "[RUN] Tag UMI (XM) from read2: 1-10"
  TagBamWithReadSequenceExtended \
    INPUT="$ubam" \
    OUTPUT="$bam_xm" \
    SUMMARY="${OUT_DIR}/${sample}.xm.summary.txt" \
    BASE_RANGE="1-10" \
    BASE_QUALITY=10 \
    BARCODED_READ=2 \
    DISCARD_READ=false \
    TAG_NAME=XM \
    NUM_BASES_BELOW_QUALITY=1
else
  echo "[SKIP] XM-tagged BAM exists: $bam_xm"
fi

########################################
# 3) Tag Cell barcode (XC)
# CB1: 87-94
# CB2: 49-56
# CB3: 11-18
########################################

if [ ! -s "$bam_xc" ]; then
  echo "[RUN] Tag Cell barcode (XC) from read2: 87-94,49-56,11-18"
   TagBamWithReadSequenceExtended \
    INPUT="$bam_xm" \
    OUTPUT="$bam_xc" \
    SUMMARY="${OUT_DIR}/${sample}.xc.summary.txt" \
    BASE_RANGE=87-94:49-56:11-18 \
    BASE_QUALITY=10 \
    BARCODED_READ=2 \
    DISCARD_READ=true \
    TAG_NAME=XC \
    NUM_BASES_BELOW_QUALITY=1
else
  echo "[SKIP] XM+XC BAM exists: $bam_xc"
fi

########################################
# 4) Filter reads failing barcode quality
########################################

if [ ! -s "$bam_filt" ]; then
  echo "[RUN] FilterBam reject XQ"
  FilterBam INPUT="$bam_xc" OUTPUT="$bam_filt" TAG_REJECT=XQ
else
  echo "[SKIP] Filtered BAM exists: $bam_filt"
fi

echo "DONE:"
ls -lh "$ubam" "$bam_filt"

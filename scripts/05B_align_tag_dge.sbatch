#!/usr/bin/env bash
#SBATCH -J align_dge
#SBATCH -p medium
#SBATCH -c 16
#SBATCH --mem=80G
#SBATCH -t 5:00:00
#SBATCH -o ../logs/%x_%A_%a.out
#SBATCH -e ../logs/%x_%A_%a.err
#SBATCH --array=1-4

set -euo pipefail

source ../scripts/env_conda.sh

unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS

export TMPDIR="$(readlink -f ../tmp)"
mkdir -p "$TMPDIR"
echo "TMPDIR=$TMPDIR"
df -h "$TMPDIR" || true

SAMPLES=../scripts/samples.tsv
REF=../ref

STAR_DIR="${REF}/star_index"
REF_FASTA="${REF}/genome.fna"
REFFLAT="${REF}/dropseq/genes.refFlat"

IN_DIR=../results/05_tagged
OUT_BASE=../results/06_align_count
mkdir -p "$OUT_BASE"

line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES")
sample=$(echo "$line" | cut -f1)

echo "Sample: $sample"

in_bam="${IN_DIR}/${sample}.xm_xc.filtered.bam"
S_OUT="${OUT_BASE}/${sample}"
mkdir -p "$S_OUT"

echo "Input BAM:"
ls -lh "$in_bam"

echo "Reference check:"
ls -lh "$REF_FASTA" "$REFFLAT"
ls -lh "$STAR_DIR" | head -n 10 || true

# 1) Trim TSO
trim_tso="${S_OUT}/${sample}.trim_tso.bam"
if [ ! -s "$trim_tso" ]; then
  echo "[RUN] TrimStartingSequence"
  TrimStartingSequence \
    INPUT="$in_bam" \
    OUTPUT="$trim_tso" \
    OUTPUT_SUMMARY="${S_OUT}/${sample}.tso.summary.txt" \
    SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
    MISMATCHES=0 \
    NUM_BASES=5
else
  echo "[SKIP] trim_tso exists"
fi

# 2) Trim polyA
trim_polya="${S_OUT}/${sample}.trim_polya.bam"
if [ ! -s "$trim_polya" ]; then
  echo "[RUN] PolyATrimmer"
  PolyATrimmer \
    INPUT="$trim_tso" \
    OUTPUT="$trim_polya" \
    OUTPUT_SUMMARY="${S_OUT}/${sample}.polya.summary.txt" \
    MISMATCHES=0 \
    NUM_BASES=6
else
  echo "[SKIP] trim_polya exists"
fi

# 3) BAM -> FASTQ
cdna_fastq="${S_OUT}/${sample}.cdna.fastq.gz"

if [ ! -s "$cdna_fastq" ]; then
  echo "[RUN] picard SamToFastq"

  picard SamToFastq \
    I="$trim_polya" \
    F="${S_OUT}/${sample}.cdna.fastq" \
    TMP_DIR="$TMPDIR" \
    VALIDATION_STRINGENCY=SILENT

  if [ -f "${S_OUT}/${sample}.cdna.fastq" ]; then
    gzip -f "${S_OUT}/${sample}.cdna.fastq"
  fi
else
  echo "[SKIP] cdna_fastq already exists"
fi

# 4) STAR alignment
star_prefix="${S_OUT}/star_"
aligned_sam="${S_OUT}/star_Aligned.out.sam"
if [ ! -s "$aligned_sam" ]; then
  echo "[RUN] STAR"
  STAR \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --genomeDir "$STAR_DIR" \
    --readFilesIn "$cdna_fastq" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$star_prefix" \
    --outSAMtype SAM \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS nM
else
  echo "[SKIP] aligned_sam exists"
fi

# 5) Sort aligned SAM by queryname
aligned_qname="${S_OUT}/${sample}.aligned.queryname.bam"
if [ ! -s "$aligned_qname" ]; then
  echo "[RUN] picard SortSam"
  picard SortSam \
    I="$aligned_sam" \
    O="$aligned_qname" \
    SORT_ORDER=queryname \
    TMP_DIR="$TMPDIR"
else
  echo "[SKIP] aligned_qname exists"
fi

# 6) Merge aligned BAM with original tags (XM/XC)
merged_bam="${S_OUT}/${sample}.merged.bam"
if [ ! -s "$merged_bam" ]; then
  echo "[RUN] picard MergeBamAlignment"
  picard MergeBamAlignment \
    REFERENCE_SEQUENCE="$REF_FASTA" \
    UNMAPPED_BAM="$trim_polya" \
    ALIGNED_BAM="$aligned_qname" \
    OUTPUT="$merged_bam" \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false \
    TMP_DIR="$TMPDIR"
else
  echo "[SKIP] merged_bam exists"
fi

# 7) Tag reads with gene function
final_bam="${S_OUT}/${sample}.gene_function_tagged.bam"
if [ ! -s "$final_bam" ]; then
  echo "[RUN] TagReadWithGeneFunction"
  TagReadWithGeneFunction \
    INPUT="$merged_bam" \
    OUTPUT="$final_bam" \
    ANNOTATIONS_FILE="$REFFLAT"
else
  echo "[SKIP] final_bam exists"
fi

# 8) DigitalExpression
dge="${S_OUT}/${sample}.DGE_with_introns.txt.gz"
dge_sum="${S_OUT}/${sample}.DGE.summary.txt"
if [ ! -s "$dge" ]; then
  echo "[RUN] DigitalExpression"
  DigitalExpression \
    INPUT="$final_bam" \
    OUTPUT="$dge" \
    SUMMARY="$dge_sum" \
    READ_MQ=10 \
    EDIT_DISTANCE=1 \
    MIN_NUM_GENES_PER_CELL=100 \
    LOCUS_FUNCTION_LIST=CODING \
    LOCUS_FUNCTION_LIST=UTR \
    LOCUS_FUNCTION_LIST=INTRONIC
else
  echo "[SKIP] DGE exists"
fi

echo "DONE:"
ls -lh "$trim_tso" "$trim_polya" "$cdna_fastq" || true
ls -lh "$aligned_sam" "$aligned_qname" "$merged_bam" "$final_bam" "$dge" "$dge_sum" || true

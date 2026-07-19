#!/usr/bin/env bash
set -euo pipefail

IN_FILE=../results/08_relaxed_mapping/mapping_summary.tsv
OUT_FILE=../results/08_relaxed_mapping/mapping_average_summary.txt

if [ ! -s "$IN_FILE" ]; then
    echo "ERROR: mapping summary file not found: $IN_FILE"
    exit 1
fi

awk -F'\t' '
NR==1 { next }  # skip header

{
    gsub(/%/, "", $3)
    gsub(/%/, "", $4)
    gsub(/%/, "", $5)
    gsub(/%/, "", $6)
    gsub(/%/, "", $7)

    n++
    uniq += $3
    multi += $4
    too_many += $5
    too_short += $6
    other += $7
}

END {
    if (n == 0) {
        print "ERROR: No sample rows found in mapping summary."
        exit 1
    }

    avg_uniq = uniq / n
    avg_multi = multi / n
    avg_too_many = too_many / n
    avg_too_short = too_short / n
    avg_other = other / n

    printf "Average mapping statistics across %d scRNA-seq datasets mapped to the NCBI genome:\n", n
    printf "\n"
    printf "Average uniquely mapped reads: %.2f%%\n", avg_uniq
    printf "Average multi-mapped reads: %.2f%%\n", avg_multi
    printf "Average reads mapped to too many loci: %.2f%%\n", avg_too_many
    printf "Average unmapped reads (too short): %.2f%%\n", avg_too_short
    printf "Average unmapped reads (other): %.2f%%\n", avg_other
    printf "\n"
    printf "Collaborator-ready summary:\n"
    printf "I mapped all four scRNA-seq datasets to the NCBI genome using STAR. The average uniquely mapped read percentage across the four datasets was %.2f%%, with an average of %.2f%% multi-mapped reads.\n", avg_uniq, avg_multi
}' "$IN_FILE" | tee "$OUT_FILE"

echo
echo "Saved summary to:"
echo "$OUT_FILE"

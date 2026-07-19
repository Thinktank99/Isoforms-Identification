#!/usr/bin/env bash
set -euo pipefail

TABLE=../results/08_relaxed_mapping/mapping_summary_table.tsv
OUT=../results/08_relaxed_mapping/mapping_statistics.txt

awk -F'\t' '
NR>1 {

    gsub(/%/,"",$3)
    gsub(/%/,"",$6)

    n++
    uniq_sum += $3
    total_sum += $6

    uniq_sq += ($3*$3)
    total_sq += ($6*$6)

}

END {

    uniq_mean = uniq_sum/n
    total_mean = total_sum/n

    uniq_sd = sqrt((uniq_sq/n) - (uniq_mean*uniq_mean))
    total_sd = sqrt((total_sq/n) - (total_mean*total_mean))

    printf "Number of samples: %d\n\n", n

    printf "Uniquely mapped reads:\n"
    printf "Mean: %.2f%%\n", uniq_mean
    printf "SD:   %.2f%%\n\n", uniq_sd

    printf "Total mapped reads:\n"
    printf "Mean: %.2f%%\n", total_mean
    printf "SD:   %.2f%%\n\n", total_sd

    printf "Suggested report sentence:\n"
    printf "Across %d scRNA-seq datasets mapped to the NCBI genome, an average of %.2f%% (± %.2f%% SD) of reads mapped uniquely and %.2f%% (± %.2f%% SD) mapped overall.\n", n, uniq_mean, uniq_sd, total_mean, total_sd

}
' "$TABLE" | tee "$OUT"

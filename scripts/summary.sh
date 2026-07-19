#!/usr/bin/env bash
set -euo pipefail

IN_DIR=../results/08_relaxed_mapping
OUT_FILE="${IN_DIR}/mapping_summary.tsv"

echo -e "Sample\tInput_reads\tUniquely_mapped_percent\tMulti_mapped_percent\tToo_many_loci_percent\tUnmapped_too_short_percent\tUnmapped_other_percent" > "$OUT_FILE"

for f in "${IN_DIR}"/*_Log.final.out; do
    sample=$(basename "$f" _Log.final.out)

    input_reads=$(grep "Number of input reads" "$f" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
    uniq=$(grep "Uniquely mapped reads %" "$f" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
    multi=$(grep "% of reads mapped to multiple loci" "$f" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
    too_many=$(grep "% of reads mapped to too many loci" "$f" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
    too_short=$(grep "% of reads unmapped: too short" "$f" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
    other=$(grep "% of reads unmapped: other" "$f" | awk -F'|' '{gsub(/ /,"",$2); print $2}')

    echo -e "${sample}\t${input_reads}\t${uniq}\t${multi}\t${too_many}\t${too_short}\t${other}" >> "$OUT_FILE"
done

echo "Summary written to:"
echo "$OUT_FILE"
cat "$OUT_FILE"

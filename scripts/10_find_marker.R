suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 10_find_markers.R <input_rds> <output_tsv>")
}

input_rds <- args[1]
output_tsv <- args[2]

obj <- readRDS(input_rds)

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

markers <- markers[markers$p_val_adj <= 1e-5, ]

write.table(markers, file = output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

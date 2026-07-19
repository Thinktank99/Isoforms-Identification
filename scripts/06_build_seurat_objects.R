options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("Seurat")
install.packages("Matrix")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 06_build_seurat_objects.R <input_dir> <output_dir>")
}

input_dir <- args[1]
output_dir <- args[2]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dge_files <- list.files(input_dir, pattern = "DGE_with_introns\\.txt\\.gz$", recursive = TRUE, full.names = TRUE)

if (length(dge_files) == 0) {
  stop("No DGE files found.")
}

read_dge <- function(f) {
  message("Reading: ", f)
  x <- read.delim(gzfile(f), header = TRUE, row.names = 1, check.names = FALSE)
  x <- as.matrix(x)
  x <- Matrix(x, sparse = TRUE)
  return(x)
}

for (f in dge_files) {
  sample_name <- basename(dirname(f))
  mat <- read_dge(f)

  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 0,
    min.features = 0
  )

  obj$sample <- sample_name

  saveRDS(obj, file = file.path(output_dir, paste0(sample_name, ".seurat.raw.rds")))
  message("Saved: ", file.path(output_dir, paste0(sample_name, ".seurat.raw.rds")))
}

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 09_integrate_cluster.R <input_dir> <output_rds>")
}

input_dir <- args[1]
output_rds <- args[2]

files <- list.files(input_dir, pattern = "\\.qc_sct\\.rds$", full.names = TRUE)
if (length(files) < 2) {
  stop("Need at least 2 QC/SCT objects for integration.")
}

objs <- lapply(files, readRDS)

features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
objs <- PrepSCTIntegration(object.list = objs, anchor.features = features)

objs <- lapply(objs, function(x) RunPCA(x, features = features, npcs = 50, verbose = FALSE))

anchors <- FindIntegrationAnchors(
  object.list = objs,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  dims = 1:50
)

combined <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:50
)

combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
combined <- FindNeighbors(combined, dims = 1:50, k.param = 100)
combined <- FindClusters(combined, resolution = 1.2)

combined <- RunUMAP(
  combined,
  dims = 1:50,
  n.neighbors = 100,
  min.dist = 0.3,
  metric = "correlation",
  seed.use = 42
)

saveRDS(combined, output_rds)

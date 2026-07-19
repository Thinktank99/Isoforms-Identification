suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(glmGamPoi)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 08_qc_doubletfinder_sct.R <input_rds> <output_rds>")
}

input_rds <- args[1]
output_rds <- args[2]

obj <- readRDS(input_rds)

# Mito pattern may need adjustment depending on your annotation
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")

# If stage metadata exists, use stage-specific UMI minimums
umi_min <- 500
if ("stage" %in% colnames(obj@meta.data)) {
  st <- unique(obj$stage)
  if (length(st) == 1) {
    if (st == "stage7") umi_min <- 650
    if (st == "stage8.1") umi_min <- 700
    if (st == "stage9.1") umi_min <- 500
  }
}

# gene filtering: keep genes detected in at least 20 cells
counts <- GetAssayData(obj, slot = "counts")
keep_genes <- Matrix::rowSums(counts > 0) >= 20
obj <- subset(obj, features = rownames(obj)[keep_genes])

# cell filtering
obj <- subset(
  obj,
  subset =
    nFeature_RNA >= 400 &
    nFeature_RNA <= 1800 &
    nCount_RNA >= umi_min &
    nCount_RNA <= 4500 &
    percent.mt <= 1
)

# Normalize enough for DoubletFinder pre-processing
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)

# DoubletFinder parameter sweep
sweep.res <- paramSweep_v3(obj, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

best.pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
if (is.na(best.pk)) best.pk <- 0.005

nExp <- round(ncol(obj) * 0.05)

obj <- doubletFinder_v3(
  obj,
  PCs = 1:20,
  pN = 0.25,
  pK = best.pk,
  nExp = nExp,
  reuse.pANN = FALSE,
  sct = FALSE
)

df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
obj <- subset(obj, subset = obj@meta.data[[df_col]] == "Singlet")

# SCTransform exactly for downstream integration
obj <- SCTransform(
  obj,
  method = "glmGamPoi",
  vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"),
  verbose = TRUE
)

saveRDS(obj, output_rds)

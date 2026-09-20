library(Seurat)
library(tidyverse)
library(glmGamPoi) # The glmGamPoi package substantially improves speed and is used by default if installed
library(future) # Called from: getGlobalsAndPackages(expr, envir = envir, globals = globals)
options(future.globals.maxSize = 100000 * 1024^5)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

seurat_preprocess <- function(seu){
  # Normalization
  seu <- SCTransform(seu, assay = "RNA", 
                     variable.features.n = 3000,
                     return.only.var.genes = FALSE) 
  
  # Run PCA
  seu <- RunPCA(seu, assay = "SCT", verbose = F)
  # Neighbours 
  seu <- FindNeighbors(seu, reduction = "pca")
  # UMAP
  seu <- RunUMAP(seu, dims = 1:50)
  # Clusters 
  seu <- FindClusters(seu, resolution = c(seq(0.1, 1.1, 0.2), 1.3, 1.5, 2.0), alg = 1, graph.name = "SCT_snn")
  seu
}

inFile <- "20250509_output/Aggr_31Samples_Sep_2024_no_MT_RB.rds"
outDir <- "20250509_output/"

# 1. Load object
seu <- readRDS(inFile)

# Start timer for preprocessing
start_time_cca = Sys.time()
# 2. Pre-process object
seu <- seurat_preprocess(seu)
# Measure time for benchmark
end_time_cca = Sys.time()
# Time difference
measured_time_cca = end_time_cca - start_time_cca 

# 3. Save pre-processed object
saveRDS(seu, file = paste0(outDir, 'Aggr_31Samples_Sep_2024_pre-processed.rds'))

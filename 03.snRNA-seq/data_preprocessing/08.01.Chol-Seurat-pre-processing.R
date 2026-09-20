# Quick pre-processing of a Seurat object. 
library(Seurat)
library(tidyverse)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

seurat_preprocess <- function(seu){
  # Normalizarion
  seu <- SCTransform(seu, assay = "RNA", 
                     variable.features.n = 3000,
                     return.only.var.genes = FALSE) # variable.features.n set to default
  
  # Run PCA
  seu <- RunPCA(seu, assay = "SCT", verbose = F)
  return(seu)
}


inFile <- '20250509_output/Aggr_Sep2024_filt_All_harmony_th=0.rds'
outDir <- '20250509_output/Subset_Chol/'
dir.create(outDir, recursive = T)

# 1. Load object
seu <- readRDS(inFile)
m <- seu@meta.data
table(seu@meta.data$cell.annotation)

Idents(seu) <- m$cell.annotation

chol.seu <- subset(seu, idents = 'Cholangiocytes')
DefaultAssay(chol.seu) <- 'RNA'
chol.seu@assays$SCT <- NULL
colnames(chol.seu@meta.data)
chol.seu@meta.data <- chol.seu@meta.data[, -grep('SCT', colnames(chol.seu@meta.data))]

chol.seu <- seurat_preprocess(chol.seu)

saveRDS(chol.seu, file = paste0(outDir, 'Aggr_31Samples_Sep_2024_Chol_pre-processed.rds'))

# Quick pre-processing of a Seurat object. 
library(Seurat)
library(future) # Called from: getGlobalsAndPackages(expr, envir = envir, globals = globals)
options(future.globals.maxSize = 100000 * 1024^5)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

seurat_preprocess <- function(seu){
  # Normalizarion
  seu <- SCTransform(seu, assay = "RNA", 
                     variable.features.n = 3000,
                     return.only.var.genes = FALSE) # variable.features.n set to default
  
  # Scale Data (including all features)
  # seu <- ScaleData(seu, assay = "SCT", features=NULL, do.scale = TRUE, do.center = TRUE)
  
  # Run PCA
  seu <- RunPCA(seu, assay = "SCT", verbose = F)
  return(seu)
}


inFile <- '20250509_output/Aggr_Sep2024_filt_All_harmony_th=0.rds'
outDir <- '20250509_output/Subset_Hep/'
dir.create(outDir, recursive = T)

# 1. Load object
seu <- readRDS(inFile)
m <- seu@meta.data
table(seu@meta.data$cell.annotation)

# remove the 122 unknown celltype cells
Idents(seu) <- m$cell.annotation

hep.seu <- subset(seu, idents = 'Hepatocytes')
DefaultAssay(hep.seu) <- 'RNA'
hep.seu@assays$SCT <- NULL
colnames(hep.seu@meta.data)
hep.seu@meta.data <- hep.seu@meta.data[, -grep('SCT', colnames(hep.seu@meta.data))]

hep.seu <- seurat_preprocess(hep.seu)

saveRDS(hep.seu, file = paste0(outDir, 'Aggr_31Samples_Sep_2024_Hep_pre-processed.rds'))

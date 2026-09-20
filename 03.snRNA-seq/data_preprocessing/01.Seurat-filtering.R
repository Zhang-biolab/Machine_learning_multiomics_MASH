library(Seurat)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")


## Perform filtering of Seurat object or object list in terms of: nCount, nFeature, percent.mito, percent.ribo
## 1. Load data
inFile <- "20250509_output/Aggr_31Samples_Sep_2024_AddMetadata.rds"
outDir <- "20250509_output/"

## 1. Filter cells based on mitochondrial and ribosomal gene content
subset.cells <- function(seu, thr.mt, thr.rp, min.thr.nFeat.RNA, min.thr.nCount.RNA, max.thr.nCount.RNA){
  ## Filter Cells based on their mito and ribo content
  
  # Checks
  if(!(all(c("percent.mt.RNA", "percent.rp.RNA", "nFeature_RNA", "nCount_RNA") %in% names(seu@meta.data)))) stop ("Not all parameters present in Seurat object meta.data")
  ## Replace NAs with a value, if not they get filtered out
  seu@meta.data$percent.mt.RNA [ is.na(seu@meta.data$percent.mt.RNA) ] <- 0
  seu@meta.data$percent.rp.RNA [ is.na(seu@meta.data$percent.rp.RNA) ] <- 0
  
  print(paste0("Dimensions PRE-filtering: ", dim(seu)))
  
  seu <- subset(seu, percent.mt.RNA < thr.mt & percent.rp.RNA < thr.rp &
                  nFeature_RNA > min.thr.nFeat.RNA & nCount_RNA > min.thr.nCount.RNA & nCount_RNA < max.thr.nCount.RNA)
  
  print(paste0("Dimensions POST-filtering: ", dim(seu)))
  
  seu
}

# 0. params
mito_thres = 10
ribo_thres = 10
nFeat_min = 800
nCount_min = 1000
nCount_max = 50000

# 1. Load object
seu <- readRDS(inFile)
# calc rna and mt %
seu[["percent.mt.RNA"]] <- PercentageFeatureSet(seu, assay="RNA",  pattern = "^MT-")
seu[["percent.rp.RNA"]] <- PercentageFeatureSet(seu, assay = "RNA", pattern = "^RPS") + PercentageFeatureSet(seu, assay = "RNA", pattern = "^RPL")

# 1. Filter cells
seu <-  subset.cells(seu, 
                     thr.mt = mito_thres, 
                     thr.rp = ribo_thres, 
                     min.thr.nFeat.RNA = nFeat_min, 
                     min.thr.nCount.RNA = nCount_min, 
                     max.thr.nCount.RNA = nCount_max)

# 2. Save pre-processed object
saveRDS(seu, file = paste0(outDir, 'Aggr_31Samples_Sep_2024_filtered.rds'))

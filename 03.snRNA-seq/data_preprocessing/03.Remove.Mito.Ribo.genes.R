library(Seurat)
library(ggplot2)
library(cowplot)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

inFile <- "20250509_output/Aggr_31Samples_Sep_2024_filtered.rds"
outDir <- "20250509_output/"

## 1. Remove mitochondrial and ribosomal genes from cells
remove.ribo.mito.genes <- function(seu){
  ## Filter Cells based on their mito and ribo content
  
  print(paste0("Dimensions PRE-filtering: ", dim(seu)))
  
  seu <- seu[!grepl("^MT-", rownames(seu)) & !grepl("^RPS", rownames(seu)) & !grepl("^RPL", rownames(seu)), ]
  
  print(paste0("Dimensions POST-filtering: ", dim(seu)))
  
  seu
}

## 1. Load data
seu <- readRDS(inFile)
# 1. Filter mito ribo genes
seu <- remove.ribo.mito.genes(seu)
# [1] "Dimensions PRE-filtering: 30407" "Dimensions PRE-filtering: 50605"
# [1] "Dimensions POST-filtering: 30407" "Dimensions POST-filtering: 50605"
# 2. Save pre-processed object
saveRDS(seu, file = paste0(outDir, 'Aggr_31Samples_Sep_2024_no_MT_RB.rds'))

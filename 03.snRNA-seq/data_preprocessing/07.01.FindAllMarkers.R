library(Seurat)
library(data.table)
library(Matrix)
library(lme4)
library(dplyr)
library(stringr)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

inFile <- "20250509_output/Aggr_Sep2024_filt_All_harmony_th=0.rds"
outDir <- '20250509_output/DEGsRes/'
dir.create(outDir, recursive = T)


# summary plots of the harmony reduction
seu <- readRDS(inFile)

# do DGE of all cell types -------------------------
# DE parameters:
Idents(seu) <- 'cell.annotation'
DefaultAssay(seu) <- "SCT"

all.markers <- FindAllMarkers(seu,
                              assay = "SCT", 
                              logfc.threshold=0,
                              test.use='wilcox',
                              min.pct=0.1,
                              only.pos = FALSE,
                              base=2) # base = which to use for logs
fwrite(file=paste0(outDir, 'CelltypeSpecific_DE_table.csv'), 
       all.markers)
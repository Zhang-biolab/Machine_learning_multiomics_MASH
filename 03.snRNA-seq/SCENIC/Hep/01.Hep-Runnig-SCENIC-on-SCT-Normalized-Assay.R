library(Seurat)
library(tidyverse)
library(foreach)
library(RcisTarget)
library(doParallel)
library(SCopeLoomR)
library(AUCell)
library(doRNG)
library(GENIE3)
library(SCENIC)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")
set.seed(123)

print(Sys.time())
seuDir <- '20250509_output/Subset_Hep/'
cisTarget.dbDir <-  paste0(wd, '20250509_output/TF-regulon-analysis/1.SCENIC/cisTarget_databases/')
outDir <- '20250509_output/TF-regulon-analysis/Hep/1.SCENIC.SCTAssay.Hep.allCells/'
dir.create(outDir, recursive = T)


Hep.seu <- readRDS(paste0(seuDir, 'Aggr_Sep2024_filt_Hep_harmony_th=0.rds'))

Idents(Hep.seu) <- paste(Hep.seu@meta.data$Patient.ID, Hep.seu@meta.data$SCT_snn_harmony_t.0.0.05, sep = "_")
Hep.seu.downs <- subset(Hep.seu, downsample = 200) # Downsample cells

# ran SCENIC on an SCT-normalized assay
exprMat <- as.matrix(Hep.seu.downs@assays$SCT@data) # data: log1p(counts)=log(counts+1), 

cellInfo <- Hep.seu.downs@meta.data[, colnames(Hep.seu.downs@meta.data) %in% c('Patient.ID', 'Fibrosis.stage', 'cell.annotation', 'SCT_snn_harmony_t.0.0.05')]
colnames(cellInfo) <- c('Patient.ID',  'Fibrosis.stage', 'CellType', 'Cluster')
head(cellInfo)
table(cellInfo$CellType)

setwd(outDir)
dir.create("int/")
saveRDS(cellInfo, file="int/cellInfo.Rds")

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

mydbDir <- cisTarget.dbDir 
mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")

# S4 method for ScenicOptions (Object to store SCENIC settings)
scenicOptions <- initializeScenic(org = "hgnc", 
                                  dbDir = mydbDir, 
                                  nCores = 10,
                                  dbs = mydbs)

scenicOptions@inputDatasetInfo$cellInfo<-"int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# GENIE3
genesKept <- geneFiltering(exprMat, 
                           scenicOptions,
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat),
                           minSamples = ncol(exprMat) * 0.01)

exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions) 

# Optional: add log (if it is not logged/normalized already)
runGenie3(exprMat_filtered, scenicOptions) 

# SCENIC steps
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) 
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 

exprMat_all <- as.matrix(Hep.seu@assays$SCT@data)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_all)
# scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
print(Sys.time())

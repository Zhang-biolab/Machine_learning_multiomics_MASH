library(Seurat)
library(Matrix)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

inFile <- '20250509_output/Aggr_Sep2024_filt_All_harmony_th=0.rds'
outDir <- '20250509_output/Cell-Cell-Interaction/1.CellPhoneDB/AllCells/'
dir.create(outDir, recursive = T)

all.seu <- readRDS(inFile)

meta_data <- data.frame(cell=rownames(all.seu@meta.data), cell_type=all.seu@meta.data[, 'cell.annotation'])
meta_data[is.na(meta_data)] <- "Unkown"

write.table(meta_data, 
          file = paste0(outDir, 'AllCells_Mild_Advanced.celltype.CellphoneDB.metadata.txt'),
          row.names = FALSE, sep = '\t', quote = FALSE)

sct.norm.counts.matrix <- as.data.frame(all.seu@assays$SCT@data)
identical(colnames(sct.norm.counts.matrix), meta_data$cell)

sct.norm.counts.matrix <- cbind(rownames(sct.norm.counts.matrix), sct.norm.counts.matrix)
colnames(sct.norm.counts.matrix)[1] <- 'GeneSymbol'
fwrite(sct.norm.counts.matrix, 
       file = paste0(outDir, 'AllCells_Mild_Advanced.celltype.CellphoneDB.count.txt'), row.names = F, quote = FALSE, sep = '\t', nThread=8)

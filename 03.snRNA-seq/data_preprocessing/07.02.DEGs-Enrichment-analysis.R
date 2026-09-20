library(Seurat)
library(ggplot2)
library(data.table)
library(Matrix)
library(lme4)
library(cowplot)
library(dplyr)
library(stringr)
library(enrichplot)
library(gprofiler2)
setwd("~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/")

EnrichR.gseapy <- '~/projects/NASH_Fibrosis_datamining/Christopher_human_NASH_snRNA_datamining/scripts/figure_scripts/EnrichR_gseapy.py'

outDir <- '20250509_output/DEGsRes/'

all.markers <- read.csv(file=paste0(wd, outDir, 'CelltypeSpecific_DE_table.csv'), header = T, sep = ',') 

#1. Enrichment analysis for up regulated genes in cell specific way
all.markers.up <- subset(all.markers, p_val < 0.05 & all.markers$avg_log2FC > 0.25)

pathwaysDir <- paste0(wd, '20250509_output/DEGsRes/CelltypeSpecific_up_pathways')
dir.create(pathwaysDir, recursive = T)

for(i in unique(all.markers.up$cluster)){
  degs <- data.frame(gene = all.markers.up[all.markers.up$cluster == i, 'gene'])
  degsDir <- paste0(pathwaysDir, '/', i)
  dir.create(degsDir, recursive = T)
  write.table(degs, file = paste0(degsDir, '/gene'), row.names = FALSE, sep = '\t', quote = FALSE)
  print(dim(degs))
  system(paste0('cd ', degsDir))
  setwd(degsDir)
  system(paste0('~/anaconda3/bin/python ', EnrichR.gseapy))
}

#2. Enrichment analysis for down regulated genes in cell specific way
all.markers.down <- subset(all.markers, p_val < 0.05 & all.markers$avg_log2FC < -0.25)

pathwaysDir <- paste0(wd, '20250509_output/DEGsRes/CelltypeSpecific_down_pathways')
dir.create(pathwaysDir, recursive = T)
for(i in unique(all.markers.down$cluster)){
  degs <- data.frame(gene = all.markers.down[all.markers.down$cluster == i, 'gene'])
  degsDir <- paste0(pathwaysDir, '/', i)
  dir.create(degsDir, recursive = T)
  write.table(degs, file = paste0(degsDir, '/gene'), row.names = FALSE, sep = '\t', quote = FALSE)
  print(dim(degs))
  system(paste0('cd ', degsDir))
  setwd(degsDir)
  system(paste0('~/anaconda3/bin/python ', EnrichR.gseapy))
}


#3. Enrichment analysis for up regulated genes in advanced fibrosis
all.markers <- read.csv(file=paste0(wd, outDir, 'AllCelltype_DE_Advanced_vs_Mild_table.csv'), header = T, sep = ',') 

all.markers.up <- subset(all.markers, p_val < 0.05 & all.markers$avg_log2FC > 0.25)

pathwaysDir <- paste0(wd, '20250509_output/DEGsRes/CelltypeSpecific_Advanced_up_pathways')
dir.create(pathwaysDir, recursive = T)

for(i in unique(all.markers.up$cell.type)){
  degs <- data.frame(gene = all.markers.up[all.markers.up$cell.type == i, 'gene'])
  degsDir <- paste0(pathwaysDir, '/', i)
  dir.create(degsDir, recursive = T)
  write.table(degs, file = paste0(degsDir, '/gene'), row.names = FALSE, sep = '\t', quote = FALSE)
  print(dim(degs))
  system(paste0('cd ', degsDir))
  setwd(degsDir)
  system(paste0('~/anaconda3/bin/python ', EnrichR.gseapy))
}

#4. Enrichment analysis for down regulated genes in advanced fibrosis
all.markers.down <- subset(all.markers, p_val < 0.05 & all.markers$avg_log2FC < -0.25)

pathwaysDir <- paste0(wd, '20250509_output/DEGsRes/CelltypeSpecific_Advanced_down_pathways')
dir.create(pathwaysDir, recursive = T)
for(i in unique(all.markers.down$cell.type)){
  degs <- data.frame(gene = all.markers.down[all.markers.down$cell.type == i, 'gene'])
  degsDir <- paste0(pathwaysDir, '/', i)
  dir.create(degsDir, recursive = T)
  write.table(degs, file = paste0(degsDir, '/gene'), row.names = FALSE, sep = '\t', quote = FALSE)
  print(dim(degs))
  system(paste0('cd ', degsDir))
  setwd(degsDir)
  system(paste0('~/anaconda3/bin/python ', EnrichR.gseapy))
}